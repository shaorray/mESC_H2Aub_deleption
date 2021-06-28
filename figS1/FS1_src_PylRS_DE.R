setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
library(DESeq2)

# get gene specific DE in PylRS nAA treatment

filenames <- sort(list.files('../data/kallisto_output_E14_PylRS/', full.names = T)) # count with combined reference GENCODE vM20 and ncRNA annotation
sampleNewName <- gsub(".*/", "\\2", filenames)

count_table <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                  as = 'matrix', what = "est_counts")

colnames(count_table) <- sampleNewName

txCounts <- count_table[grepl("^ENS", rownames(count_table)), ] %>%
  keepOneTx(rowname_gene_id = F)

txCounts <- rbind(txCounts, count_table[!grepl("^ENS|^chrS", rownames(count_table)), ])

dds <- DESeqDataSetFromMatrix(round(as.matrix(txCounts)),
                              colData = data.frame(condition = gsub("_R(1|2)$", "\\2", sampleNewName) ),
                              design = ~ condition)
dds <- DESeq(dds)

res_CpK1 <- results(dds, contrast = c("condition", "E14_cl1_CPK", "E14_cl1_ctrl"))
res_CpK2 <- results(dds, contrast = c("condition", "E14_cl2_CPK", "E14_cl2_ctrl"))

CpK_DE_idx <- which((abs(res_CpK1$log2FoldChange) > 1.5 & res_CpK1$padj < 0.05) |
                      (abs(res_CpK2$log2FoldChange) > 1.5 & res_CpK2$padj < 0.05))

CpK_DE_tx <- rownames(txCounts)[CpK_DE_idx]

# save results
Gene_input <- importRanges("/mnt/0E471D453D8EE463/genomeDir/GENCODE/gencode.mouse.v1.annotation.gtf")
Gene_input$gene_id <- gsub("\\..*", "", Gene_input$gene_id)
Gene_input$transcript_id <- gsub("\\..*", "", Gene_input$transcript_id)

CpK_DE.gr <- Gene_input[Gene_input$transcript_id %in% 
                          Gene_input[Gene_input$transcript_id %in% CpK_DE_tx]$gene_id]
mcols(CpK_DE.gr) <- mcols(CpK_DE.gr)[, c("gene_id", "transcript_id")] 

CpK_DE2.gr <- CpK_DE_tx[!grepl("ENS", CpK_DE_tx)] %>%
  gsub(".*_(chr.*)", "\\1", .) %>%
  sapply(., strsplit, "_") %>% 
  Reduce(rbind, .) %>% 
  `colnames<-`(c("Seqname", "Start", "End", "Strand")) %>%
  as.data.frame() %>%
makeGRangesFromDataFrame()

mcols(CpK_DE2.gr) <- data.frame(gene_id = CpK_DE_tx[!grepl("ENS", CpK_DE_tx)],
                                transcript_id = NA)

CpK_DE.gr <- c(CpK_DE.gr, CpK_DE2.gr)
rm(CpK_DE2.gr)
names(CpK_DE.gr) <- NULL

idx <- which(rownames(res_CpK1) %in% CpK_DE.gr$gene_id)
CpK_DE.gr$CpK_log2FoldChange <- (res_CpK1[CpK_DE.gr$gene_id, "log2FoldChange"] +
                                   res_CpK2[CpK_DE.gr$gene_id, "log2FoldChange"]) / 2
CpK_DE.gr$CpK_padj <- (res_CpK1[CpK_DE.gr$gene_id, "padj"] +
                         res_CpK2[CpK_DE.gr$gene_id, "padj"]) / 2

export.gff3(CpK_DE.gr, "../data/PylRS_CpK_DE_genes_mm9.gtf")


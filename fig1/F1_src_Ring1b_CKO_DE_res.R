res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                       keys = as.character(GSE132753_DESeq2_res$Gene_name),
                       keytype = "GENENAME",
                       columns = "GENEID")

#-------------------- Ring1b CKO assay, Tamburri et al., GSE134053 --------------------------#
# provided results
GSE134053_DE_res <- read.table("../data/Ring1b_CKO_DE_res/Tamburri/B1_GSE134053_Ring1_OHT-vs-Ring1_ETA_diffexp_log2fc1.5_pval0.05.tsv",
                  header = 1)

GSE134053_DE_res$gene_id <- res$GENEID[match(GSE134053_DE_res$Geneid, res$GENENAME)]
GSE134053_DE_res <- GSE134053_DE_res[!is.na(GSE134053_DE_res$gene_id) &
                                       !duplicated(GSE134053_DE_res$gene_id), ]
rownames(GSE134053_DE_res) <- GSE134053_DE_res$gene_id

# rerun kallisto
filenames <- sort(list.files('/mnt/0E471D453D8EE463/GEO_Ring1b_RNA_seq/output', full.names = T)) # count with combined reference GENCODE vM20 and ncRNA annotation
sampleNewName <- gsub(".*/", "\\2", filenames)

count_table <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                  as = 'matrix', what = "est_counts")

colnames(count_table) <- sampleNewName

txCounts <- count_table[grepl("^ENS", rownames(count_table)), ] %>%
  keepOneTx(rowname_gene_id = T)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(round(as.matrix(txCounts)),
                              colData = data.frame(condition = gsub("_(rep1|rep2)$", "\\2", sampleNewName) ),
                              design = ~ condition)
dds <- DESeq(dds)

# res_Eed_WT_ETA <- results(dds, contrast = c("condition", "Eed_ETA", "WT_ETA"))
# res_I53S_WT_ETA <- results(dds, contrast = c("condition", "I53S_ETA", "WT_ETA"))
# res_Ring1_WT_ETA <- results(dds, contrast = c("condition", "Ring1_ETA", "WT_ETA"))

# res_Eed_OHT <- results(dds, contrast = c("condition", "Eed_OHT", "Eed_ETA"))
res_I53S_OHT <- results(dds, contrast = c("condition", "I53S_OHT", "I53S_ETA"))
res_Ring1_OHT <- results(dds, contrast = c("condition", "Ring1_OHT", "Ring1_ETA"))

# res_Eed_WT_OHT <- results(dds, contrast = c("condition", "Eed_OHT", "WT_OHT"))
res_I53S_WT_OHT <- results(dds, contrast = c("condition", "I53S_OHT", "WT_OHT"))
res_Ring1_WT_OHT <- results(dds, contrast = c("condition", "Ring1_OHT", "WT_OHT"))

GSE134053_DE_res.list <- list(res_I53S_OHT, res_Ring1_OHT, res_I53S_WT_OHT, res_Ring1_WT_OHT)
names(GSE134053_DE_res.list) <- c("I53S_OHT", "Ring1b_OHT", "I53S_WT_OHT", "Ring1b_WT_OHT")

rm(res_I53S_OHT, res_Ring1_OHT, res_I53S_WT_OHT, res_Ring1_WT_OHT)

#-------------------- Ring1b CPM assay, Blackledge et al., GSE132753 --------------------------#
GSE132753_DESeq2_res <- read.table("../data/Ring1b_CKO_DE_res/GSE132753_DESeq2_Output_PRC1CPM_PRC1CKO_RNAseq.txt",
                                   header = 1)

GSE132753_DESeq2_res$gene_id <- res$GENEID[match(GSE132753_DESeq2_res$Gene_name, res$GENENAME)]
GSE132753_DESeq2_res <- GSE132753_DESeq2_res[!is.na(GSE132753_DESeq2_res$gene_id) &
                                               !duplicated(GSE132753_DESeq2_res$gene_id), ]
rownames(GSE132753_DESeq2_res) <- GSE132753_DESeq2_res$gene_id

# GSE132753_DESeq2_res <- GSE132753_DESeq2_res[, c(2, 9, 8, 11, 10)]

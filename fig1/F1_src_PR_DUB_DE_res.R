# -------------------------- PR-DUB OE assay, Chen et al. 2021 -------------------------- #

# Chen_res_table <- read.table("../data/GSE153530_FPKM_all.xlsx", sep = "\t")

filenames <- sort(list.files('../data/kallisto_output_PR_DUB', full.names = T)) # count with combined reference GENCODE vM20 and ncRNA annotation
sampleNewName <- gsub(".*/", "\\2", filenames)


# get read counts for DESeq2
count_table_RNA_rc <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                         as = 'matrix', what = "est_counts")
colnames(count_table_RNA_rc) <- sampleNewName

geneRC_RNA <- count_table_RNA_rc[grepl("^ENS", rownames(count_table_RNA_rc)), ] %>%
  keepOneTx(rowname_gene_id = T, is_gene_sum = T)


library(DESeq2)
# genes read count
dds_geneRC <- DESeqDataSetFromMatrix(round(as.matrix(geneRC_RNA)),
                                     colData = data.frame(condition = gsub("_Rep.*", "", colnames(geneRC_RNA))),
                                     design = ~ condition)
dds_geneRC <- DESeq(dds_geneRC)

res.gene.list <- list(two_cell = results(dds_geneRC, contrast = c("condition", "RNAseq_2cell_PR-DUB_(OE)", "RNAseq_2cell_PR-DUB_(CTR)")),
                      four_cell = results(dds_geneRC, contrast = c("condition", "RNAseq_4cell_PR-DUB_(OE)", "RNAseq_4cell_PR-DUB_(CTR)"))  )


# etract DEGs
PR_DUB_DE_gene.list <- with(res.gene.list,
                     list("two_cell_up" = rownames(two_cell)[two_cell$baseMean > 5 & two_cell$log2FoldChange > 1 & two_cell$padj < 0.05],
                          "two_cell_down" = rownames(two_cell)[two_cell$baseMean > 5 & two_cell$log2FoldChange < (-1) & two_cell$padj < 0.05],
                          "four_cell_up" = rownames(four_cell)[four_cell$baseMean > 5 & four_cell$log2FoldChange > 1 & four_cell$padj < 0.05],
                          "four_cell_down" = rownames(four_cell)[four_cell$baseMean > 5 & four_cell$log2FoldChange < (-1) & four_cell$padj < 0.05] 
                          ))

PR_DUB_DE_gene.list <- lapply(PR_DUB_DE_gene.list, function(x) x[!is.na(x)])

saveRDS(PR_DUB_DE_gene.list, "data/PR_DUB_DE_gene.list.RData")


# intersection
intersect_two_sets(derepressed_gene_ids, unique(c(PR_DUB_DE_gene.list$two_cell_up, PR_DUB_DE_gene.list$four_cell_up)))



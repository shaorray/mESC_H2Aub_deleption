# -------------------------- BAP1 KO assay, Fursova et al. 2021 -------------------------- #

DESeq_res_files <- list.files("../data/GEO_BAP1_CKO_res/2021_Fursova/", pattern = "GEO.txt",
                              full.names = TRUE)

Fursova_res_table <- read.table(DESeq_res_files[1], header = 1)[, c("Gene_name", "log2FoldChange", "padj")]

res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                       keys = as.character(na.omit(Fursova_res_table$Gene_name)),
                       keytype = "GENENAME",
                       columns = "GENEID")
Fursova_res_table$gene_id <- res$GENEID[match(Fursova_res_table$Gene_name, res$GENENAME)]

Fursova_res_list <- list("up" = Fursova_res_table$gene_id[which(Fursova_res_table$log2FoldChange > 1 & Fursova_res_table$padj < 0.05)],
                         "down" = Fursova_res_table$gene_id[which(Fursova_res_table$log2FoldChange < -1 & Fursova_res_table$padj < 0.05)])

# -------------------------- BAP1 KO assay, Conway et al. 2021 -------------------------- #

Conway_res_table <- read.table("../data/GEO_BAP1_CKO_res/2021_Conway/GSE162739_Bap1KO_EV_DMSO-vs-WT_EV_DMSO_diffexp_log2fc1_pval0.05_fpkm0.tsv",
                               header = 1)
Conway_res_table$ENSEMBL <- gsub("\\..*", "", Conway_res_table$ENSEMBL)

Conway_res_list <- list("up" = Conway_res_table$ENSEMBL[which(Conway_res_table$log2FoldChange > 1 & Conway_res_table$padj < 0.05)],
                        "down" = Conway_res_table$ENSEMBL[which(Conway_res_table$log2FoldChange < -1 & Conway_res_table$padj < 0.05)])

# -------------------------- BAP1 KO assay, Kolovos et al. 2020 ------------------------- #

# rerun kallisto
filenames <- sort(list.files('../data/GEO_BAP1_CKO_res/2020_Kolovos/RNA_seq_kallisto_output', full.names = T)) # count with combined reference GENCODE vM20 and ncRNA annotation
sampleNewName <- gsub(".*/", "\\2", filenames)

count_table <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                  as = 'matrix', what = "est_counts")
colnames(count_table) <- sampleNewName

txCounts <- count_table[grepl("^ENS", rownames(count_table)), ] %>%
  keepOneTx(rowname_gene_id = T, is_gene_sum = T)

dds <- DESeqDataSetFromMatrix(round(as.matrix(txCounts)),
                              colData = data.frame(condition = gsub("_(rep1|rep2).*", "", sampleNewName) ),
                              design = ~ condition)
dds <- DESeq(dds)

res_BAP1_KO_WT <- results(dds, contrast = c("condition", "RNA-seq_in_BAP1_KO_mESCs", "RNA-seq_in_wt_mESCs"))
res_BAP1_rescued_WT <- results(dds, contrast = c("condition", "RNA-seq_in_BAP1_KO_rescued_with_wt_BAP1_mESCs", "RNA-seq_in_wt_mESCs"))

res_BAP1_KO_WT$gene_id <- rownames(res_BAP1_KO_WT)
Kolovos_res_table <- res_BAP1_KO_WT

Kolovos_res_list <- list("up" = Kolovos_res_table$gene_id[which(Kolovos_res_table$log2FoldChange > 1 & Kolovos_res_table$padj < 0.05)],
                         "down" = Kolovos_res_table$gene_id[which(Kolovos_res_table$log2FoldChange < -1 & Kolovos_res_table$padj < 0.05)])

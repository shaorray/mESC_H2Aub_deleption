res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                       keys = as.character(GSE132753_DESeq2_res$Gene_name),
                       keytype = "GENENAME",
                       columns = "GENEID")

#--------------------------------- Ring1b CKO assay, Tamburri et al., GSE134053 ----------------------------------#
# provided results
GSE134053_DE_res <- read.table("../data/GEO_Ring1b_CKO_res/Tamburri/B1_GSE134053_Ring1_OHT-vs-Ring1_ETA_diffexp_log2fc1.5_pval0.05.tsv",
                  header = 1) 
# Ring1_OHT-vs-Ring1_ETA: Parental (OHT vs -OHT) CKO
# WT_OHT-vs-WT_ETA: Ring1b WT (OHT vs -OHT) rescue
# Ring1_OHT-vs-WT_OHT: Ring1b CKO vs rescue

# WT_OHT: +OHT; [Ring1b WT]; Ring1b CKO rescue;
# WT_ETA: -OHT; [Ring1b WT]

# Ring1_OHT: +OHT; []; Ring1b CKO; 
# Ring1_ETA: -OHT; []

# I53S_OHT: +OHT; [I53S]; 
# I53S_ETA: -OHT; []

GSE134053_DE_res$gene_id <- res$GENEID[match(GSE134053_DE_res$Geneid, res$GENENAME)]
GSE134053_DE_res <- GSE134053_DE_res[!is.na(GSE134053_DE_res$gene_id) &
                                       !duplicated(GSE134053_DE_res$gene_id), ]
rownames(GSE134053_DE_res) <- GSE134053_DE_res$gene_id

# rerun kallisto
if (TRUE) {
  
  filenames <- sort(list.files('/mnt/0E471D453D8EE463/GEO_Ring1b_RNA_seq/output', full.names = T)) # count with combined reference GENCODE vM20 and ncRNA annotation
  filenames <- filenames[!grepl("E14_", filenames) & grepl("ETA|OHT", filenames)]
  sampleNewName <- gsub(".*/", "\\2", filenames)
  
  count_table <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                    as = 'matrix', what = "est_counts")
  colnames(count_table) <- sampleNewName
  
  txCounts <- count_table[grepl("^ENS", rownames(count_table)), ] %>%
    keepOneTx(rowname_gene_id = T, is_gene_sum = T)
  
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
  res_Ring1_OHT <- results(dds, contrast = c("condition", "Ring1_OHT", "Ring1_ETA")) # Ring1b CKO
  
  # res_Eed_WT_OHT <- results(dds, contrast = c("condition", "Eed_OHT", "WT_OHT"))
  res_I53S_WT_OHT <- results(dds, contrast = c("condition", "I53S_OHT", "WT_OHT"))
  res_Ring1_WT_OHT <- results(dds, contrast = c("condition", "Ring1_OHT", "WT_OHT"))
  
  GSE134053_DE_res.list <- list(res_I53S_OHT, res_Ring1_OHT, res_I53S_WT_OHT, res_Ring1_WT_OHT)
  names(GSE134053_DE_res.list) <- c("I53S_OHT", "Ring1b_OHT", "I53S_WT_OHT", "Ring1b_WT_OHT")
  
  rm(res_I53S_OHT, res_Ring1_OHT, res_I53S_WT_OHT, res_Ring1_WT_OHT)
}

#--------------------------------- Ring1b CPM assay, Blackledge et al., GSE132753 --------------------------------#
GSE132753_DESeq2_res <- read.table("../data/GEO_Ring1b_CKO_res/GSE132753_DESeq2_Output_PRC1CPM_PRC1CKO_RNAseq.txt",
                                   header = 1)

GSE132753_DESeq2_res$gene_id <- res$GENEID[match(GSE132753_DESeq2_res$Gene_name, res$GENENAME)]
GSE132753_DESeq2_res <- GSE132753_DESeq2_res[!is.na(GSE132753_DESeq2_res$gene_id) &
                                               !duplicated(GSE132753_DESeq2_res$gene_id), ]
rownames(GSE132753_DESeq2_res) <- GSE132753_DESeq2_res$gene_id

# GSE132753_DESeq2_res <- GSE132753_DESeq2_res[, c(2, 9, 8, 11, 10)]

# ----------------------------- Ring1b cKO, Dobrinic et al., GSE159399 ------------------------------------------ #
GSE159399_DESeq2_res_2h <- read.table("../data/GEO_Ring1b_CKO_res/Dobrinic/GSE159399_RING1AKO.RING1BAID_spikenormalised_DESeq2_NucRNAseq_IAA_2h_vs_UNT.txt",
                                   header = 1)
GSE159399_DESeq2_res_4h <- read.table("../data/GEO_Ring1b_CKO_res/Dobrinic/GSE159399_RING1AKO.RING1BAID_spikenormalised_DESeq2_NucRNAseq_IAA_4h_vs_UNT.txt",
                                      header = 1)
GSE159399_DESeq2_res_8h <- read.table("../data/GEO_Ring1b_CKO_res/Dobrinic/GSE159399_RING1AKO.RING1BAID_spikenormalised_DESeq2_NucRNAseq_IAA_8h_vs_UNT.txt",
                                      header = 1)
GSE159399_DESeq2_res_24h <- read.table("../data/GEO_Ring1b_CKO_res/Dobrinic/GSE159399_RING1AKO.RING1BAID_spikenormalised_DESeq2_NucRNAseq_IAA_24h_vs_UNT.txt",
                                      header = 1)

res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                       keys = as.character(na.omit(GSE159399_DESeq2_res_2h$geneSymbol)),
                       keytype = "GENENAME",
                       columns = "GENEID")

GSE159399_DE_tab <- rbind(cbind("2h",  GSE159399_DESeq2_res_2h$geneSymbol[which(GSE159399_DESeq2_res_2h$log2FoldChange > 1.5 & 
                                                                            GSE159399_DESeq2_res_2h$padj < 0.05)] %>% as.character()),
                          cbind("4h", GSE159399_DESeq2_res_2h$geneSymbol[which(GSE159399_DESeq2_res_4h$log2FoldChange > 1.5 & 
                                                                           GSE159399_DESeq2_res_4h$padj < 0.05)] %>% as.character()),
                          cbind("8h", GSE159399_DESeq2_res_2h$geneSymbol[which(GSE159399_DESeq2_res_8h$log2FoldChange > 1.5 & 
                                                                           GSE159399_DESeq2_res_8h$padj < 0.05)] %>% as.character())
                          )
GSE159399_DE_tab <- GSE159399_DE_tab %>% as.data.frame() %>% `colnames<-`(., c('time', "gene_name"))
GSE159399_DE_tab$gene_id <- res$GENEID[match(GSE159399_DE_tab$gene_name, res$GENENAME)]


# ------------------------------------------- Ring1b CKO up-regulated genes ------------------------------------ #
# Klose
GSE132753_DE_genes_CKO <- rownames(GSE132753_DESeq2_res)[which(GSE132753_DESeq2_res$LFC_PRC1CKO > 1.5 & 
                                                                 GSE132753_DESeq2_res$padj_PRC1CKO < 0.05)]
# Pasini
GSE134053_DE_genes_CKO <- rownames(GSE134053_DE_res.list$I53S_OHT)[which(GSE134053_DE_res.list$Ring1b_OHT$log2FoldChange > 1.5 &
                                                                           GSE134053_DE_res.list$Ring1b_OHT$padj < 0.05) ]

Ring1b_CKO_genes <- unique(c(GSE132753_DE_genes_CKO, GSE134053_DE_genes_CKO))

intersect_two_sets(derepressed_gene_ids, GSE132753_DE_genes_CKO)
intersect_two_sets(derepressed_gene_ids, GSE134053_DE_genes_CKO)
intersect_two_sets(derepressed_gene_ids, unique(GSE159399_DE_tab$gene_id))


# ------------------------------------------- Ring1b CPM up-regulated genes ------------------------------------ #
# Klose
GSE132753_DE_genes_CPM <- rownames(GSE132753_DESeq2_res)[which(GSE132753_DESeq2_res$LFC_PRC1CPM > 1.5 & 
                                                                 GSE132753_DESeq2_res$padj_PRC1CPM < 0.05)]
# Pasini
GSE134053_DE_genes_CPM <- rownames(GSE134053_DE_res.list$I53S_OHT)[which(GSE134053_DE_res.list$I53S_OHT$log2FoldChange > 1.5 &
                                                                           GSE134053_DE_res.list$I53S_OHT$padj < 0.05) ]


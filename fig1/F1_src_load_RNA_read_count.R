
# mm9 genes
if (FALSE) {
  gene.gr <- importRanges("/mnt/0E471D453D8EE463/genomeDir/GENCODE/gencode.mouse.v1.annotation.gtf")
  gene.gr <- gene.gr[gene.gr$type == "gene" & gene.gr$gene_type %in% c("protein_coding", "lincRNA")]
  
  gene.gr$gene_id <- gsub("\\..*$", "", gene.gr$gene_id)
  gene.gr$transcript_id <- gsub("\\..*$", "", gene.gr$transcript_id)

  gene.gr <- gene.gr[!is.na(gene.gr$gene_id)]
  names(gene.gr) <- gene.gr$gene_id
  
  saveRDS(gene.gr, "../data/gene.gr.RData")
}

# ------------------------------------------------------------------------------------------- #
# --------------------------------------- RNA-seq ------------------------------------------- #
# ------------------------------------------------------------------------------------------- #

filenames <- sort(list.files('../data/kallisto_output_BGI/', full.names = TRUE)) # count with combined reference GENCODE vM20 and ncRNA annotation
sampleNewName <- gsub(".*/", "\\2", filenames)

.design <- c("P0_R1", "P0_R2", "P0_R3",
             "P12_R1", "P12_R2", "P12_R3",
             "P12_C12_R1", "P12_C12_R2", "P12_C12_R3",
             "P12_C24_R1", "P12_C24_R2", "P12_C24_R3",
             "P12_C36_R1", "P12_C36_R2", "P12_C36_R3")
names(.design) <- 
  sampleNewName[order(gsub("(^.).*", "\\1", sampleNewName), as.numeric(gsub("^.", "", sampleNewName)))]
  
# get read counts for DESeq2
count_table_RNA_rc <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                         as = 'matrix', what = "est_counts")
colnames(count_table_RNA_rc) <- sampleNewName

txRC_RNA <- count_table_RNA_rc[grepl("^ENS", rownames(count_table_RNA_rc)), names(.design)] 

geneRC_RNA <- count_table_RNA_rc[grepl("^ENS", rownames(count_table_RNA_rc)), names(.design)] %>%
  keepOneTx(rowname_gene_id = T, is_gene_sum = TRUE)

flyRC_RNA <- count_table_RNA_rc[grep("^FBtr", rownames(count_table_RNA_rc)), ]


if (TRUE) {
  txRC_RNA.norm <- sweep(txRC_RNA, 2, size_factor_cal(flyRC_RNA[rowMeans(flyRC_RNA) > 1, ]), "/")
  data.frame(size = size_factor_cal(txRC_RNA.norm[rowMeans(txRC_RNA.norm) > 1, ]),
             sample = gsub("_R.", "", .design))[grep("^P(0|1)", .design), ] %>% 
    dplyr::group_by(sample) %>% 
    dplyr::summarise(mean = mean(size)) %>% 
    dplyr::mutate(mean = mean / mean[1], 
                  sample = factor(sample, levels = sample)) %>% 
    ggplot(aes(x = sample, y = mean, fill = sample)) + 
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = c("grey25", "salmon", "salmon2", "coral2", "brown1", "darksalmon")) +
    xlab("") + ylab("Relative size factors") +
    ggtitle("Total RNA-seq sizes (n = 14032)") +
    ggpubr::theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  ggsave("FigS1_RNAseq_sample_sizes.png", path = "../figS1/figs", 
         width = 4, height = 4)
}

if (FALSE) {
  count_table_RNA_tpm <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                           as = 'matrix', what = "tpm")
  colnames(count_table_RNA_tpm) <- sampleNewName
  geneRPK_RNA <- count_table_RNA_tpm[grepl("^ENS", rownames(count_table_RNA_tpm)), names(.design)] %>%
    keepOneTx(rowname_gene_id = TRUE, is_gene_sum = TRUE)
  flyRPK_RNA <- count_table_RNA_tpm[grep("^FBtr", rownames(count_table_RNA_tpm)), ]
  geneRPK_RNA.norm <- sweep(geneRPK_RNA, 2, size_factor_cal(flyRPK_RNA[rowMeans(flyRPK_RNA) > 1, ]), "/")
  colnames(geneRPK_RNA.norm) <- unname(.design)
  saveRDS(geneRPK_RNA.norm, "../fig6/data/geneRPK_RNA.norm.rds")
}

# kallisto est read count DE ---------------------------------------------------------------------------- #
library(DESeq2)

if (FALSE) {
  # transcripts read count
  dds_txRC <- DESeqDataSetFromMatrix(round(as.matrix(txRC_RNA)),
                                     colData = data.frame(condition = gsub("_R.", "", .design)),
                                     design = ~ condition)
  sizeFactors(dds_txRC) <- size_factor_cal(flyRC_RNA[rowMeans(flyRC_RNA) > 1, ])
  dds_txRC <- DESeq(dds_txRC)
  
  res.tx.list <- list(P12 = results(dds_txRC, contrast = c("condition", "P12", "P0")),
                      P12_C12 = results(dds_txRC, contrast = c("condition", "P12_C12", "P0")),
                      P12_C24 = results(dds_txRC, contrast = c("condition", "P12_C24", "P0")),
                      P12_C36 = results(dds_txRC, contrast = c("condition", "P12_C36", "P0")))
  saveRDS(res.tx.list, "data/RNAseq_res.tx.list.RData")
}

# genes read count
dds_geneRC <- DESeqDataSetFromMatrix(round(as.matrix(geneRC_RNA)),
                                     colData = data.frame(condition = gsub("_R.", "", .design)),
                                     design = ~ condition)
sizeFactors(dds_geneRC) <- SizeFactorCal(flyRPK_RNA[rowMeans(flyRPK_RNA) > 1, ])
dds_geneRC <- DESeq(dds_geneRC)

res.gene.list <- list(P12 = results(dds_geneRC, contrast = c("condition", "P12", "P0")),
                      P12_C12 = results(dds_geneRC, contrast = c("condition", "P12_C12", "P0")),
                      P12_C24 = results(dds_geneRC, contrast = c("condition", "P12_C24", "P0")),
                      P12_C36 = results(dds_geneRC, contrast = c("condition", "P12_C36", "P0"))
)

saveRDS(res.gene.list, "data/RNAseq_res.gene.list.RData")

# list of DE genes ------------------------------------------------------------------------- #
# kallisto transcripts

filter_DE <- function(res, min_baseMean, DE_threshold, DE_direction, padj_threshold) {
  if (DE_direction == "up") {
    res$baseMean > min_baseMean & res$log2FoldChange > DE_threshold & res$padj < padj_threshold
  } else {
    res$baseMean > min_baseMean & res$log2FoldChange < -DE_threshold & res$padj < padj_threshold
  }
}

min_baseMean <- 5
DE_threshold <- 1.5

DE_gene.list <- with(res.gene.list,
                     list("P12_up" = rownames(P12)[filter_DE(P12, min_baseMean, DE_threshold, "up", 0.05)],
                     "P12_down" = rownames(P12)[filter_DE(P12, min_baseMean, DE_threshold, "down", 0.05)],
                     "P12_C12_up" = rownames(P12_C12)[filter_DE(P12_C12, min_baseMean, DE_threshold, "up", 0.05)],
                     "P12_C12_down" = rownames(P12_C12)[filter_DE(P12, min_baseMean, DE_threshold, "down", 0.05)],
                     "P12_C24_up" = rownames(P12_C24)[filter_DE(P12_C24, min_baseMean, DE_threshold, "up", 0.05)],
                     "P12_C24_down" = rownames(P12_C24)[filter_DE(P12_C24, min_baseMean, DE_threshold, "down", 0.05)],
                     "P12_C36_up" = rownames(P12_C36)[filter_DE(P12_C36, min_baseMean, DE_threshold, "up", 0.05)],
                     "P12_C36_down" = rownames(P12_C36)[filter_DE(P12_C36, min_baseMean, DE_threshold, "down", 0.05)]))

DE_gene.list <- lapply(DE_gene.list, function(x) x[!is.na(x)])
lengths(DE_gene.list)

saveRDS(DE_gene.list, "../fig1/data/DE_gene.list.RData")


# log2FC shrinkage ------------------------------------------------------------------------- #

res.gene.shrink.list <- list(P12 = lfcShrink(dds_geneRC, contrast = c("condition", "P12", "P0")),
                      P12_C12 = lfcShrink(dds_geneRC, contrast = c("condition", "P12_C12", "P0")),
                      P12_C24 = lfcShrink(dds_geneRC, contrast = c("condition", "P12_C24", "P0")),
                      P12_C36 = lfcShrink(dds_geneRC, contrast = c("condition", "P12_C36", "P0"))
)

DE_gene.shrink.list <- with(res.gene.shrink.list,
                     list("P12_up" = rownames(P12)[filter_DE(P12, min_baseMean, DE_threshold, "up", 0.05)],
                          "P12_down" = rownames(P12)[filter_DE(P12, min_baseMean, DE_threshold, "down", 0.05)],
                          "P12_C12_up" = rownames(P12_C12)[filter_DE(P12_C12, min_baseMean, DE_threshold, "up", 0.05)],
                          "P12_C12_down" = rownames(P12_C12)[filter_DE(P12, min_baseMean, DE_threshold, "down", 0.05)],
                          "P12_C24_up" = rownames(P12_C24)[filter_DE(P12_C24, min_baseMean, DE_threshold, "up", 0.05)],
                          "P12_C24_down" = rownames(P12_C24)[filter_DE(P12_C24, min_baseMean, DE_threshold, "down", 0.05)],
                          "P12_C36_up" = rownames(P12_C36)[filter_DE(P12_C36, min_baseMean, DE_threshold, "up", 0.05)],
                          "P12_C36_down" = rownames(P12_C36)[filter_DE(P12_C36, min_baseMean, DE_threshold, "down", 0.05)]))

DE_gene.shrink.list <- lapply(DE_gene.shrink.list, function(x) x[!is.na(x)])
lengths(DE_gene.shrink.list)

# ---------------------------------- bam read counts DE ------------------------------------ #
if (FALSE) {
  Gene_input <- importRanges("/mnt/0E471D453D8EE463/genomeDir/GENCODE/gencode.mouse.v1.annotation.gtf")
  Gene_input$gene_id <- gsub("\\..*", "", Gene_input$gene_id)
  exon.gr <- Gene_input[Gene_input$type == "exon" & 
                          Gene_input$gene_type %in% c("protein_coding", "lincRNA"), ]
  bam_files <- list.files("/mnt/E0767589767560E8/UPPMAX/BGI/mESC/bam_mm9", "P.*bam$", full.names = T)
  
  require(Rsubread)
  exon_RC_RNA <- featureCounts(bam_files,
                               annot.ext = data.frame(GeneID = exon.gr$gene_id,
                                                      Chr = seqnames(exon.gr),
                                                      Start = start(exon.gr),
                                                      End = end(exon.gr),
                                                      Strand = strand(exon.gr)),
                               isPairedEnd = TRUE, nthreads = 8)$counts 
  
  dds_exonRC <- DESeqDataSetFromMatrix(round(as.matrix(exon_RC_RNA)),
                                       colData = data.frame(condition = gsub("(P.).*", "\\1", colnames(exon_RC_RNA))),
                                       design = ~ condition)
  # sizeFactors(dds_exonRC) <- SizeFactorCal(flyRPK_RNA[rowMeans(flyRPK_RNA) > 1, ])
  dds_exonRC <- DESeq(dds_exonRC)
  
  res.exon.list <- list(P12 = results(dds_exonRC, contrast = c("condition", "P2", "P1")),
                        P12_C12 = results(dds_exonRC, contrast = c("condition", "P3", "P1")),
                        P12_C24 = results(dds_exonRC, contrast = c("condition", "P4", "P1")),
                        P12_C36 = results(dds_exonRC, contrast = c("condition", "P5", "P1")))
  saveRDS(res.exon.list, "../fig1/data/RNAseq_res.exon.list.RData")
  
  

  gencode.gr <- Gene_input[Gene_input$type == "gene" & 
                             Gene_input$gene_type %in% c("protein_coding", "lincRNA"), ]
  
  gene_RC_RNA <- featureCounts(bam_files,
                               annot.ext = data.frame(GeneID = gencode.gr$gene_id,
                                                      Chr = seqnames(gencode.gr),
                                                      Start = start(gencode.gr),
                                                      End = end(gencode.gr),
                                                      Strand = strand(gencode.gr)),
                               isPairedEnd = TRUE, nthreads = 8)$counts 
  
}

# --------------------------------------- not in use -------------------------------------------- #

if (FALSE) {
  # bam -> featureCounts gencode read count DE
  dds_gencode <- DESeqDataSetFromMatrix(round(as.matrix(count_table_rnaseq.norm)),
                                        colData = data.frame(condition = gsub("_R.", "", .design)),
                                        design = ~ condition)
  sizeFactors(dds_gencode) <- 1 # read counts are spike-in normalized
  dds_gencode <- DESeq(dds_gencode)
  
  res_P12 <- results(dds_gencode, contrast = c("condition", "P12", "P0"))
  res_P12_C12 <- results(dds_gencode, contrast = c("condition", "P12_C12", "P0"))
  res_P12_C24 <- results(dds_gencode, contrast = c("condition", "P12_C24", "P0"))
  res_P12_C36 <- results(dds_gencode, contrast = c("condition", "P12_C36", "P0"))
  
  res_gencode.list <- list(P12 = res_P12,
                           P12_C12 = res_P12_C12,
                           P12_C24 = res_P12_C24,
                           P12_C36 = res_P12_C36)
  saveRDS(res_gencode.list, "data/RNAseq.res_gencode.list.RData")
  
  # list of DE genes -------------------------------------------------------------------------
  # gencode exons
  DE_gencode.list <- with(res_gencode.list,
                          list("P12_up" = rownames(res_P12)[res_P12$log2FoldChange > 1 & res_P12$padj < 0.05],
                               "P12_down" = rownames(res_P12)[res_P12$log2FoldChange < (-1) & res_P12$padj < 0.05],
                               "P12_C12_up" = rownames(res_P12_C12)[res_P12_C12$log2FoldChange > 1 & res_P12_C12$padj < 0.05],
                               "P12_C12_down" = rownames(res_P12_C12)[res_P12_C12$log2FoldChange < (-1) & res_P12_C12$padj < 0.05],
                               "P12_C24_up" = rownames(res_P12_C24)[res_P12_C24$log2FoldChange > 1 & res_P12_C24$padj < 0.05],
                               "P12_C24_down" = rownames(res_P12_C24)[res_P12_C24$log2FoldChange < (-1) & res_P12_C24$padj < 0.05],
                               "P12_C36_up" = rownames(res_P12_C36)[res_P12_C36$log2FoldChange > 1 & res_P12_C36$padj < 0.05],
                               "P12_C36_down" = rownames(res_P12_C36)[res_P12_C36$log2FoldChange < (-1) & res_P12_C36$padj < 0.05]))
  
  DE_gencode.list <- lapply(DE_gencode.list, function(x) x[!is.na(x)])
  
  saveRDS(DE_gencode.list, "data/DE_gencode.list.RData")
  
  library(ggpubr)
  res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                         keys = as.character(rownames(res_P12)),
                         keytype = "GENEID",
                         columns = "GENENAME")
  
  g1 <- ggpubr::ggmaplot(res_P12[res$GENEID, ], genenames = res$GENENAME, 
                         top = 50, select.top.method = "fc", ggtheme = theme_pubclean(), 
                         fc = 1.5, fdr =  0.05, size = 0.5, ylim = c(-15, 15), 
                         main = "P12 vs P0")
  
  g2 <- ggpubr::ggmaplot(res_P12_C12[res$GENEID, ], genenames = res$GENENAME, 
                         top = 50, select.top.method = "fc", ggtheme = theme_pubclean(),
                         fc = 1.5, fdr =  0.05, size = 0.5, ylim = c(-15, 15), 
                         main = "P12_C12 vs P0")
  
  g3 <- ggpubr::ggmaplot(res_P12_C24[res$GENEID, ], genenames = res$GENENAME, 
                         top = 50, select.top.method = "fc", ggtheme = theme_pubclean(),
                         fc = 1.5, fdr =  0.05, size = 0.5, ylim = c(-15, 15), 
                         main = "P12_C24 vs P0")
  
  g4 <- ggpubr::ggmaplot(res_P12_C36[res$GENEID, ], genenames = res$GENENAME, 
                         top = 50, select.top.method = "fc", ggtheme = theme_pubclean(),
                         fc = 1.5, fdr =  0.05, size = 0.5, ylim = c(-15, 15), 
                         main = "P12_C36 vs P0")
  
  ggsave(grid.arrange(g1, g2, g3, g4, nrow = 1),
         filename = "RNAseq_gencode_MAplot_BAP1_pulse_chase.png",
         width = 16, height = 5, path = "figs")
}



if (FALSE) {
  # load ChIP-seq data
  bw_files <- list.files("/mnt/E0767589767560E8/UPPMAX/Ezh2i_BAP1_ChIP",
                         ".bw$", full.names = T)
  sample_names <- gsub(".*ChIP\\/(.*)", "\\1", bw_files)
  
  sample_table <- data.frame(Target = gsub("(_BAP1|_Ezh2i|-NTD|_NT|_Trp).*", "", sample_names),
                             Treatment = gsub("_RC.fltd.bw", "", sample_names))
  
  source("../util/getCoverage.R")
  count_table_ChIP <- .countBW(bw_files = bw_files, 
                               intervals = promoters(gene.gr, upstream = 2000, downstream = 2000), 
                               fast = T)
  count_table_ChIP <- count_table_ChIP[match(rownames(count_table_ttseq),
                                             res$GENEID[match(rownames(count_table_ChIP), 
                                                              res$ENTREZID)]), ]
  
  sample_ChIP_sizes <- t(matrix(size_factor_cal(count_table_ChIP), nrow = 9))
  colnames(sample_ChIP_sizes) <- sample_table$Treatment[1:9]
  rownames(sample_ChIP_sizes) <- unique(sample_table$Target)
  
  pheatmap::pheatmap(sample_ChIP_sizes[, -8] / sample_ChIP_sizes[, 8],
                     cluster_rows = F, cluster_cols = F,
                     main = "Relative changes to control")
}

if (FALSE) {
  # correlate ChIP changes
  log2FC_list <- list(log2(count_table_ttseq.norm[, c(3,2,6,5,4,7,8)] / count_table_ttseq.norm[, 1]))
  for (i in unique(sample_table$Target)) {
    tmp_counts <- count_table_ChIP[, sample_table$Target == i]
    tmp_FC <- (log2(tmp_counts[, c(3,1,4,2,5,6,7)] / tmp_counts[, 8]))
    log2FC_list <- c(log2FC_list, list(tmp_FC))
  }
  
  Treatment = c("BAP1_6h", "BAP1_12h", "BAP1_6h_Ezh2_1d", "BAP1_12h_Ezh2_1d", 
                "Ezh2_1d", "Ezh2_2d", "Ezh2_7d")
  
  g_list <- list()
  for (i in 1:8) {
    for (j in 1:7) {
      dat <- data.frame(Tx_log2FC = log2FC_list[[1]][, j],
                        ChIP_log2FC = log2FC_list[[i+1]][, j])
      dat <- dat[apply(dat, 1, function(x) all(!is.infinite(x) & !is.na(x))), ]
      
      g_list <- c(g_list, 
                  list(ggplot(dat, aes(x = ChIP_log2FC, y = Tx_log2FC)) +
                         geom_hex(bins = 50) +
                         scale_fill_gradient(low="grey75", high=colors_9[i]) +
                         ggtitle(paste(unique(sample_table$Target)[i], Treatment[j])) +
                         theme_minimal() +
                         theme(legend.position = "none")))
    }
  }
  
  ggsave(do.call(grid.arrange, c(g_list, nrow = 4)),
         filename = "cor_all_Tx_ChIP.png",
         width = 40, height = 10, path = "figs")
  
  # estimate kinetics by bins
  sample_bw_bins <- binBwBam(file_names = bw_files)
  colnames(sample_bw_bins) <- paste(sample_table$Target, sample_table$Treatment, sep = "_")
  
  estimate_kinetics <- function(dat, point) {
    point = cbind(1, point)
    dat = log(dat)
    (solve(t(point) %*% point) %*% t(point) %*% t(dat))[2, ]
  }
  
  bin_H2Aub_kinetics <- estimate_kinetics(sample_bw_bins[, c(8, 3, 1)], c(0, 6, 12))
  bin_H3K27m3_kinetics <- estimate_kinetics(sample_bw_bins[, c(17, 14:15)], c(0, 1, 2))
  
  # estimate kinetics by genes
  gene_H2Aub_kinetics <- estimate_kinetics(count_table_ChIP[, c(8, 3, 1)], c(0, 6, 12))
  gene_H3K27m3_kinetics <- estimate_kinetics(count_table_ChIP[, c(17, 14:15)], c(0, 1, 2))
}


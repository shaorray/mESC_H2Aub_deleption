source("background_calibration.R")

# load annotations
gene.gr <- readRDS("../mESC_Ezh2i_BAP1/data/gene.gr.RData")

# house keeping genes
res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                       keys = gsub("\\..*", "", as.character(gene.gr$gene_id)),
                       keytype = "GENEID",
                       columns = "GENENAME")
houseKeeping.genes <- res$GENEID[match(c("Rer1", "Rpl13a", "Hprt1", "Rpl27", "Tbp", "Gapdh"), res$GENENAME)]
houseKeeping.genes <- houseKeeping.genes[!is.na(houseKeeping.genes)]

# highly expressed genes
total_RPK <- readRDS("../mESC_Ezh2i_BAP1/fig1/data/txRPK_RNA.norm.RData")

high_expr.genes <- rownames(total_RPK)[head(order(rowMeans(total_RPK), decreasing = T), 100)]
high_expr.genes <- high_expr.genes[high_expr.genes %in% gsub("\\..*", "", as.character(gene.gr$gene_id))]


# load file paths
chip.path <- "bam/Rybp_P0_rep1.mm9.bam"
input.path <- "bam/sfINPUT_P0_rep1.mm9.bam"
input.path <- "../../bam/Rybp_P0_rep2.mm9.bam"


chip.paths <- list.files("/mnt/E0767589767560E8/UPPMAX/PHILIP_LOF3_20210205R/bam",
                         pattern = "Ring1b.*rep.*bam$", full.names = T)

input.paths <- list.files("/mnt/E0767589767560E8/UPPMAX/PHILIP_LOF3_20210205R/bam",
                          pattern = "INPUT.*rep.*bam$", full.names = T)


res.list <- list()
for (i in seq_along(chip.paths)) {
  res.list <- c(res.list, 
                list(NCIS.run(chip.path = chip.paths[i], 
                               input.path = input.paths[i], 
                               max.binsize = 10000,
                               quant = 0.75)))
}

source("../mESC_Ezh2i_BAP1/util/utils.R")
source("../mESC_Ezh2i_BAP1/util/getCoverage.R")

bam_files_2 <- list.files(path = "/mnt/E0767589767560E8/UPPMAX/PHILIP_LOF3_20210205R/bam",
                          pattern = "(Ring1b|INPUT).*rep.*mm9.bam$",
                          full.names = T)

ChIP_TSS_B2_BAM_mat <- .countBam(bam_files = bam_files_2,
                                 intervals = promoters(gene.gr, upstream = 2000, downstream = 2000),
                                 stranded = F, 
                                 paired.end = 'ignore')
colnames(ChIP_TSS_B2_BAM_mat) <- gsub(".mm9.bam", "", colnames(ChIP_TSS_B2_BAM_mat))

input_sf <- SizeFactorCal(ChIP_TSS_B2_BAM_mat[, grep("INPUT", colnames(ChIP_TSS_B2_BAM_mat))])
Ring1b_input_norm <- sweep(ChIP_TSS_B2_BAM_mat[, grep("Ring1b", colnames(ChIP_TSS_B2_BAM_mat))], 
                           2, input_sf, "/")

background_sf <- unlist(lapply(res.list, function(x) x$est))

Ring1b_back_norm <- sweep(ChIP_TSS_B2_BAM_mat[, grep("Ring1b", colnames(ChIP_TSS_B2_BAM_mat))], 
                          2, background_sf * input_sf, "/")


PcG_enriched_genes <- readRDS("../mESC_Ezh2i_BAP1/PcG_enriched_genes.rds")

library(dplyr)
library(ggplot2)

g1 <- Ring1b_input_norm %>%
  reshape::melt() %>%
  dplyr::mutate(condition = gsub("Ring1b_", "", variable)) %>%
  dplyr::mutate(.time = gsub("_rep.*", "", condition)) %>%
  ggplot(aes(x = condition, y = log1p(value), fill = .time)) +
  geom_boxplot(outlier.colour = NA) + ylim(c(2, 6)) +
  xlab("") + ylab("log1p TSS reads") + 
  ggtitle("All genes (n = 23640) Input normalization") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

g2 <- Ring1b_back_norm %>%
  reshape::melt() %>%
  dplyr::mutate(condition = gsub("Ring1b_", "", variable)) %>%
  dplyr::mutate(.time = gsub("_rep.*", "", condition)) %>%
  ggplot(aes(x = condition, y = log1p(value), fill = .time)) +
  geom_boxplot(outlier.colour = NA) + ylim(c(2, 6)) +
  xlab("") + ylab("log1p TSS reads") + 
  ggtitle("All genes (n = 23640) Background normalization") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

g3 <- Ring1b_input_norm[PcG_enriched_genes, ] %>%
  reshape::melt() %>%
  dplyr::mutate(condition = gsub("Ring1b_", "", variable)) %>%
  dplyr::mutate(.time = gsub("_rep.*", "", condition)) %>%
  ggplot(aes(x = condition, y = log1p(value), fill = .time)) +
  geom_boxplot(outlier.colour = NA) + ylim(c(2, 8)) +
  xlab("") + ylab("log1p TSS reads") + 
  ggtitle("PcG enriched genes (n = 861) Input normalization") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

g4 <- Ring1b_back_norm[PcG_enriched_genes, ] %>%
  reshape::melt() %>%
  dplyr::mutate(condition = gsub("Ring1b_", "", variable)) %>%
  dplyr::mutate(.time = gsub("_rep.*", "", condition)) %>%
  ggplot(aes(x = condition, y = log1p(value), fill = .time)) +
  geom_boxplot(outlier.colour = NA) + ylim(c(2, 8)) +
  xlab("") + ylab("log1p TSS reads") + 
  ggtitle("PcG enriched genes (n = 861) Background normalization") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(grid.arrange(g1, g2, g3, g4, nrow = 2),
       filename = "BAP1_PC_B2_Ring1b_normalization.png",
       path = "./figs", 
       height = 8, width = 12)

# ------------------------------------------------------------------------------------
chip.paths <- list.files("/mnt/E0767589767560E8/UPPMAX/PHILIP_LOF3_20210205R/bam",
                         pattern = "H3K27.*rep.*bam$", full.names = T)

H3K27me3.res.list <- list()
for (i in seq_along(chip.paths)) {
  H3K27me3.res.list <- c(H3K27me3.res.list, 
                         list(NCIS.run(chip.path = chip.paths[i], 
                                        input.path = input.paths[i], 
                                        max.binsize = 10000,
                                        quant = 0.75)))
}


H3K27me3_TSS_B2_BAM_mat <- .countBam(bam_files = list.files(path = "/mnt/E0767589767560E8/UPPMAX/PHILIP_LOF3_20210205R/bam",
                                                        pattern = "H3K27me3.*rep.*mm9.bam$",
                                                        full.names = T),
                                 intervals = promoters(gene.gr, upstream = 2000, downstream = 2000),
                                 stranded = F, 
                                 paired.end = 'ignore')
colnames(H3K27me3_TSS_B2_BAM_mat) <- gsub(".mm9.bam", "", colnames(H3K27me3_TSS_B2_BAM_mat))


H3K27me3_input_norm <- sweep(H3K27me3_TSS_B2_BAM_mat, 2, input_sf, "/")

background_sf <- unlist(lapply(H3K27me3.res.list, function(x) x$est))

H3K27me3_back_norm <- sweep(H3K27me3_TSS_B2_BAM_mat, 2, background_sf * input_sf, "/")


g1 <- H3K27me3_input_norm %>%
  reshape::melt() %>%
  dplyr::mutate(condition = gsub("H3K27me3_", "", variable)) %>%
  dplyr::mutate(.time = gsub("_rep.*", "", condition)) %>%
  ggplot(aes(x = condition, y = log1p(value), fill = .time)) +
  geom_boxplot(outlier.colour = NA) + ylim(c(1, 5)) +
  xlab("") + ylab("log1p TSS reads") + 
  ggtitle("All genes (n = 23640) Input normalization") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

g2 <- H3K27me3_back_norm %>%
  reshape::melt() %>%
  dplyr::mutate(condition = gsub("H3K27me3_", "", variable)) %>%
  dplyr::mutate(.time = gsub("_rep.*", "", condition)) %>%
  ggplot(aes(x = condition, y = log1p(value), fill = .time)) +
  geom_boxplot(outlier.colour = NA) + ylim(c(1, 5)) +
  xlab("") + ylab("log1p TSS reads") + 
  ggtitle("All genes (n = 23640) Background normalization") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

g3 <- H3K27me3_input_norm[PcG_enriched_genes, ] %>%
  reshape::melt() %>%
  dplyr::mutate(condition = gsub("H3K27me3_", "", variable)) %>%
  dplyr::mutate(.time = gsub("_rep.*", "", condition)) %>%
  ggplot(aes(x = condition, y = log1p(value), fill = .time)) +
  geom_boxplot(outlier.colour = NA) + ylim(c(2, 6)) +
  xlab("") + ylab("log1p TSS reads") + 
  ggtitle("PcG enriched genes (n = 861) Input normalization") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

g4 <- H3K27me3_back_norm[PcG_enriched_genes, ] %>%
  reshape::melt() %>%
  dplyr::mutate(condition = gsub("H3K27me3_", "", variable)) %>%
  dplyr::mutate(.time = gsub("_rep.*", "", condition)) %>%
  ggplot(aes(x = condition, y = log1p(value), fill = .time)) +
  geom_boxplot(outlier.colour = NA) + ylim(c(2, 6)) +
  xlab("") + ylab("log1p TSS reads") + 
  ggtitle("PcG enriched genes (n = 861) Background normalization") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(grid.arrange(g1, g2, g3, g4, nrow = 2),
       filename = "BAP1_PC_B2_H3K27me3_normalization.png",
       path = "./figs", 
       height = 8, width = 12)


# --------------------------------------------------------------------------
chip.paths <- list.files("/mnt/E0767589767560E8/UPPMAX/PHILIP_LOF3_20210205R/bam",
                         pattern = "H2Aub.*rep.*bam$", full.names = T)

H2Aub.res.list <- list()
for (i in seq_along(chip.paths)) {
  H2Aub.res.list <- c(H2Aub.res.list, 
                         list(NCIS.run(chip.path = chip.paths[i], 
                                        input.path = input.paths[i], 
                                        min.binsize = 200,
                                        max.binsize = 1000,
                                        quant = 0.85)))
}


H2Aub_TSS_B2_BAM_mat <- .countBam(bam_files = list.files(path = "/mnt/E0767589767560E8/UPPMAX/PHILIP_LOF3_20210205R/bam",
                                                            pattern = "H2Aub.*rep.*mm9.bam$",
                                                            full.names = T),
                                     intervals = promoters(gene.gr, upstream = 2000, downstream = 2000),
                                     stranded = F, 
                                     paired.end = 'ignore')
colnames(H2Aub_TSS_B2_BAM_mat) <- gsub(".mm9.bam", "", colnames(H2Aub_TSS_B2_BAM_mat))


H2Aub_input_norm <- sweep(H2Aub_TSS_B2_BAM_mat, 2, input_sf, "/")

background_sf <- unlist(lapply(H2Aub.res.list, function(x) x$est))

H2Aub_back_norm <- sweep(H2Aub_TSS_B2_BAM_mat, 2, background_sf * input_sf, "/")


g1 <- H2Aub_input_norm %>%
  reshape::melt() %>%
  dplyr::mutate(condition = gsub("H2Aub_", "", variable)) %>%
  dplyr::mutate(.time = gsub("_rep.*", "", condition)) %>%
  ggplot(aes(x = condition, y = log1p(value), fill = .time)) +
  geom_boxplot(outlier.colour = NA) + ylim(c(1, 5)) +
  xlab("") + ylab("log1p TSS reads") + 
  ggtitle("All genes (n = 23640) Input normalization") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

g2 <- H2Aub_back_norm %>%
  reshape::melt() %>%
  dplyr::mutate(condition = gsub("H2Aub_", "", variable)) %>%
  dplyr::mutate(.time = gsub("_rep.*", "", condition)) %>%
  ggplot(aes(x = condition, y = log1p(value), fill = .time)) +
  geom_boxplot(outlier.colour = NA) + ylim(c(1, 5)) +
  xlab("") + ylab("log1p TSS reads") + 
  ggtitle("All genes (n = 23640) Background normalization") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

g3 <- H2Aub_input_norm[PcG_enriched_genes, ] %>%
  reshape::melt() %>%
  dplyr::mutate(condition = gsub("H2Aub_", "", variable)) %>%
  dplyr::mutate(.time = gsub("_rep.*", "", condition)) %>%
  ggplot(aes(x = condition, y = log1p(value), fill = .time)) +
  geom_boxplot(outlier.colour = NA) + ylim(c(2, 5.5)) +
  xlab("") + ylab("log1p TSS reads") + 
  ggtitle("PcG enriched genes (n = 861) Input normalization") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

g4 <- H2Aub_back_norm[PcG_enriched_genes, ] %>%
  reshape::melt() %>%
  dplyr::mutate(condition = gsub("H2Aub_", "", variable)) %>%
  dplyr::mutate(.time = gsub("_rep.*", "", condition)) %>%
  ggplot(aes(x = condition, y = log1p(value), fill = .time)) +
  geom_boxplot(outlier.colour = NA) + ylim(c(2, 5.5)) +
  xlab("") + ylab("log1p TSS reads") + 
  ggtitle("PcG enriched genes (n = 861) Background normalization") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(grid.arrange(g1, g2, g3, g4, nrow = 2),
       filename = "BAP1_PC_B2_H2Aub_normalization.png",
       path = "./figs", 
       height = 8, width = 12)

# ------------------------------------------------------------------------------------
# CALIBRATION data set
chip.paths <- list.files("/mnt/0E471D453D8EE463/GEO_calibration/", recursive = T,
                         pattern = ".*bam$", full.names = T)

input.paths <- chip.paths[grep("IN", chip.paths)]
chip.paths <- chip.paths[grep("H3K27me3", chip.paths)]

H3K27me3.res.list <- list()
for (i in seq_along(chip.paths)) {
  H3K27me3.res.list <- c(H3K27me3.res.list, 
                         list(NCIS.run(chip.path = chip.paths[i], 
                                        input.path = input.paths[i], 
                                        which.bg = NULL, 
                                        min.binsize = 1000,
                                        max.binsize = 2000,
                                        quant = 0.75)))
}

which.bg <- sortSeqlevels(H3K27me3.res.list[[1]]$bg.bin)
for (i in seq_along(H3K27me3.res.list)[-1]) {
  bg.tmp <- sortSeqlevels(H3K27me3.res.list[[i]]$bg.bin)
  seqlengths(bg.tmp) <- seqlengths(which.bg) <- pmax(seqlengths(bg.tmp), seqlengths(which.bg))
  which.bg <- GenomicRanges::intersect(which.bg, bg.tmp)
}

H3K27me3.res.list2 <- list()
for (i in seq_along(chip.paths)) {
  H3K27me3.res.list2 <- c(H3K27me3.res.list2, 
                         list(NCIS.run(chip.path = chip.paths[i], 
                                        input.path = input.paths[i], 
                                        which.bg = which.bg + 1000,
                                        min.binsize = 1000,
                                        max.binsize = 2000,
                                        quant = 0.75)))
}

# check the used background bin size
unlist(lapply(H3K27me3.res.list, function(x) sum(width(x$empty.bin)) / 1e6))
unlist(lapply(H3K27me3.res.list, function(x) sum(width(x$bg.bin)) / 1e6))

unlist(lapply(H3K27me3.res.list2, function(x) sum(width(x$bg.bin)) / 1e6))


H3K27me3_TSS_BAM_mat <- .countBam(bam_files = chip.paths,
                                     intervals = promoters(gene.gr, upstream = 2000, downstream = 2000),
                                     stranded = F, 
                                     paired.end = 'ignore')

INPUT_TSS_BAM_mat <- .countBam(bam_files = input.paths,
                               intervals = promoters(gene.gr, upstream = 2000, downstream = 2000),
                               stranded = F, 
                               paired.end = 'ignore')
colnames(H3K27me3_TSS_BAM_mat) <- colnames(INPUT_TSS_BAM_mat) <- gsub(".mm9.bam", "", colnames(H3K27me3_TSS_BAM_mat))
# H3K27me3_TSS_BAM_mat <- H3K27me3_TSS_BAM_mat[rowSums(H3K27me3_TSS_BAM_mat == 0) < 4, ]

input_sf <- SizeFactorCal(INPUT_TSS_BAM_mat)
H3K27me3_input_norm <- sweep(H3K27me3_TSS_BAM_mat, 2, input_sf, "/")

background_sf <- unlist(lapply(H3K27me3.res.list2, function(x) x$est))
H3K27me3_back_norm <- sweep(H3K27me3_TSS_BAM_mat, 2, background_sf * input_sf, "/")

colnames(H3K27me3_input_norm) <- colnames(H3K27me3_back_norm) <- paste0((gsub(".*\\/ESC_(H3K27me3.*)\\/bam.*", "\\1", chip.paths)), "_", 1:2)

dat <- data.frame(ratio = gsub(".*_.*_(.*)_.", "\\1", colnames(H3K27me3_input_norm)),
                  bg_sf = background_sf,
                  input_sf = input_sf,
                  
                  input_norm_rc = colSums(H3K27me3_input_norm),
                  input_norm_sf = SizeFactorCal(H3K27me3_input_norm),
                  bg_norm_rc = colSums(H3K27me3_back_norm),
                  bg_norm_sf = SizeFactorCal(H3K27me3_back_norm),
                  
                  input_norm_rc_PcG = colSums(H3K27me3_input_norm[PcG_enriched_genes, ]),
                  input_norm_sf_PcG = SizeFactorCal(H3K27me3_input_norm[PcG_enriched_genes, ]),
                  bg_norm_rc_PcG = colSums(H3K27me3_back_norm[PcG_enriched_genes, ]),
                  bg_norm_sf_PcG = SizeFactorCal(H3K27me3_back_norm[PcG_enriched_genes, ]),
                  
                  input_norm_rc_hk = colSums(H3K27me3_input_norm[houseKeeping.genes, ]),
                  input_norm_sf_hk = SizeFactorCal(H3K27me3_input_norm[houseKeeping.genes, ]),
                  bg_norm_rc_hk = colSums(H3K27me3_back_norm[houseKeeping.genes, ]),
                  bg_norm_sf_hk = SizeFactorCal(H3K27me3_back_norm[houseKeeping.genes, ]),
                  
                  input_norm_rc_he = colSums(H3K27me3_input_norm[high_expr.genes, ]),
                  input_norm_sf_he = SizeFactorCal(H3K27me3_input_norm[high_expr.genes, ]),
                  bg_norm_rc_he = colSums(H3K27me3_back_norm[high_expr.genes, ]),
                  bg_norm_sf_he = SizeFactorCal(H3K27me3_back_norm[high_expr.genes, ])
                  )

dat <- dat[!grepl("H3K27me3", dat$ratio), ]
dat$ratio <- as.numeric(as.character(dat$ratio))

library(ggrepel)


# --------- all genes read counts ---------
g1 <- ggplot(dat, aes(x = ratio, y = input_norm_rc / 1e6, label = as.character(ratio))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggrepel::geom_label_repel() +
  xlab("H3K27me3 ratio") + ylab("H3K27me3 read counts (million)") + ggtitle("[Input normalization] All genes (n=23640)") +
  theme_bw()

g2 <- ggplot(dat, aes(x = ratio, y = bg_norm_rc / 1e6, label = as.character(ratio))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggrepel::geom_label_repel() +
  xlab("H3K27me3 ratio") + ylab("H3K27me3 read counts (million)") + ggtitle("[Background normalization] All genes (n=23640)") +
  theme_bw()


# --------- PcG genes read counts ---------
g3 <- ggplot(dat, aes(x = ratio, y = input_norm_rc_PcG / 1e6, label = as.character(ratio))) +
  geom_point(color = "blue3") +
  geom_smooth(method = "lm") +
  ggrepel::geom_label_repel() +
  xlab("H3K27me3 ratio") + ylab("H3K27me3 read counts (million)") + ggtitle("[Input normalization] PcG genes (n=861)") +
  theme_bw()

g4 <- ggplot(dat, aes(x = ratio, y = bg_norm_rc_PcG / 1e6, label = as.character(ratio))) +
  geom_point(color = "blue3") +
  geom_smooth(method = "lm") +
  ggrepel::geom_label_repel() +
  xlab("H3K27me3 ratio") + ylab("H3K27me3 read counts (million)") + ggtitle("[Background normalization] PcG genes (n=861)") +
  theme_bw()



# --------- house keeping genes should be unchanged ---------
g5 <- ggplot(dat, aes(x = ratio, y = input_norm_rc_hk / 1e3, label = as.character(ratio))) +
  geom_point(color = "red3") +
  geom_smooth(method = "lm") +
  ggrepel::geom_label_repel() +
  xlab("H3K27me3 ratio") + ylab("H3K27me3 read counts (k)") + ggtitle("[Input normalization] Housekeeping genes (n=5)") +
  theme_bw()

g6 <- ggplot(dat, aes(x = ratio, y = bg_norm_rc_hk / 1e3, label = as.character(ratio))) +
  geom_point(color = "red3") +
  geom_smooth(method = "lm") +
  ggrepel::geom_label_repel() +
  xlab("H3K27me3 ratio") + ylab("H3K27me3 read counts (k)") + ggtitle("[Background normalization] Housekeeping genes (n=5)") +
  theme_bw()


# --------- highly expressed genes should be unchanged ---------
g7 <- ggplot(dat, aes(x = ratio, y = bg_norm_rc / 1e3 * ratio, label = as.character(ratio))) +
  geom_point(color = "green3") +
  geom_smooth(method = "lm") +
  ggrepel::geom_label_repel() +
  xlab("H3K27me3 ratio") + ylab("H3K27me3 read counts (k)") + ggtitle("[Input normalization] Highly expressed genes (n=84)") +
  theme_bw()

g8 <- ggplot(dat, aes(x = ratio, y = bg_norm_rc_he / 1e3, label = as.character(ratio))) +
  geom_point(color = "green3") +
  geom_smooth(method = "lm") +
  ggrepel::geom_label_repel() +
  xlab("H3K27me3 ratio") + ylab("H3K27me3 read counts (k)") + ggtitle("[Background normalization] Highly expressed genes (n=84)") +
  theme_bw()


# --------- adjusted background normalization ---------
g9 <- ggplot(dat, aes(x = ratio, y = bg_norm_rc / 1e6 * ratio, label = as.character(ratio))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggrepel::geom_label_repel() +
  xlab("H3K27me3 ratio") + ylab("H3K27me3 read counts (million)") + ggtitle("[Adj. Background normalization] All genes (n=23640)") +
  theme_bw()

g10 <- ggplot(dat, aes(x = ratio, y = bg_norm_rc_PcG / 1e6 * ratio, label = as.character(ratio))) +
  geom_point(color = "blue3") +
  geom_smooth(method = "lm") +
  ggrepel::geom_label_repel() +
  xlab("H3K27me3 ratio") + ylab("H3K27me3 read counts (million)") + ggtitle("[Adj. Background normalization] PcG genes (n=861)") +
  theme_bw()

g11 <- ggplot(dat, aes(x = ratio, y = bg_norm_rc_hk / 1e3 * ratio, label = as.character(ratio))) +
  geom_point(color = "red3") +
  geom_smooth(method = "lm") +
  ggrepel::geom_label_repel() +
  xlab("H3K27me3 ratio") + ylab("H3K27me3 read counts (k)") + ggtitle("[Adj. Background normalization] Housekeeping genes (n=5)") +
  theme_bw()

g12 <- ggplot(dat, aes(x = ratio, y = bg_norm_rc_he / 1e3 * ratio, label = as.character(ratio))) +
  geom_point(color = "green3") +
  geom_smooth(method = "lm") +
  ggrepel::geom_label_repel() +
  xlab("H3K27me3 ratio") + ylab("H3K27me3 read counts (k)") + ggtitle("[Adj. Background normalization] Highly expressed genes (n=84)") +
  theme_bw()


ggsave(grid.arrange(g1, g3, g5, g7, g2, g4, g6, g8, g9, g10, g11, g12, ncol = 4),
       filename = "Calibration_H3K27me3_normalization.png",
       path = "./figs", 
       height = 12, width = 16)

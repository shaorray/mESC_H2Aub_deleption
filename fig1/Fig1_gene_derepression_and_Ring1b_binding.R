# Rui Shao 2021 May
# Figure 1
# BAP1 pulse chase RNA-seq derepression

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
source("../util/getCoverage.R")

# -------------------read results from 'load_read_count.R'--------------------------- #
# 
txRC_RNA <- readRDS("data/txRC_RNA_seq.RData")
DE_gene.list <- readRDS("../data/DE_gene.list.RData")
gene.gr <- readRDS("../data/gene.gr.RData")
gene.gr <- gene.gr[seqnames(gene.gr) %in% paste0("chr", c(1:100, "X", "Y"))]
res.list <- readRDS("data/RNAseq_res.list.RData")

blacklist.gr <- importRanges("../data/blacklist_peak_mm9.bed")
CpK_DE.gr <- importRanges("../data/PylRS_CpK_DE_genes_mm9.gff3")

use_gene_ids <- intersect.Vector(gene.gr$gene_id, rownames(log2FC_mat))
use_gene_ids <- use_gene_ids[!use_gene_ids %in% CpK_DE.gr$gene_id[grep("ENS", CpK_DE.gr$gene_id)]]

# ----------------------------------------------------------------------------------- #
# rlog transformation
convert_rlog <- function(mat) {
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(round(as.matrix(mat)),
                                colData = data.frame(condition = colnames(mat)),
                                design = ~ condition)
  rld <- rlog(dds)
  rld_val <- assay(rld)
  return(rld_val)
}

# ----------------------------------------------------------------------------------- #
g1.1 <- plot_ma(baseMean = res.list$P12$baseMean,
        log2FC = res.list$P12$log2FoldChange, 
        pval = res.list$P12$padj, 
        p_col = 1, cut_baseMean = 0, 
        title = "P12")

g1.2 <- plot_ma(baseMean = res.list$P12_C12$baseMean,
        log2FC = res.list$P12_C12$log2FoldChange, 
        pval = res.list$P12_C12$padj, 
        p_col = 1, cut_baseMean = 0, 
        title = "P12C12")

g1.3 <- plot_ma(baseMean = res.list$P12_C24$baseMean,
        log2FC = res.list$P12_C24$log2FoldChange, 
        pval = res.list$P12_C24$padj, 
        p_col = 1, cut_baseMean = 0, 
        title = "P12C24")

g1.4 <- plot_ma(baseMean = res.list$P12_C36$baseMean,
        log2FC = res.list$P12_C36$log2FoldChange, 
        pval = res.list$P12_C36$padj, 
        p_col = 1, cut_baseMean = 0, 
        title = "P12C36")

ggsave(plot = grid.arrange(g1.1, g1.2, g1.3, g1.4, ncol = 4),
       filename = paste0("Fig1_MAplot_RNAseq_PC_DE.png"), 
       path = "figs/",
       device = "png", width = 16, height = 4)


# -------------------------- cluster log2FC from DESeq2 ----------------------------- #
order_cls <- function(cls, val) {
  factor(cls, levels = order(aggregate(val, list(cls), mean)[, 2], decreasing = T)) %>%
    as.numeric()
}
 
log2FC_mat <- sapply(res.list, function(x) x$log2FoldChange)
rownames(log2FC_mat) <- rownames(res.list[[1]])
log2FC_mat <- log2FC_mat[!is.na(rowSums(log2FC_mat)), ] %>% trim_quantile()
log2FC_mat <- log2FC_mat[use_gene_ids, ]
log2FC_mat <- sweep(log2FC_mat, 2, colMedians(log2FC_mat), "-")
colnames(log2FC_mat) <- c("P12", "P12_C12", "P12_C24", "P12_C36")
saveRDS(log2FC_mat, "data/log2FC_mat.RData")

# kmeans
set.seed(1)
log2FC_cls <- kmeans(log2FC_mat, centers = 4)$cluster
log2FC_cls <- order_cls(log2FC_cls, rowMeans(log2FC_mat))
names(log2FC_cls) <- rownames(log2FC_mat)

genes_derepressed <- intersect.Vector(rownames(log2FC_mat)[log2FC_cls == 1], use_gene_ids)
genes_W_derepressed <- intersect.Vector(rownames(log2FC_mat)[log2FC_cls == 2], use_gene_ids)
genes_unchanged <- intersect.Vector(rownames(log2FC_mat)[log2FC_cls == 3], use_gene_ids)
genes_repressed <- intersect.Vector(rownames(log2FC_mat)[log2FC_cls == 4], use_gene_ids)

# plot PCA
txRC_RNA_cmb <- sapply(unique(gsub("_R.*", "", .design)), 
                       function(sample)
                         rowMeans(txRC_RNA.norm[, names(.design[grep(paste0(sample, "_R"), .design)])]))
txRC_RNA.cmb <- txRC_RNA.cmb[, c("P0", "P12", "P12_C12", "P12_C24", "P12_C36")]
txRC_RNA.cmb.rlog <- convert_rlog(txRC_RNA.cmb)

plot_pca <- function(mat, gene_ids) {
  mat <- mat[intersect.Vector(rownames(mat), gene_ids), ]
  pca <- prcomp(t(mat))
  pca$sdev <- pca$sdev^2
  pc_var <- pca$sdev[1:2] / sum(pca$sdev) * 100
  pc <- pca$x[, 1:2] %>% as.data.frame()
  grid_dat <- `colnames<-`(cbind(pc[-5, ], pc[c(2:5), ]), c("x1", "y1", "x2", "y2"))
  pc$Sample <- rownames(pc)
  
  ggplot(pc, aes(x = PC2, y = PC1, color = Sample)) +
    geom_segment(data = grid_dat,
                 mapping = aes(x = y1, y = x1, xend = y2, yend = x2),
                 arrow = arrow(type="open", angle=20, length = unit(x = 14, 'pt')),
                 size = 0.6, lty = 1,
                 color = add.alpha("black", 0.5)) + 
    geom_point(size = 4) +
    scale_color_viridis_d(begin = 0.04, end = 0.84) +
    ylab(paste0("PC1 ", round(pc_var[1]), "% variance")) +
    xlab(paste0("PC2 ", round(pc_var[2]), "% variance")) +
    theme_setting
}

plot_pca(txRC_RNA.cmb.rlog, rownames(txRC_RNA.cmb.rlog))
ggsave(filename = paste0("Fig1_RNAseq_PCA_derepressed.png"), 
       path = "figs/",
       device = "png", width = 4.5, height = 5)

# overlap with PcG target peaks
Ring1b_Ezh2_peaks <- importRanges("../data/GSM3100852_Ring1b_EZH2_Peaks.mm9.bed")
gene_Ring1b <- gene.gr$gene_id[findOverlaps(gene.gr, Ring1b_Ezh2_peaks[Ring1b_Ezh2_peaks$name == "Ring1b"]) %>% queryHits()]
gene_Ezh2 <- gene.gr$gene_id[findOverlaps(gene.gr, Ring1b_Ezh2_peaks[Ring1b_Ezh2_peaks$name == "Ezh2"]) %>% queryHits()]
gene_Both <- gene.gr$gene_id[findOverlaps(gene.gr, Ring1b_Ezh2_peaks[Ring1b_Ezh2_peaks$name == "both"]) %>% queryHits()]


use_gene.gr <- gene.gr[use_gene_ids]
log2FC_cls_use <- log2FC_cls[use_gene_ids]
use_gene.gr$Type <- c("Derepressed", "Weakly\nDerepressed", "Unchanged", "Repressed")[log2FC_cls_use]

use_gene.gr$PcG <- "None"
use_gene.gr$PcG[queryHits(findOverlaps(use_gene.gr, Ring1b_Ezh2_peaks))] <- "Target"

PcG_counts <- mcols(use_gene.gr) %>% as.data.frame() %>%
  dplyr::count(Type, PcG, sort = TRUE)

PcG_counts$Type <- factor(PcG_counts$Type, 
                          levels = c("Derepressed", "Weakly\nDerepressed", "Unchanged", "Repressed"))
PcG_counts$PcG <- factor(PcG_counts$PcG,
                         levels = c("Target", "None"))

library("ggforce")
dat_ggforce <- PcG_counts  %>%
  ggforce::gather_set_data(1:2) %>%
  arrange(x, Type, desc(PcG))

ggplot(dat_ggforce1, aes(x = x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = Type), alpha = 0.8, axis.width = 0.2,
                     n=100, strength = 0.5) +
  geom_parallel_sets_axes(axis.width = 0.25, fill = "gray90",  size = 0) +
  geom_parallel_sets_labels(colour = 'gray35', size = 0, angle = 0, fontface="bold") +
  scale_fill_manual(values  = viridis(4, option = "E", direction = -1)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x  = element_blank()
  )

ggsave(filename = "FigS1_DE_group_PcG_overlaps1.png", 
       path = "figs", width = 2.5, height = 5)

# plot log2FC heatmap
cbind(log2FC_mat, cls = log2FC_cls) %>% 
  as.data.frame() %>%
  arrange(cls, rowMeans(log2FC_mat[, 3:4])) %>%
  cbind(id = seq_len(nrow(log2FC_mat))) %>%
# log2FC_mat[order(log2FC_cls, rowMeans(log2FC_mat[, 3:4])), ] %>%
  `rownames<-`(NULL) %>% 
  reshape::melt(id.vars = c("cls", "id")) %>% 
  ggplot(aes(x = variable, y = id, fill = value)) +
  geom_tile(width = 0.99) +
  facet_grid(cls ~ ., scales = "free", space = "free") +
  labs(title = 'Genes (n = 16734)',  x = '', y = '') +
  scale_fill_gradientn(name = "Log2FC", 
                       colours = c('blue', "white", 'red'), 
                       values = (c(min(log2FC_mat), 0, max(log2FC_mat)) - min(log2FC_mat)) /
                         diff(range(log2FC_mat)) ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.x = element_text(size = 11, hjust = 1, angle = 45),
        axis.text.y = element_blank(), 
        strip.text.y = element_blank(), 
        legend.direction = "horizontal", 
        legend.key.size = unit(13, "pt"), 
        plot.margin = unit(c(1, 1, 0, 1), "lines")) 

ggsave(filename = paste0("Fig1_RNAseq_log2FC_class_Hetmap.png"), 
       path = "figs/",
       device = "png", width = 5, height = 5)

# ----------------------------------------------------------------------------------- #
# plot TSS heatmap by pulsing time
source("F1_src_gene_tss_coverage_heatmap.R")

# -----------------------------TSS occupancy values---------------------------------- #
bw_files_1 <- list.files(path = "/mnt/E0767589767560E8/UPPMAX/Ezh2i_BAP1_ChIP",
                              pattern = ".*bw$", full.names = T)
bw_files_1 <- bw_files_1[!grepl("_Ezh2i|-Ezh2i|Trp", bw_files_1)]
bw_files_2 <- list.files(path = "/mnt/E0767589767560E8/UPPMAX/BAP1_PC",
                              pattern = "_P.*bw$", full.names = T)

all_ChIP_TSS_mat <- 
  .countBW(bw_files = c(rev(bw_files_1), bw_files_2),
           intervals = promoters(gene.gr[use_gene_ids], upstream = 1000, downstream = 1000),
           blacklist = blacklist.gr,
           fast = F) #%>% trim_quantile()

rownames(all_ChIP_TSS_mat) <- use_gene_ids
colnames(all_ChIP_TSS_mat) <-
  gsub(".*(ChIP/|BAP1_PC/)(.*).fltd.bw",
       "\\2",
       c(rev(bw_files_1), bw_files_2))
saveRDS(all_ChIP_TSS_mat, "data/all_ChIP_TSS_mat.RData")

Ring1b_PC_ChIP_TSS_mat <- all_ChIP_TSS_mat[, grep("Ring1b", colnames(all_ChIP_TSS_mat))]
Ezh2_PC_ChIP_TSS_mat <- all_ChIP_TSS_mat[, grep("Ezh2", colnames(all_ChIP_TSS_mat))]
H2AUb_PC_ChIP_TSS_mat <- all_ChIP_TSS_mat[, grep("H2Aub", colnames(all_ChIP_TSS_mat))]
Pol2_PC_ChIP_TSS_mat <- all_ChIP_TSS_mat[, grep("Pol2_", colnames(all_ChIP_TSS_mat))]

H3K27m3_PC_ChIP_TSS_mat <- all_ChIP_TSS_mat[, grep("H3K27m3", colnames(all_ChIP_TSS_mat))]


# Ring1b enriched genes
idx <- rowMeans(Ring1b_PC_ChIP_TSS_mat) > quantile(rowMeans(Ring1b_PC_ChIP_TSS_mat), 0.9)


# --------------------------------Plot ecdf by DE groups------------------------------- #
df <- data.frame(all_ChIP_TSS_mat[, grep("NT|P0_ALL", colnames(all_ChIP_TSS_mat))], 
            groups = c("Derepressed", "Weakly Derepressed", "Unchanged", "Weakly Repressed")[log2FC_cls])
df$groups <- factor(df$groups, c("Derepressed", "Weakly Derepressed", "Unchanged", "Weakly Repressed"))

g3.1 <- ggplot(df, aes(log2(H2Aub_NT_RC), colour = groups)) +
  stat_ecdf(lwd = 1) +
  theme_classic() + 
  scale_color_viridis_d(option = "E", direction = -1) +
  xlim(c(-3, 3)) +
  xlab("log2 TSS density") + ylab("ECDF") +
  ggtitle("H2Aub  ") +
  theme_setting +
  theme(legend.position = "none")

g3.2 <- ggplot(df, aes(log2(H3K27m3_NT_RC), colour = groups)) +
  stat_ecdf(lwd = 1) +
  theme_classic() + 
  scale_color_viridis_d(option = "E", direction = -1) +
  xlim(c(-3, 3)) +
  xlab("log2 TSS density") + ylab("ECDF") +
  ggtitle("H3K27me3  ") +
  theme_setting +
  theme(legend.position = "none")

g3.3 <- ggplot(df, aes(log2(Cbx7_P0_ALL), colour = groups)) +
  stat_ecdf(lwd = 1) +
  theme_classic() + 
  scale_color_viridis_d(option = "E", direction = -1) +
  xlim(c(-3, 4)) +
  xlab("log2 TSS density") + ylab("ECDF") +
  ggtitle("Cbx7  ") +
  theme_setting +
  theme(legend.position = "none")

g3.4 <- ggplot(df, aes(log2(H3K4m3_NT_RC), colour = groups)) +
  stat_ecdf(lwd = 1) +
  theme_classic() + 
  scale_color_viridis_d(option = "E", direction = -1) +
  xlim(c(-3, 3)) +
  xlab("log2 TSS density") + ylab("ECDF") +
  ggtitle("H3K4me3  ") +
  theme_setting +
  theme(legend.position = "none")

g3.5 <- ggplot(df, aes(log2(Pol2.S2p_NT_RC), colour = groups)) +
  stat_ecdf(lwd = 1) +
  theme_classic() + 
  scale_color_viridis_d(option = "E", direction = -1) +
  xlim(c(-4, 4)) +
  xlab("log2 TSS density") + ylab("ECDF") +
  ggtitle("Pol2-S2p  ") +
  theme_setting +
  theme(legend.position = "none")

g3.6 <- ggplot(df, aes(log2(Ring1b_P0_ALL), colour = groups)) +
  stat_ecdf(lwd = 1) +
  theme_classic() + 
  scale_color_viridis_d(option = "E", direction = -1) +
  xlim(c(-3, 6)) +
  xlab("log2 TSS density") + ylab("ECDF") +
  ggtitle("Ring1b  ") +
  theme_setting +
  theme(legend.position = "none")

g3.7 <- ggplot(df, aes(log2(Ezh2_NT_RC), colour = groups)) +
  stat_ecdf(lwd = 1) +
  theme_classic() + 
  scale_color_viridis_d(option = "E", direction = -1) +
  xlim(c(-3, 4)) +
  xlab("log2 TSS density") + ylab("ECDF") +
  ggtitle("Ezh2  ") +
  theme_setting +
  theme(legend.position = "none")

g3.8 <- ggplot(df, aes(log2(Rybp_P0_ALL), colour = groups)) +
  stat_ecdf(lwd = 1) +
  theme_classic() + 
  scale_color_viridis_d(option = "E", direction = -1) +
  xlim(c(-3, 4)) +
  xlab("log2 TSS density") + ylab("ECDF") +
  ggtitle("Rybp  ") +
  theme_setting +
  theme(legend.position = "none")

g3.9 <- ggplot(df, aes(log2(Pol2.NTD_NT_RC), colour = groups)) +
  stat_ecdf(lwd = 1) +
  theme_classic() + 
  scale_color_viridis_d(option = "E", direction = -1) +
  xlim(c(-4, 4)) +
  xlab("log2 TSS density") + ylab("ECDF") +
  ggtitle("Pol2  ") +
  theme_setting +
  theme(legend.position = "none")

g3.10 <- ggplot(df, aes(log2(Pol2.S5p_NT_RC), colour = groups)) +
  stat_ecdf(lwd = 1) +
  theme_classic() + 
  scale_color_viridis_d(option = "E", direction = -1) +
  xlim(c(-4, 4)) +
  xlab("log2 TSS density") + ylab("ECDF") +
  ggtitle("Pol2-S5p  ") +
  theme_setting +
  theme(legend.position = "none")

ggsave(plot = grid.arrange(g3.1, g3.2, g3.3, g3.4, g3.5, 
                           g3.6, g3.7, g3.8, g3.9, g3.10, ncol = 5),
       filename = paste0("Fig1_ECDF_ChIP_DE_groups.png"), 
       path = "figs/",
       device = "png", width = 15, height = 6)

# ---------------------------------GO analysis--------------------------------- #
source("../util/go.R")
enrichGeneSets(genes_derepressed, ontology = "BP")


# -----------------------------------PRC1 types-------------------------------------- #
Ring1b_binding_genes <-
  rownames(all_ChIP_TSS_mat)[kink_index(all_ChIP_TSS_mat[, "Ring1b_P0_ALL"], 
                                        method = "weight")]

Cbx7_binding_genes <-
  rownames(all_ChIP_TSS_mat)[kink_index(all_ChIP_TSS_mat[, "Cbx7_P0_ALL"], 
                                        method = "weight")]
Rybp_binding_genes <-
  rownames(all_ChIP_TSS_mat)[kink_index(all_ChIP_TSS_mat[, "Rybp_P0_ALL"], 
                                        method = "weight")]

Rsf1_peaks <- importRanges("../data/Rsf1_peaks.narrowPeak")
Rsf1_peaks <- Rsf1_peaks[Rsf1_peaks$qValue > 10]
Rsf1_binding_genes <- names(gene.gr)[findOverlaps(gene.gr, Rsf1_peaks) %>% countQueryHits() > 0]


PRC1_type_DE <- log2FC_mat[unique(c(Ring1b_binding_genes, 
                                    Cbx7_binding_genes,
                                    Rybp_binding_genes)), ] %>% 
  `colnames<-`(c("P12", "P12C12", "P12C24", "P12C36")) %>%
  reshape::melt()

PRC1_type_DE$PRC1 <- ifelse(PRC1_type_DE$X1 %in% intersect.Vector(Ring1b_binding_genes, 
                                                                  intersect.Vector(Cbx7_binding_genes, 
                                                   Rybp_binding_genes)), 
                            "111", ifelse(!PRC1_type_DE$X1 %in% c(Cbx7_binding_genes, Rybp_binding_genes),
                                           "100", 
                                          ifelse(!PRC1_type_DE$X1 %in% c(Cbx7_binding_genes, Ring1b_binding_genes),
                                                 "001",
                                                 ifelse(!PRC1_type_DE$X1 %in% c(Rybp_binding_genes, Ring1b_binding_genes),
                                                        "010",
                                                        ifelse(!PRC1_type_DE$X1 %in% Rybp_binding_genes,
                                                               "110",
                                                               ifelse(!PRC1_type_DE$X1 %in% Cbx7_binding_genes,
                                                                      "101",
                                                                      "011"
                                                                      )
                                                              )
                                                        )
                                                 )
                                          )
                            )
PRC1_type_DE$PRC1[PRC1_type_DE$PRC1 %in% c("001", "010", "011")] <- "Cbx7 / Rybp\nn=901" # 901
PRC1_type_DE$PRC1[PRC1_type_DE$PRC1 == "100"] <- "Ring1b only\nn=343" # 343
PRC1_type_DE$PRC1[PRC1_type_DE$PRC1 == "110"] <- "Canonical\nn=351" # 351
PRC1_type_DE$PRC1[PRC1_type_DE$PRC1 == "101"] <- "Non-canonical\nn=211" # 211
PRC1_type_DE$PRC1[PRC1_type_DE$PRC1 == "111"] <- "Both\nn=240" # 240
PRC1_type_DE$PRC1 <- factor(PRC1_type_DE$PRC1,
                            levels = c("Cbx7 / Rybp\nn=901", 
                                       "Ring1b only\nn=343",
                                       "Non-canonical\nn=211",
                                       "Canonical\nn=351", 
                                       "Both\nn=240"))

ggplot(PRC1_type_DE, aes(x = PRC1, y = value, fill = X2)) +
  geom_hline(yintercept = 0, lty = 1, lwd = 1, col = "red") +
  geom_boxplot(outlier.size = 0, notch = T) +
  scale_fill_viridis_d(begin = 0.2, end = 0.84) +
  xlab("") + ylab("RNA-seq log2FC") +
  labs(fill = "Pulse-chase") +
  ggpubr::theme_pubclean() +
  theme(axis.text = element_text(face = "bold", size = 12), 
        axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 13))

ggsave(filename = paste0("Fig1_RNA_DE_PRC1_groups.png"), 
       path = "figs/",
       device = "png", width = 5, height = 5)







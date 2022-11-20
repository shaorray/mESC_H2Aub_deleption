# Rui Shao 2022 Oct
# Figure S10
# double depletion - recovery

# to test if H2Aub and H3K27me3 could restore to the untreated condition after depletion 

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
source("../util/getCoverage.R")

# -------------------------------------------------------------------------------------------------

# load data

# native_files <- list.files("/mnt/E0767589767560E8/UPPMAX/LOF_20210121/bw", ".bw", full.names = TRUE)
# native_files <- native_files[!grepl("supersmooth|P1_", native_files)]
# 
# native_conds <- data.frame("Mark" = gsub(".*P2_(.*)_Ezh2.*", "\\1", native_files),
#                            "Treat" = gsub(".*_(Ezh2.*).fltd.*", "\\1", native_files))
# native_conds$Mark <- gsub("H3K27m.*", "H3K27me3", native_conds$Mark)


cross_linked_files <- list.files("/mnt/E0767589767560E8/UPPMAX/LOF3_20210205/bw", ".bw", full.names = TRUE)
cross_linked_files <- cross_linked_files[!grepl("supersmooth|_R", cross_linked_files)]

cross_linked_conds <- data.frame("Mark" = gsub(".*bw/(.*)_(Ezh2|P).*", "\\1", cross_linked_files),
                                 "Treat" = gsub(".*_(Ezh2.*|P.*)_ALL.fltd.*", "\\1", cross_linked_files))
cross_linked_conds$Treat <- gsub("Ezh2i7DP", "Ezh2i7D_P", cross_linked_conds$Treat)


# TSS coverage

# TSS_mat_native <- .countBW(native_files, 
#                            intervals = promoters(gene.gr[PcG_enriched_genes], upstream = 5000, downstream = 5000),
#                            fast = TRUE)
TSS_mat_cross_linked <- .countBW(cross_linked_files, 
                                 intervals = promoters(gene.gr[PcG_enriched_genes], upstream = 5000, downstream = 5000),
                                 fast = TRUE)
TSS_mat_cross_linked_all <- .countBW(cross_linked_files,
                                     intervals = promoters(gene.gr, upstream = 5000, downstream = 5000),
                                     fast = TRUE)

# TSS_cov_native <- convert_coverage(file_names = native_files,
#                                    intervals = promoters(gene.gr[PcG_enriched_genes], upstream = 5000, downstream = 5000),
#                                    bin_width = 100,
#                                    new_len = 100)
TSS_cov_cross_linked <- convert_coverage(file_names = cross_linked_files,
                                         intervals = promoters(gene.gr[PcG_enriched_genes], upstream = 5000, downstream = 5000),
                                         bin_width = 100,
                                         new_len = 100)


# -------------------------------------------------------------------------------------------------
# plot PCA
plot_recovery_pca <- function(x, cond_names) {
  
  idx <- match(paste(x, names(cond_names)), paste(cross_linked_conds$Mark, cross_linked_conds$Treat))
  pca <- prcomp(log1p(t(trim_quantile(TSS_mat_cross_linked[, idx], 0.99))))
  
  pca$sdev <- pca$sdev^2
  pc_var <- pca$sdev[1:2] / sum(pca$sdev) * 100
  
  dat <- data.frame(
    PC1 = pca$x[, 1] ,
    PC2 = pca$x[, 2] ,
    Sample = cond_names
  )
  
  grid_dat <- data.frame(x1 = dat$PC1[-nrow(dat)], 
                         y1 = dat$PC2[-nrow(dat)], 
                         x2 = dat$PC1[-1], 
                         y2 = dat$PC2[-1], 
                         Sample = dat$Sample[-1],
                         Conditions = c(rep("Depletion", 2), rep("Recovery", 3)))
  
  ggplot(dat, aes(x = PC1, y = PC2)) +
    geom_point(size = 3) +
    geom_segment(data = grid_dat,
                 mapping = aes(x = x1, y = y1,
                               xend = x2, yend = y2,
                               group = Sample,
                               color = Conditions),
                 arrow = arrow(type = "open", angle = 20, length = unit(x = 10, 'pt')),
                 size = 2, alpha = 0.5) + 
    ggforce::geom_mark_ellipse(aes(label = Sample, fill = Sample), 
                               show.legend = FALSE, 
                               expand = unit(0, "mm"),
                               label.fontsize = 12, 
                               label.buffer = unit(0.1, 'mm'),
                               label.fill = NA,
                               con.type = "none") +
    xlab(paste0("PC1 (",round(pc_var[1]), "% variance)" )) + 
    ylab(paste0("PC2 (",round(pc_var[2]), "% variance)" )) +
    ggtitle(x) +
    xlim(range(dat$PC1) * 1.3) + ylim(range(dat$PC2) * 1.3) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_text(hjust = 0, vjust = 0, size = 13), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggsave(paste0("Fig4_PCA_ChIP_depletion_recovery_", x, ".png"), 
         device = "png", path = "../fig4/figs", width = 4, height = 4)
}

cond_names <- c("P0", "Ezh2i (7d)", "BAP1", "1d", "2d", "7d")
names(cond_names) <- c("P0", "Ezh2i7D", "Ezh2i7D_P12C12", "Ezh2i7D_P12C1D", "Ezh2i7D_P12C2D", "Ezh2i7D_P12C7D")

for (x in unique(cross_linked_conds$Mark)) {
  plot_recovery_pca(x, cond_names)
}


# TSS density
plot_recovery_line("H3K27me3")
plot_recovery_line("H2Aub")
plot_recovery_line("Ezh2")
plot_recovery_line("Ring1b")
plot_recovery_line("Cbx7")
plot_recovery_line("Rybp")
plot_recovery_line("Pol2")

plot_recovery_line <- function(x) {
  
  idx <- match(paste(x, names(cond_names)), paste(cross_linked_conds$Mark, cross_linked_conds$Treat))
  
  cond_names2 <- c("P0", "Ezh2i\n (7d)", "BAP1", "1d", "2d", "7d")
  
  dat <- TSS_mat_cross_linked[, idx] %>%
    trim_quantile(0.99) %>%
    "*"(10) %>% 
    log1p() %>% 
    t() %>% 
    scale() %>%
    t() %>% 
    `colnames<-`(., cond_names2) 
  
  dat_mean <- colMeans(dat) %>% as.matrix() %>% 
    reshape::melt() %>% 
    dplyr::mutate(X1 = factor(X1, levels = cond_names2),
                  X2 = "mu")
  dat_mean$Conditions <- c(rep("Depletion", 2), rep("Recovery", 4))
  
  dat <- dat %>% 
    reshape::melt() %>% 
    dplyr::mutate(X2 = factor(X2, levels = cond_names2))
  
  ggplot(dat, aes(x = X2, y = value, group = X1)) +
    geom_line(size = 0.2, color = "grey60", alpha = 0.1) +
    geom_line(data = dat_mean, aes(x = X1, y = value, group = X2, color = Conditions), 
              size = 2, alpha = 0.5) +
    scale_x_discrete(expand = c(0.05, 0.05)) +
    xlab("") + ylab("Z-scores") + ggtitle(x) +
    theme_minimal() +
    theme(axis.line = element_line(),
          axis.ticks = element_line(), 
          axis.text = element_text(size = 12),
          panel.grid.major = element_blank())
  
  ggsave(paste0("Fig4_line_plot_ChIP_depletion_recovery_", x, ".png"), 
         device = "png", path = "../fig4/figs", width = 5, height = 3.5)
}


# TSS coverage

plot_TSS_coverage <- function(x) {
  idx <- match(paste(x, names(cond_names)), paste(cross_linked_conds$Mark, cross_linked_conds$Treat))
  
  dat <- data.frame(pos = rep(c(1:100), 6),
                    Coverage = c(apply(TSS_cov_cross_linked[, , idx], 3, quantile_mean)),
                    Condition = rep(cond_names, each = 100),
                    Period = c(rep("Depletion", each = 300), rep("Recovery", each = 300)))
  dat$Condition <- factor(dat$Condition, levels = cond_names)
  
  ggplot(dat, aes(x = pos, y = Coverage, group = Condition, color = Condition)) +
    geom_line(size = 1) +
    facet_grid(. ~ Condition) +
    scale_x_continuous(name = "", breaks = c(1, 50, 100), labels = c("  -5 kb", "TSS", "5 kb  ")) +
    ylab("Mean RPGC") + ggtitle(x) +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid = element_blank(), 
          axis.line = element_line(),
          axis.ticks = element_line(), 
          axis.text = element_text(size = 10))
}


g_list <- sapply(c("H3K27me3", "H2Aub", "Ezh2", "Ring1b"), function(x) list(plot_TSS_coverage(x)))


ggsave(plot = cowplot::plot_grid(plotlist = g_list,
                                 ncol = 1),
       filename = "Fig4_TSS_coverage_ChIP_depletion_recovery2.png", 
       path = "../fig4/figs",
       device = "png", width = 12, height = 8)


# recovery rate with simple linear regression

idx <- match(paste("H2Aub", names(cond_names)), paste(cross_linked_conds$Mark, cross_linked_conds$Treat))
fit <- lm(t(as.matrix(TSS_mat_cross_linked_all)[, idx][, -(1:2)]) ~ c(0, 1, 2, 7))
H2Aub_slope_terms <- coef(fit)[2, ]

idx <- match(paste("H3K27me3", names(cond_names)), paste(cross_linked_conds$Mark, cross_linked_conds$Treat))
fit <- lm(t(as.matrix(TSS_mat_cross_linked_all)[, idx][, -1]) ~ c(0, 1, 2, 3, 8))
H3K27me3_slope_terms <- coef(fit)[2, ]


PCL2_peak <- importRanges("../data/peaks/2019_Hojfeldt_WT_mESC_PCL2_ChIP_peaks.narrowPeak")
JARID2_peak <- importRanges("../data/peaks/2019_Hojfeldt_WT_mESC_JARID2_ChIP_peaks.narrowPeak")

PCL2_idx <- countQueryHits(findOverlaps(promoters(gene.gr[Ezh2_binding_genes]), PCL2_peak)) > 0
JARID2_idx <- countQueryHits(findOverlaps(promoters(gene.gr[Ezh2_binding_genes]), JARID2_peak)) > 0

dat <- rbind(data.frame(Type = "PRC2.1 (PCL2)\n (n=2528)",
                        Slope = H3K27me3_slope_terms[names(gene.gr) %in% Ezh2_binding_genes][PCL2_idx & !JARID2_idx]),
             data.frame(Type = "PRC2.2 (JARID2)\n (n=1113)",
                        Slope = H3K27me3_slope_terms[names(gene.gr) %in% Ezh2_binding_genes][JARID2_idx]))

library(ggpubr)

g1 <- ggplot(dat, aes(x = Type, y = Slope, fill = Type)) +
  geom_boxplot(outlier.color = NA) +
  coord_cartesian(ylim = c(-0.1, 0.3)) +
  stat_compare_means(aes(label = ..p.signif..), na.rm = TRUE,
                     method = "t.test", label.y = 0.3, paired = FALSE) +
  scale_fill_manual(values = colors_20[c(1,2)]) +
  xlab("\nEzh2 binding genes (n=8363)") + ylab("Linear slope (0-7d)") + ggtitle("H3K27me3 recovery") +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        legend.position = "none", 
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))


dat <- rbind(data.frame(Type = "PRC2.1 (PCL2)\n (n=2528)",
                        Slope = H2Aub_slope_terms[names(gene.gr) %in% Ezh2_binding_genes][PCL2_idx & !JARID2_idx]),
             data.frame(Type = "PRC2.2 (JARID2)\n (n=1113)",
                        Slope = H2Aub_slope_terms[names(gene.gr) %in% Ezh2_binding_genes][JARID2_idx]))

g2 <- ggplot(dat, aes(x = Type, y = Slope, fill = Type)) +
  geom_boxplot(outlier.color = NA) +
  coord_cartesian(ylim = c(-0.1, 0.3)) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test", label.y = 0.3, paired = FALSE) +
  scale_fill_manual(values = colors_20[c(1,2)]) +
  xlab("\nEzh2 binding genes (n=8363)") + ylab("Linear slope (1-7d)") + ggtitle("H2Aub recovery") +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        legend.position = "none", 
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

ggsave(grid.arrange(g1, g2, nrow = 1), 
       filename = "Fig4_TSS_PRC2_recovery_slope.png", 
       path = "../fig4/figs",
       device = "png", width = 7, height = 3.5)


x1 = intersect(Cbx7_binding_genes, Ring1b_binding_genes)
x1 = x1[!x1 %in% Rybp_binding_genes]

x2 = intersect(Rybp_binding_genes, Ring1b_binding_genes)
x2 = x2[!x2 %in% Cbx7_binding_genes]


dat <- rbind(data.frame(Type = "vPR1 (Rybp)\n (n=1710)",
                        Slope = H3K27me3_slope_terms[names(gene.gr) %in% x2]),
             data.frame(Type = "cPRC1 (Cbx7)\n (n=1360)",
                        Slope = H3K27me3_slope_terms[names(gene.gr) %in% x1])
             )

g1 <- ggplot(dat, aes(x = Type, y = Slope, fill = Type)) +
  geom_boxplot(outlier.color = NA) +
  coord_cartesian(ylim = c(-0.1, 0.3)) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test", label.y = 0.3, paired = FALSE) +
  scale_fill_manual(values = colors_20[c(19, 6)]) +
  xlab("\nRing1b binding genes (n=8277)") + ylab("Linear slope (0-7d)") + ggtitle("H3K27me3 recovery") +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        legend.position = "none", 
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))



dat <- rbind(data.frame(Type = "vPR1 (Rybp)\n (n=1710)",
                        Slope = H2Aub_slope_terms[names(gene.gr) %in% x2]),
             data.frame(Type = "cPRC1 (Cbx7)\n (n=1360)",
                        Slope = H2Aub_slope_terms[names(gene.gr) %in% x1])
             )

g2 <- ggplot(dat, aes(x = Type, y = Slope, fill = Type)) +
  geom_boxplot(outlier.color = NA) +
  coord_cartesian(ylim = c(-0.1, 0.3)) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test", label.y = 0.3, paired = FALSE) +
  scale_fill_manual(values = colors_20[c(19, 6)]) +
  xlab("\nRing1b binding genes (n=8277)") + ylab("Linear slope (1-7d)") + ggtitle("H2Aub recovery") +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        legend.position = "none", 
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))


ggsave(grid.arrange(g1, g2, nrow = 1), 
       filename = "Fig4_TSS_PRC1_recovery_slope.png", 
       path = "../fig4/figs",
       device = "png", width = 7, height = 3.5)

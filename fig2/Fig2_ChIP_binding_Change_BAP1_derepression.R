# Rui Shao 2021 May
# Figure 2
# Ring1b binding and BAP1 pulse derepression

# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# source("../util/utils.R")
# source("../util/getCoverage.R")

library(ggplot2)
library(ggpubr)
library(gridExtra)


# sample_colors <- c("H2Aub" = colors_20[16], 
#                    "H3K27me3" = colors_20[6],
#                    "H3K27ac" = colors_20[3],
#                    "H3K4me3" =  colors_20[2],
#                    "Ring1b" =  colors_20[15],
#                    "Ezh2" =  colors_20[5],
#                    "Pol II" =  colors_20[10],
#                    "RNAPII-NTD" =  colors_20[10],
#                    "Pol II-S5p" =  colors_20[8], 
#                    "Pol II-S2p" =  colors_20[11],
#                    "Cbx7" = colors_9[1],
#                    "Rybp" = colors_9[3],
#                    "TT-seq" = colors_9[9])

# list of polycomb domains (Kundu et al. 2017) (mm9)
# HoxA         chr6 51526356 52761443
# HoxD         chr2 73683519 75068833
# Pax6         chr2 105290036 105753314 (no derepression)
# Pax1/Nkx2.2  chr2 146637019 147259205
# Igf2bp3      chr6 48985949 49491383

# polycomb_domain_Kundu <- GRanges(seqnames = paste0("chr", c(6, 2, 2, 2, 6)),
#                                  IRanges(start = c(51526356, 73683519, 105290036, 146637019, 48985949),
#                                          end = c(52761443, 75068833, 105753314, 147259205, 49491383)))

# gene_intersect <- H2Aub_repressed_genes[H2Aub_repressed_genes %in% 
#                                           intersect.Vector(rownames(ChIP_TSS_B1_mat),
#                                                            rownames(ChIP_TSS_B2_mat))] # 1903 genes



# -----------------------------------------------------------------------------------------

plot_boxplot <- function(target, gene_ids, .color = 1) {
  
  mat <- ChIP_TSS_B2_mat[match(gene_ids, rownames(ChIP_TSS_B2_mat)),
                         grep(target, colnames(ChIP_TSS_B2_mat))]
  mat <- `colnames<-`(log2(mat + 0.1), gsub(".*_(P.*)", "\\1", colnames(mat)))
  medians <- matrixStats::colMedians(as.matrix(mat), na.rm = T)
  # cmps <- lapply(colnames(mat)[-1], function(x) c(colnames(mat)[1], x))
  fold_change_tab <- data.frame(x = colnames(mat),
                                y = medians[1] + 3,
                                label = round(2^(medians - medians[1]), 2))
  
  mat %>% 
    as.data.frame() %>%
    reshape::melt() %>%
    ggplot(aes(x = variable, y = value)) +
    geom_hline(yintercept = medians[1], alpha = 0.1, size = 1) +
    stat_boxplot(geom = "errorbar", width = 0.1) +
    geom_boxplot(color = .color, outlier.color = NA, width = 0.7) +
    # stat_compare_means(comparisons = cmps, size = 0) +
    geom_text(data = fold_change_tab, aes(x = x, y = y, label = label)) +
    ylim(c(medians[1] - 3, medians[1] + 3)) +
    xlab("") + ylab(expression("log"[2] ~ "TSS RPGC")) + ggtitle(target) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
          panel.grid = element_blank(), 
          legend.position = "none", 
          plot.title = element_text(size = 14, hjust = 0.05))
}


plot_boxplot(target = "Ring1b", 
             gene_ids = PcG_enriched_genes, 
             .color = sample_colors["Ring1b"])

# Ring1b enriched genes
ggsave(plot = do.call(grid.arrange,
                      c(lapply(c("H2Aub", "Ring1b", "Cbx7", "Rybp", "H3K27me3", "Ezh2"), 
                               function(x) plot_boxplot(target = x, 
                                                        gene_ids = Ring1b_enriched_genes, 
                                                        .color = sample_colors[x]) ), ncol = 3)),
       "figs/Fig2_boxplot_batch2_Ring1b_enriched_genes.png", height = 5, width = 7)

# PcG enriched genes
ggsave(plot = do.call(grid.arrange,
                      c(lapply(c("H2Aub", "Ring1b", "Cbx7", "Rybp", "H3K27me3", "Ezh2"), 
                               function(x) plot_boxplot(target = x, 
                                                        gene_ids = PcG_enriched_genes, 
                                                        .color = sample_colors[x]) ), ncol = 3)),
       "figs/Fig2_boxplot_batch2_PcG_enriched_genes.png", height = 5, width = 7)

# all genes
ggsave(plot = do.call(grid.arrange,
                      c(lapply(c("H2Aub", "Ring1b", "Cbx7", "Rybp", "H3K27me3", "Ezh2"), 
                               function(x) plot_boxplot(target = x, 
                                                        gene_ids = gene.gr$gene_id, 
                                                        .color = sample_colors[x]) ), ncol = 3)),
       "figs/Fig2_boxplot_batch2_all_genes.png", height = 5, width = 7)



# ----------------------------------- log2 fold-changes correlation ---------------------------------- #
# Ring1b enriched genes 

P12_log2FC <- cbind(ChIP_TSS_B1_mat_log2FC[, grepl("BAP1-12$", colnames(ChIP_TSS_B1_mat_log2FC)) &
                                             !grepl("IN|S2p|S5p", colnames(ChIP_TSS_B1_mat_log2FC))],
                    ChIP_TSS_B2_mat_log2FC[, grep("(Cbx7|Rybp)_P12$", colnames(ChIP_TSS_B2_mat_log2FC))])
P12_log2FC <- P12_log2FC[complete.cases(P12_log2FC), ] # 23437 genes

dev.off()
pdf("figs/Fig2_corrplot_ChIP_log2FC_PcG_enriched.pdf", width = 5, height = 4.5)
r1 <- P12_log2FC[intersect(PcG_enriched_genes, rownames(P12_log2FC)), ] %>% # 1480 genes
  as.data.frame() %>%
  dplyr::filter(is.finite(rowSums(.)) & complete.cases(.)) %>%
  `colnames<-`(., unique(gsub("(_.*)", "\\2", colnames(P12_log2FC)))) %>%
  cor()
r1[matrix(rep(seq_len(nrow(r1)), 2), ncol = 2)] <- NA
pheatmap::pheatmap(r1, 
                   display_numbers = TRUE, 
                   number_color = "white",
                   na_col = NA, 
                   border_color = "white",
                   color = add.alpha(colorRampPalette(rev(viridis(10, option = "B")))(50), 0.7),
                   breaks = seq(0, 0.5, length.out = 50))
dev.off()


# -------------------------------- log2FC by PcG taget and bivalent TSS ------------------------------------ #

P12_log2FC_bivalent <- as.data.frame(P12_log2FC[intersect(bivalent_genes_chromHMM, rownames(P12_log2FC)),
                                                c("H3K27me3_BAP1-12", "H2Aub_BAP1-12", "Ring1b_BAP1-12", "Ezh2_BAP1-12", "Cbx7_P12")])

P12_log2FC_bivalent$PcG_enriched <- ifelse(rownames(P12_log2FC_bivalent) %in% PcG_enriched_genes, 
                                              "PcG enriched (n=1354)", "Bivalent TSS (n=4404)")
P12_log2FC_bivalent$Expression <- ifelse(rownames(P12_log2FC_bivalent) %in% derepressed_gene_ids,
                                 "Up", "Stable")

P12_log2FC_bivalent %>%
  reshape::melt(id.vars = c("PcG_enriched", "Expression")) %>%
  dplyr::mutate(variable = factor(gsub("_P12", "", gsub("_BAP1-12", "", variable)),
                                  levels = rev(c("Ring1b", "Cbx7", "Ezh2", "H2Aub", "H3K27me3"))),
                PcG_enriched = factor(PcG_enriched, 
                                         levels = rev(c("PcG enriched (n=1354)",
                                                    "Bivalent TSS (n=4404)"))),
                Expression = factor(Expression,
                                    levels = rev(c("Up", "Stable")))) %>% # "Up \n(n=1292)", "Unchanged \n(n=2183)"
  ggplot(aes(x = variable, y = value, color = Expression)) +
  geom_hline(yintercept = 0, alpha = 0.1, size = 1) +
  stat_boxplot(geom = "errorbar", position = position_dodge(0.7), width = 0.2) +
  geom_boxplot(position = position_dodge(), outlier.color = NA,  width = 0.7) +
  ggpubr::stat_compare_means(aes(label = ..p.signif..), method = "t.test") +
  facet_grid(. ~ PcG_enriched) +
  scale_color_manual(name = "Expression", values = rev(c("red3", "darkblue"))) +
  ylim(c(-2.3, 1.5)) +
  xlab("") + ylab((expression("P12 TSS log"[2] ~ "fold-change"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.ticks = element_line(),
        axis.title = element_text(size = 12),
        legend.key.height = unit(1.2, 'cm'),
        panel.grid = element_blank(), 
        strip.text = element_text(size = 13))

ggsave("figs/Fig2_boxplot_bivalent_TSS_cross_Ring1b_DE2.png", width = 6.5, height = 4)

# --------------------------------- Ezh2 change on sub-units sites ------------------------------------- #
PCL2_peak <- importRanges("../data/peaks/2019_Hojfeldt_WT_mESC_PCL2_ChIP_peaks.narrowPeak")
JARID2_peak <- importRanges("../data/peaks/2019_Hojfeldt_WT_mESC_JARID2_ChIP_peaks.narrowPeak")

dat_Ezh2_change <- ChIP_TSS_B2_mat_log2FC[Ezh2_binding_genes, grep("^Ezh2", colnames(ChIP_TSS_B2_mat_log2FC))] %>% as.data.frame()
dat_Ezh2_change$type <- "Other \n(n=1204)"

PCL2_idx <- countQueryHits(findOverlaps(promoters(gene.gr[Ezh2_binding_genes]), PCL2_peak)) > 0
JARID2_idx <- countQueryHits(findOverlaps(promoters(gene.gr[Ezh2_binding_genes]), JARID2_peak)) > 0

dat_Ezh2_change$type[PCL2_idx] <- "PCL2 \n(n=1169)"
dat_Ezh2_change$type[JARID2_idx] <- "JARID2 \n(n=933)"

table(dat_Ezh2_change$type)
# JARID2  Other   PCL2 
#  933   1204   1169

dat_Ezh2_change$type <- factor(dat_Ezh2_change$type, levels = c( "Other \n(n=1204)", "JARID2 \n(n=933)", "PCL2 \n(n=1169)"))

cmp <- list(c("JARID2 \n(n=933)", "PCL2 \n(n=1169)"), c("JARID2 \n(n=933)", "Other \n(n=1204)") )

ggplot(dat_Ezh2_change, aes(x = type, y = `Ezh2_P12C36`, color = type)) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(aes(label = ..p.signif..),
                     comparisons = cmp, method = "t.test",
                     label.y = c(1.0, 1.6)) +
  ylim(c(-2.2, 2)) +
  xlab("") + ylab("Ezh2 TSS occupancy log2FC at P12") + ggtitle("Ezh2 binding genes (n=3306)") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none")

ggsave("../figS2/figs/FigS2_boxplot_Ezh2_log2FC_Jarid2_PCL_peaks.png", width = 4, height = 4)

# --------------------------------- PRC2 change on sub-units sites ------------------------------------- #
dat_Ring1b_change <- ChIP_TSS_B2_mat_log2FC[Ring1b_binding_genes, grep("^Ring1b", colnames(ChIP_TSS_B2_mat_log2FC))] %>% 
  as.data.frame()
dat_Ring1b_change$type <- "Other \n(n=1875)"
dat_Ring1b_change$type[Ring1b_binding_genes %in% Rybp_binding_genes] <- "Rybp \n(n=603)"
dat_Ring1b_change$type[Ring1b_binding_genes %in% Cbx7_binding_genes] <- "Cbx7 \n(n=695)"

table(dat_Ring1b_change$type)
# Other  Rybp  Cbx7 
# 1875   603   695 

dat_Ring1b_change$type <- factor(dat_Ring1b_change$type, 
                                 levels = c( "Other \n(n=1875)", "Cbx7 \n(n=695)", "Rybp \n(n=603)"))

cmp <- list(c("Rybp \n(n=603)", "Cbx7 \n(n=695)"), c("Other \n(n=1875)", "Cbx7 \n(n=695)"))

ggplot(dat_Ring1b_change, aes(x = type, y = `Ring1b_P12C36`, color = type)) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_boxplot(outlier.colour = NA) +
  stat_compare_means(aes(label = ..p.signif..),
                     comparisons = cmp, method = "t.test",
                     label.y = c(0.6, 1.1)) +
  ylim(c(-2.2, 1.6)) +
  xlab("") + ylab("Ring1b TSS occupancy log2FC at P12") + ggtitle("Ring1b binding genes (n=3173)") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none")

ggsave("../figS2/figs/FigS2_boxplot_Ring1b_log2FC_Cbx7_Rybp_binding.png", width = 4, height = 4)


# ----------------------------------------------------------------------------------------- #
# Response Time

get_response_index <- function(x) x[1] / x[2] # this simple formula can effectively report the order of changes

dat_tss <- ChIP_TSS_B1_mat_log2FC[intersect.Vector(rownames(ChIP_TSS_B1_mat_log2FC), PcG_binding_genes), 
                                      grepl("BAP1-6$|BAP1-12$", colnames(ChIP_TSS_B1_mat_log2FC)) &
                                    !grepl("IN", colnames(ChIP_TSS_B1_mat_log2FC))] 

response_index_mat <- data.frame(Type = "TT-seq",
                                 RI = cbind(TTseq_res.list$res_E0B6$log2FoldChange,
                                            TTseq_res.list$res_E0B12$log2FoldChange
                                 )[match(rownames(dat_tss), rownames(TTseq_res.list$res_E0B6)), ] %>%
                                   apply(., 1, get_response_index) %>% trim_quantile(),
                                 Pos = "TSS")

response_index_mat <- rbind(response_index_mat,
                            sapply(c("H2Aub", "H3K27me3", "H3K27ac", "Pol II-S2p", "Pol II-S5p", "Ring1b", "Ezh2", "Pol II"),
                                   function(i) {
                                     dat_tss[, rev(grep(i, colnames(dat_tss)))] %>% 
                                       apply(., 1, get_response_index) %>% 
                                       trim_quantile(0.9)
                                   }) %>% reshape::melt() %>% 
                              `colnames<-`(., c("gene_id", "Type", "RI")) %>% 
                              dplyr::mutate(gene_id = NULL, Pos = "TSS")
) %>% dplyr::filter(complete.cases(.) & is.finite(RI)) 

response_index_mat$Type <- factor(response_index_mat$Type,
                                  levels = response_index_mat %>%
                                    dplyr::group_by(Type) %>% 
                                    dplyr::summarise(median = median(RI, na.rm = TRUE)) %>%
                                    dplyr::arrange(median) %>% 
                                    "$"("Type") %>% 
                                    as.character() )


ggplot(response_index_mat, aes(x = RI, y = Type, fill = Type)) +
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = -Inf, xmax = 0),
            fill = "grey97") +
  geom_rect(aes(ymin = -Inf, ymax = Inf, xmin = 1, xmax = Inf),
            fill = "grey97") +
  geom_vline(xintercept =  c(0, 1), alpha = 0.2, size = 0.5, lty = 2) +
  geom_vline(xintercept =  0.5, alpha = 0.2, size = 0.5) +
  ggridges::geom_density_ridges(rel_min_height = 0.005, 
                                quantile_lines = TRUE,
                                quantiles = 2, 
                                lty = 1) +
  scale_x_continuous(name = "Response Index",
                     breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(-0.5, 1.5)) +
  scale_fill_manual(values = add.alpha(sample_colors[levels(response_index_mat$Type)], 0.7)) +
  ylab("") + 
  theme_minimal() +
  theme(axis.text = element_text(size = 13), 
        axis.title = element_text(size = 14),
        legend.position = "none",
        panel.grid = element_blank()
  )

ggsave(filename = paste0("../fig2/figs/Fig2_ridgeplot_response_index_PcG_binding_genes.png"), 
       device = "png", width = 6, height = 4.5)

# ---------------------------------------- Triptolide ------------------------------------------- #
# supplementary plot

dat_Trp <-
  ChIP_TSS_B1_mat_log2FC[intersect.Vector(rownames(ChIP_TSS_B1_mat_log2FC), PcG_binding_genes), 
                             grep("Trp", colnames(ChIP_TSS_B1_mat_log2FC))] %>% 
  as.data.frame() %>%
  dplyr::select(starts_with(c("H2Aub", "Ring1b", "H3K27me3", "Ezh2", "H3K4me3"))) %>%
  as.matrix() %>%
  reshape::melt() %>% cbind(., Pos = "TSS") %>% 
  dplyr::mutate(X2 = factor(gsub("_Trp-9", "", X2), 
                            levels = c("H3K4me3", "H2Aub", "Ring1b", "H3K27me3", "Ezh2")))

t_tests = dat_Trp %>%
  group_by(X2) %>%
  summarise(P = t.test(value, mu = 0)$p.value,
            Label = paste("p <", formatC(P, format = "e", digits = 2)),
            Sig = ifelse(P < 0.05, "*", "ns"))


g1 <- ggplot(dat_Trp, aes(x = X2, y = value, color = X2)) +
  geom_hline(yintercept =  0, color = "red", alpha = 0.5, size = 1) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_boxplot(outlier.colour = NA) +
  geom_text(data = t_tests, aes(label = Label, y = 1.5), size = 3.5) +
  scale_color_manual(values = sample_colors[levels(dat_Trp$X2)]) +
  xlab("") + ylab("TSS log2FC (Triptolide 9h)\n") + ggtitle("PcG binding genes (n=2093)") +
  scale_y_continuous(breaks = c(-2, -1, 0, 1), labels = c(-2, -1, 0, 1), limits = c(-2, 1.5)) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 15), 
        legend.position = "none",
        panel.grid = element_blank())

dat_Trp2 <-
  ChIP_TSS_B1_mat_log2FC[intersect.Vector(rownames(ChIP_TSS_B1_mat_log2FC), Ring1b_enriched_genes), 
                             grep("Pol II.*Trp", colnames(ChIP_TSS_B1_mat_log2FC))] %>% 
  reshape::melt() %>% cbind(., Pos = "TSS") %>% 
  dplyr::mutate(X2 = factor(gsub("_Trp-9", "", X2), 
                            levels = c("Pol II-NTD", "Pol II-S2p", "Pol II-S5p")))

t_tests2 = dat_Trp2 %>%
  group_by(X2) %>%
  summarise(P = t.test(value, mu = 0)$p.value,
            Label = paste("p <", formatC(P, format = "e", digits = 0)),
            Sig = ifelse(P < 0.05, "*", "ns"))

g2 <- 
  ggplot(dat_Trp2, aes(x = X2, y = value, color = X2)) +
  geom_hline(yintercept =  0, color = "red", alpha = 0.5, size = 1) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_boxplot(outlier.colour = NA) +
  geom_text(data = t_tests2, aes(label = Label, y = 3.7), size = 3.5) +
  scale_color_manual(values = sample_colors[levels(dat_Trp$X2)]) +
  xlab("") + ylab(" \n") + ggtitle(" ") +
  scale_y_continuous(breaks = c(-4, -2, 0, 2), labels = c(-4, -2, 0, 2), limits = c(-5, 3.7)) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 15), 
        legend.position = "none",
        panel.grid = element_blank())

ggsave(grid.arrange(g1, g2, nrow = 1, widths = c(5, 3)),
       filename = "../figS2/figs/FigS2_boxplot_ChIP_log2FC_Trp.png", 
       device = "png", width = 9, height = 5)


# -------------------------------------- Ring1b change vs TT-seq change ------------------------------------------ #
tmp_gene_ids_1 <- gene.gr$gene_id %>% 
  # intersect.Vector(., rownames(log2FC_mat)) %>% 
  intersect.Vector(., rownames(ChIP_TSS_B1_mat_log2FC)) %>% 
  intersect.Vector(., rownames(ChIP_TSS_B2_mat_log2FC))

dat_1 <- data.frame(Ring1b_mean = log1p(rowMeans(ChIP_TSS_B1_mat[tmp_gene_ids_1, c("Ring1b_NT", "Ring1b_BAP1-6", "Ring1b_BAP1-12")])),
                    Ring1b_log2FC_P12 = ChIP_TSS_B1_mat_log2FC[tmp_gene_ids_1, "Ring1b_BAP1-12"],
                    
                    H2Aub_mean = log1p(rowMeans(ChIP_TSS_B1_mat[tmp_gene_ids_1, c("H2Aub_NT", "H2Aub_BAP1-6", "H2Aub_BAP1-12")])),
                    H2Aub_log2FC_P12 = ChIP_TSS_B1_mat_log2FC[tmp_gene_ids_1, "H2Aub_BAP1-12"],
                    
                    RNA_log2FC_P12 = res.gene.list$P12[tmp_gene_ids_1, "log2FoldChange"],
                    
                    TT_log2FC = TTseq_res.list$res_E0B12[tmp_gene_ids_1, "log2FoldChange"])


dat_1 <- dat_1[complete.cases(dat_1) & is.finite(rowSums(dat_1)), ]
dat_1_Ring1b <- dat_1[intersect(Ring1b_enriched_genes, rownames(dat_1)), ]
# dat_1 <- dat_1[!rownames(dat_1) %in% Ring1b_enriched_genes, ]

# Ring1b density ~ Ring1b log2FC
g1 <- ggplot(dat_1, 
       aes(x = Ring1b_mean, y = Ring1b_log2FC_P12,
           color = get_dens(X = Ring1b_mean, Y = Ring1b_log2FC_P12)) 
       ) +
  geom_vline(xintercept = 0, lty = 1, alpha = 0.2) +
  geom_hline(yintercept = 0, lty = 1, alpha = 0.2) +
  geom_point(size = 1, alpha = 0.5) + 
  geom_point(data = dat_1_Ring1b,
             mapping = aes(x = Ring1b_mean, y = Ring1b_log2FC_P12), color = "blue", alpha = 0.2, pch = 4) +
  scale_color_viridis(option = "B", direction = -1, begin = 0.2, end = 0.9) +
  xlab("Ring1b TSS density (log2)") + ylab("Ring1b TSS log2FC at P12") + 
  annotate(geom = "text", x = -Inf, y = Inf, hjust = 0, vjust = 1,
           label = paste("All genes\n", 
                         "r =", round(cor(dat_1$Ring1b_mean, dat_1$Ring1b_log2FC_P12, use = "complete.obs"), 2),
                         "\nn =", nrow(dat_1))) +
  annotate(geom = "text", x = Inf, y = Inf, hjust = 1, vjust = 1, color = "blue",
           label = paste("Ring1b enriched genes\n", 
                         "r =", round(cor(dat_1_Ring1b$Ring1b_mean, dat_1_Ring1b$Ring1b_log2FC_P12, use = "complete.obs"), 2),
                         "\nn =", nrow(dat_1_Ring1b))) +
  ggpubr::theme_pubclean() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5))


# TT-seq log2FC ~ Ring1b log2FC
g2 <- ggplot(dat_1, 
       aes(x = TT_log2FC, y = Ring1b_log2FC_P12,
           color = get_dens(X = TT_log2FC, Y = Ring1b_log2FC_P12)) 
) +
  geom_vline(xintercept = 0, lty = 1, alpha = 0.2) +
  geom_hline(yintercept = 0, lty = 1, alpha = 0.2) +
  geom_point(size = 1, alpha = 0.5) + 
  geom_point(data = dat_1_Ring1b,
             mapping = aes(x = TT_log2FC, y = Ring1b_log2FC_P12), color = "blue", alpha = 0.2, pch = 4) +
  scale_color_viridis(option = "B", direction = -1, begin = 0.2, end = 0.9) +
  xlab("TT-seq log2FC at P12") + ylab("Ring1b TSS log2FC at P12") + 
  annotate(geom = "text", x = -Inf, y = Inf, hjust = 0, vjust = 1,
           label = paste("All genes\n", 
                         "r =", round(cor(dat_1$TT_log2FC, dat_1$Ring1b_log2FC_P12, use = "complete.obs"), 2),
                         "\nn =", nrow(dat_1))) +
  annotate(geom = "text", x = Inf, y = Inf, hjust = 1, vjust = 1, color = "blue",
           label = paste("Ring1b enriched genes\n", 
                         "r =", round(cor(dat_1_Ring1b$TT_log2FC, dat_1_Ring1b$Ring1b_log2FC_P12, use = "complete.obs"), 2),
                         "\nn =", nrow(dat_1_Ring1b))) +
  ggpubr::theme_pubclean() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5))


# H2Aub density ~ H2Aub log2FC
g3 <- ggplot(dat_1, 
       aes(x = H2Aub_mean, y = H2Aub_log2FC_P12,
           color = get_dens(X = Ring1b_mean, Y = H2Aub_log2FC_P12)) 
) +
  geom_vline(xintercept = 0, lty = 1, alpha = 0.2) +
  geom_hline(yintercept = 0, lty = 1, alpha = 0.2) +
  geom_point(size = 1, alpha = 0.5) + 
  geom_point(data = dat_1_Ring1b,
             mapping = aes(x = H2Aub_mean, y = H2Aub_log2FC_P12), color = "blue", alpha = 0.2, pch = 4) +
  scale_color_viridis(option = "B", direction = -1, begin = 0.2, end = 0.9) +
  xlab("H2Aub TSS density (log2)") + ylab("H2Aub TSS log2FC at P12") + 
  annotate(geom = "text", x = -Inf, y = Inf, hjust = 0, vjust = 1,
           label = paste("All genes\n", 
                         "r =", round(cor(dat_1$H2Aub_mean, dat_1$H2Aub_log2FC_P12, use = "complete.obs"), 2),
                         "\nn =", nrow(dat_1))) +
  annotate(geom = "text", x = Inf, y = Inf, hjust = 1, vjust = 1, color = "blue",
           label = paste("Ring1b enriched genes\n", 
                         "r =", round(cor(dat_1_Ring1b$H2Aub_mean, dat_1_Ring1b$H2Aub_log2FC_P12, use = "complete.obs"), 2),
                         "\nn =", nrow(dat_1_Ring1b))) +
  ggpubr::theme_pubclean() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5))


# TT-seq log2FC ~ H2Aub log2FC
g4 <- ggplot(dat_1, 
       aes(x = TT_log2FC, y = H2Aub_log2FC_P12,
           color = get_dens(X = TT_log2FC, Y = H2Aub_log2FC_P12)) 
) +
  geom_vline(xintercept = 0, lty = 1, alpha = 0.2) +
  geom_hline(yintercept = 0, lty = 1, alpha = 0.2) +
  geom_point(size = 1, alpha = 0.5) + 
  geom_point(data = dat_1_Ring1b,
             mapping = aes(x = TT_log2FC, y = H2Aub_log2FC_P12), color = "blue", alpha = 0.2, pch = 4) +
  scale_color_viridis(option = "B", direction = -1, begin = 0.2, end = 0.9) +
  xlab("TT-seq log2FC at P12") + ylab("H2Aub TSS log2FC at P12") + 
  annotate(geom = "text", x = -Inf, y = Inf, hjust = 0, vjust = 1,
           label = paste("All genes\n", 
                         "r =", round(cor(dat_1$TT_log2FC, dat_1$H2Aub_log2FC_P12, use = "complete.obs"), 2),
                         "\nn =", nrow(dat_1))) +
  annotate(geom = "text", x = Inf, y = Inf, hjust = 1, vjust = 1, color = "blue",
           label = paste("Ring1b enriched genes\n", 
                         "r =", round(cor(dat_1_Ring1b$TT_log2FC, dat_1_Ring1b$H2Aub_log2FC_P12, use = "complete.obs"), 2),
                         "\nn =", nrow(dat_1_Ring1b))) +
  ggpubr::theme_pubclean() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5))


ggsave(grid.arrange(g1, g2, ncol = 2),
       filename = "FigS2_scatter_Ring1b_ChIP_TT_log2FC_correlation.png", 
       path = "../figS2/figs",
       device = "png", width = 8, height = 3)

ggsave(grid.arrange(g3, g4, ncol = 2),
       filename = "FigS2_scatter_H2Aub_ChIP_TT_log2FC_correlation.png", 
       path = "../figS2/figs",
       device = "png", width = 8, height = 3)



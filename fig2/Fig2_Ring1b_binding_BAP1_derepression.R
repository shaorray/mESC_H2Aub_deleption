# Rui Shao 2021 May
# Figure 2
# Ring1b binding and BAP1 pulse derepression

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")

all_ChIP_TSS_mat <- readRDS("../fig1/data/all_ChIP_TSS_mat.RData")

# input normalization
all_ChIP_rerun_TSS_mat <-
  .countBW(bw_files = list.files(path = "/mnt/E0767589767560E8/UPPMAX/PHILIP_LOF3_20210205R",
                                 pattern = "unscale.*bw$", full.names = T),
           intervals = promoters(gene.gr[use_gene_ids], upstream = 1000, downstream = 1000),
           blacklist = blacklist.gr,
           fast = F) 
colnames(all_ChIP_rerun_TSS_mat) <- gsub(".mm9.*", "", colnames(all_ChIP_rerun_TSS_mat))
all_ChIP_rerun_TSS_mat <- 
  all_ChIP_rerun_TSS_mat[, !grepl("sfI", colnames(all_ChIP_rerun_TSS_mat))] /
  all_ChIP_rerun_TSS_mat[, rep(which(grepl("sfI", colnames(all_ChIP_rerun_TSS_mat))), 7)]

all_ChIP_rerun_TSS_mat_cmb <- sapply(unique(gsub("_rep.", "", colnames(all_ChIP_rerun_TSS_mat))), 
                                     function(x) 
                                       rowMeans(all_ChIP_rerun_TSS_mat[, paste0(x, "_rep", 1:2)]))

# ------------------------------ChIP occupancy changes------------------------------ #
# z-score tiles of median occupancy during BAP1 pulse-chase
dat <- all_ChIP_TSS_mat[, grep("_RC", colnames(all_ChIP_TSS_mat))] %>% log()

brewer_cols <- c("Spectral", "RdYlGn",  "RdYlBu",  "RdGy",  "RdBu", "PuOr", "PRGn", "PiYG", "BrBG")

sample_names <- c("H2Aub", "Ring1b", "H3K27m3", "Ezh2","Pol2-NTD", "Pol2-S5p", "Pol2-S2p", "H3K4m3")

g_list1 <- g_list2 <- list()
for (i in seq_along(sample_names)) {
  tmp <- dat[, grep(sample_names[i], colnames(dat))] %>% as.matrix()
  keep.idx <- apply(tmp, 1, function(x) !any(is.na(x) | is.infinite(x)))
  tmp_box <- cbind(tmp, log2FC_cls) %>% 
    as.data.frame() %>% dplyr::filter(keep.idx) %>% 
    reshape::melt(id.vars = "log2FC_cls") %>% 
    dplyr::mutate(log2FC_cls = factor(paste0("C", log2FC_cls)))
  
  ggplot(tmp_box, aes(x = log2FC_cls, y = value, fill = variable)) +
    geom_boxplot(notch = T, outlier.size = 0) +
    xlab("") + ylab("log TSS Occupancy") + ggtitle(sample_names[i]) +
    ylim(quantile(tmp_box$value, c(0.0002, 0.9998))) +
    scale_fill_viridis_d(begin = 0.1, end = 0.4, option = "C") +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          strip.text.y = element_blank(), 
          legend.position = "none") -> g1
  
  g_list1 %c=% list(g1)
  
  # # z-score tile
  # tmp <- (tmp - mean(tmp[keep.idx, ])) / sd(tmp[keep.idx, ])
  # tmp_cls_median <- apply(tmp[keep.idx, ], 2, 
  #                         function(x) aggregate(x, list(log2FC_cls[keep.idx]), "mean")[, 2])
  # 
  # tmp_cls_sd <- apply(tmp, 2, function(x) aggregate(x, list(log2FC_cls), "sd")[, 2])
  # 
  # rownames(tmp_cls_median) <- c("Derepressed", "Weakly Derepressed", "Unchanged", "Weakly Repressed")
  # reshape::melt(tmp_cls_median) %>% 
  #   dplyr::mutate(X1 = factor(X1, levels = rev(c("Derepressed", "Weakly Derepressed",
  #                                                "Unchanged", "Weakly Repressed")))) %>%
  #   dplyr::mutate(X2 = factor(rep(c("P0", "P6", "P12"), each = 4), 
  #                             levels = c("P0", "P6", "P12"))) %>%
  #   ggplot(aes(x = X2, y = X1, fill = value)) +
  #   geom_tile(width = 0.99) +
  #   scale_fill_gradientn(name = "Z-scores",
  #                        colours = c("grey50", "white", ifelse(i < 5, "red", "darkgreen")), #RColorBrewer::brewer.pal(n = 10, brewer_cols[ifelse(i < 5, 4, 9)]),
  #                        values = (c(min(tmp_cls_median), 0, max(tmp_cls_median)) - min(tmp_cls_median)) /
  #                          diff(range(tmp_cls_median))
  #   ) +
  #   xlab("") + ylab("") + ggtitle(sample_names[i]) +
  #   theme_minimal() +
  #   theme(panel.grid = element_blank(), 
  #         axis.text.x = element_text(size = 12),
  #         axis.text.y = element_blank()) -> g2
  # g_list2 %c=% list(g2)
}

# ggsave(plot = do.call(grid.arrange, 
#                       c(g_list1, 
#                         ncol = 4)),
#        filename = paste0("Fig2_ChIP_TSS_occupancy_DE_groups.png"), 
#        path = "figs/",
#        device = "png", width = 11, height = 5)

# ggsave(plot = do.call(grid.arrange, 
#                       c(g_list2, 
#                         ncol = 4)),
#        filename = paste0("Fig2_ChIP_Z_scores_DE_groups.png"), 
#        path = "figs",
#        device = "png", width = 11, height = 5)

# batch2
dat <- all_ChIP_TSS_mat[, grep("Cbx7|Rybp", colnames(all_ChIP_TSS_mat))] %>% log()
dat <- dat[, -grep("C12|C36", colnames(dat))]
sample_names <- c("Cbx7", "Rybp")

for (i in seq_along(sample_names)) {
  tmp <- dat[, grepl(sample_names[i], colnames(dat))] %>% as.matrix()
  keep.idx <- apply(tmp, 1, function(x) !any(is.na(x) | is.infinite(x)))
  tmp_box <- cbind(tmp, log2FC_cls) %>% 
    as.data.frame() %>% dplyr::filter(keep.idx) %>% 
    reshape::melt(id.vars = "log2FC_cls") %>% 
    dplyr::mutate(log2FC_cls = factor(paste0("C", log2FC_cls)))
  
  ggplot(tmp_box, aes(x = log2FC_cls, y = value, fill = variable)) +
    geom_boxplot(notch = T, outlier.size = 0) +
    xlab("") + ylab("log TSS Occupancy") + ggtitle(sample_names[i]) +
    ylim(quantile(tmp_box$value, c(0.0002, 0.9998))) +
    scale_fill_viridis_d(begin = 0.1, end = 0.7, option = "C") +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          strip.text.y = element_blank(), 
          legend.position = "none") -> g
  
  g_list1 %c=% list(g)
}

ggsave(plot = do.call(grid.arrange, 
                      c(g_list1[c(1:4,9,5:8,10)], 
                        ncol = 5)),
       filename = paste0("Fig2_ChIP_TSS_occupancy_DE_groups2.png"), 
       path = "figs/",
       device = "png", width = 11, height = 5)


# -----------------calculate the rate of Ring1b/Ezh2 binding changes--------------------- #
# in batch 1 (P6 ~ P12) and batch 2 (P12 ~ P12C36)
source("F2_src_ChIP_change_rates.R")
Ring1b_rates_tss <- get_change_rates("Ring1b") %>% dplyr::filter(complete.cases(.))
Ezh2_rates_tss <- get_change_rates("Ezh2")
H3K27me3_rates_tss <- get_change_rates("H3K27me3")
Cbx7_rates_tss <- get_change_rates("Cbx7")
Rybp_rates_tss <- get_change_rates("Rybp")

Ring1b_binding_genes <-
  rownames(all_ChIP_TSS_mat) %>%
  "["(kink_index(all_ChIP_TSS_mat[, "Ring1b_P0_ALL"],
                 method = "poisson")) %>% 
  intersect.Vector(rownames(Ring1b_rates_tss)) %>%
  intersect.Vector(use_gene_ids)

g1.1 <- 
  ggplot(Ring1b_rates_tss[Ring1b_binding_genes, ],
         aes(x = NT, y = Rate, 
             color = get_dens(NT, Rate))) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  annotate(geom = "text", 
           x = 5, y = 0.05,
           hjust = "right", vjust = "top",
           label = paste("r = ", round(with(Ring1b_rates_tss[Ring1b_binding_genes, ], 
                                            cor(NT, Rate)),
                                       3))) +
  xlab("Initial level 1 (P0)") + 
  ylab("Rate 1 (P0-P6-P12)") +
  ggtitle("Ring1b target") +
  theme_setting +
  theme(legend.position = "none")

g1.2 <- 
  ggplot(Ring1b_rates_tss[Ring1b_binding_genes, c("P0", "P12_Rate")] %>% 
           filter(complete.cases(.)),
         aes(x = P0, y = P12_Rate, 
             color = get_dens(P0, P12_Rate))) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  annotate(geom = "text", 
           x = 4, y = 0.05,
           hjust = "right", vjust = "top",
           label = paste("r = ", 
                         round(with(Ring1b_rates_tss[Ring1b_binding_genes, 
                                                     c("P0", "P12_Rate")] %>%
                                      filter(complete.cases(.)), 
                                    cor(P0, P12_Rate)),
                                       3))) +
  xlab("Initial level 2 (P0)") + 
  ylab("Rate 2 (P0-P12-P12C12-P12C24)") +
  ggtitle("Ring1b target") +
  theme_setting +
  theme(legend.position = "none")

g1.3 <- ggplot(Ring1b_rates_tss[Ring1b_binding_genes, c("Rate", "P12_Rate")],
       aes(x = Rate, y = P12_Rate, 
           color = get_dens(Rate, P12_Rate))) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  annotate(geom = "text", 
           x = 0.07, y = 0.05,
           hjust = "right", vjust = "top",
           label = paste0(paste("r = ", round(with(Ring1b_rates_tss[Ring1b_binding_genes, ], 
                                                   cor(Rate, P12_Rate)),
                                              3)),
                          "\nn = ", length(Ring1b_binding_genes))) +
  xlab("Rate 1 (P0-P6-P12)") + 
  ylab("Rate 2 (P0-P12)") +
  ggtitle("Ring1b target") +
  theme_setting +
  theme(legend.position = "none")
  
ggsave(plot = grid.arrange(g1.1, g1.2, g1.3, ncol = 3),
       filename = paste0("Fig2_point_Ring1b_binding_change_rate.png"), 
       path = "figs/",
       device = "png", width = 12, height = 4) 

pheatmap::pheatmap(all_ChIP_rerun_TSS_mat_cmb[Ring1b_binding_genes,
                                              grep("Ring1b", colnames(all_ChIP_rerun_TSS_mat_cmb))] %>% log,
                   show_rownames = F, cluster_cols = F)

# ------------------------------cluster Ring1b occupancy------------------------------ #
set.seed(1)
Ring1b_PC_ChIP_TSS_mat <- 
  all_ChIP_rerun_TSS_mat_cmb[Ring1b_binding_genes,
                             grep("Ring1b", colnames(all_ChIP_rerun_TSS_mat_cmb))] %>%
  log() %>%
  trim_quantile()

Ring1b_cls <- kmeans(t(scale(t(Ring1b_PC_ChIP_TSS_mat))), centers = 3)
Ring1b_cls_M <- mclust::Mclust(t(scale(t(Ring1b_PC_ChIP_TSS_mat))), G = 5)

Ring1b_cls$cluster <- order_cls(Ring1b_cls$cluster, rowMeans(Ring1b_PC_ChIP_TSS_mat[, 1:3]))
# Ring1b_cls_M$classification <- order_cls(Ring1b_cls_M$classification,
#                                          rowMeans(Ring1b_PC_ChIP_TSS_mat[, 1:3]))
Ring1b_cls_M$classification <- factor(Ring1b_cls_M$classification, levels = c(2, 5, 4, 3, 1)) %>% 
  as.numeric()

ord_val <- Ring1b_PC_ChIP_TSS_mat[, 1]
ord_val[Ring1b_cls_M$classification == 2] <- 
  Ring1b_PC_ChIP_TSS_mat[Ring1b_cls_M$classification == 2, 2]
ord_val[Ring1b_cls_M$classification == 3] <- 
  Ring1b_PC_ChIP_TSS_mat[Ring1b_cls_M$classification == 3, 3]
ord_val[Ring1b_cls_M$classification == 4] <- 
  Ring1b_PC_ChIP_TSS_mat[Ring1b_cls_M$classification == 4, 4]
ord_val[Ring1b_cls_M$classification == 5] <- 
  Ring1b_PC_ChIP_TSS_mat[Ring1b_cls_M$classification == 5, 5]

# Ezh2_PC_ChIP_TSS_mat <- all_ChIP_TSS_mat[, grep("Ezh2", colnames(all_ChIP_TSS_mat))] %>%
#   trim_quantile()
# 
# Ezh2_cls <- kmeans(Ezh2_PC_ChIP_TSS_mat[, 1:3], centers = 3)
# 
# Ezh2_cls$cluster <- order_cls(Ezh2_cls$cluster, rowMeans(Ezh2_PC_ChIP_TSS_mat[, 1:3]))


g2.1 <- cbind(Ring1b_PC_ChIP_TSS_mat, Ring1b_cls_M$classification) %>% 
  as.data.frame() %>% `colnames<-`(c("P0", "P12", "P12C12", "P12C24", "P12C36", "cls")) %>%
  arrange(cls, ord_val) %>%
  cbind(id = seq_len(nrow(Ring1b_PC_ChIP_TSS_mat))) %>%
  `rownames<-`(NULL) %>% 
  reshape::melt(id.vars = c("cls", "id")) %>% 
  ggplot(aes(x = variable, y = id, fill = value)) +
  geom_tile(width = 0.99) +
  facet_grid(cls ~ ., scales = "free", space = "free") +
  labs(title = 'Ring1b binding genes (n = 1180)',  x = '', y = '') +
  scale_fill_gradientn(name = "Ring1b", 
                       colours = c("deepskyblue", "white", 'brown4', 'grey15') ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_blank(), 
        strip.text.y = element_blank(), 
        legend.position = "left",
        legend.key.size = unit(13, "pt"), 
        plot.margin = unit(c(1, 1, 0, 1), "lines")) 

g2.2 <- (table(Ring1b_cls_M$classification, log2FC_cls[Ring1b_binding_genes]) /
  c(table(Ring1b_cls_M$classification))) %>%
  as.data.frame() %>%
  dplyr::mutate(Var1 = factor(Var1, levels = 5:1)) %>%
  dplyr::mutate(Var2 = factor(Var2, levels = c(1,2,4,3))) %>%
  dplyr::mutate(Freq = Freq * ifelse(Var2 %in% c(3, 4), -1, 1)) %>%
  ggplot(aes(x = Freq, y = Var1, fill = Var2)) +
  geom_bar(stat="identity", width = 0.6) +
  scale_x_continuous(breaks = c(-0.4, 0, 0.4), labels = c(0.4, 0, 0.4)) +
  scale_fill_manual(values = `names<-`(viridis(n = 4, option = "E", direction = -1)[c(1,2,4,3)], c(1,2,4,3)), 
                    breaks = c(1,2,3,4),
                    labels = c("Derepressed", "W.Derepressed", "Unchanged", "Repressed")) +
  xlab("\nFrequency") + ylab("Ring1b binding class\n") + labs(fill = "Gene expression") +
  ggpubr::theme_pubclean() +
  theme(legend.position = "right")

ggsave(plot = grid.arrange(g2.1, g2.2, ncol = 2),
       filename = paste0("Fig2_Ring1b_groups_log2FC_groups.png"), 
       path = "figs/",
       device = "png", width = 9, height = 7)

cbind(rowMeans(log2FC_mat[Ring1b_binding_genes, ]),
      rowMeans(Ring1b_PC_ChIP_TSS_mat[, c(3:5)]) - rowMeans(Ring1b_PC_ChIP_TSS_mat[, c(1,1)]),
      Ring1b_cls_M$classification) %>%
  as.data.frame() %>% `colnames<-`(c("RNA-seq log2FC", "Ring1b logFC", "Ring1b change")) %>%
  ggplot(aes(x = `RNA-seq log2FC`, y = `Ring1b logFC`, color = as.factor(`Ring1b change`))) +
  geom_point() +
  labs(color = "Ring1b binding class") +
  ggpubr::theme_pubclean()



library("ggforce")
PcG_binding_cls <- data.frame(Ring1b = as.character(Ring1b_cls$cluster),
                              Ezh2 = as.character(Ezh2_cls$cluster)) %>%
  dplyr::count(Ring1b, Ezh2, sort = TRUE) %>%
  ggforce::gather_set_data(1:2) %>%
  arrange(x, Ezh2, desc(Ring1b))

g2.3 <- ggplot(PcG_binding_cls, aes(x = rev(x), id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = Ring1b), alpha = 0.8, axis.width = 0.2,
                     n=100, strength = 0.5) +
  geom_parallel_sets_axes(axis.width = 0.25, fill = "gray90",  size = 0) +
  geom_parallel_sets_labels(colour = 'gray35', size = 0, angle = 0, fontface="bold") +
  scale_fill_manual(values  = viridis(4, end = 0.7, option = "A", direction = -1)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x  = element_blank(),
    plot.margin = unit(c(0.3, 0, 1.5, 0), "lines")
  )



# DE class ~ Ring1b class
DE_binding_cls <- data.frame(Ring1b = as.character(Ring1b_cls$cluster),
                            DE_cls = as.character(log2FC_cls)) %>%
  dplyr::count(Ring1b, DE_cls, sort = TRUE) %>%
  ggforce::gather_set_data(1:2) %>%
  arrange(x, DE_cls, desc(Ring1b))

ggplot(DE_binding_cls, aes(x = (x), id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = Ring1b), alpha = 0.8, axis.width = 0.2,
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
    # axis.text.x = element_blank(),
    axis.title.x  = element_blank(),
    plot.margin = unit(c(0.3, 0, 1.5, 0), "lines")
  )

ggsave(filename = paste0("Fig2_DE_cls_Ring1b_groups.png"), 
       path = "figs/",
       device = "png", width = 3, height = 7)

# cluster Ring1b binding and rate
library(mclust)
mcls_res <- Mclust(Ring1b_rates_tss[Ring1b_binding_genes, 
                                    grep("P.*Rate", colnames(Ring1b_rates_tss))], G = 3)

ggplot(Ring1b_rates_tss[Ring1b_binding_genes, ],
       aes(x = NT, y = Rate, 
           color = factor(mcls_res$classification))) +
  geom_point(size = 2) +
  scale_color_viridis_d() +
  xlab("Initial level 1 (P0)") + 
  ylab("Rate 1 (P0-P6-P12)") +
  ggtitle("Ring1b target") +
  theme_setting +
  theme(legend.position = "none")

data.frame(gene_ids = Ring1b_binding_genes, 
           cls = mcls_res$classification, 
           DE_group = log2FC_cls[Ring1b_binding_genes]) %>%
  dplyr::count(cls, DE_group, sort = TRUE) %>%
  ggforce::gather_set_data(1:2) %>%
  arrange(x, cls, desc(DE_group)) %>%
  ggplot(aes(x = x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = DE_group), alpha = 0.8, axis.width = 0.2,
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
    # axis.text.x = element_blank(),
    axis.title.x  = element_blank(),
    plot.margin = unit(c(0.3, 0, 1.5, 0), "lines")
  )


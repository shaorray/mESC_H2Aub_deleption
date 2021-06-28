setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
source("../util/getCoverage.R")

Ring1b_PC_ChIP_TSS_mat1 <- 
  .countBW(bw_files = list.files(path = "/mnt/E0767589767560E8/UPPMAX/BAP1_PC/",
                                 pattern = "^Ring1b_P", full.names = T),
           intervals = promoters(gene.gr[use_gene_ids], upstream = 1000, downstream = 1000),
           blacklist = blacklist.gr,
           fast = F)

Ring1b_PC_ChIP_TSS_mat2 <- 
  .countBW(bw_files = list.files(path = "/mnt/E0767589767560E8/UPPMAX/PHILIP_LOF3_20210205R/",
                                 pattern = "Ring1b", full.names = T),
           intervals = promoters(gene.gr[use_gene_ids], upstream = 1000, downstream = 1000),
           blacklist = blacklist.gr,
           fast = F)

all_PC_ChIP_TSS_mat <- 
  .countBW(bw_files = list.files(path = "/mnt/E0767589767560E8/UPPMAX/PHILIP_LOF3_20210205R/",
                                 pattern = ".bw$", full.names = T),
           intervals = promoters(gene.gr[use_gene_ids], upstream = 1000, downstream = 1000),
           blacklist = blacklist.gr,
           fast = F)
colnames(all_PC_ChIP_TSS_mat) <- gsub("_pooled.mm9.scaled.bw", "", colnames(all_PC_ChIP_TSS_mat))

# compare
Ring1b_binding_genes <-
  rownames(all_PC_ChIP_TSS_mat)[kink_index(all_PC_ChIP_TSS_mat[, "Ring1b_P0"], 
                                        method = "weight")]

g1 <- data.frame(x = Ring1b_PC_ChIP_TSS_mat1[Ring1b_binding_genes, 1],
           y = Ring1b_PC_ChIP_TSS_mat2[Ring1b_binding_genes, 1]) %>%
  ggplot(aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  xlab("P0 Previous") + ylab("P0 Rerun") +
  ggtitle("Ring1b TSS +/-1kb (Ring1b enriched genes)") +
  theme_classic()

g2 <- data.frame(x = Ring1b_PC_ChIP_TSS_mat1[Ring1b_binding_genes, 2],
           y = Ring1b_PC_ChIP_TSS_mat2[Ring1b_binding_genes, 2]) %>%
  ggplot(aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  xlab("P12 Previous") + ylab("P12 Rerun") +
  ggtitle("Ring1b TSS +/-1kb (Ring1b enriched genes)") +
  theme_classic()

g3 <- data.frame(x = Ring1b_PC_ChIP_TSS_mat1[Ring1b_binding_genes, 3],
                 y = Ring1b_PC_ChIP_TSS_mat2[Ring1b_binding_genes, 3]) %>%
  ggplot(aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  xlab("P12C12 Previous") + ylab("P12 Rerun") +
  ggtitle("Ring1bC12 TSS +/-1kb (Ring1b enriched genes)") +
  theme_classic()

g4 <- data.frame(x = Ring1b_PC_ChIP_TSS_mat1[Ring1b_binding_genes, 4],
           y = Ring1b_PC_ChIP_TSS_mat2[Ring1b_binding_genes, 4]) %>%
  ggplot(aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  xlab("P12C24 Previous") + ylab("P12 Rerun") +
  ggtitle("Ring1bC24 TSS +/-1kb (Ring1b enriched genes)") +
  theme_classic()

g5 <- data.frame(x = Ring1b_PC_ChIP_TSS_mat1[Ring1b_binding_genes, 5],
                 y = Ring1b_PC_ChIP_TSS_mat2[Ring1b_binding_genes, 5]) %>%
  ggplot(aes(x = x, y = y)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  xlab("P12C36 Previous") + ylab("P12 Rerun") +
  ggtitle("Ring1bC36 TSS +/-1kb (Ring1b enriched genes)") +
  theme_classic()

g6 <- rbind(Ring1b_PC_ChIP_TSS_mat1[Ring1b_binding_genes, ] %>%
  `colnames<-`(c("P0", "P12", "P12C12", "P12C24", "P12C36")) %>% 
  reshape::melt() %>% add_column(run = "Previous"),
  Ring1b_PC_ChIP_TSS_mat2[Ring1b_binding_genes, ] %>%
    `colnames<-`(c("P0", "P12", "P12C12", "P12C24", "P12C36")) %>% 
    reshape::melt() %>% add_column(run = "Rerun")) %>%
  ggplot(aes(x = variable, y = value, color = run)) +
  geom_boxplot(outlier.size = 0) +
  scale_y_log10() +
  xlab("") + ylab("Ring1b TSS +/-1kb occupancy") + ggtitle("Ring1b enriched genes") +
  theme_classic()

ggsave(plot = grid.arrange(g1, g2, g3, g4, g5, g6, 
                           nrow = 3),
       filename = paste0("test_Ring1b_rerun.png"), 
       path = "figs/",
       device = "png", width = 8, height = 12)

# heatmap
cov_mat_Ring1b <- convertCoverage(file_names = list.files(path = "/mnt/E0767589767560E8/UPPMAX/PHILIP_LOF3_20210205R/",
                                                          pattern = "Ring1b", full.names = T), 
                                 intervals = promoters(gene.gr[use_gene_ids], 5000, 5000),
                                 bin_width = 50, 
                                 new_len = 100)
if (T) {
  bg_val <- 1
  P0_max <-  quantile(log2(cov_mat_Ring1b[, ,1] + bg_val), 0.996)
  row_order <- order(rowMeans(log2(cov_mat_Ring1b[, ,1] + bg_val)))
  
  g2.1 <- plot_TSS_heatmap(log2(cov_mat_Ring1b[row_order, ,1] + bg_val), 
                           .col = "darkgreen",
                           .max = P0_max, 
                           .min = log2(bg_val))
  g2.2 <- plot_TSS_heatmap(log2(cov_mat_Ring1b[row_order, ,2] + bg_val), 
                           .col = "darkgreen",
                           .max = P0_max, 
                           .min = log2(bg_val))
  g2.3 <- plot_TSS_heatmap(log2(cov_mat_Ring1b[row_order, ,3] + bg_val),
                           .col = "darkgreen",
                           .max = P0_max, 
                           .min = log2(bg_val))
  g2.4 <- plot_TSS_heatmap(log2(cov_mat_Ring1b[row_order, ,4] + bg_val),
                           .col = "darkgreen",
                           .max = P0_max, 
                           .min = log2(bg_val))
  g2.5 <- plot_TSS_heatmap(log2(cov_mat_Ring1b[row_order, ,5] + bg_val), 
                           .col = "darkgreen",
                           .max = P0_max, 
                           .min = log2(bg_val), .show_legend = T)
  
  ggsave(plot = grid.arrange(g2.1, g2.2, g2.3, g2.4, g2.5,
                             nrow = 1, widths = c(3, 3, 3, 3, 4.2)),
         filename = paste0("test_TSS_Hetmap_Ring1b_PC.png"), 
         path = "figs/",
         device = "png", width = 10, height = 4)
}






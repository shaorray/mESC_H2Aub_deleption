# ----------------------------------------------------------------------------------- #
bw_files_H2Aub <- list.files(path = "/mnt/E0767589767560E8/UPPMAX/BAP1_PC",
                             pattern = "H2Aub_P.*bw$", full.names = T)

cov_mat_H2Aub <- convertCoverage(file_names = bw_files_H2Aub, 
                                 intervals = promoters(gene.gr[use_gene_ids], 5000, 5000),
                                 bin_width = 50, 
                                 new_len = 100)

plot_TSS_heatmap <- function(mat, 
                             .title = "", .col = "black", 
                             .max = 4, .min = log2(0.01), 
                             .show_legend = F) 
{
  
  `rownames<-`(mat, NULL) %>% 
    trim_quantile() %>%
    reshape::melt() %>% 
    ggplot(aes(X2, X1, fill=value)) +
    geom_tile() +
    scale_x_continuous(expand = c(0,0),
                       breaks = c(1, 51, 100), 
                       labels = c('-5 kb', 'TSS', '+5 kb')) +
    scale_y_continuous(expand = c(0,0), 
                       breaks = c(1, nrow(mat))) +
    scale_fill_gradient(name = "",
                        low = "white",
                        high = .col,
                        limits = c(.min, .max), 
                        breaks = c(.min, .max),
                        labels = c("", "")) +
    labs(title = .title,  x = '', y = '') +
    # guides(fill = guide_legend(title = "")) +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, size = 1),
          axis.ticks.x = element_line(), 
          axis.text.x = element_text(size = 12),
          axis.text.y = element_blank(),
          plot.margin = unit(c(1, 1, 0, 1), "lines")) -> g
  
  if (.show_legend) {
    return(g)
  } else {
    return(g + theme(legend.position = "none"))
  }
}

plot_TSS_coverage <- function(mat, .col = 1, .min, .max) {
  data.frame(x = seq_len(ncol(mat)), y = colMeans(mat, na.rm = T)) %>%
    ggplot(aes(x = x, y = y)) +
    geom_line(size = 1, color = .col) + 
    scale_y_continuous(limits = c(.min, .max)) +
    xlab("") + ylab("log2 Coverage") +
    theme_setting +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}


if (T) {
  bg_val <- 1
  P0_max <-  quantile(log2(cov_mat_H2Aub[, ,1] + bg_val), 0.996)
  row_order <- order(rowMeans(log2(cov_mat_H2Aub[, ,1] + bg_val)))
  
  c_min <- apply(log2(cov_mat_H2Aub + bg_val), 3, function(x) min(colMeans(x))) %>% min()
  c_max <- apply(log2(cov_mat_H2Aub + bg_val), 3, function(x) max(colMeans(x))) %>% max()
  
  # g2.1c <- plot_TSS_coverage(log2(cov_mat_H2Aub[row_order, ,1] + bg_val), 
  #                      .max = c_max, 
  #                      .min = c_min)
  # g2.2c <- plot_TSS_coverage(log2(cov_mat_H2Aub[row_order, ,2] + bg_val), 
  #                            .max = c_max, 
  #                            .min = c_min)
  # g2.3c <- plot_TSS_coverage(log2(cov_mat_H2Aub[row_order, ,3] + bg_val),
  #                            .max = c_max, 
  #                            .min = c_min)
  # g2.4c <- plot_TSS_coverage(log2(cov_mat_H2Aub[row_order, ,4] + bg_val),
  #                            .max = c_max, 
  #                            .min = c_min)
  # g2.5c <- plot_TSS_coverage(log2(cov_mat_H2Aub[row_order, ,5] + bg_val), 
  #                            .max = c_max, 
  #                            .min = c_min)
  
  g2.1 <- plot_TSS_heatmap(log2(cov_mat_H2Aub[row_order, ,1] + bg_val), 
                           .max = P0_max, 
                           .min = log2(bg_val))
  g2.2 <- plot_TSS_heatmap(log2(cov_mat_H2Aub[row_order, ,2] + bg_val), 
                           .max = P0_max, 
                           .min = log2(bg_val))
  g2.3 <- plot_TSS_heatmap(log2(cov_mat_H2Aub[row_order, ,3] + bg_val),
                           .max = P0_max, 
                           .min = log2(bg_val))
  g2.4 <- plot_TSS_heatmap(log2(cov_mat_H2Aub[row_order, ,4] + bg_val),
                           .max = P0_max, 
                           .min = log2(bg_val))
  g2.5 <- plot_TSS_heatmap(log2(cov_mat_H2Aub[row_order, ,5] + bg_val), 
                           .max = P0_max, 
                           .min = log2(bg_val), .show_legend = T)
  
  ggsave(plot = grid.arrange(g2.1, g2.2, g2.3, g2.4, g2.5,
                             nrow = 1, widths = c(3, 3, 3, 3, 4.2)),
         filename = paste0("Fig1_TSS_Hetmap_H2Aub_PC.png"), 
         path = "figs/",
         device = "png", width = 10, height = 4)
}

# ----------------------------------------------------------------------------------- #
bw_files_Ring1b <- list.files(path = "/mnt/E0767589767560E8/UPPMAX/BAP1_PC",
                              pattern = "Ring1b_P.*bw$", full.names = T)

cov_mat_Ring1b <- convertCoverage(file_names = bw_files_Ring1b, 
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
         filename = paste0("Fig1_TSS_Hetmap_Ring1b_PC.png"), 
         path = "figs/",
         device = "png", width = 10, height = 4)
}

# ----------------------------------------------------------------------------------- #
bw_files_Ezh2 <- list.files(path = "/mnt/E0767589767560E8/UPPMAX/BAP1_PC",
                            pattern = "Ezh2_P.*bw$", full.names = T)

cov_mat_Ezh2 <- convertCoverage(file_names = bw_files_Ezh2, 
                                intervals = promoters(gene.gr[use_gene_ids], 5000, 5000),
                                bin_width = 50, 
                                new_len = 100)

if (T) {
  bg_val <- 1
  P0_max <-  quantile(log2(cov_mat_Ezh2[, ,1] + bg_val), 0.996)
  row_order <- order(rowMeans(log2(cov_mat_Ezh2[, ,1] + bg_val)))
  
  g2.1 <- plot_TSS_heatmap(log2(cov_mat_Ezh2[row_order, ,1] + bg_val), 
                           .col = "darkred",
                           .max = P0_max, 
                           .min = log2(bg_val))
  g2.2 <- plot_TSS_heatmap(log2(cov_mat_Ezh2[row_order, ,2] + bg_val), 
                           .col = "darkred",
                           .max = P0_max, 
                           .min = log2(bg_val))
  g2.3 <- plot_TSS_heatmap(log2(cov_mat_Ezh2[row_order, ,3] + bg_val),
                           .col = "darkred",
                           .max = P0_max, 
                           .min = log2(bg_val))
  g2.4 <- plot_TSS_heatmap(log2(cov_mat_Ezh2[row_order, ,4] + bg_val),
                           .col = "darkred",
                           .max = P0_max, 
                           .min = log2(bg_val))
  g2.5 <- plot_TSS_heatmap(log2(cov_mat_Ezh2[row_order, ,5] + bg_val), 
                           .col = "darkred",
                           .max = P0_max, 
                           .min = log2(bg_val), .show_legend = T)
  
  ggsave(plot = grid.arrange(g2.1, g2.2, g2.3, g2.4, g2.5,
                             nrow = 1, widths = c(3, 3, 3, 3, 4.2)),
         filename = paste0("Fig1_TSS_Hetmap_Ezh2_PC.png"), 
         path = "figs/",
         device = "png", width = 10, height = 4)
}

# ----------------------------------------------------------------------------------- #
bw_files_H3K27me3 <- list.files(path = "/mnt/E0767589767560E8/UPPMAX/BAP1_PC",
                                pattern = "H3K27me3_P.*bw$", full.names = T)

cov_mat_H3K27me3 <- convertCoverage(file_names = bw_files_H3K27me3, 
                                    intervals = promoters(gene.gr[use_gene_ids], 5000, 5000),
                                    bin_width = 50, 
                                    new_len = 100)

if (T) {
  bg_val <- 1
  P0_max <-  quantile(log2(cov_mat_H3K27me3[, ,1] + bg_val), 0.996)
  row_order <- order(rowMeans(log2(cov_mat_H3K27me3[, ,1] + bg_val)))
  
  g2.1 <- plot_TSS_heatmap(log2(cov_mat_H3K27me3[row_order, ,1] + bg_val), 
                           .col = "midnightblue",
                           .max = P0_max, 
                           .min = log2(bg_val))
  g2.2 <- plot_TSS_heatmap(log2(cov_mat_H3K27me3[row_order, ,2] + bg_val), 
                           .col = "midnightblue",
                           .max = P0_max, 
                           .min = log2(bg_val))
  g2.3 <- plot_TSS_heatmap(log2(cov_mat_H3K27me3[row_order, ,3] + bg_val),
                           .col = "midnightblue",
                           .max = P0_max, 
                           .min = log2(bg_val))
  g2.4 <- plot_TSS_heatmap(log2(cov_mat_H3K27me3[row_order, ,4] + bg_val),
                           .col = "midnightblue",
                           .max = P0_max, 
                           .min = log2(bg_val))
  g2.5 <- plot_TSS_heatmap(log2(cov_mat_H3K27me3[row_order, ,5] + bg_val), 
                           .col = "midnightblue",
                           .max = P0_max, 
                           .min = log2(bg_val), .show_legend = T)
  
  ggsave(plot = grid.arrange(g2.1, g2.2, g2.3, g2.4, g2.5,
                             nrow = 1, widths = c(3, 3, 3, 3, 4.2)),
         filename = paste0("Fig1_TSS_Hetmap_H3K27me3_PC.png"), 
         path = "figs/",
         device = "png", width = 10, height = 4)
}

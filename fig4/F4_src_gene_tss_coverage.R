# ----------------------------------------------------------------------------------- #
if (TRUE) {
  # list file path
  bw_files_H2Aub <- list.files(path = "../data/bw/batch1/input_norm",
                               pattern = "H2Aub_(NT|BAP1-12).bw$", full.names = T)[2:1]
  bw_files_H3K4me3 <- list.files(path = "../data/bw/batch1/input_norm",
                                 pattern = "H3K4m3_(NT|BAP1-12).bw$", full.names = T)[2:1]
  bw_files_H3K27me3 <- list.files(path = "../data/bw/batch1/input_norm",
                                  pattern = "H3K27m3_(NT|BAP1-12).bw$", full.names = T)[2:1]
  
  bw_files_Pol2 <- list.files(path = "../data/bw/batch1/input_norm",
                              pattern = "Pol2-NTD_(NT|BAP1-12).bw$", full.names = T)[2:1]
  bw_files_Pol2S5p <- list.files(path = "../data/bw/batch1/input_norm",
                                 pattern = "Pol2-S5p_(NT|BAP1-12).bw$", full.names = T)[2:1]
  
  bw_files_Ring1b <- list.files(path = "../data/bw/batch1/input_norm",
                                pattern = "Ring1b_(NT|BAP1-12).bw$", full.names = T)[2:1]
  bw_files_Ezh2 <- list.files(path = "../data/bw/batch1/input_norm",
                              pattern = "Ezh2_(NT|BAP1-12).bw$", full.names = T)[2:1]
  
  derepressed_gene.gr <- gene.gr[gene.gr$gene_id %in% H2Aub_repressed_genes]
  
  # read in TSS coverage from bw files
  cov_mat_H2Aub <- list(get_cov_matrix(intervals = promoters(derepressed_gene.gr, 3000, 3000), 
                                       input_file = bw_files_H2Aub[1], 
                                       bin_width = 50,
                                       new_len = 100, 
                                       .df = 25),
                        get_cov_matrix(intervals = promoters(derepressed_gene.gr, 3000, 3000), 
                                       input_file = bw_files_H2Aub[2], 
                                       bin_width = 50,
                                       new_len = 100, 
                                       .df = 25) )
  
  cov_mat_H3K4me3 <- list(get_cov_matrix(intervals = promoters(derepressed_gene.gr, 3000, 3000), 
                                         input_file = bw_files_H3K4me3[1], 
                                         bin_width = 50,
                                         new_len = 100, 
                                         .df = 25),
                          get_cov_matrix(intervals = promoters(derepressed_gene.gr, 3000, 3000), 
                                         input_file = bw_files_H3K4me3[2], 
                                         bin_width = 50,
                                         new_len = 100, 
                                         .df = 25) )
  
  cov_mat_H3K27me3 <- list(get_cov_matrix(intervals = promoters(derepressed_gene.gr, 3000, 3000), 
                                          input_file = bw_files_H3K27me3[1], 
                                          bin_width = 50,
                                          new_len = 100, 
                                          .df = 25),
                           get_cov_matrix(intervals = promoters(derepressed_gene.gr, 3000, 3000), 
                                          input_file = bw_files_H3K27me3[2], 
                                          bin_width = 50,
                                          new_len = 100, 
                                          .df = 25) )
  
  cov_mat_Pol2 <- list(get_cov_matrix(intervals = promoters(derepressed_gene.gr, 3000, 3000), 
                                      input_file = bw_files_Pol2[1], 
                                      bin_width = 50,
                                      new_len = 100, 
                                      .df = 25),
                       get_cov_matrix(intervals = promoters(derepressed_gene.gr, 3000, 3000), 
                                      input_file = bw_files_Pol2[2], 
                                      bin_width = 50,
                                      new_len = 100, 
                                      .df = 25) )
  
  cov_mat_Pol2S5p <- list(get_cov_matrix(intervals = promoters(derepressed_gene.gr, 3000, 3000), 
                                         input_file = bw_files_Pol2S5p[1], 
                                         bin_width = 50,
                                         new_len = 100, 
                                         .df = 25),
                          get_cov_matrix(intervals = promoters(derepressed_gene.gr, 3000, 3000), 
                                         input_file = bw_files_Pol2S5p[2], 
                                         bin_width = 50,
                                         new_len = 100, 
                                         .df = 25) )
  
  cov_mat_Ring1b <- list(get_cov_matrix(intervals = promoters(derepressed_gene.gr, 3000, 3000), 
                                        input_file = bw_files_Ring1b[1], 
                                        bin_width = 50,
                                        new_len = 100, 
                                        .df = 25),
                         get_cov_matrix(intervals = promoters(derepressed_gene.gr, 3000, 3000), 
                                        input_file = bw_files_Ring1b[2], 
                                        bin_width = 50,
                                        new_len = 100, 
                                        .df = 25) )
  
  cov_mat_Ezh2 <- list(get_cov_matrix(intervals = promoters(derepressed_gene.gr, 3000, 3000), 
                                      input_file = bw_files_Ezh2[1], 
                                      bin_width = 50,
                                      new_len = 100, 
                                      .df = 25),
                       get_cov_matrix(intervals = promoters(derepressed_gene.gr, 3000, 3000), 
                                      input_file = bw_files_Ezh2[2], 
                                      bin_width = 50,
                                      new_len = 100, 
                                      .df = 25) )
}

# -------------------------------------------------------------------------------------- #
if (FALSE) {
  # set order and group by derepression groups
  
  derepressed_gene_tab$is_H2Aub_repressed <- derepressed_gene_tab$gene_id %in% derepressed_gene.gr$gene_id
  
  row_order <- NULL
  for (i in rev(unique(derepressed_gene_tab[, 2]))) {
    tmp_idx <- which(derepressed_gene_tab[derepressed_gene_tab$is_H2Aub_repressed, 2] == i)
    row_order %c=% tmp_idx[cov_mat_H2Aub[[1]][tmp_idx, ] %>%
                             trim_quantile() %>%
                             log1p() %>%
                             rowMeans(., na.rm = T) %>%
                             order(.)]
  }
  .breaks = derepressed_gene_tab[derepressed_gene_tab$is_H2Aub_repressed, 2] %>%
    table() %>% unname() %>% rev() %>% cumsum()
  
  
  
  plot_heatmap_condition_pair <-  function(cov_mat, row_order,
                                           .col = "red", bg_val = 1,
                                           .breaks = c(1, 10, 20),
                                           file_name = "H2Aub" ) {
    require(ggbreak)
    .max <-  max(c(log2(trim_quantile(cov_mat[[1]]) + bg_val),
                   log2(trim_quantile(cov_mat[[1]]) + bg_val)))
    .min <- 0
    
    g2.1 <- plot_TSS_heatmap(mat = log2(cov_mat[[1]][row_order, ] + bg_val), 
                             .max = .max, 
                             .min = .min, 
                             .col = .col) +
      scale_y_break(c(.breaks[1], .breaks[1])) +
      scale_y_break(c(.breaks[2], .breaks[2])) +
      scale_y_break(c(.breaks[3], .breaks[3]))
    
    ggsave(plot = g2.1,
           filename = paste0("../fig4/figs/Fig4_tmp_TSS_Heatmap_", file_name, "_1.png"), 
           device = "png", width = 1.6, height = 8)
    
    g2.2 <- plot_TSS_heatmap(mat = log2(cov_mat[[2]][row_order, ] + bg_val), 
                             .max = .max, 
                             .min = .min, 
                             .col = .col) +
      scale_y_break(c(.breaks[1], .breaks[1])) +
      scale_y_break(c(.breaks[2], .breaks[2])) +
      scale_y_break(c(.breaks[3], .breaks[3]))
    
    ggsave(plot = g2.2,
           filename = paste0("../fig4/figs/Fig4_tmp_TSS_Hetmap_", file_name, "_2.png"), 
           device = "png", width = 1.6, height = 8)
  }
  
  plot_heatmap_condition_pair(cov_mat_Ring1b, row_order = row_order, .col = sample_colors["Ring1b"],
                              bg_val = 0.001, .breaks = .breaks, file_name = "Ring1b")
  
  plot_heatmap_condition_pair(cov_mat_Ezh2, row_order = row_order, .col = sample_colors["Ezh2"],
                              bg_val = 0.001, .breaks = .breaks, file_name = "Ezh2")
  
  plot_heatmap_condition_pair(cov_mat_H2Aub, row_order = row_order, .col = sample_colors["H2Aub"],
                              bg_val = 0.001, .breaks = .breaks, file_name = "H2Aub")
  
  plot_heatmap_condition_pair(cov_mat_H3K27me3, row_order = row_order, .col = sample_colors["H3K27me3"],
                              bg_val = 0.001, .breaks = .breaks, file_name = "H3K27me3")
  
  plot_heatmap_condition_pair(cov_mat_H3K4me3, row_order = row_order, .col = sample_colors["H3K4me3"],
                              bg_val = 0.001, .breaks = .breaks, file_name = "H3K4me3")
  
  plot_heatmap_condition_pair(cov_mat_Pol2, row_order = row_order, .col = sample_colors["Pol II"],
                              bg_val = 0.001, .breaks = .breaks, file_name = "Pol II")
}


# -------------------------------------------------------------------------------------- #
if (TRUE) {

  H2Aub_repressed_genes_gr_fig4 <- sort(gene.gr[H2Aub_repressed_genes_fig4])
  
  
  cov_mat_H2Aub_1 <- get_cov_matrix(input_file = bw_files_H2Aub[1], 
                                    intervals = promoters(gene.gr[H2Aub_repressed_genes_fig4], 15000, 15000),
                                    bin_width = 50,
                                    new_len = 100, 
                                    .df = 25)
  
  cov_mat_H3K4me3_1 <- get_cov_matrix(input_file = bw_files_H3K4me3[1], 
                                      intervals = promoters(gene.gr[H2Aub_repressed_genes_fig4], 15000, 15000), 
                                         bin_width = 50,
                                         new_len = 100, 
                                         .df = 25)
                          
  cov_mat_H3K27ac_1 <- get_cov_matrix(input_file = "../data/bw/batch1/input_norm/H3K27ac_NT.bw", 
                                       intervals = promoters(gene.gr[H2Aub_repressed_genes_fig4], 15000, 15000),
                                       bin_width = 50,
                                       new_len = 100, 
                                       .df = 25)
  
  cov_mat_H3K27me3_1 <- get_cov_matrix(input_file = bw_files_H3K27me3[1], 
                                        intervals = promoters(gene.gr[H2Aub_repressed_genes_fig4], 15000, 15000),
                                       bin_width = 50,
                                       new_len = 100, 
                                       .df = 25)
  
  cov_mat_Ring1b_1 <- get_cov_matrix(input_file = bw_files_Ring1b[1], 
                                      intervals = promoters(gene.gr[H2Aub_repressed_genes_fig4], 15000, 15000),
                                     bin_width = 50,
                                     new_len = 100, 
                                     .df = 25)
  
  cov_mat_Ezh2_1 <- get_cov_matrix(input_file = bw_files_Ezh2[1], 
                                    intervals = promoters(gene.gr[H2Aub_repressed_genes_fig4], 15000, 15000),
                                   bin_width = 50,
                                   new_len = 100, 
                                   .df = 25)
  
  cov_mat_Pol2_1 <- get_cov_matrix(input_file = bw_files_Pol2[1], 
                                   intervals = promoters(gene.gr[H2Aub_repressed_genes_fig4], 15000, 15000),
                                   bin_width = 50,
                                   new_len = 100, 
                                   .df = 25)
  
  cov_mat_Pol2_2 <- get_cov_matrix(input_file = "../data/bw/batch2/input_norm/Pol2_P0_pool.bw", 
                                   intervals = promoters(gene.gr[H2Aub_repressed_genes_fig4], 15000, 15000),
                                   bin_width = 50,
                                   new_len = 100, 
                                   .df = 25)
  
  
  cov_mat_Pol2S5p_1 <- get_cov_matrix(input_file = bw_files_Pol2S5p[1], 
                                     intervals = promoters(gene.gr[H2Aub_repressed_genes_fig4], 15000, 15000),
                                     bin_width = 50,
                                     new_len = 100, 
                                     .df = 25)
  
  
  cov_mat_Cbx7_2 <- get_cov_matrix(input_file = "../data/bw/batch2/input_norm/Cbx7_P0_pool.bw", 
                                    intervals = promoters(gene.gr[H2Aub_repressed_genes_fig4], 15000, 15000),
                                   bin_width = 50,
                                   new_len = 100, 
                                   .df = 25)
  
  cov_mat_Rybp_2 <- get_cov_matrix(input_file = "../data/bw/batch2/input_norm/Rybp_P0_pool.bw", 
                                    intervals = promoters(gene.gr[H2Aub_repressed_genes_fig4], 15000, 15000),
                                   bin_width = 50,
                                   new_len = 100, 
                                   .df = 25)

}
# -------------------------------------------------------------------------------------- #

plot_TSS_heatmap <- function(mat, 
                             .title = "", .col = "black", 
                             .max = 4, .min = 0, 
                             .show_legend = F) 
{
  mat[mat < 0] <- 0
  `rownames<-`(mat, seq_len(nrow(mat))) %>% 
    trim_quantile() %>%
    reshape::melt() %>% 
    ggplot(aes(x = X2, y = X1, fill = value)) +
    geom_tile() +
    scale_x_continuous(expand = c(0,0),
                       breaks = c(1, 51, 100), 
                       labels = c('', '', '')
                       # labels = c('-3 kb', 'TSS', '+3 kb')
                       ) +
    # scale_y_cut(breaks = .breaks) +
    # scale_y_break(breaks = .breaks) +
    scale_y_continuous(expand = c(0,0)) +
    
    scale_fill_gradientn(name = "", 
                         colours = c(add.alpha("white", 0.9), .col, add.alpha("black", 0.9)),
                         limits = c(.min, .max), 
                         breaks = c(.min, (.min + .max) / 2, .max),
                         labels = round(c(0, .max / 2, .max), 1),
                         guide = guide_colorbar(label = TRUE,
                                                draw.ulim = TRUE, 
                                                draw.llim = TRUE,
                                                frame.colour = "black",
                                                ticks = TRUE, 
                                                nbin = 10,
                                                label.position = "bottom",
                                                barwidth = 6,
                                                barheight = 1, 
                                                direction = 'horizontal')) +
    
    labs(title = .title,  x = '', y = '') +
    
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, size = 1),
          axis.ticks.x = element_line(), 
          axis.text.y = element_blank(),
          plot.margin = unit(c(0, 1, 0, 1), "lines"),
          legend.position = "bottom", 
          # legend.margin = margin(unit(c(0, 0, 0, 1), "lines")),
          legend.justification = "center")
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

# wrapped function
plot_one_factor_heatmap <- function(mat, .col = "black") {
  mat %>%
    as.matrix(ncol = 1) %>%
    `rownames<-`(NULL) %>% 
    reshape::melt() %>% 
    ggplot(aes(x = X2, y = X1, fill = value)) +
    geom_tile() +
    labs(title = '',  x = '', y = '') +
    scale_fill_manual(values = c("white", .col)) +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # plot.margin = margin(0,1,0, 1, "lines"),
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.text.y = element_blank(), 
          strip.text.y = element_blank(), 
          legend.position = "bottom",
          legend.direction = "horizontal", 
          legend.key.size = unit(13, "pt"),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "white")) +
    scale_y_break(c(.breaks[1], .breaks[1])) +
    scale_y_break(c(.breaks[2], .breaks[2])) 
}

# c_min <- apply(cov_mat_H2Aub, 3, function(x) min(colMeans(x))) %>% min()
# c_max <- apply(cov_mat_H2Aub, 3, function(x) max(colMeans(x))) %>% max()
# 
# g2.1c <- plot_TSS_coverage(cov_mat_H2Aub[, ,1],
#                            .max = c_max,
#                            .min = c_min)


# -------------------------------------------------------------------------------------- #
# histogram of control RPKs
# baseMeans_de <- res.gene.list$P12[derepressed_gene_ids, "baseMean"]






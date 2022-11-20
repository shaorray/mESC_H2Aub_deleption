

# ---------------------------------------- plot coverage ---------------------------------------- #
get_fig1_tss_coverage_array <- function(target = "Ring1b",
                                        gene.gr = gene.gr,
                                        tmp_gene_ids = PcG_enriched_genes) {
  
  target_name <- ifelse(target == "Pol II", "Pol2-NTD", 
                        ifelse(target == "H3K27me3", "H3K27m3", 
                               ifelse(target == "H3K4me3", "H3K4m3", target)))
  
  target_color <- unname(sample_colors[target])
  # append TSS coverage
  # with bw files
  tmp_bw_files <- list.files("../data/bw/batch1/input_norm",
                             paste0(target_name, "_(NT|BAP1-6\\.|BAP1-12\\.)"), full.names = TRUE)
  tmp_bw_files <- rev(tmp_bw_files)
  
  tmp_1_array <- convert_coverage(file_names = tmp_bw_files, 
                                  intervals = promoters(gene.gr[tmp_gene_ids],
                                                        upstream = 2120,
                                                        downstream = 2120), 
                                  new_len = 56, bin_width = 100)
  
  return(tmp_1_array / b1_RPGC_scale)
}


plot_fig1_tss_coverage <- function(target = "Ring1b",
                                   gene.gr = gene.gr,
                                   tmp_gene_ids = PcG_enriched_genes) {
  
  tmp_1_array <- get_fig1_tss_coverage_array(target = target,
                                             gene.gr = gene.gr,
                                             tmp_gene_ids = tmp_gene_ids)
  
  tmp_1_dat <- data.frame(pos = rep(seq_len(50), 3),
                          cov = c(colMedians(log1p(tmp_1_array[, 4:53, 1])),
                                  colMedians(log1p(tmp_1_array[, 4:53, 2])),
                                  colMedians(log1p(tmp_1_array[, 4:53, 3]))),
                          Time = rep(c("P0", "P6", "P12"), each = 50))
  
  tmp_1_dat$Time <- factor(tmp_1_dat$Time, levels = c("P0", "P6", "P12"))
  
  ggsave(ggplot(tmp_1_dat, aes(x = pos, y = cov, group = Time, color = Time)) +
           geom_line(lwd = 2) +
           # geom_smooth(method = "gam") +
           scale_color_manual(values = c(colors_20[16], colors_9[c(6, 3)])) +
           scale_x_continuous(name = "",
                              breaks = c(1, 25, 50),
                              labels = c("-2 kb", "TSS", "2 kb")) +
           scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                              limits = c(0, max(tmp_1_dat$cov) * 1.1) ) +
           ylab("Median occupancy (log1p)") + ggtitle(target) +
           theme_setting +
           theme(legend.position = "bottom",
                 legend.justification = "center"),
         filename = paste0("Fig1_TSS_coverage_", target, "_Batch1.png"), 
         path = "../fig1/figs",
         device = "png", width = 3, height = 3.7)
}

for (i in c("H2Aub", "Ring1b", "H3K27me3", "Ezh2", "H3K4me3", "Pol II")) {
  plot_fig1_tss_coverage(target = i, 
                         gene.gr = gene.gr,
                         tmp_gene_ids = PcG_enriched_genes)
}

# -------------------------------------- plot heatmap ---------------------------------------- #

plot_TSS_heatmap <- function(mat, .max, .min, .order, .col = "black") {
  mat[mat >= .max] <- .max
  mat[mat <= .min] <- .min
  mat <- mat[.order, ]
  
  mat %>%
  `rownames<-`(seq_len(nrow(mat))) %>% 
    `colnames<-`(seq_len(ncol(mat))) %>% 
    # trim_quantile(0.99) %>%
    reshape::melt() %>% 
    ggplot(aes(x = X2, y = X1, fill = value)) +
    geom_tile() +
    geom_vline(xintercept = 26, lty = 2, cex = 0.1) +
    scale_x_continuous(expand = c(0,0),
                       breaks = c(1, 26, 50), 
                       labels = c('', '', '')
                       # labels = c('-3 kb', 'TSS', '+3 kb')
    ) +
    scale_y_continuous(expand = c(0,0)) +
    
    scale_fill_gradientn(name = "", 
                         colours = c(add.alpha("white", 0.7), add.alpha(.col, 0.7), .col),
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
    
    labs(title = '',  x = '', y = '') +
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, size = 1),
          axis.ticks.x = element_line(), 
          axis.text.y = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          legend.position = "none", 
          # legend.margin = margin(unit(c(0, 0, 0, 0), "lines")),
          legend.justification = "center")
}


plot_fig1_tss_heatmap <- function(target = "Ring1b",
                                  gene.gr = gene.gr,
                                  tmp_gene_ids_1 = PcG_enriched_genes,
                                  tmp_gene_ids_2 = highly_expressed_genes) {
  
  tmp_array <- get_fig1_tss_coverage_array(target = target,
                                           gene.gr = gene.gr,
                                           tmp_gene_ids = c(tmp_gene_ids_1, tmp_gene_ids_2))
  
  g_list <- list()
  
  .max <- quantile(tmp_array, .99)
  .order_1 <- order(rowMedians(tmp_array[tmp_gene_ids_1, 13:42, 1], na.rm = TRUE))
  for(i in 1:3) {
    g_list <- c(g_list, 
                list(plot_TSS_heatmap(tmp_array[tmp_gene_ids_1, 4:53, i], 
                                      .max = .max,
                                      .min = 2,
                                      .order = .order_1,
                                      .col = sample_colors[target])))
  }
  
  .order_2 <- order(rowMedians(tmp_array[tmp_gene_ids_2, 13:42, 1], na.rm = TRUE))
  for(i in 1:3) {
    g_list <- c(g_list, 
                list(plot_TSS_heatmap(tmp_array[tmp_gene_ids_2, 4:53, i], 
                                      .max = .max,
                                      .min = 2,
                                      .order = .order_2,
                                      .col = sample_colors[target])))
  }
  
  g_list
}




tmp_1_array <- get_fig1_tss_coverage_array(target = "Pol II",
                                           gene.gr = gene.gr,
                                           tmp_gene_ids = PcG_enriched_genes)
.order_1 <- order(rowMeans(tmp_1_array[, 24:32, 1], na.rm = TRUE))

highly_expressed_genes_2 <- intersect(highly_expressed_genes, names(gene.gr)) # 1123
tmp_2_array <- get_fig1_tss_coverage_array(target = "Pol II",
                                           gene.gr = gene.gr,
                                           tmp_gene_ids = highly_expressed_genes_2)
.order_2 <- order(rowMeans(tmp_2_array[, 24:32, 1], na.rm = TRUE))



for (i in c("H2Aub", "Ring1b", "H3K27me3", "Ezh2", "H3K4me3", "Pol II")) {
  g_list <- plot_fig1_tss_heatmap(target = i, 
                                  gene.gr = gene.gr, 
                                  tmp_gene_ids_1 = PcG_enriched_genes,
                                  tmp_gene_ids_2 = highly_expressed_genes_2)
  
  ggsave(filename = paste0("figs/Fig1_", i, "_TSS_heatmap.png"),
         plot = do.call(grid.arrange, args = c(g_list, nrow = 1)), 
         width = 12, height = 2)
}



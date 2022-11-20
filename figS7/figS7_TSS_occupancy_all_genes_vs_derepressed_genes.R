# --------------------------------- show occupancy scatter plots ------------------------------------- #
# 

plot_scatter_2 <- function(X, Y,
                           .offset = 10,
                           row_names,
                           target_ids = NA,
                           .xlab = "", 
                           .ylab = "", 
                           xlim = NULL, 
                           ylim = NULL) {
  
  X <- log2(X + .offset) - log2(.offset)
  Y <- log2(Y + .offset) - log2(.offset)
  
  r1 <- round(cor(X[row_names %in% target_ids], Y[row_names %in% target_ids]), 2)
  r2 <- round(cor(X, Y), 2)
  
  cbind(X, Y) %>% 
    `rownames<-`(., row_names) %>%
    `colnames<-`(., c("x", "y")) %>% 
    as.data.frame() %>%
    dplyr::filter(complete.cases(.) & is.finite(rowSums(.)) & !is.na(rowSums(.))) -> dat
  
  
  if (!is.null(xlim) & !is.null(ylim)) 
    dat <- dplyr::filter(dat, x > min(xlim) & x < max(xlim) & y > min(ylim) & y < max(ylim))
  
  dat$highlight <- rownames(dat) %in% target_ids
  
  g <- ggplot(dat, aes(x = x, y = y)) +
    geom_point(cex = 2, pch = 15, color = add.alpha("grey20", 0.1)) +
    annotate("text", x = -Inf, y = Inf,
             hjust = 0, vjust = 1.5,
             label = paste0(" n = ", nrow(dat), "\n  r = ", r2)) +
    annotate("text", x = Inf, y = Inf,
             hjust = 1.3, vjust = 1.5,
             label = paste0("  n = ", sum(dat$highlight), "\nr = ", r1),
             color = 'red') +
    scale_x_continuous(name = .xlab, limits = xlim) +
    scale_y_continuous(name = .ylab, limits = ylim) +
    scale_fill_manual(values = add.alpha(c("black", "red"), 0.5)) +
    scale_color_gradientn(colours = c("white", "grey50", "grey30", "grey15", "black")) +
    theme_setting +
    theme(legend.position = "none")
  
  g + geom_point(data = dat[dat$highlight, ], aes(x = x, y = y, color = get_dens(x, y)),
                 size = 1, pch = 19) +
    scale_color_gradientn(colours = add.alpha(rev(colors_n[1:6]), 0.2))
}


tmp <- ChIP_TSS_B1_mat[, grep("_NT", colnames(ChIP_TSS_B1_mat))]
tmp <- tmp[complete.cases(tmp), ]

# plot
g1 <- plot_scatter_2(X = tmp$H3K27m3_NT,
                     Y = tmp$Ring1b_NT, 
                     row_names = rownames(tmp),
                     xlim = c(0, 3), ylim = c(0, 4), 
                     .xlab = "H3K27me3 P0 (log2 density)", 
                     .ylab = "Ring1b P0 (log2 density)",
                     target_ids = H2Aub_repressed_genes)

g2 <- plot_scatter_2(X = tmp$H3K27m3_NT,
                     Y = tmp$Ezh2_NT, 
                     row_names = rownames(tmp),
                     xlim = c(0, 3), ylim = c(0, 4), 
                     .xlab = "H3K27me3 P0 (log2 density)", 
                     .ylab = "Ezh2 P0 (log2 density)",
                     target_ids = H2Aub_repressed_genes)

g3 <- plot_scatter_2(X = tmp$H3K27m3_NT,
                     Y = tmp$H3K4m3_NT, 
                     row_names = rownames(tmp),
                     xlim = c(0, 3), ylim = c(0, 3), 
                     .xlab = "H3K27me3 P0 (log2 density)", 
                     .ylab = "H3K4me3 P0 (log2 density)",
                     target_ids = H2Aub_repressed_genes)

g4 <- plot_scatter_2(X = tmp$H3K27m3_NT,
                     Y = tmp$`Pol2-NTD_NT`, 
                     row_names = rownames(tmp),
                     xlim = c(0, 3), ylim = c(0, 4), 
                     .xlab = "H3K27me3 P0 (log2 density)", 
                     .ylab = "Pol II P0 (log2 density)",
                     target_ids = H2Aub_repressed_genes)


g5 <- plot_scatter_2(X = tmp$Ring1b_NT,
                     Y = tmp$H2Aub_NT, 
                     row_names = rownames(tmp),
                     xlim = c(0, 4), ylim = c(0, 4), 
                     .xlab = "Ring1b P0 (log2 density)", 
                     .ylab = "H2Aub P0 (log2 density)",
                     target_ids = H2Aub_repressed_genes)

g6 <- plot_scatter_2(X = tmp$Ring1b_NT,
                     Y = tmp$H3K27ac_NT , 
                     row_names = rownames(tmp),
                     xlim = c(0, 4), ylim = c(0, 1.5), 
                     .xlab = "Ring1b P0 (log2 density)", 
                     .ylab = "H3K27ac P0 (log2 density)",
                     target_ids = H2Aub_repressed_genes)


g7 <- plot_scatter_2(X = tmp$Ring1b_NT,
                     Y = tmp$H3K4m3_NT, 
                     row_names = rownames(tmp),
                     xlim = c(0, 6), ylim = c(0, 4), 
                     .xlab = "Ring1b P0 (log2 density)", 
                     .ylab = "H3K4me3 P0 (log2 density)",
                     target_ids = H2Aub_repressed_genes)

g8 <- plot_scatter_2(X = tmp$Ring1b_NT,
                     Y = tmp$`Pol2-NTD_NT`, 
                     row_names = rownames(tmp),
                     xlim = c(0, 4), ylim = c(0, 3), 
                     .xlab = "Ring1b P0 (log2 density)", 
                     .ylab = "Pol II P0 (log2 density)",
                     target_ids = H2Aub_repressed_genes)

ggsave(plot = grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, ncol = 4),
       filename = paste0("FigS7_Scatter_H3K27me3_Ring1b_binding_All_vs_Derepressed.png"), 
       path = "../figS7/figs/",
       device = "png", width = 16, height = 8) 

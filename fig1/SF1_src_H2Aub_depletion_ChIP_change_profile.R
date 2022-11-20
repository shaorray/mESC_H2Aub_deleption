

# ----------------------------------- Supplement plot --------------------------------------- #

# reproducibility 
if (FALSE) {
  # H2Aub log2FC comparison
  dat <- data.frame(x = ChIP_TSS_B1_mat_log2FC[, "H2Aub_BAP1-12"], 
                    y = ChIP_TSS_B1_mat_log2FC[, "H2Aub_R1bKO"]) %>% 
    `rownames<-`(., rownames(ChIP_TSS_B1_mat_log2FC)) %>%
    dplyr::filter(complete.cases(.) & is.finite(rowSums(.)) & !is.na(rowSums(.)))
  
  
  g <- ggplot(dat, aes(x = x, y = y, color = get_dens(x, y))) +
    geom_point(cex = 0.5) +
    annotate("text", x = -Inf, y = Inf,
             hjust = -0.5, vjust = 1.2,
             label = paste0("All genes\n", "r = ", round(cor(dat$x, dat$y), 3), "\nn = ", nrow(dat))) +
    scale_x_continuous(name = "TSS H2Aub log2FC [BAP1 P12 batch 1]") +
    scale_y_continuous(name = "TSS H2Aub log2FC [Ring1b CKO]") +
    scale_color_gradientn(colours = c("grey80", "grey50", "grey30", "grey15", "grey10", "black")) +
    theme_setting +
    theme(legend.position = "none")
  
  dat_highlight <- dat[rownames(dat) %in% PcG_enriched_genes, ] %>% 
    dplyr::filter(complete.cases(.) &
                    !is.infinite(x) & 
                    !is.infinite(y))
  
  g + geom_point(data = dat_highlight, aes(x = x, y = y),
                 size = 0.2, pch = 19, color = add.alpha("red", 0.5)) +
    annotate("text", x = Inf, y = Inf,
             hjust = 1.2, vjust = 1.2,
             label = paste0("PcG enriched genes\n", 
                            "r = ", round(cor(dat_highlight$x, dat_highlight$y), 3), "\n",
                            "n = ", nrow(dat_highlight)),
             color = "red")
  ggsave("FigS1_scatter_plot_H2Aub_log2FC_BAP1_P12_batch1_Ring1b_KO.png", path = "../figS1/figs", 
         device = "png", width = 5, height = 5)
  
  
  # Ezh2
  dat <- data.frame(x = ChIP_TSS_B1_mat_log2FC[, "Ezh2_BAP1-12"], 
                    y = ChIP_TSS_B1_mat_log2FC[, "Ezh2_R1bKO"]) %>% 
    `rownames<-`(., rownames(ChIP_TSS_B1_mat_log2FC)) %>%
    dplyr::filter(complete.cases(.) & is.finite(rowSums(.)) & !is.na(rowSums(.)))
  
  
  g <- ggplot(dat, aes(x = x, y = y, color = get_dens(x, y))) +
    geom_point(cex = 0.5) +
    annotate("text", x = -Inf, y = Inf,
             hjust = -0.5, vjust = 1.2,
             label = paste0("All genes\n", "r = ", round(cor(dat$x, dat$y), 3), "\nn = ", nrow(dat))) +
    scale_x_continuous(name = "TSS Ezh2 log2FC [BAP1 P12 batch 1]") +
    scale_y_continuous(name = "TSS Ezh2 log2FC [Ring1b CKO]") +
    scale_color_gradientn(colours = c("grey80", "grey50", "grey30", "grey15", "grey10", "black")) +
    theme_setting +
    theme(legend.position = "none")
  
  dat_highlight <- dat[rownames(dat) %in% PcG_enriched_genes, ] %>% 
    dplyr::filter(complete.cases(.) &
                    !is.infinite(x) & 
                    !is.infinite(y))
  
  g + geom_point(data = dat_highlight, aes(x = x, y = y),
                 size = 0.2, pch = 19, color = add.alpha("red", 0.5)) +
    annotate("text", x = Inf, y = Inf,
             hjust = 1.2, vjust = 1.2,
             label = paste0("PcG enriched genes\n", 
                            "r = ", round(cor(dat_highlight$x, dat_highlight$y), 3), "\n",
                            "n = ", nrow(dat_highlight)),
             color = "red")
  ggsave("FigS1_scatter_plot_Ezh2_log2FC_BAP1_P12_batch1_Ring1b_KO.png", path = "../figS1/figs", 
         device = "png", width = 5, height = 5)

  
  # Ezh2
  dat <- data.frame(x = ChIP_TSS_B1_mat_log2FC[, "Ezh2_BAP1-12"], 
                    y = ChIP_TSS_B1_mat_log2FC[, "Ezh2_R1bKO"]) %>% 
    `rownames<-`(., rownames(ChIP_TSS_B1_mat_log2FC)) %>%
    dplyr::filter(complete.cases(.) & is.finite(rowSums(.)) & !is.na(rowSums(.)))
  
  
  g <- ggplot(dat, aes(x = x, y = y, color = get_dens(x, y))) +
    geom_point(cex = 0.5) +
    annotate("text", x = -Inf, y = Inf,
             hjust = -0.5, vjust = 1.2,
             label = paste0("All genes\n", "r = ", round(cor(dat$x, dat$y), 3), "\nn = ", nrow(dat))) +
    scale_x_continuous(name = "TSS Ezh2 log2FC [BAP1 P12 batch 1]") +
    scale_y_continuous(name = "TSS Ezh2 log2FC [Ring1b CKO]") +
    scale_color_gradientn(colours = c("grey80", "grey50", "grey30", "grey15", "grey10", "black")) +
    theme_setting +
    theme(legend.position = "none")
  
  dat_highlight <- dat[rownames(dat) %in% PcG_enriched_genes, ] %>% 
    dplyr::filter(complete.cases(.) &
                    !is.infinite(x) & 
                    !is.infinite(y))
  
  g + geom_point(data = dat_highlight, aes(x = x, y = y),
                 size = 0.2, pch = 19, color = add.alpha("red", 0.5)) +
    annotate("text", x = Inf, y = Inf,
             hjust = 1.2, vjust = 1.2,
             label = paste0("PcG enriched genes\n", 
                            "r = ", round(cor(dat_highlight$x, dat_highlight$y), 3), "\n",
                            "n = ", nrow(dat_highlight)),
             color = "red")
  ggsave("FigS1_scatter_plot_Ezh2_log2FC_BAP1_P12_batch1_Ring1b_KO.png", path = "../figS1/figs", 
         device = "png", width = 5, height = 5)
  
  # H3K27me3 log2FC, BAP1 vs Ring1b CKO
  dat <- data.frame(x = ChIP_TSS_B1_mat_log2FC[, "H3K27me3_BAP1-12"], 
                    y = ChIP_TSS_B1_mat_log2FC[, "H3K27me3_R1bKO"]) %>% 
    `rownames<-`(., rownames(ChIP_TSS_B1_mat_log2FC)) %>%
    dplyr::filter(complete.cases(.) & is.finite(rowSums(.)) & !is.na(rowSums(.)))
  
  g <- ggplot(dat, aes(x = x, y = y, color = get_dens(x, y))) +
    geom_point(cex = 0.5) +
    annotate("text", x = -Inf, y = Inf,
             hjust = -0.5, vjust = 1.2,
             label = paste0("All genes\n", "r = ", round(cor(dat$x, dat$y), 3), "\nn = ", nrow(dat))) +
    scale_x_continuous(name = "TSS H3K27me3 log2FC [BAP1 P12 batch 1]") +
    scale_y_continuous(name = "TSS H3K27me3 log2FC [Ring1b CKO]") +
    scale_color_gradientn(colours = c("grey80", "grey50", "grey30", "grey15", "grey10", "black")) +
    theme_setting +
    theme(legend.position = "none")
  
  dat_highlight <- dat[rownames(dat) %in% PcG_enriched_genes, ] %>% 
    dplyr::filter(complete.cases(.) &
                    !is.infinite(x) & 
                    !is.infinite(y))
  
  g + geom_point(data = dat_highlight, aes(x = x, y = y),
                 size = 0.2, pch = 19, color = add.alpha("red", 0.5)) +
    annotate("text", x = Inf, y = Inf,
             hjust = 1.2, vjust = 1.2,
             label = paste0("PcG enriched genes\n", 
                            "r = ", round(cor(dat_highlight$x, dat_highlight$y), 3), "\n",
                            "n = ", nrow(dat_highlight)),
             color = "red")
  
  ggsave("FigS1_scatter_plot_H3K27me3_log2FC_BAP1_batch1_Ring1b_CKO.png", path = "../figS1/figs", 
         device = "png", width = 5, height = 5)
  
  # Pol II log2FC, BAP1 vs Ring1b CKO
  dat <- data.frame(x = ChIP_TSS_B1_mat_log2FC[, "Pol II-NTD_BAP1-12"], 
                    y = ChIP_TSS_B1_mat_log2FC[, "Pol II-NTD_R1bKO"]) %>% 
    `rownames<-`(., rownames(ChIP_TSS_B1_mat_log2FC)) %>%
    dplyr::filter(complete.cases(.) & is.finite(rowSums(.)) & !is.na(rowSums(.)))
  
  g <- ggplot(dat, aes(x = x, y = y, color = get_dens(x, y))) +
    geom_point(cex = 0.5) +
    annotate("text", x = -Inf, y = Inf,
             hjust = -0.5, vjust = 1.2,
             label = paste0("All genes\n", "r = ", round(cor(dat$x, dat$y), 3), "\nn = ", nrow(dat))) +
    scale_x_continuous(name = "TSS Pol II log2FC [BAP1 P12 batch 1]") +
    scale_y_continuous(name = "TSS Pol II log2FC [Ring1b CKO]") +
    scale_color_gradientn(colours = c("grey80", "grey50", "grey30", "grey15", "grey10", "black")) +
    theme_setting +
    theme(legend.position = "none")
  
  dat_highlight <- dat[rownames(dat) %in% PcG_enriched_genes, ] %>% 
    dplyr::filter(complete.cases(.) &
                    !is.infinite(x) & 
                    !is.infinite(y))
  
  g + geom_point(data = dat_highlight, aes(x = x, y = y),
                 size = 0.2, pch = 19, color = add.alpha("red", 0.5)) +
    annotate("text", x = Inf, y = Inf,
             hjust = 1.2, vjust = 1.2,
             label = paste0("PcG enriched genes\n", 
                            "r = ", round(cor(dat_highlight$x, dat_highlight$y), 3), "\n",
                            "n = ", nrow(dat_highlight)),
             color = "red")
  
  ggsave("FigS1_scatter_plot_Pol II_log2FC_BAP1_batch1_Ring1b_CKO.png", path = "../figS1/figs", 
         device = "png", width = 5, height = 5)
  
  # Ring1b log2FC, BAP1 vs Ring1b CKO
  dat <- data.frame(x = ChIP_TSS_B1_mat_log2FC[, "Ring1b_BAP1-12"], 
                    y = ChIP_TSS_B1_mat_log2FC[, "Ring1b_R1bKO"]) %>% 
    `rownames<-`(., rownames(ChIP_TSS_B1_mat_log2FC)) %>%
    dplyr::filter(complete.cases(.) & is.finite(rowSums(.)) & !is.na(rowSums(.)))
  
  g <- ggplot(dat, aes(x = x, y = y, color = get_dens(x, y))) +
    geom_point(cex = 0.5) +
    annotate("text", x = -Inf, y = Inf,
             hjust = -0.5, vjust = 1.2,
             label = paste0("All genes\n", "r = ", round(cor(dat$x, dat$y), 3), "\nn = ", nrow(dat))) +
    scale_x_continuous(name = "TSS Ring1b log2FC [BAP1 P12 batch 1]") +
    scale_y_continuous(name = "TSS Ring1b log2FC [Ring1b CKO]") +
    scale_color_gradientn(colours = c("grey80", "grey50", "grey30", "grey15", "grey10", "black")) +
    theme_setting +
    theme(legend.position = "none")
  
  dat_highlight <- dat[rownames(dat) %in% PcG_enriched_genes, ] %>% 
    dplyr::filter(complete.cases(.) &
                    !is.infinite(x) & 
                    !is.infinite(y))
  
  g + geom_point(data = dat_highlight, aes(x = x, y = y),
                 size = 0.2, pch = 19, color = add.alpha("red", 0.5)) +
    annotate("text", x = Inf, y = Inf,
             hjust = 1.2, vjust = 1.2,
             label = paste0("PcG enriched genes\n", 
                            "r = ", round(cor(dat_highlight$x, dat_highlight$y), 3), "\n",
                            "n = ", nrow(dat_highlight)),
             color = "red")
  
  ggsave("FigS1_scatter_plot_Ring1b_log2FC_BAP1_batch1_Ring1b_CKO.png", path = "../figS1/figs", 
         device = "png", width = 5, height = 5)
  
  # correlation barplot
  dat <- NULL
  for (i in c("H2Aub", "H3K27me3", "Ezh2", "Ring1b", "H3K4me3", "H3K27ac", "Pol II-NTD")) {
    dat <- rbind(dat,
                 c(i,
                   single_variance_explained(ChIP_TSS_B1_mat_log2FC[, paste0(i, "_BAP1-12")],
                                             ChIP_TSS_B1_mat_log2FC[, paste0(i, "_R1bKO")], is.cor = T),
                   single_variance_explained(ChIP_TSS_B1_mat_log2FC[PcG_enriched_genes, paste0(i, "_BAP1-12")],
                                             ChIP_TSS_B1_mat_log2FC[PcG_enriched_genes, paste0(i, "_R1bKO")], is.cor = T)
                 )
    )
    
  }
  dat <- as.data.frame(dat)
  colnames(dat) <- c("ChIP", "All_gene", "PcG_gene")
  dat <- data.frame(ChIP = rep(dat$ChIP, 2),
                    Corr = c(as.character(dat$All_gene), as.character(dat$PcG_gene)),
                    Group = c(rep("All genes", 7), rep("PcG enriched genes", 7)))
  dat$ChIP <- factor(dat$ChIP, c("H2Aub", "H3K27me3", "Ezh2", "Ring1b", "H3K4me3", "H3K27ac", "Pol II-NTD"))
  dat$Corr <- as.numeric(as.character(dat$Corr))
  
  ggplot(dat, aes(x = ChIP, y = Corr, fill = Group)) +
    geom_hline(yintercept = c(0.1, 0.3), lty = 2) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = c("grey50", "red2")) +
    xlab("") + ylab("Log2FC Pearson's correlation") + ggtitle("BAP1 P12 ~ Ring1b KO") +
    theme_setting +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("FigS1_barplot_ChIP_log2FC_BAP1_P12_Ring1b_KO.png", path = "../figS1/figs", 
         device = "png", width = 6, height = 5)
  
  # H2Aub change reproducibility
  
  # Batch 1 vs Batch at P0
  dat <- data.frame(x = ChIP_TSS_B1_mat[, "H2Aub_NT"], 
                    y = ChIP_TSS_B2_mat[rownames(ChIP_TSS_B1_mat), "H2Aub_P0"]) %>% 
    "*"(10) %>% log1p() %>% trim_quantile(0.9995) %>% 
    `colnames<-`(., c("x", "y")) %>%
    `rownames<-`(., rownames(ChIP_TSS_B1_mat)) %>% 
    as.data.frame() %>%
    dplyr::filter(complete.cases(.) & is.finite(rowSums(.)) & !is.na(rowSums(.)))
  
  g <- ggplot(dat, aes(x = x, y = y, color = get_dens(x, y))) +
    geom_point(cex = 0.5) +
    annotate("text", x = -Inf, y = Inf,
             hjust = -0.5, vjust = 1.2,
             label = paste0("All genes\n", "r = ", round(cor(dat$x, dat$y), 3), "\nn = ", nrow(dat))) +
    scale_x_continuous(name = "TSS H2Aub untreated [batch 1]") +
    scale_y_continuous(name = "TSS H2Aub untreated [batch 2]") +
    scale_color_gradientn(colours = c("grey80", "grey50", "grey30", "grey15", "grey10", "black")) +
    theme_setting +
    theme(legend.position = "none")
  
  dat_highlight <- dat[rownames(dat) %in% PcG_enriched_genes, ] %>% 
    dplyr::filter(complete.cases(.) &
                    !is.infinite(x) & 
                    !is.infinite(y))
  
  g + geom_point(data = dat_highlight, aes(x = x, y = y),
                 size = 0.2, pch = 19, color = add.alpha("red", 0.5)) +
    annotate("text", x = Inf, y = Inf,
             hjust = 1.2, vjust = 1.2,
             label = paste0("PcG enriched genes\n", 
                            "r = ", round(cor(dat_highlight$x, dat_highlight$y), 3), "\n",
                            "n = ", nrow(dat_highlight)),
             color = "red")
  ggsave("FigS1_scatter_plot_H2Aub_batch1_vs_batch2_NT.png", path = "../figS1/figs", 
         device = "png", width = 5, height = 5)
  
}
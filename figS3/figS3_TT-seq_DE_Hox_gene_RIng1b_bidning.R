# --------------------------------------------------------------------------------------------------- #
library(ggrepel)

Hox_genes <- gene.gr[grep("Hox", gene.gr$gene_name)]
Hox_gene_ids <- intersect.Vector(Hox_genes$gene_id,
                                 rownames(ChIP_TSS_B1_mat_log2FC))

dat_p <- cbind(TTseq_res.list$res_E0B12[Hox_gene_ids, "log2FoldChange"],
               ChIP_TSS_B1_mat_log2FC[Hox_gene_ids, "Ring1b_BAP1-12"]) %>%
  as.data.frame() %>% `colnames<-`(c("RNA_log2FC", "ChIP_log2FC")) 
dat_p$name <- Hox_genes$gene_name[match(Hox_gene_ids, Hox_genes$gene_id)]


H2Aub_repressed_genes_tmp <- intersect.Vector(intersect.Vector(H2Aub_repressed_genes,
                                                               rownames(TTseq_res.list$res_E0B12)),
                                              rownames(ChIP_TSS_B1_mat_log2FC))



g_list <- list()
for (i in c("Ring1b", "Ezh2", "H3K4me3", "H2Aub", "H3K27me3", "Pol II-NTD")) {
  dat <- cbind(TTseq_res.list$res_E0B12[H2Aub_repressed_genes_tmp, "log2FoldChange"],
               ChIP_TSS_B1_mat_log2FC[H2Aub_repressed_genes_tmp, paste0(i, "_BAP1-12")]) %>%
    as.data.frame() %>% `colnames<-`(c("RNA_log2FC", "ChIP_log2FC")) %>% 
    dplyr::filter(complete.cases(.) & is.finite(rowSums(.)))
  
  dat_p$ChIP_log2FC <- ChIP_TSS_B1_mat_log2FC[Hox_gene_ids, paste0(i, "_BAP1-12")]
  idx <- abs(dat_p$ChIP_log2FC) < 0.5 & dat_p$RNA_log2FC > 1.5
  
  dat_p_i2 <- dat_p[!idx, ]
  
  g <- ggplot(dat, aes(x = RNA_log2FC, y = ChIP_log2FC,
                       color = get_dens(X = RNA_log2FC, Y = ChIP_log2FC)),
              size = 0.5,
              alpha = 0.2) +
    geom_vline(xintercept = 0, lty = 1, alpha = 0.5) +
    geom_hline(yintercept = 0, lty = 1, alpha = 0.5) +
    geom_point(size = 1) + 
    scale_color_gradientn(colours = rev(colors_n)) +
    xlab("TT-seq log2FC [P12]") + ylab("ChIP log2FC [P12]") + ggtitle(i) +
    annotate(geom = "text", x = -Inf, y = Inf, hjust = 0, vjust = 1,
             label = paste("r =", round(cor(dat$RNA_log2FC, dat$ChIP_log2FC), 2),
                           "\nn =", nrow(dat))) +
    ggpubr::theme_pubclean() +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5))
  if (i != "H2Aub") {
    dat_p_i <- dat_p[idx, ]
    g <- g + geom_point(data = dat_p_i, 
                        mapping = aes(x = RNA_log2FC, y = ChIP_log2FC), 
                        color = "red2", size = 2, pch = 15) +
      geom_text_repel(data = dat_p_i, 
                      aes(x = RNA_log2FC, y = ChIP_log2FC,
                          label = name,
                          segment.shape = 0, 
                          color = 0.15),
                      parse = TRUE,
                      force = 10,
                      point.padding = 1e-06,
                      nudge_x = 6,
                      nudge_y = 1,
                      direction = "y",
                      hjust = 1, 
                      vjust = 0.5, 
                      box.padding = 0.5, 
                      segment.size = 0.2,
                      segment.ncp = 1,
                      segment.square = T)
  }
  
  g <- g + geom_point(data = dat_p_i2, 
                      mapping = aes(x = RNA_log2FC, y = ChIP_log2FC), 
                      color = "red2", size = 2, pch = 0)
  
  g_list <- c(g_list, list(g))
}

ggsave(plot = do.call(what = grid.arrange, c(g_list, ncol = 3)),
       filename = paste0("FigS3_scatter_RNA_ChIP_LFC_correlation2.png"), 
       path = "../figS3/figs/",
       device = "png", width = 12, height = 8)

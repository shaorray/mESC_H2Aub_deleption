
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

TT_seq_txTPM_use_genes <- TT_seq_txTPM_norm_cmb[intersect.Vector(use_gene_ids,
                                                                 rownames(TT_seq_txTPM_norm_cmb)),
                                                c(1,3,2,4, 7,8, 5)]


Pol_gene_den <- .countBW(bw_files = list.files("/mnt/E0767589767560E8/UPPMAX/Ezh2i_BAP1_ChIP",
                                               "Pol.*bw$", full.names = T), 
                         intervals = gene.gr[use_gene_ids], fast = F)

Pol_gene_den <- Pol_gene_den[, !grepl("Trp", colnames(Pol_gene_den))]

for (i in c("NTD", "S2p", "S5p")) {
  Pol_gene_den_tmp <- 
    Pol_gene_den[rownames(TT_seq_txTPM_use_genes),
                 grep(i, colnames(Pol_gene_den))][, c(8, 3, 1, 5, 6, 7, 2)]
  
  speed_gene_tmp <- log10(TT_seq_txTPM_use_genes / Pol_gene_den_tmp)
  speed_gene_tmp <- speed_gene_tmp[is.finite(rowSums(speed_gene_tmp)) & 
                                     !is.na(rowSums(speed_gene_tmp)), ]
  colnames(speed_gene_tmp) <- c("P0", "P6", "P12", "Ezh2i_1", "Ezh2i_2", "Ezh2i_7", "Ezh2i_1_P12")
  
  speed_gene_tmp$DE_group <- log2FC_cls[rownames(speed_gene_tmp)]
  
  speed_table_tmp <- reshape::melt(speed_gene_tmp, id.vars = "DE_group")
  
  ggplot(speed_table_tmp, aes(x = as.factor(DE_group), y = value, fill = as.factor(variable))) +
    geom_boxplot(notch = T, outlier.size = 0) +
    xlab("") + ylab("log Est velocity") + ggtitle(i) +
    ylim(quantile(speed_table_tmp$value, c(0.002, 0.998))) +
    scale_fill_viridis_d(begin = 0.1, end = 0.9, option = "C") +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          strip.text.y = element_blank(), 
          legend.position = "none")
  ggsave(filename = paste0("Fig4_est_velocity_", i, "_DE_cls.png"), 
         path = "figs/",
         device = "png", width = 4.5, height = 5)
}

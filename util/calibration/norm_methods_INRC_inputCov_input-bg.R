# 2022 Jun 16
# compare INRC, input coverage, input-background norm

library(ggplot2)

PcG_enriched_genes <- readRDS("../PcG_enriched_genes.rds")

bw_files_2 <- list.files(path = "../bw_batch2_input_norm_merge",
                         pattern = "*bw$", full.names = T)
ChIP_TSS_B2_mat_input_only <- .countBW(bw_files = bw_files_2,
                             intervals = gene.gr,
                             blacklist = blacklist.gr, fast = F)
colnames(ChIP_TSS_B2_mat_input_only) <- colnames(ChIP_TSS_B2_mat)
rownames(ChIP_TSS_B2_mat_input_only) <- rownames(ChIP_TSS_B2_mat)

colnames(ChIP_TSS_B2_BAM_mat_norm) <- gsub("Pol II", "RNAPII-NTD", colnames(ChIP_TSS_B2_BAM_mat_norm))


plot_boxplot <- function(mat, gene_ids = NULL, gene_set = "All genes", 
                         target = "Ring1b", norm_method = "INRC", .ylim = c(0, 4)) {
  if (!is.null(gene_ids)) mat <- mat[intersect(gene_ids, rownames(mat)), ]
  mat %>%
    select(grep(target, colnames(mat))) %>%
    reshape::melt() %>%
    dplyr::mutate(condition = gsub(paste0(target, "_"), "", variable)) %>%
    ggplot(aes(x = condition, y = log1p(value), color = condition)) +
    geom_boxplot(outlier.colour = NA) + ylim(.ylim) +
    xlab("") + ylab(paste(target, "log1p TSS occupancy")) + 
    ggtitle(paste0(gene_set, " (n = ", nrow(mat), ") ", norm_method)) +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_boxplot_helper <- function(gene_ids = NULL, gene_set = "All genes", target, .ylim = c(0, 2)) {
  g1 <- plot_boxplot(ChIP_TSS_B2_BAM_mat_norm, gene_ids = gene_ids, gene_set = gene_set,
                     target = target, norm_method = "INRC", .ylim = .ylim)
  g2 <- plot_boxplot(log1p(ChIP_TSS_B2_mat_input_only), gene_ids = gene_ids, gene_set = gene_set,
                     target = target, norm_method = "Input bin", .ylim = .ylim)
  g3 <- plot_boxplot(log1p(ChIP_TSS_B2_mat), gene_ids = gene_ids, gene_set = gene_set,
                     target = target, norm_method = "Input-bg", .ylim = .ylim)
  grid.arrange(g1, g2, g3, nrow = 1)
}

ggsave("bg_norm/figs/Ring1b.png",
       grid.arrange(
         plot_boxplot_helper(gene_ids = NULL, gene_set = "All genes", target = "Ring1b", .ylim = c(0, 1.5)),
         plot_boxplot_helper(gene_ids = PcG_enriched_genes, gene_set = "PcG genes", target = "Ring1b", .ylim = c(0, 2)),
         nrow = 2), width = 10, height = 6
)

ggsave("bg_norm/figs/Rybp.png",
       grid.arrange(
         plot_boxplot_helper(gene_ids = NULL, gene_set = "All genes", target = "Rybp", .ylim = c(0, 1.5)),
         plot_boxplot_helper(gene_ids = PcG_enriched_genes, gene_set = "PcG genes", target = "Rybp", .ylim = c(0, 2)),
         nrow = 2), width = 10, height = 6
)

ggsave("bg_norm/figs/Cbx7.png",
       grid.arrange(
         plot_boxplot_helper(gene_ids = NULL, gene_set = "All genes", target = "Cbx7", .ylim = c(0, 1.5)),
         plot_boxplot_helper(gene_ids = PcG_enriched_genes, gene_set = "PcG genes", target = "Cbx7", .ylim = c(0, 2)),
         nrow = 2), width = 10, height = 6
)

ggsave("bg_norm/figs/H3K27me3.png",
       grid.arrange(
         plot_boxplot_helper(gene_ids = NULL, gene_set = "All genes", target = "H3K27me3", .ylim = c(0, 1.5)),
         plot_boxplot_helper(gene_ids = PcG_enriched_genes, gene_set = "PcG genes", target = "H3K27me3", .ylim = c(0, 2)),
         nrow = 2), width = 10, height = 6
)

ggsave("bg_norm/figs/H2Aub.png",
       grid.arrange(
         plot_boxplot_helper(gene_ids = NULL, gene_set = "All genes", target = "H2Aub", .ylim = c(0, 1)),
         plot_boxplot_helper(gene_ids = PcG_enriched_genes, gene_set = "PcG genes", target = "H2Aub", .ylim = c(0, 1.5)),
         nrow = 2), width = 10, height = 6
)

ggsave("bg_norm/figs/RNAPII-NTD.png",
       grid.arrange(
         plot_boxplot_helper(gene_ids = NULL, gene_set = "All genes", target = "RNAPII-NTD", .ylim = c(0, 1.5)),
         plot_boxplot_helper(gene_ids = PcG_enriched_genes, gene_set = "PcG genes", target = "RNAPII-NTD", .ylim = c(0, 1.5)),
         nrow = 2), width = 10, height = 6
)

system("convert bg_norm/figs/*png -append bg_norm/figs/norm_method_comparison_long.png")

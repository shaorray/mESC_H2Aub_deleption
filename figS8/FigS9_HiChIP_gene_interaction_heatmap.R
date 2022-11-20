
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")

gene.mm10.gr <- importRanges("/mnt/0E471D453D8EE463/genomeDir/Mus_musculus.GRCm38.79.gtf")
gene.mm10.gr <- gene.mm10.gr[gene.mm10.gr$type == "gene" & gene.mm10.gr$gene_biotype %in% c("protein_coding", "lincRNA")]


E14_H3K4me3_pairs <- readRDS("data/E14_H3K4me3_pairs.rds")
RING_H3K4me3_pairs <- readRDS("data/RING_H3K4me3_pairs.rds")
EED_H3K4me3_pairs <- readRDS("data/EED_H3K4me3_pairs.rds")



# type 1: by locus

# PcG enriched genes
E14_PcG_mat = foreach(gene_id = PcG_enriched_genes, .combine = "+") %dopar% {
  get_gene_interaction_matrix(gene_id, E14_H3K4me3_pairs)
}

RING_PcG_mat = foreach(gene_id = PcG_enriched_genes, .combine = "+") %dopar% {
  get_gene_interaction_matrix(gene_id, RING_H3K4me3_pairs)
}

EED_PcG_mat = foreach(gene_id = PcG_enriched_genes, .combine = "+") %dopar% {
  get_gene_interaction_matrix(gene_id, EED_H3K4me3_pairs)
}

g1 <- plot_gene_interaction_matrix(E14_PcG_mat) + ggtitle("WT")
g2 <- plot_gene_interaction_matrix(RING_PcG_mat) + ggtitle("PRC1 KO")
g3 <- plot_gene_interaction_matrix(EED_PcG_mat) + ggtitle("PRC2 KO")

ggsave(plot = grid.arrange(g1, g2, g3, nrow = 1),
       filename = "figs/FigS8_HIChIP_PcG_enriched_gene_interaction_heatmap.png", width = 15, height = 3)


# highly active genes

foreach(gene_id = intersect.Vector(highly_expressed_genes, gene.mm10.gr$gene_id),
                            .combine = "+") %dopar% {
  get_gene_interaction_matrix(gene_id, E14_H3K4me3_pairs)
} %>% plot_gene_interaction_matrix() + ggtitle("WT") -> g1

foreach(gene_id = intersect.Vector(highly_expressed_genes, gene.mm10.gr$gene_id),
                             .combine = "+") %dopar% {
  get_gene_interaction_matrix(gene_id, RING_H3K4me3_pairs)
} %>% plot_gene_interaction_matrix() + ggtitle("PRC1 KO") -> g2

foreach(gene_id = intersect.Vector(highly_expressed_genes, gene.mm10.gr$gene_id),
                            .combine = "+") %dopar% {
  get_gene_interaction_matrix(gene_id, EED_H3K4me3_pairs)
} %>% plot_gene_interaction_matrix() + ggtitle("PRC2 KO") -> g3

ggsave(plot = grid.arrange(g1, g2, g3, nrow = 1),
       filename = "figs/FigS8_HIChIP_highly_expressed_gene_interaction_heatmap.png", width = 15, height = 3)


# bivalent genes

foreach(gene_id = intersect.Vector(Bivalent_genes_chromHMM, gene.mm10.gr$gene_id),
        .combine = "+") %dopar% {
          get_gene_interaction_matrix(gene_id, E14_H3K4me3_pairs)
        } %>% plot_gene_interaction_matrix() + ggtitle("WT") -> g1

foreach(gene_id = intersect.Vector(Bivalent_genes_chromHMM, gene.mm10.gr$gene_id),
        .combine = "+") %dopar% {
          get_gene_interaction_matrix(gene_id, RING_H3K4me3_pairs)
        } %>% plot_gene_interaction_matrix() + ggtitle("PRC1 KO") -> g2

foreach(gene_id = intersect.Vector(Bivalent_genes_chromHMM, gene.mm10.gr$gene_id),
        .combine = "+") %dopar% {
          get_gene_interaction_matrix(gene_id, EED_H3K4me3_pairs)
        } %>% plot_gene_interaction_matrix() + ggtitle("PRC2 KO") -> g3

ggsave(plot = grid.arrange(g1, g2, g3, nrow = 1),
       filename = "figs/FigS8_HIChIP_bivalent_gene_interaction_heatmap.png", width = 15, height = 3)


# H2Aub repressed genes

foreach(gene_id = intersect.Vector(H2Aub_repressed_genes, gene.mm10.gr$gene_id),
        .combine = "+") %dopar% {
          get_gene_interaction_matrix(gene_id, E14_H3K4me3_pairs)
        } %>% plot_gene_interaction_matrix() + ggtitle("WT") -> g1

foreach(gene_id = intersect.Vector(H2Aub_repressed_genes, gene.mm10.gr$gene_id),
        .combine = "+") %dopar% {
          get_gene_interaction_matrix(gene_id, RING_H3K4me3_pairs)
        } %>% plot_gene_interaction_matrix() + ggtitle("PRC1 KO") -> g2

foreach(gene_id = intersect.Vector(H2Aub_repressed_genes, gene.mm10.gr$gene_id),
        .combine = "+") %dopar% {
          get_gene_interaction_matrix(gene_id, EED_H3K4me3_pairs)
        } %>% plot_gene_interaction_matrix() + ggtitle("PRC2 KO") -> g3

ggsave(plot = grid.arrange(g1, g2, g3, nrow = 1),
       filename = "figs/FigS8_HIChIP_H2Aub_repressed_gene_interaction_heatmap.png", width = 15, height = 3)


# functions -------------------------------------------------------------------------------------

plot_gene_interaction_matrix <- function(mat, 
                                         flank_bin_num = 25,
                                         body_bin_num = 80,
                                         dist_bin_num = 50) {
  
  ggplot(reshape::melt(mat / sum(mat)), 
         aes(x = X2, y = X1, fill = value)) +
    geom_tile() +
    geom_hline(yintercept = dist_bin_num * log10(200) / 8, lty = 2, alpha = 0.5) +
    annotate(geom = "text", x = flank_bin_num + body_bin_num, y = dist_bin_num * log10(200) / 8, 
             hjust = 0, vjust = 1.2, label = "200 bp") +
    scale_fill_viridis(name = "Frequency", option = "B", direction = -1) +
    scale_x_continuous(name = "", 
                       breaks = c(1, flank_bin_num, flank_bin_num + body_bin_num, flank_bin_num * 2 + body_bin_num), 
                       labels = c("-5 kb", "TSS", "TES", "5 kb")) +
    scale_y_continuous(name = "Interaction distance (bp) \n", 
                       breaks = seq(0, dist_bin_num, length.out = 5), 
                       labels = formatC(10^(seq(0, dist_bin_num, length.out = 5) / dist_bin_num * 8), format = "e", digits = 0) ) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.ticks.x = element_line(),
          legend.position = "none")
  
}


get_gene_interaction_matrix <- function(gene_id, 
                                        gene_pair_table,
                                        flank_bin_num = 25,
                                        body_bin_num = 80,
                                        dist_bin_num = 50) {
  
  gene_strand = as.character(strand(gene.mm10.gr[gene.mm10.gr$gene_id == gene_id]))
  gene_start = start(flank(gene.mm10.gr[gene.mm10.gr$gene_id == gene_id], width = 1, start = TRUE))
  gene_length = width(gene.mm10.gr[gene.mm10.gr$gene_id == gene_id])
  
  
  idx_left = which(gene_pair_table$left_gene_id == gene_id)  
  idx_right = which(gene_pair_table$right_gene_id == gene_id)
  
  if (length(c(idx_left, idx_right)) == 0) 
    return(matrix(0, nrow = dist_bin_num, ncol = body_bin_num + 2 * flank_bin_num))
  
  tmp_pairs = rbind(gene_pair_table[idx_left, ],
                    `colnames<-`(gene_pair_table[idx_right, c(4:6, 1:3, 8, 7)],
                                 colnames(gene_pair_table)))
  
  tmp_pairs$dist = abs(tmp_pairs$V2 - tmp_pairs$V5)
  
  # keep same chromosome
  tmp_pairs = tmp_pairs[with(tmp_pairs, as.character(V1) == as.character(V4) & dist < 10^8), ]
  
  
  if (gene_strand == "+") {
    tmp_pairs[, c(2, 5)] = tmp_pairs[, c(2, 5)] - gene_start
  } else {
    tmp_pairs[, c(2, 5)] = gene_start - tmp_pairs[, c(2, 5)]
  }
  
  # slice gene interval (upstream, gene body, downstream)
  # upstream
  up_pairs = tmp_pairs[tmp_pairs$V2 < 0, ]
  up_cuts = cut(up_pairs$V2, breaks = seq(-5000, 0, length.out = flank_bin_num + 1))
  up_mat = matrix(0, nrow = dist_bin_num, ncol = flank_bin_num)
  for (i in seq_len(flank_bin_num)) {
    bin_i = which(up_cuts == levels(up_cuts)[i])
    if (length(bin_i) > 0)
      up_mat[, i] = 
        unname(
          table(
            cut(up_pairs[bin_i, "dist"],
                breaks = 10^seq(0, 8, length.out = dist_bin_num + 1))
          )
        ) 
  }
  # gene body
  mid_pairs = tmp_pairs[tmp_pairs$V2 > 0 & tmp_pairs$V2 < gene_length, ]
  mid_cuts = cut(mid_pairs$V2, breaks = seq(0, gene_length, length.out = body_bin_num + 1))
  mid_mat = matrix(0, nrow = dist_bin_num, ncol = body_bin_num)
  for (i in seq_len(body_bin_num)) {
    bin_i = which(mid_cuts == levels(mid_cuts)[i])
    if (length(bin_i) > 0)
      mid_mat[, i] = 
        unname(
          table(
            cut(mid_pairs[bin_i, "dist"],
                breaks = 10^seq(0, 8, length.out = dist_bin_num + 1)
            )
          )
        ) 
  }
  mid_mat = mid_mat / ((gene_length / body_bin_num) / (5000 / flank_bin_num))
  # downstream
  down_pairs = tmp_pairs[tmp_pairs$V2 > gene_length & tmp_pairs$V2 < gene_length + 5000, ]
  down_cuts = cut(down_pairs$V2, breaks = seq(gene_length, gene_length + 5000, length.out = dist_bin_num + 1))
  down_mat = matrix(0, nrow = dist_bin_num, ncol = flank_bin_num)
  for (i in seq_len(flank_bin_num)) {
    bin_i = which(down_cuts == levels(down_cuts)[i])
    if (length(bin_i) > 0)
      down_mat[, i] = 
        unname(
          table(
            cut(down_pairs[bin_i, "dist"],
                breaks = 10^seq(0, 8, length.out = dist_bin_num + 1)
            )
          )
        ) 
  }
  
  out_mat = cbind(up_mat, mid_mat, down_mat)
  if (any(rowSums(out_mat) > 1e3)) {
    message(paste("Discard gene:", gene_id, ", due to large counts"))
    return(0)
  } else {
    return(out_mat)
  }
}




# type 2: triangle




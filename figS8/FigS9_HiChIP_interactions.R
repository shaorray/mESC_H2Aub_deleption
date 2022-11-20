setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")

# load processed TSS pairs -------------------------------------------------
if (TRUE) {
  E14_H3K4me3_pairs <- readRDS("data/E14_H3K4me3_pairs.rds")
  RING_H3K4me3_pairs <- readRDS("data/RING_H3K4me3_pairs.rds")
  EED_H3K4me3_pairs <- readRDS("data/EED_H3K4me3_pairs.rds")
  E14_H3K27me3_pairs <- readRDS("data/E14_H3K27me3_pairs.rds")
}

# plot interaction -----------------------------------------------------------------

H2Aub_repressed_genes_tbl <- cbind("gene_id" = H2Aub_repressed_genes, # from RNA-seq de-repressed genes
                                   "gene_name" = gene.mm10.gr$gene_name[match(H2Aub_repressed_genes, gene.mm10.gr$gene_id)])

annotation_list <- list("ATAC" = ATAC.mm10.gr,
                        "Bivalent state" = ChromHMM_mm10[ChromHMM_mm10$name == "BivalentChromatin"],
                        "TU" = intergenic_cmb.mm10[intergenic_cmb.mm10$is.TU])

save_plot_wrapper(H2Aub_repressed_genes_tbl[grep("Hox", H2Aub_repressed_genes_tbl[, 2]), 2])

# summarize interaction changes ----------------------------------------------------

E14_H3K4me3_pair_count <- table(na.omit(c(E14_H3K4me3_pairs$left_gene_id, E14_H3K4me3_pairs$right_gene_id)))
RING_H3K4me3_pair_count <- table(na.omit(c(RING_H3K4me3_pairs$left_gene_id, RING_H3K4me3_pairs$right_gene_id)))
EED_H3K4me3_pair_count <- table(na.omit(c(EED_H3K4me3_pairs$left_gene_id, EED_H3K4me3_pairs$right_gene_id)))

pair_count_tab <- data.frame(gene_id = unique(c(names(E14_H3K4me3_pair_count), 
                                                names(RING_H3K4me3_pair_count),
                                                names(EED_H3K4me3_pair_count))))

pair_count_tab$E14_pair_count <- as.numeric(E14_H3K4me3_pair_count[match(pair_count_tab$gene_id,
                                                              names(E14_H3K4me3_pair_count) )])
pair_count_tab$RING_pair_count <- as.numeric(RING_H3K4me3_pair_count[match(pair_count_tab$gene_id,
                                                              names(RING_H3K4me3_pair_count) )])
pair_count_tab$EED_pair_count <- as.numeric(EED_H3K4me3_pair_count[match(pair_count_tab$gene_id,
                                                              names(EED_H3K4me3_pair_count) )])

smoothScatter(log1p(pair_count_tab$E14_pair_count), log1p(pair_count_tab$RING_pair_count) )
with(pair_count_tab[pair_count_tab$gene_id %in% H2Aub_repressed_genes, ], 
     points(log1p(E14_pair_count), log1p(RING_pair_count), col = 2, pch = 19, cex = 0.1))

# functions ------------------------------------------------------------------------
save_plot_wrapper <- function(gene_names) {
  sapply(gene_names, function(gene_name) {
    plot_interaction_profile(gene_name = gene_name, 
                             all_gene_gr = gene.mm10.gr, 
                             pairs_table_list = list(E14_H3K4me3_pairs,
                                                     RING_H3K4me3_pairs,
                                                     EED_H3K4me3_pairs),
                             pairs_table_titles = c("WT H3K4me3 HiChIP", 
                                                    "RING_KO H3K4me3 HiChIP",
                                                    "EED_KO H3K4me3 HiChIP"),
                             bin_size = 1e4,
                             cut_min_int = 1,
                             distance_limits = c(-2e5, 2e5), 
                             annotation_list = annotation_list)
    ggsave(paste0("figs/FigS8_HiChIP_interaction_", gene_name, ".png"))
  })
}


plot_interaction_profile <- function(gene_name,
                                     all_gene_gr,
                                     pairs_table_list, 
                                     pairs_table_titles,
                                     bin_size = 1e4,
                                     distance_limits = c(-2e5, 2e5),
                                     cut_min_int = 0,
                                     annotation_list = NULL, 
                                     anno_box_height = 0.2) {
  
  # get gene range
  gene_gr <- all_gene_gr[which(all_gene_gr$gene_name == gene_name)]
  
  # plot range
  view_interval <- gene_gr
  start(view_interval) <- start(view_interval) + distance_limits[1]
  end(view_interval) <- end(view_interval) + distance_limits[2]
  strand(view_interval) <- "*"
  
  g_list <- list()
  for(i in seq_along(pairs_table_list)) {
    g_list <- c(g_list,
                list(plot_arch(table_of_pairs = pairs_table_list[[i]], 
                               pairs_title = pairs_table_titles[i],
                               tmp_gene_gr = gene_gr, 
                               view_interval = view_interval,
                               bin_size = bin_size,
                               cut_min_int = cut_min_int)))
  }
  
  g_list <- c(g_list,
              list(plot_anno_interval(view_interval = view_interval,
                                      gene_gr = gene_gr, 
                                      gene_name = gene_name,
                                      annotation_list = annotation_list, 
                                      anno_box_height = anno_box_height)))
  
  do.call(cowplot::plot_grid, c(g_list, align = "v", ncol = 1))
  # cowplot::plot_grid(g1, g2, align = "v", nrow = 2, rel_heights = c(2/3, 1/3))
}


plot_arch <- function(table_of_pairs, 
                      pairs_title = "Condition1",
                      tmp_gene_gr,
                      gene_name, 
                      view_interval,
                      bin_size = 2e4,
                      cut_min_int = 4) {
  # get gene interactions
  tmp_gene_id <- tmp_gene_gr$gene_id
  chr_tmp <- as.character(seqnames(tmp_gene_gr))
  tss_tmp <- start(flank(tmp_gene_gr, width = 1, start = TRUE))
  
  
  idx <- which(table_of_pairs$V1 == chr_tmp &
                 table_of_pairs$V4 == chr_tmp &
                 (table_of_pairs$left_gene_id == tmp_gene_id |
                    table_of_pairs$right_gene_id == tmp_gene_id))
  
  tmp_interactions <- table_of_pairs[idx, ]
  
  # sort the pairs
  tmp_gr <- rbind(`colnames<-`(tmp_interactions[which(is.na(tmp_interactions$left_gene_id)), c(1,2,2,3)],
                               c("seqname", "start", "end", "strand")),
                  `colnames<-`(tmp_interactions[which(is.na(tmp_interactions$right_gene_id)), c(4,5,5,6)],
                               c("seqname", "start", "end", "strand"))) 
  # tmp_gr$strand <- "*"
  tmp_gr <- sort(GenomicRanges::makeGRangesFromDataFrame(tmp_gr), ignore.strand = TRUE)
  
  tmp_gr$distance <- start(tmp_gr) - tss_tmp
  
  # filter pairs
  tmp_gr <- tmp_gr[start(tmp_gr) > start(view_interval) & start(tmp_gr) < end(view_interval)]
  tmp_gr <- tmp_gr[countQueryHits(findOverlaps(tmp_gr, tmp_gene_gr + bin_size, ignore.strand = T)) == 0]
  
  # bin pairs
  bin_int <- findInterval(tmp_gr$distance,
                          seq(start(view_interval) - tss_tmp, end(view_interval) - tss_tmp, by = bin_size))
  
  bin_tbl <- table(bin_int)
  bin_tbl <- bin_tbl[bin_tbl > cut_min_int]
  
  # position x
  tss_bin <- (tss_tmp - start(view_interval)) / bin_size
  # bin_tbl <- bin_tbl[names(bin_tbl) != tss_bin]
  
  interaction_x <- as.numeric(names(bin_tbl))
  
  dat <- data.frame(x = c(rep(tss_bin, length(bin_tbl)), 
                          (tss_bin + interaction_x) / 2 - sign(interaction_x - tss_bin) * abs(interaction_x - tss_bin) * 0.41,
                          (tss_bin + interaction_x) / 2 - sign(interaction_x - tss_bin) * abs(interaction_x - tss_bin) * 0.2,
                          (tss_bin + interaction_x) / 2, 
                          (tss_bin + interaction_x) / 2 + sign(interaction_x - tss_bin) * abs(interaction_x - tss_bin) * 0.2,
                          (tss_bin + interaction_x) / 2 + sign(interaction_x - tss_bin) * abs(interaction_x - tss_bin) * 0.41,
                          interaction_x), 
                    y = c(rep(0, length(bin_tbl)), 
                          unname(bin_tbl) * 0.3,
                          unname(bin_tbl) * 0.91,
                          unname(bin_tbl),
                          unname(bin_tbl) * 0.91,
                          unname(bin_tbl) * 0.3,
                          rep(0, length(bin_tbl))),
                    group = rep(seq_along(bin_tbl), 7),
                    pairs = rep(bin_tbl / max(bin_tbl), 7)
                    )
  
  ggplot(data = dat, 
  aes(x = x, y = y, group = group, color = pairs)) +
    geom_line(size = 0.3,
              stat = "smooth",
              method = NULL,
              se = FALSE,
              lineend = "round") +
    annotate(geom = "text", x = -Inf, y = Inf, label = pairs_title, vjust = 1, hjust = -0.1) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, width(view_interval) / bin_size)) +
    scale_color_viridis_c(direction = -1, begin = 0.1, end = 0.9) +
    labs(y = "Contact frequency", color = "") +
    theme_minimal() +
    theme(legend.position = "none", 
          panel.grid = element_blank(), 
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.y = element_line(size = 0.2),
          axis.ticks.y = element_line(),
          plot.margin = unit(c(0, 2, 2, 2), "pt"))
  
}



plot_anno_interval <- function(view_interval,
                               gene_gr, 
                               gene_name,
                               annotation_list = NULL, 
                               anno_box_height = 0.1) {
  
  tss_gr <- flank(gene_gr, width = 1, start = T)
  
  annotation_list <- suppressWarnings(sapply(annotation_list, function(x) x[countQueryHits(findOverlaps(x, view_interval)) > 0]))
  annotation_list <- c(annotation_list, "Gene" = list(gene_gr))
  
  anno_tile_param <- data.frame(starts = unlist(sapply(annotation_list, start)),
                                ends = unlist(sapply(annotation_list, end)),
                                ys = rep(seq_along(annotation_list), times = lengths(annotation_list)),
                                anno = rep(names(annotation_list), times = lengths(annotation_list)))
  
  num_of_annos <- length(annotation_list)
  
  ggplot(anno_tile_param, aes(xmin = starts, xmax = ends, 
                              ymin = ys - anno_box_height, ymax = ys + anno_box_height,
                              fill = anno)) +
    geom_rect() +
    annotate(geom = "text",
             x = (start(gene_gr) + end(gene_gr)) / 2,
             y = num_of_annos,
             hjust = 0.5,
             vjust = 1.7,
             cex = 4,
             label = gene_name) +
    # add TSS arrow
    geom_segment(aes(x = start(tss_gr), 
                     y = num_of_annos + .5,
                     xend = start(tss_gr) + ifelse(strand(tss_gr) == "+", 1, -1) * width(view_interval) / 100, 
                     yend = num_of_annos + .5), 
                 size = 0.3,
                 arrow = arrow(length = unit(2, "pt"))) +
    geom_segment(aes(x = start(tss_gr),
                     y = num_of_annos,
                     xend = start(tss_gr),
                     yend = num_of_annos + .5),
                 size = 0.3) +
    
    scale_x_continuous(name = "",
                       expand = c(0, 0), 
                       limits = c(start(view_interval), end(view_interval)), 
                       labels = scales::comma_format()) +
    scale_y_continuous(name = "", 
                       breaks = unique(anno_tile_param$ys),
                       labels = unique(anno_tile_param$anno), 
                       limits = c(0, max(anno_tile_param$ys) + 1)) +
    theme_minimal() +
    theme(legend.position = "none", 
          panel.grid = element_blank(), 
          axis.ticks = element_line(), 
          axis.line.x = element_line(),
          axis.line.y = element_line(size = 0.2),
          axis.text.y = element_text(face = "bold"),
          axis.title.x = element_blank(),
          plot.margin = unit(c(2, 2, 0, 2), "pt"))
  
}




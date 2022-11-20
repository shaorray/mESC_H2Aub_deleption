# Rui Shao Nov 2021

# ---------------------------------------- occupancy changes ---------------------------------------- #

# all genes
conds <- c("BAP1_P6","BAP1_P12", "Ezh2i_1d", "Ezh2i_2d", "Ezh2i_7d", 
           "Ezh2i_1d_BAP1_6h", "Ezh2i_1d_BAP1_12h", 
           "Ring1b_CKO")
sample_names <- c("H2Aub", "H3K27me3", "Ring1b", "Ezh2", "H3K4me3", "H3K27ac", "Pol II-NTD", "Pol II-S5p")

dat_tile <- NULL
for (i in sample_names) {
  idx <- grepl(paste0("^", i, "_"), colnames(ChIP_TSS_B1_mat_log2FC)) & !grepl("Trp", colnames(ChIP_TSS_B1_mat_log2FC))
  print(sum(idx))
  dat_tile <- rbind(dat_tile, colMedians(ChIP_TSS_B1_mat_log2FC[, which(idx)], na.rm = TRUE))
}

colnames(dat_tile) <- c("Ezh2i_1d_BAP1_12h", "BAP1_P12", "Ezh2i_1d_BAP1_6h", "BAP1_P6",
                        "Ezh2i_1d", "Ezh2i_2d", "Ezh2i_7d", "Ring1b_CKO")
rownames(dat_tile) <- sample_names


dat_tile %>% t() %>%
  reshape::melt() %>%
  dplyr::mutate(X2 = factor(X2, levels = rev(sample_names)), 
                X1 = factor(X1, levels = conds)) %>%
  ggplot(aes(x = X1, y = X2, fill = value)) +
  geom_tile(width = 0.9, height = 0.9) +
  geom_text(aes(label = round(value, 2)), color = "white") +
  scale_fill_gradientn(name = expression('log'[2]~FC), 
                       colours = c('blue', "white", add.alpha('red', 0.2)), 
                       values = (c(min(dat_tile), 0, max(dat_tile)) - min(dat_tile)) / diff(range(dat_tile)),
                       guide = guide_colorbar(label = TRUE,
                                              draw.ulim = TRUE, 
                                              draw.llim = TRUE,
                                              frame.colour = "black", 
                                              frame.linewidth = 1,
                                              ticks = TRUE,
                                              ticks.colour = "black", 
                                              ticks.linewidth = 1)) +
  xlab("") + ylab("") + ggtitle("Median TSS changes (all genes, n=23640)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank(), 
        axis.ticks = element_line()) 

ggsave(filename = "Fig1_Heatmap_Samples_norm_Batch_1_all_Genes_ChIP_TSS_log2FC.png",
       device = "png", path = "../fig1/figs", width = 6, height = 6)


# H2Aub repressed genes
tmp_genes <- intersect(rownames(ChIP_TSS_B1_mat_log2FC), H2Aub_repressed_genes)

dat_tile <- NULL
for (i in sample_names) {
  idx <- grepl(paste0("^", i, "_"), colnames(ChIP_TSS_B1_mat_log2FC)) & !grepl("Trp", colnames(ChIP_TSS_B1_mat_log2FC))
  print(sum(idx))
  dat_tile <- rbind(dat_tile, colMedians(ChIP_TSS_B1_mat_log2FC[tmp_genes, which(idx)], na.rm = TRUE))
}

colnames(dat_tile) <- c("Ezh2i_1d_BAP1_12h", "BAP1_P12", "Ezh2i_1d_BAP1_6h", "BAP1_P6",
                        "Ezh2i_1d", "Ezh2i_2d", "Ezh2i_7d", "Ring1b_CKO")
rownames(dat_tile) <- sample_names


dat_tile %>% t() %>%
  reshape::melt() %>%
  dplyr::mutate(X2 = factor(X2, levels = rev(sample_names)), 
                X1 = factor(X1, levels = conds)) %>%
  ggplot(aes(x = X1, y = X2, fill = value)) +
  geom_tile(width = 0.9, height = 0.9) +
  geom_text(aes(label = round(value, 2)), color = "white") +
  scale_fill_gradientn(name = expression('log'[2]~FC), 
                       colours = c('blue', "white", add.alpha('red', 0.2)), 
                       values = (c(min(dat_tile), 0, max(dat_tile)) - min(dat_tile)) / diff(range(dat_tile)),
                       guide = guide_colorbar(label = TRUE,
                                              draw.ulim = TRUE, 
                                              draw.llim = TRUE,
                                              frame.colour = "black", 
                                              frame.linewidth = 1,
                                              ticks = TRUE,
                                              ticks.colour = "black", 
                                              ticks.linewidth = 1)) +
  xlab("") + ylab("") + ggtitle("Median TSS changes (H2Aub repressed, n=1893)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank(), 
        axis.ticks = element_line()) 

ggsave(filename = "Fig1_Heatmap_Samples_norm_Batch_1_H2Aub_repressed_Genes_ChIP_TSS_log2FC.png",
       device = "png", path = "../fig1/figs", width = 6, height = 6)

# ---------------------------------- ChIP log2FC cross chromHMM states# ---------------------------------- #

ChromHMM <- importRanges("../data/mESC_E14_12_dense.annotated.mm10.bed")
ChromHMM$name <- gsub(".*_", "\\2", ChromHMM$name)

ChromHMM_mm9 <- liftOver(ChromHMM, import.chain("../data/liftOver_chains/mm10ToMm9.over.chain")) %>% unlist()

# batch 1
bw_files_1 <-list.files(path = "../data/bw/batch1/input_norm",
                        pattern = ".*bw$",
                        full.names = T)
  
bw_files_1 <- bw_files_1[!grepl("_Ezh2i|-Ezh2i|Trp|2D|7D|R1bKO|IN", bw_files_1)]


ChIP_ChromHMM_B1_mat <- .countBW(bw_files = bw_files_1,
                                 intervals = ChromHMM_mm9,
                                 blacklist = blacklist.gr,
                                 fast = FALSE) 

col_names <- gsub(".*input_norm/(.*).bw", "\\1", bw_files_1)
col_names[grep("H3K4m3", col_names)] <- gsub("H3K4m3", "H3K4me3", col_names[grep("H3K4m3", col_names)])
col_names[grep("Pol2-NTD", col_names)] <- gsub("Pol2-NTD", "Pol II-NTD", col_names[grep("Pol2-NTD", col_names)])
col_names[grep("Pol2-S5p", col_names)] <- gsub("Pol2-S5p", "Pol II-S5p", col_names[grep("Pol2-S5p", col_names)])
col_names[grep("Pol2-S2p", col_names)] <- gsub("Pol2-S2p", "Pol II-S2p", col_names[grep("Pol2-S2p", col_names)])
col_names[grep("H3K27m3", col_names)] <- gsub("H3K27m3", "H3K27me3", col_names[grep("H3K27m3", col_names)])
col_names[grep("Trp-9", col_names)] <- gsub("Trp-9", "Triptolide", col_names[grep("Trp-9", col_names)])
col_names[grep("BAP1-6", col_names)] <- gsub("BAP1-6", "P6", col_names[grep("BAP1-6", col_names)])
col_names[grep("BAP1-12", col_names)] <- gsub("BAP1-12", "P12", col_names[grep("BAP1-12", col_names)])
col_names[grep("Ezh2i-1D", col_names)] <- gsub("Ezh2i-1D", "Ezh2i", col_names[grep("Ezh2i-1D", col_names)])
col_names[grep("_NT", col_names)] <- gsub("_NT", "_P0", col_names[grep("_NT", col_names)])

colnames(ChIP_ChromHMM_B1_mat) <- col_names
ChIP_ChromHMM_B1_mat <- ChIP_ChromHMM_B1_mat / b1_RPGC_scale

# aggregate by chromHMM states
ChIP_ChromHMM_B1_mat_median <- apply(ChIP_ChromHMM_B1_mat,
                                     MARGIN = 2, 
                                     function(x) aggregate(x, list(ChromHMM_mm9$name), median, na.rm = T)[, 2])

# for (i in unique(gsub("_.*", "", colnames(ChIP_ChromHMM_B1_mat_median))) ) { # max-min normalization
#   idx <- grep(paste0(i, "_"), colnames(ChIP_ChromHMM_B1_mat_median))
# 
#   tmp <- ChIP_ChromHMM_B1_mat_median[, idx]
#   ChIP_ChromHMM_B1_mat_median[, idx] <- (tmp - min(tmp)) / (max(tmp) - min(tmp))
# }

# to log2 Fold-change
ChIP_ChromHMM_B1_mat_cmb <- apply(ChIP_ChromHMM_B1_mat, 
                                  2, function(x) aggregate(x, list(ChromHMM_mm9$name), sum, na.rm = T)[, 2])
rownames(ChIP_ChromHMM_B1_mat_cmb) <- rownames(ChIP_ChromHMM_B1_mat_median) <- sort(unique(ChromHMM_mm9$name))

#
ChIP_ChromHMM_B1_mat_cmb_norm_lfc <- NULL
for (i in unique(gsub("_.*", "", colnames(ChIP_ChromHMM_B1_mat_cmb))) ) {
  idx <- grep(paste0("^", i, "_"), colnames(ChIP_ChromHMM_B1_mat_cmb))
  print(length(idx))
  ChIP_ChromHMM_B1_mat_cmb_norm_lfc <- cbind(ChIP_ChromHMM_B1_mat_cmb_norm_lfc,
                                                 log2(ChIP_ChromHMM_B1_mat_cmb[, idx] / 
                                                        ChIP_ChromHMM_B1_mat_cmb[, idx[3]]))
}

dat <- reshape::melt(ChIP_ChromHMM_B1_mat_cmb_norm_lfc[, !grepl("Pol II-S|H3K27ac", colnames(ChIP_ChromHMM_B1_mat_cmb_norm_lfc))])
dat$target <- factor(gsub("(.*)_.*", "\\1", dat$X2), 
                     levels = c("H2Aub", "H3K27me3", "Ring1b", "Ezh2",  "H3K4me3", "Pol II-NTD"))
dat$condition <- factor(gsub(".*_(.*)", "\\1", dat$X2), levels = rev(c("P0", "P6", "P12")))

dat$Occupancy <- reshape::melt(ChIP_ChromHMM_B1_mat_median[, !grepl("Pol II-S|H3K27ac", colnames(ChIP_ChromHMM_B1_mat_median))])$value

if (FALSE) {
  dat$Occupancy[dat$condition != "P0"] <- NA
  
  ggplot(dat, aes(x = X1, y = condition, fill = value)) +
    geom_tile(height = 0.95, width = 0.95) +
    geom_point(aes(size = Occupancy, alpha = Occupancy), 
               color = "blue4", pch = 19) +
    facet_grid(target~.) +
    scale_fill_gradientn(name = "log2FC", 
                         colours = c('blue', 'white', 'grey99', 'grey95'), 
                         values = (c(min(dat$value), 0, 0.1, max(dat$value)) - min(dat$value)) /
                           diff(range(dat$value)),
                         guide = guide_colorbar(label = TRUE,
                                                draw.ulim = TRUE, 
                                                draw.llim = TRUE,
                                                frame.colour = "black", 
                                                frame.linewidth = 1,
                                                ticks = TRUE,
                                                ticks.colour = "black", 
                                                ticks.linewidth = 1)) +
    xlab("") + ylab("") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12.5, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12),
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, size = 1), 
          strip.text = element_text(size = 12) )
  
  ggsave(filename = paste0("Fig1_heatmap_chromHMM_log2FC_sample_b1_bw.png"), 
         path = "../fig1/figs/",
         device = "png", width = 6, height = 10)
}

ggplot(dat, aes(x = X1, y = condition, fill = Occupancy)) +
  geom_tile(height = 0.9, width = 0.9) +
  facet_grid(target~.) +
  scale_fill_gradientn(name = "Occupancy", colours = brewer.pal(9, "Reds")) +
  xlab("") + ylab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12.5, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_line(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, size = 1), 
        strip.text = element_text(size = 12) )

ggsave(filename = paste0("Fig1_heatmap_chromHMM_RPGC_sample_b1_bw_tile.png"), 
       path = "../fig1/figs/",
       device = "png", width = 5, height = 7.5)


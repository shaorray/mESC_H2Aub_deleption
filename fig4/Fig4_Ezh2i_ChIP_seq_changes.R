# Rui Shao 2022 Jun

# evaluate H3K27me3's role in PcG mediated repression

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
source("../util/getCoverage.R")

# -------------------------------- compare ChIP TSS binding -----------------------------------#
# plot tss heatmap

sample_names <- c("H2Aub", "H3K27me3", "Ring1b", "Ezh2", "H3K4me3", "H3K27ac", "Pol II")
dat_tile <- NULL
for (i in sample_names) {
  idx <- grepl(paste0("^", i, "_"), colnames(ChIP_TSS_B1_BAM_all_log2FC)) & 
    !grepl("NT|BAP1-6-Ezh2i", colnames(ChIP_TSS_B1_BAM_all_log2FC))
  dat_tile <- rbind(dat_tile, colMedians(ChIP_TSS_B1_BAM_all_log2FC[H2Aub_repressed_genes, idx]))
}
colnames(dat_tile) <- c("BAP1_6h", "BAP1_12h", "Ezh2i_1d", "Ezh2i_2d", "Ezh2i_7d", 
                        "Ezh2i_1d_BAP1_12h", "Ring1b_KO", "Triptolide_9h")
rownames(dat_tile) <- sample_names
# dat_tile <- dat_tile[, !grepl("Ezh2i|Trp", colnames(dat_tile))]

dat_tile %>% t() %>%
  reshape::melt() %>%
  dplyr::mutate(X2 = factor(X2, levels = rev(sample_names)), 
                X1 = factor(X1, levels = (colnames(dat_tile)))) %>%
  ggplot(aes(x = X1, y = X2, fill = value)) +
  geom_tile(width = 0.9, height = 0.9) +
  geom_text(aes(label = round(value, 2)), color = "white") +
  scale_fill_gradientn(name = "log2FC", 
                       colours = c('blue', "white", 'red'), 
                       values = (c(min(dat_tile), 0, max(dat_tile)) - min(dat_tile)) /
                         diff(range(dat_tile)),
                       guide = guide_colorbar(label = TRUE,
                                              draw.ulim = TRUE, 
                                              draw.llim = TRUE,
                                              frame.colour = "black", 
                                              frame.linewidth = 1,
                                              ticks = TRUE,
                                              ticks.colour = "black", 
                                              ticks.linewidth = 1)) +
  xlab("") + ylab("") + ggtitle("H2Aub repressed genes (n = 1903)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank(), 
        axis.ticks = element_line())

ggsave(filename = "Fig4_Heatmap_Samples_norm_Batch_1_all_Genes_ChIP_TSS_log2FC.pdf",
       device = "pdf", path = "../fig4/figs", width = 7, height = 6)

# -------------------------------------- plot chromHMM ----------------------------------------#
ChromHMM <- importRanges("/mnt/0E471D453D8EE463/genomeDir/chromHMM/mESC_E14_12_dense.annotated.mm10.bed")
ChromHMM$name <- gsub(".*_", "\\2", ChromHMM$name)
ChromHMM_mm9 <- liftOver(ChromHMM, import.chain("../data/liftOver_chains/mm10ToMm9.over.chain")) %>% unlist()

# BAM files
bam_files_1_tmp <-
  list.files(path = "/mnt/E0767589767560E8/UPPMAX/PHILIP_LOF_20201122/bam",
             pattern = "(_BAP1|Ezh2i|_NT).*RC.mm9.fltd.bam$",
             full.names = T)

ChIP_ChromHMM_B1_BAM_tmp <- .countBam(bam_files = bam_files_1_tmp,
                                      intervals = ChromHMM_mm9,
                                      stranded = F, 
                                      paired.end = 'ignore') 
colnames(ChIP_ChromHMM_B1_BAM_tmp) <- gsub("_RC.mm9.fltd.bam", "", colnames(ChIP_ChromHMM_B1_BAM_tmp))

# aggregate by chromHMM states
ChIP_ChromHMM_B1_mat_cmb2 <- apply(ChIP_ChromHMM_B1_BAM_tmp, 2, 
                                   function(x) aggregate(x, list(ChromHMM_mm9$name), sum)[, 2])

rownames(ChIP_ChromHMM_B1_mat_cmb2) <- sort(unique(ChromHMM_mm9$name))

sf <- size_factor_cal(ChIP_ChromHMM_B1_mat_cmb2[, grep("^IN", colnames(ChIP_ChromHMM_B1_mat_cmb2))])

ChIP_ChromHMM_B1_mat_cmb_norm <- NULL
for (i in unique(gsub("_.*", "", colnames(ChIP_ChromHMM_B1_mat_cmb2)))) {
  ChIP_ChromHMM_B1_mat_cmb_norm <- 
    cbind(ChIP_ChromHMM_B1_mat_cmb_norm,
          sweep(
            ChIP_ChromHMM_B1_mat_cmb2[, grep(paste0("^", i, "_"), colnames(ChIP_ChromHMM_B1_mat_cmb2))], 
            MARGIN = 2, STATS = sf, FUN = "/"))
}
ChIP_ChromHMM_B1_mat_cmb_norm <- ChIP_ChromHMM_B1_mat_cmb_norm[, -grep("^IN", colnames(ChIP_ChromHMM_B1_mat_cmb_norm))]

ChIP_ChromHMM_B1_mat_cmb_norm_lfc <- NULL
for (i in unique(gsub("_.*", "", colnames(ChIP_ChromHMM_B1_mat_cmb_norm))) ) {
  idx <- grep(paste0("^", i, "_"), colnames(ChIP_ChromHMM_B1_mat_cmb_norm))
  print(length(idx))
  ChIP_ChromHMM_B1_mat_cmb_norm_lfc <- cbind(ChIP_ChromHMM_B1_mat_cmb_norm_lfc,
                                             log2(ChIP_ChromHMM_B1_mat_cmb_norm[, idx] / 
                                                    ChIP_ChromHMM_B1_mat_cmb_norm[, idx[8]]))
}
ChIP_ChromHMM_B1_mat_cmb_norm_lfc <- ChIP_ChromHMM_B1_mat_cmb_norm_lfc[, -grep("NT$", colnames(ChIP_ChromHMM_B1_mat_cmb_norm_lfc))]


dat <- reshape::melt(ChIP_ChromHMM_B1_mat_cmb_norm_lfc)
dat$target <- gsub("(.*)_.*", "\\1", dat$X2)
dat$target <- gsub("H3K27m3", "H3K27me3", dat$target)
dat$target <- gsub("H3K4m3", "H3K4me3", dat$target)
dat$target <- gsub("Pol2", "Pol II", dat$target)
dat$target <- gsub("Pol II-NTD", "Pol II", dat$target)

dat$target <- factor(dat$target, 
                     levels = c("H2Aub", "H3K27me3", "Ring1b", "Ezh2",  "H3K27ac", "H3K4me3",
                                "Pol II", "Pol II-S2p", "Pol II-S5p"))
dat$condition <- factor(gsub(".*_(.*)", "\\1", dat$X2), 
                        levels = rev(c("Ezh2i-1D", "Ezh2i-2D", "Ezh2i-7D",
                                       "BAP1-6", "BAP1-12",
                                       "BAP1-6-Ezh2i", "BAP1-12-Ezh2i")))

ggplot(dat, aes(x = X1, y = condition, fill = value)) +
  geom_tile(height = 0.95, width = 0.95) +
  facet_grid(target~.) +
  scale_fill_gradientn(name = "log2FC", 
                       colours = c('blue', "white", 'red'), 
                       values = (c(min(dat$value), 0, max(dat$value)) - min(dat$value)) /
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

ggsave(filename = paste0("Fig4_heatmap_chromHMM_log2FC_all_sample_b1_bw.pdf"), 
       path = "../fig4/figs/",
       device = "pdf", width = 6, height = 15)
ggsave(filename = paste0("fig4_heatmap_chromHMM_log2FC_all_sample_b1_bam.png"), 
       path = "../fig4/figs/",
       device = "png", width = 6, height = 15)

# --------------------------------------- plot ridge ------------------------------------------#

plot_ridge <- function(mat_log2FC, 
                       .sample, 
                       .row_name = c("P6", "P12"), 
                       .range = c(-5, 5)) {
  idx <- grep(paste0("^", .sample), colnames(mat_log2FC))[seq_along(.row_name) + 1]
  tmp <- mat_log2FC[, idx]
  p_vals <- apply(tmp, 2, function(x) get_sig_symbols(t.test(x)$p.val) )
  
  # .breaks <- seq(ceiling(min(tmp)), floor(max(tmp)))
  
  tmp %>% 
    as.matrix() %>%
    `colnames<-`(., .row_name) %>%
    reshape::melt() %>%
    dplyr::mutate(X2 = factor(X2, levels = rev(.row_name))) %>%
    ggplot(aes(x = value, y = X2, group = X2)) +
    geom_vline(xintercept = 0, color = "grey80", size = 2) +
    ggridges::geom_density_ridges(rel_min_height = 0.01, 
                                  quantile_lines = TRUE,
                                  quantiles = 2, lty = 1,
                                  fill = add.alpha(sample_colors[.sample], 0.8)) +
    annotate(geom = "text", 
             x = Inf, y = seq_along(.row_name),
             hjust = 1, label = rev(p_vals)) +
    # scale_x_continuous(breaks = .breaks, labels = 2^.breaks) +
    scale_fill_manual(values = unname(rep(sample_colors[1], 3))) +
    xlim(.range) +
    xlab("") + ylab("") + ggtitle(.sample) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text = element_text(size = 13))
}


ggsave(filename = "Fig4_Ridge_plot_Samples_norm_Batch_1_Derepressed_Genes.pdf", 
       do.call(grid.arrange, c(lapply(c("H2Aub", "H3K27me3", "Ring1b", "Ezh2",
                                        "H3K4me3", "H3K27ac", "Pol II",
                                        "Pol II-S2p", "Pol II-S5p"), 
                                      plot_ridge, 
                                      mat_log2FC = ChIP_TSS_B1_BAM_all_log2FC[derepressed_gene_ids, 
                                                                              !grepl("BAP1-6-Ezh2i|Trp", colnames(ChIP_TSS_B1_BAM_all_log2FC))], 
                                      .row_name = c("P6", "P12", 
                                                    "Ezh2i_1d", "Ezh2i_2d", "Ezh2i_7d", 
                                                    "Ezh2i_1d_P12",
                                                    "Ring1b_KO")),
                               nrow = 1)),
       device = "pdf", path = "../fig4/figs", width = 27, height = 6)



# --------------------------------------- plot fingerprint ------------------------------------------#


plot_fingerprint <- function(x) {
  
  idx <- match(paste0(x, c("_NT", "_BAP1-12", "_Ezh2i-1D", "_BAP1-12-Ezh2i")),
               colnames(ChIP_TSS_B1_mat))
  
  dat <- ChIP_TSS_B1_mat[, idx] %>% "*"(25) %>% "+"(1) %>% log10() 
  colnames(dat) <- c("P0", "P12", "Ezh2i", "Ezh2i_P12")
  
  q_995 <- quantile(dat, 0.995, na.rm = TRUE)
  
  # use the order of column 1 as reference
  dat <- as.data.frame(dat[order(dat[, 1], decreasing = TRUE), ])
  dat$x_order <- seq_len(nrow(dat))
  
  
  reshape::melt(dat, id.vars = "x_order") %>%
    # `colnames<-`(., c("x_order", "Sample", "value"))
    ggplot(aes(x = x_order, y = value, group = variable, color = variable)) +
    geom_point(cex = 0.1, alpha = 0.1) + 
    geom_smooth(se = FALSE, alpha = 0.8) +
    xlab("All gene by P0 rank (n=23640)") + ylab("log10(TSS occupancy + 1)") +
    ggtitle(x) +
    ylim(c(0, q_995)) +
    scale_x_continuous(breaks = c(1, nrow(dat)), labels = c(1, nrow(dat))) +
    scale_color_manual(values = c("#EB4037", "#342F59", "#4A7B7B", "#B88828")) +
    scale_linetype_manual(values = c(1, 1, 2, 2)) +
    theme_minimal() +
    theme(axis.line = element_line(),
          axis.ticks = element_line(), 
          axis.text = element_text(size = 12),
          panel.grid = element_blank())
    
    
    ggsave(paste0("Fig4_fingerprint_plot_all_gene_TSS_", x, ".png"), 
           device = "png", path = "../fig4/figs", width = 4, height = 3)
  
}


plot_fingerprint("Ring1b")
plot_fingerprint("Ezh2")
plot_fingerprint("H3K27m3")
plot_fingerprint("H2Aub")

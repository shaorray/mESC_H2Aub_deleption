# Rui Shao 2021 May
# Figure 3
# DE genes comparison

# ----------------------------------------------------------------------------------- #
dat <- data.frame(gene_id = Reduce(c, DE_gene.list[grep("up", names(DE_gene.list))]),
                  Time = sapply(c("P12", "P12_C12", "P12_C24", "P12_C36"), 
                                function(i) {
                                  .time = lengths(DE_gene.list[grep("up", names(DE_gene.list))])
                                  .time = .time[paste0(i, "_up")]
                                  rep(i, .time)
                                }) %>% Reduce("c", .)
)
dat$Time <- factor(dat$Time, c("P12", "P12_C12", "P12_C24", "P12_C36"))
dat$`Derepressed genes` <- "BAP1_specific"
dat$`Derepressed genes`[dat$gene_id %in% Ring1b_CKO_genes] <- "PRC1_CKO"
dat <- dat[dat$gene_id %ni% CpK_DE.gr$gene_id, ]
# dat$`Derepressed genes`[dat$gene_id %in% c(GSE132753_DE_genes_CPM, GSE134053_DE_genes_CPM)] <- "Ring1b_CPM"

# barplot, gene counts frequency
plot_freq_types <- function(dat, .by = "Group", .levels = NA, .split = "Type", .color = 1:10, top_n = 10) {
  # Args:
  #   dat: gene_id (unique), .by (row features), .levels (row names), .split (bar split feature)
  if (is.na(.levels[1])) {
    .levels <- seq_len(as.character(nchar(dat[1, .by])))
  }
  
  # limit less frequent groups
  top_groups <- dplyr::count(dat, Group, sort = TRUE) %>% head(top_n) %>% "$"("Group")
  dat <- dat[dat$Group %in% top_groups, ]
  
  dat$Group <- factor(dat[, .by], levels = names(sort(table(dat[, .by]), decreasing = T)))
  dat$Type <- dat[, .split]
  
  # plot1
  dat_g1 <- dplyr::count(dat, Group, Type, sort = TRUE)
  dat_g1_text <- data.frame(Group = unique(dat_g1$Group),
                            Fraction = sapply(unique(dat_g1$Group), 
                                              function(x) {
                                                tmp <- dat_g1[dat_g1$Group == x, ]
                                                round(tmp[tmp$Type != "BAP1_specific", 'n'] / sum(tmp[, 'n']), 2) * 100
                                              }) %>% paste0(., "%"),
                            Type = "BAP1_specific")
  Group_n <- dplyr::count(dat, Group, sort = TRUE)
  dat_g1_text$n  <- Group_n[match(dat_g1_text$Group, Group_n$Group), "n"]
  
  g1 <- ggplot(dat_g1, aes(x = Group, y = n, fill = Type)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(data = dat_g1_text, 
              aes(x = Group, y = n, label = Fraction),
              vjust = -0.25, hjust = 0.5,
              size = 3) +
    xlab("") + ylab("Number of genes") +
    ggpubr::theme_pubclean() +
    scale_fill_manual(name = .split, values = c("grey70", .color)) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.margin = margin(c(0,1,1,1)))
  
  # plot2
  dat_tile <- NULL
  for (m in seq_along(levels(dat$Group))) {
    tmp <- cbind(x = m, y = .levels, 
                 color = ifelse(strsplit(levels(dat$Group)[m], "")[[1]] == 1, 1, 0))
    dat_tile <- rbind(dat_tile, tmp)
  }
  
  dat_tile <- as.data.frame(dat_tile) %>% 
    dplyr::mutate(x = factor(x, levels = order(as.numeric(x)))) %>%
    dplyr::mutate(y = factor(y, levels = rev(.levels)))
  
  dat_seg <- NULL
  for (i in unique(dat_tile$x)) {
    tmp <- dat_tile[dat_tile$x == i, ]
    tmp <- tmp[tmp$color == 1, ]
    if (nrow(tmp) > 1) {
      dat_seg <- rbind(dat_seg,
                       data.frame(x = tmp$x[-nrow(tmp)], y = tmp$y[-nrow(tmp)],
                                  xend = tmp$x[-1], yend = tmp$y[-1]))
    }
  }
  levels(dat_seg$x) <- levels(dat_seg$xend) <-levels(dat_tile$x)
  levels(dat_seg$y) <- levels(dat_seg$yend) <-levels(dat_tile$y)
  dat_seg$color <- 1
  
  g2 <- dat_tile %>%
    ggplot(aes(x = x, y = y, color = factor(color), group = x)) +
    geom_point(size = 4) +
    geom_segment(data = dat_seg, aes(x = x, y = y, xend = xend, yend = yend), lty = 2) +
    scale_color_manual(values = c("grey90", .color)) +
    xlab("") + ylab("") +
    ggpubr::theme_pubclean() +
    theme(panel.grid = element_blank(), 
          legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.margin = margin(c(1,1,0,1)))
  
  cowplot::plot_grid(g1, g2, ncol = 1, align = "v", rel_heights = c(1, 0.4))
}


dat2 <- dat[dat$Time == "P12", ]
for(x in c("P12_C12", "P12_C24", "P12_C36")) {
  dat2 <- rbind(dat2, dat[dat$Time == x & dat$gene_id %ni% dat2$gene_id, ])
}

dat2$Group <- paste0(as.numeric(dat2$gene_id %in% DE_gene.list$P12_up),
                     paste0(as.numeric(dat2$gene_id %in% DE_gene.list$P12_C12_up), 
                            paste0(as.numeric(dat2$gene_id %in% DE_gene.list$P12_C24_up),
                                   as.numeric(dat2$gene_id %in% DE_gene.list$P12_C36_up))))
lengths(DE_gene.list)

plot_freq_types(dat2, .by = "Group", 
                .levels = c("P12 (n=1342)", "P12_C12 (n=1711)", "P12_C24 (n=2160)", "P12_C36 (n=2653)"), 
                .split = "Derepressed genes", .color = "cyan3", top_n = 9)

ggsave(filename = paste0("Fig3_frequency_bar_derepressed_genes.pdf"), 
       path = "../fig3/figs",
       device = "pdf", width = 7, height = 5)

if (TRUE) {
  # GO
  source("../util/go.R")
  DE_gene_BAP1_specific.list <- lapply(DE_gene.list[grep("up", names(DE_gene.list))], 
                                       function(x) x[x %ni% Ring1b_CKO_genes & x %ni% CpK_DE.gr$gene_id])
  names(DE_gene_BAP1_specific.list) <- gsub("_up", "", names(DE_gene_BAP1_specific.list))
  enriched_GO_list_heatmap(DE_gene_BAP1_specific.list,
                           title = "BAP1 specific up-regulation",
                           ontology = "BP", top_n_term = 8,
                           colorset = c("white", "grey30"))
  ggsave(filename = paste0("Fig3_GO_BAP1_Specific_up_gene.png"), 
         path = "../fig3/figs",
         device = "png", width = 7, height = 4.5)
  
  DE_gene_Ring1b_CKO_Ov.list <- lapply(DE_gene.list[grep("up", names(DE_gene.list))], 
                                       function(x) x[x %in% Ring1b_CKO_genes & x %ni% CpK_DE.gr$gene_id])
  names(DE_gene_Ring1b_CKO_Ov.list) <- gsub("_up", "", names(DE_gene_Ring1b_CKO_Ov.list))
  enriched_GO_list_heatmap(DE_gene_Ring1b_CKO_Ov.list,
                           title = "Ring1b_CKO overlapped",
                           ontology = "BP", top_n_term = 8,
                           colorset = c("white", "cyan3"))
  ggsave(filename = paste0("Fig3_GO_BAP1_Ring1b_CKO_derepressed_gene.png"), 
         path = "../fig3/figs",
         device = "png", width = 6, height = 4.5)
}



# ----------------------------------------------------------------------------------- #
if (F) {
  # PRC1 CKO / CPM, jaccard index matrix
  
  get_jaccard <- function(v1, v2) length(intersect.Vector(v1, v2)) / length(unique(c(v1, v2)))
  
  BAP1_PRC1_DE_jaccard_mat <- matrix(NA, nrow = 4, ncol = 4)
  rownames(BAP1_PRC1_DE_jaccard_mat) <- c("PRC1 CPM\nBlackledge et al.", "PRC1 CPM\nTamburri et al.",
                                          "PRC1 CKO\nBlackledge et al.", "PRC1 CKO\nTamburri et al.")
  colnames(BAP1_PRC1_DE_jaccard_mat) <- c("P12", "P12_C12", "P12_C24", "P12_C36")
  
  for (j in colnames(BAP1_PRC1_DE_jaccard_mat)) {
    BAP1_PRC1_DE_jaccard_mat[, j] <- c(get_jaccard(GSE132753_DE_genes_CPM,
                                                   dat$gene_id[dat$Time == j]),
                                       get_jaccard(GSE134053_DE_genes_CPM,
                                                   dat$gene_id[dat$Time == j]),
                                       get_jaccard(GSE132753_DE_genes_CKO,
                                                   dat$gene_id[dat$Time == j]),
                                       get_jaccard(GSE134053_DE_genes_CKO,
                                                   dat$gene_id[dat$Time == j]))
  }
  
  g1 <- BAP1_PRC1_DE_jaccard_mat %>%
    reshape::melt() %>%
    dplyr::mutate(X1 = factor(X1, levels = rev(c(rownames(BAP1_PRC1_DE_jaccard_mat))))) %>%
    ggplot(aes(x = X2, y = X1, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 2))) +
    scale_fill_viridis_c(name = "Jaccard", option = "B", direction = -1,
                         limits = c(0.07, 0.31)) +
    xlab("") + ylab("") +
    theme_minimal() +
    theme_setting +
    theme(panel.border = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  # PRC1 acute depletion, jaccard index matrix
  BAP1_PRC1_DE_jaccard_mat <- matrix(NA, nrow = 3, ncol = 4)
  rownames(BAP1_PRC1_DE_jaccard_mat) <- c("IAA 2h\nDobrinic et al.",
                                          "IAA 4h\nDobrinic et al.",
                                          "IAA 8h\nDobrinic et al.")
  colnames(BAP1_PRC1_DE_jaccard_mat) <- c("P12", "P12_C12", "P12_C24", "P12_C36")
  
  
  for (i in c("2h", "4h", "8h")) {
    for (j in colnames(BAP1_PRC1_DE_jaccard_mat)) {
      BAP1_PRC1_DE_jaccard_mat[grep(i, rownames(BAP1_PRC1_DE_jaccard_mat)), j] <-
        get_jaccard(GSE159399_DE_tab$gene_id[GSE159399_DE_tab$time == i], dat$gene_id[dat$Time == j])
    }
  }
  
  g2 <- BAP1_PRC1_DE_jaccard_mat %>%
    reshape::melt() %>%
    dplyr::mutate(X1 = factor(X1, levels = rev(rownames(BAP1_PRC1_DE_jaccard_mat)))) %>%
    ggplot(aes(x = X2, y = X1, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 2))) +
    scale_fill_viridis_c(name = "Jaccard", option = "B", direction = -1,
                         limits = c(0.07, 0.31)) +
    xlab("") + ylab("") +
    theme_minimal() +
    theme_setting +
    theme(panel.border = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(cowplot::plot_grid(g1, g2, ncol = 2),
         filename = "../figS4/figs/FigS4_Jaccard_index_BAP1_PRC1_CKO_CPM_IAA.pdf",
         device = "pdf", width = 12, height = 4)
}

# ----------------------------------------------------------------------------------- #
# de-repression gene sets by ChIP TSS peak matches
if (F) {
  # histone marks
  dat2$H3K4me3 <- ifelse(dat2$gene_id %in% H3K4me3_enriched_genes, 1, 0)
  dat2$H3K27me3 <- ifelse(dat2$gene_id %in% H3K27me3_enriched_genes, 1, 0)
  dat2$H2AZ <- ifelse(dat2$gene_id %in% H2AZ_enriched_genes, 1, 0)
  dat2$H2Aub <- ifelse(dat2$gene_id %in% H2Aub_enriched_genes, 1, 0)
  dat2$H33 <- ifelse(dat2$gene_id %in% H33_enriched_genes, 1, 0)
  
  dat2$Group <- paste0(dat2$H3K4me3, paste0(dat2$H3K27me3, paste0(dat2$H2Aub, paste0(dat2$H2AZ, dat2$H33))))
  
  plot_freq_types(dat2, .by = "Group", 
                  .levels = c("H3K4me3 (n=8656)", "H3K27me3 (n=8242)", "H2Aub (n=9831)", "H2A.Z (n=5625)", "H3.3 (n=5297)"), 
                  .split = "Derepressed genes", .color = c("lightcyan4"))
  
  ggsave(filename = paste0("FigS1_frequency_bar_derepressed_genes_histone.pdf"), 
         path = "../figS1/figs",
         device = "pdf", width = 7, height = 5)
}

if (F) {
  # PcG factors
  dat2$Ring1b <- ifelse(dat2$gene_id %in% Ring1b_binding_genes, 1, 0)
  dat2$Ezh2 <- ifelse(dat2$gene_id %in% Ezh2_binding_genes, 1, 0)
  dat2$Rybp <- ifelse(dat2$gene_id %in% Rybp_binding_genes, 1, 0)
  dat2$Cbx7 <- ifelse(dat2$gene_id %in% Cbx7_binding_genes, 1, 0)
  dat2$Rsf1 <- ifelse(dat2$gene_id %in% Rsf1_binding_genes, 1, 0)
  
  dat2$Group <- paste0(dat2$Ring1b, paste0(dat2$Ezh2, paste0(dat2$Cbx7, paste0(dat2$Rybp, dat2$Rsf1))))
  
  plot_freq_types(dat2, .by = "Group", 
                  .levels = c("Ring1b (n=8277)", "Ezh2 (n=8363)", "Cbx7 (n=5963)", "Rybp (n=6241)", "Rsf1 (n=6500)"), 
                  .split = "Derepressed genes", .color = c("steelblue2"))
  
  ggsave(filename = paste0("FigS1_frequency_bar_derepressed_genes_PRC1.pdf"), 
         path = "../figS1/figs",
         device = "pdf", width = 7, height = 5)
}

# ----------------------------------------------------------------------------------- #
# public data MA plots
if (F) {
  # plot MA
  g2.1 <- plot_ma(baseMean = rowMeans(GSE132753_DESeq2_res[use_gene_ids, 14:27]),
                  log2FC = GSE132753_DESeq2_res[use_gene_ids, ]$LFC_PRC1CKO, 
                  pval = GSE132753_DESeq2_res[use_gene_ids, ]$padj_PRC1CKO, 
                  p_col = 1, cut_baseMean = 0, 
                  title = "PRC1CKO")
  
  g2.2 <- plot_ma(baseMean = rowMeans(GSE132753_DESeq2_res[use_gene_ids, 14:27]),
                  log2FC = GSE132753_DESeq2_res[use_gene_ids, ]$LFC_PRC1CPM, 
                  pval = GSE132753_DESeq2_res[use_gene_ids, ]$padj_PRC1CPM, 
                  p_col = 1, cut_baseMean = 0, 
                  title = "PRC1CPM")
  
  ggsave(plot = grid.arrange(g2.1, g2.2, ncol = 1),
         filename = paste0("FigS1_MAplot_Klose_PRC1CKO_PRC1CPM.png"), 
         path = "figs",
         device = "png", width = 3, height = 6)
  
  # Ring1_OHT vs Ring1_ETA
  g2.3 <- plot_ma(baseMean = GSE134053_DE_res$baseMean,
                  log2FC = GSE134053_DE_res$log2FoldChange, 
                  pval = GSE134053_DE_res$padj, 
                  p_col = 1, cut_baseMean = 0, 
                  title = "Ring1b_OHT_provided")
  
  g2.4 <- plot_ma(baseMean = GSE134053_DE_res.list$Ring1b_OHT$baseMean,
                  log2FC = GSE134053_DE_res.list$Ring1b_OHT$log2FoldChange, 
                  pval = GSE134053_DE_res.list$Ring1b_OHT$padj, 
                  p_col = 1, cut_baseMean = 0, 
                  title = "Ring1b_OHT_rerun")
  
  g2.5 <- plot_ma(baseMean = GSE134053_DE_res.list$Ring1b_WT_OHT$baseMean,
                  log2FC = GSE134053_DE_res.list$Ring1b_WT_OHT$log2FoldChange, 
                  pval = GSE134053_DE_res.list$Ring1b_WT_OHT$padj, 
                  p_col = 1, cut_baseMean = 0, 
                  title = "Ring1b_OHT_vs_WT_OHT_rerun")
  
  g2.6 <- plot_ma(baseMean = GSE134053_DE_res.list$I53S_OHT$baseMean,
                  log2FC = GSE134053_DE_res.list$I53S_OHT$log2FoldChange, 
                  pval = GSE134053_DE_res.list$I53S_OHT$padj, 
                  p_col = 1, cut_baseMean = 0, 
                  title = "I53S_OHT_rerun")
  
  g2.7 <- plot_ma(baseMean = GSE134053_DE_res.list$I53S_WT_OHT$baseMean,
                  log2FC = GSE134053_DE_res.list$I53S_WT_OHT$log2FoldChange, 
                  pval = GSE134053_DE_res.list$I53S_WT_OHT$padj, 
                  p_col = 1, cut_baseMean = 0, 
                  title = "I53S_OHT_vs_WT_OHT_rerun")
  
  ggsave(plot = grid.arrange(g2.4, g2.5, g2.6, g2.7, g2.3, ncol = 2),
         filename = paste0("FigS1_MAplot_Pasini_PRC1CKO_I53S.png"), 
         path = "figs",
         device = "png", width = 6, height = 9)
  
  g3.1.1 <- cbind(P12C24 = asinh(res.list$P12[use_gene_ids, "baseMean"]), 
                  Klose = asinh(rowMeans(GSE132753_DESeq2_res[use_gene_ids, 14:27]))) %>%
    as.data.frame() %>% dplyr::filter(complete.cases(.)) %>%
    ggplot(aes(x = P12C24, y = Klose, color = get_dens(X = P12C24, Y = Klose))) +
    geom_point() +
    scale_color_viridis() +
    ggtitle("baseMean") +
    ggpubr::theme_pubclean() +
    theme(legend.position = "none") 
  
  g3.1.2 <- cbind(P12C24 = trim_quantile(res.list$P12[use_gene_ids, "log2FoldChange"]), 
                  Klose = trim_quantile(GSE132753_DESeq2_res[use_gene_ids, "LFC_PRC1CKO"])) %>%
    as.data.frame() %>% dplyr::filter(complete.cases(.)) %>%
    ggplot(aes(x = P12C24, y = Klose, color = get_dens(X = P12C24, Y = Klose))) +
    geom_point() +
    scale_color_viridis() +
    ggtitle("log2FoldChange") +
    ggpubr::theme_pubclean() +
    theme(legend.position = "none") 
  
  g3.2.1 <- cbind(Pasini = asinh(GSE134053_DE_res.list$Ring1b_OHT[use_gene_ids, "baseMean"]), 
                  Klose = asinh(rowMeans(GSE132753_DESeq2_res[use_gene_ids, 14:27]))) %>%
    as.data.frame() %>% dplyr::filter(complete.cases(.)) %>%
    ggplot(aes(x = Pasini, y = Klose, color = get_dens(X = Pasini, Y = Klose))) +
    geom_point() +
    scale_color_viridis() +
    ggtitle("baseMean") +
    ggpubr::theme_pubclean() +
    theme(legend.position = "none") 
  
  g3.2.2 <- cbind(Pasini = trim_quantile(res.list$P12[use_gene_ids, "log2FoldChange"]), 
                  Klose = trim_quantile(GSE132753_DESeq2_res[use_gene_ids, "LFC_PRC1CKO"])) %>%
    as.data.frame() %>% dplyr::filter(complete.cases(.)) %>%
    ggplot(aes(x = Pasini, y = Klose, color = get_dens(X = Pasini, Y = Klose))) +
    geom_point() +
    scale_color_viridis() +
    ggtitle("log2FoldChange") +
    ggpubr::theme_pubclean() +
    theme(legend.position = "none") 
  
  g3.3.1 <- cbind(P12C24 = asinh(res.list$P12[use_gene_ids, "baseMean"]), 
                  Pasini = asinh(GSE134053_DE_res.list$Ring1b_OHT[use_gene_ids, "baseMean"])) %>%
    as.data.frame() %>% dplyr::filter(complete.cases(.)) %>%
    ggplot(aes(x = P12C24, y = Pasini, color = get_dens(X = P12C24, Y = Pasini))) +
    geom_point() +
    scale_color_viridis() +
    ggtitle("baseMean") +
    ggpubr::theme_pubclean() +
    theme(legend.position = "none") 
  
  g3.3.2 <- cbind(P12C24 = trim_quantile(res.list$P12[use_gene_ids, "log2FoldChange"]), 
                  Pasini = trim_quantile(GSE134053_DE_res.list$Ring1b_OHT[use_gene_ids, "log2FoldChange"])) %>%
    as.data.frame() %>% dplyr::filter(complete.cases(.)) %>%
    ggplot(aes(x = P12C24, y = Pasini, color = get_dens(X = P12C24, Y = Pasini))) +
    geom_point() +
    scale_color_viridis() +
    ggtitle("log2FoldChange") +
    ggpubr::theme_pubclean() +
    theme(legend.position = "none") 
  
  ggsave(plot = grid.arrange(g3.1.1, g3.2.1, g3.3.1, g3.1.2, g3.2.2, g3.3.2, ncol = 3),
         filename = paste0("FigS1_scatter_RNAseq_Klose_Pasini_bM_LFC.png"), 
         path = "figs/",
         device = "png", width = 15, height = 10)
  
}



# ------------------------------ GO comparison ---------------------------------- #

if (F) {
  # BAP1_PC and Ring1b_CKO overlapped
  enrichGeneSets(Ring1b_CKO_genes[Ring1b_CKO_genes %in% derepressed_gene_ids], 
                 method = "GO", ontology = "BP", is.GeneRatio = T, 
                 title = "Up-regulation of Ring1b_CKO in BAP1_PC")
  
  ggsave(filename = paste0("FigS4_GO_BP_BAP1_PC_and_Ring1b_CKO.png"), 
         path = "../figS4/figs/",
         device = "png", width = 6, height = 5)
  
  # BAP1_PC unique
  enrichGeneSets(derepressed_gene_ids[derepressed_gene_ids %ni% Ring1b_CKO_genes], 
                 method = "GO", ontology = "BP", is.GeneRatio = T, 
                 title = "Up-regulation of BAP1_PC not in Ring1b_CKO")
  
  ggsave(filename = paste0("FigS4_GO_BP_BAP1_PC_unique.png"), 
         path = "../figS4/figs/",
         device = "png", width = 6, height = 5)
}

# Ring1b_CKO unique
enrichGeneSets(Ring1b_CKO_genes[Ring1b_CKO_genes %ni% derepressed_gene_ids], 
               method = "GO", ontology = "BP", is.GeneRatio = TRUE, 
               title = "PRC1 CKO specific up-regulation")

ggsave(filename = paste0("FigS4_GO_BP_Ring1b_CKO_unique.png"), 
       path = "../figS4/figs/",
       device = "png", width = 8, height = 5)




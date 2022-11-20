
# ------------------------------------- functions --------------------------------- #
# rlog transformation
convert_rlog <- function(mat) {
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(round(as.matrix(mat)),
                                colData = data.frame(condition = colnames(mat)),
                                design = ~ condition)
  rld <- rlogTransformation(dds)
  rld_val <- assay(rld)
  return(rld_val)
}

over_lap_fraction <- function(subject_genes, query_genes_1, query_genes_2) {
  c(sum(subject_genes %in% query_genes_1), 
    sum(subject_genes %ni% c(query_genes_1, query_genes_2)),
    sum(subject_genes %in% query_genes_2)) / length(subject_genes) * 100
}

# ---------------------------------------- plot PCA -----------------------------------------------#
geneRC_RNA.norm <- sweep(geneRC_RNA, 2, size_factor_cal(flyRC_RNA), "/")
geneRC_RNA.norm <- sapply(unique(gsub("_R.*", "", .design)), 
                          function(sample)
                            rowMeans(geneRC_RNA.norm[, names(.design[grep(paste0(sample, "_R"), .design)])]))
geneRC_RNA.norm <- geneRC_RNA.norm[, c("P0", "P12", "P12_C12", "P12_C24", "P12_C36")]
geneRC_RNA.norm.rlog <- convert_rlog(geneRC_RNA.norm)

plot_pca <- function(mat, gene_ids) {
  mat <- mat[intersect.Vector(rownames(mat), gene_ids), ]
  pca <- prcomp(t(mat))
  pca$sdev <- pca$sdev^2
  pc_var <- pca$sdev[1:2] / sum(pca$sdev) * 100
  pc <- pca$x[, 1:2] %>% as.data.frame()
  grid_dat <- `colnames<-`(cbind(pc[-5, ], pc[c(2:5), ]), c("x1", "y1", "x2", "y2"))
  pc$Sample <- rownames(pc)
  
  ggplot(pc, aes(x = PC2, y = PC1, color = Sample)) +
    geom_segment(data = grid_dat,
                 mapping = aes(x = y1, y = x1, xend = y2, yend = x2),
                 arrow = arrow(type="open", angle=20, length = unit(x = 14, 'pt')),
                 size = 0.6, lty = 1,
                 color = add.alpha("black", 0.5)) + 
    geom_point(size = 4) +
    scale_color_viridis_d(begin = 0.04, end = 0.84) +
    ylab(paste0("PC1 ", round(pc_var[1]), "% variance")) +
    xlab(paste0("PC2 ", round(pc_var[2]), "% variance")) +
    theme_setting
}

plot_pca(log1p(geneRC_RNA.norm), rownames(geneRC_RNA.norm.rlog)[rowMeans(geneRC_RNA.norm.rlog) > 2 & rowVars(geneRC_RNA.norm.rlog) > quantile(rowVars(geneRC_RNA.norm.rlog), 0.75)])
ggsave(filename = "FigS4_RNAseq_PCA_derepressed.png", 
       path = "../figS4/figs/",
       device = "png", width = 5, height = 5)

# ----------------------------------------------------------------------------------- #
lengths(DE_gene.list)

g1.1 <- plot_ma(baseMean = res.gene.shrink.list$P12$baseMean,
                log2FC = res.gene.shrink.list$P12$log2FoldChange, 
                pval = res.gene.shrink.list$P12$padj,
                DE_log2FC = 1.5,
                p_col = 1, cut_baseMean = 5,
                up_num = 1342, down_num = 954,
                title = "P12")

g1.2 <- plot_ma(baseMean = res.gene.shrink.list$P12_C12$baseMean,
                log2FC = res.gene.shrink.list$P12_C12$log2FoldChange, 
                pval = res.gene.shrink.list$P12_C12$padj, 
                DE_log2FC = 1.5,
                p_col = 1, cut_baseMean = 5, 
                up_num = 1711, down_num = 954,
                title = "P12C12")

g1.3 <- plot_ma(baseMean = res.gene.shrink.list$P12_C24$baseMean,
                log2FC = res.gene.shrink.list$P12_C24$log2FoldChange, 
                pval = res.gene.shrink.list$P12_C24$padj, 
                DE_log2FC = 1.5,
                p_col = 1, cut_baseMean = 5, 
                up_num = 2160, down_num = 119,
                title = "P12C24")

g1.4 <- plot_ma(baseMean = res.gene.shrink.list$P12_C36$baseMean,
                log2FC = res.gene.shrink.list$P12_C36$log2FoldChange, 
                pval = res.gene.shrink.list$P12_C36$padj, 
                DE_log2FC = 1.5,
                p_col = 1, cut_baseMean = 5, 
                up_num = 2653, down_num = 78,
                title = "P12C36")

ggsave(plot = grid.arrange(g1.1, g1.2, g1.3, g1.4, ncol = 4),
       filename = paste0("Fig3_MAplot_RNAseq_PC_DE.png"), 
       path = "../fig3/figs/",
       device = "png", width = 16, height = 4)

# ------------------------------------ DE gene lists intersection ------------------------------ #

# plot DE gene frequency with target gene subset, e.g. PcG
plot_target_gene_frequency <- function(enriched_genes, .ylab) {
  dat <- data.frame(Time = rep(c("P12", "P12C12", "P12C24", "P12C36"), each = 3),
                    Freq = c(over_lap_fraction(enriched_genes, DE_gene.list$P12_C12_up, DE_gene.list$P12_down),
                             over_lap_fraction(enriched_genes, DE_gene.list$P12_C12_up, DE_gene.list$P12_C24_down),
                             over_lap_fraction(enriched_genes, DE_gene.list$P12_C24_up, DE_gene.list$P12_C24_down),
                             over_lap_fraction(enriched_genes, DE_gene.list$P12_C36_up, DE_gene.list$P12_C36_down)),
                    Expression = rep(c("Up", "Unchanged", "Down"), 4))
  dat$Expression <- factor(dat$Expression, levels = c("Up", "Unchanged", "Down"))
  
  dat_text <- dat
  dat_text$Num <- dat_text$Freq * length(enriched_genes) / 100
  dat_text$Freq <- sapply(c("P12", "P12C12", "P12C24", "P12C36"),
                          function(x) 
                            dat_text$Freq[dat_text$Time == x] %>% 
                            "/"(2) %>% 
                            "+"(rev(c(0, cumsum(dat_text$Freq[dat_text$Time == x][c(3, 2)]))))
  ) %>% c()
  dat_text$Freq[dat_text$Expression == "Down"] <- dat_text$Freq[dat_text$Expression == "Down"] + 5
  
  ggplot(dat, aes(x = Time, y = Freq, fill = Expression)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    geom_text(data = dat_text, 
              aes(x = Time, y = Freq, label = Num),
              size = 6) +
    xlab("") + ylab(.ylab) +
    scale_fill_manual(values = c("#B22332", "#A9D6EF", "#2364A7")) +
    ggpubr::theme_pubclean() +
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 16),
          legend.text=element_text(size = 14))
}


# with PcG enrichment
plot_target_gene_frequency(PcG_enriched_genes, "Frequency (%)") +
  ggtitle("PcG (Ring1b + Ezh2) enriched genes (n = 1480)")

ggsave(filename = "Fig3_barplot_DE_group_PcG_enriched_genes.png", 
       path = "../fig3/figs", device = "png", width = 6, height = 6)


# ------------------------------------ Compare with ChIP signals ------------------------------- #
# boxplot by response
ChIP_TSS_mat_P0 <- cbind(ChIP_TSS_B1_mat[, grepl("_NT", colnames(ChIP_TSS_B1_mat)) &
                                           !grepl("BAP1|H3|S2|S5", colnames(ChIP_TSS_B1_mat)) &
                                           !grepl("IN|Ezh2i", colnames(ChIP_TSS_B1_mat))],
                         ChIP_TSS_B2_mat[, grepl("(Cbx7|Rybp)_P0", colnames(ChIP_TSS_B2_mat)) ] )
colnames(ChIP_TSS_mat_P0) <- gsub("(\\W|_).*", "", colnames(ChIP_TSS_mat_P0))
colnames(ChIP_TSS_mat_P0) <- gsub("Pol2", "Pol II", colnames(ChIP_TSS_mat_P0))

# append public ChIP data
ChIP_TSS_mat_P0 <- cbind(ChIP_TSS_mat_P0,
                         "BAP1 (Kweon et al.)" = ChIP_Kweon_TSS$`2019_Kweon_WT-BAP1.mm9.bw`,
                         "FOXK1 (Kolovos et al.)" = ChIP_Kolovos_TSS[, "FOXK1"],
                         "FOXK2 (Kolovos et al.)" = ChIP_Kolovos_TSS[, "FOXK2"],
                         "HCFC1 (Conway et al.)" = ChIP_Conway_TSS[, "HCFC1"])

ChIP_TSS_mat_P0 <- `colnames<-`(trim_quantile(ChIP_TSS_mat_P0), 
                                colnames(ChIP_TSS_mat_P0))
rownames(ChIP_TSS_mat_P0) <- names(gene.gr)

ChIP_TSS_mat_P0 <- sweep(ChIP_TSS_mat_P0, 2, size_factor_cal(ChIP_TSS_mat_P0), "/")


response_table <- cbind(gene_id = names(gene.gr),
                        group = "UC")
response_table[response_table[, 1] %in% Reduce(c, DE_gene.list[grep("up", names(DE_gene.list))]), 2] <- "Up"
response_table[response_table[, 1] %in% Reduce(c, DE_gene.list[grep("down", names(DE_gene.list))]), 2] <- "Down"
rownames(response_table) <- names(gene.gr)

response_table_cmb <- cbind(response_table[rownames(ChIP_TSS_B1_mat), ], ChIP_TSS_mat_P0)

g1 <- response_table_cmb[, 2:7] %>% 
  as.data.frame() %>%
  reshape::melt(id.vars = "group") %>%
  dplyr::mutate(value = as.numeric(as.character(value)),
                group = factor(group, levels = c("Up", "UC", "Down"))) %>%
  ggplot(aes(x = group, y = log2(value + 1), fill = group)) +
  geom_boxplot(outlier.color = NA, width = 0.7) +
  facet_grid(.~variable) +
  ggpubr::stat_compare_means(comparisons = list(c("Up", "UC"),
                                                c("Up", "Down")),
                             aes(label = ..p.signif..),
                             method = "t.test", na.rm = T, vjust = 0.5) +
  scale_fill_manual(values = c("#B22332", "#A9D6EF", "#2364A7")) +
  xlab("") + ylab("TSS RPGC (log2)") + ggtitle("") +
  theme_setting +
  theme(legend.position = "none",
        strip.text = element_text(size = 13), 
        plot.margin = margin(0,0,0,0))

g2 <- response_table_cmb[, c(2, 8:12)] %>% 
  as.data.frame() %>%
  reshape::melt(id.vars = "group") %>%
  dplyr::mutate(value = as.numeric(as.character(value)),
                group = factor(group, levels = c("Up", "UC", "Down"))) %>%
  ggplot(aes(x = group, y = log2(value + 1), fill = group)) +
  geom_boxplot(outlier.color = NA, width = 0.7) +
  facet_grid(.~variable) +
  ggpubr::stat_compare_means(comparisons = list(c("Up", "UC"),
                                                c("Up", "Down")),
                             aes(label = ..p.signif..),
                             method = "t.test", na.rm = T, vjust = 0.5) +
  scale_fill_manual(values = c("#B22332", "#A9D6EF", "#2364A7")) +
  xlab("") + ylab("TSS RPGC (log2)") + ggtitle("") +
  theme_setting +
  theme(legend.position = "none",
        strip.text = element_text(size = 13), 
        plot.margin = margin(0,0,0,0))

ggsave(grid.arrange(g1, g2, nrow = 2),
       filename = "Fig3_boxplot_DE_group_ChIP_TSS_density.png", 
       path = "../fig3/figs", device = "png", width = 10, height = 6.5)

# ----------------------------- Down regulation GO  ---------------------------------- #
source("../util/go.R")
enrichGeneSets(genes_derepressed, ontology = "BP")

enrichGeneSets(unique(Reduce(c, DE_gene.list[grep("down", names(DE_gene.list))])),
               ontology = "BP", title = "BAP1 PC down-regulated genes (n=1032)", 
               colorset = rev(RColorBrewer::brewer.pal(9, "Blues")), 
               is.GeneRatio = TRUE, 
               col.limits = c(0, 0.05))

ggsave(filename = "FigS4_GO_BP_BAP1_PC_down_regulated.png",
       path = "../figS4/figs",
       device = "png", width = 10.2, height = 5)

# ----------------------------------- PRC1 types -------------------------------------- #
if (F) {
  PRC1_type_DE <- log2FC_mat[intersect.Vector(c(Ring1b_enriched_genes, 
                                                Cbx7_binding_genes,
                                                Rybp_binding_genes), 
                                              rownames(log2FC_mat)), ] %>% 
    `colnames<-`(c("P12", "P12C12", "P12C24", "P12C36")) %>%
    reshape::melt()
  
  PRC1_type_DE$PRC1 <- ifelse(PRC1_type_DE$X1 %in% intersect.Vector(Ring1b_binding_genes, 
                                                                    intersect.Vector(Cbx7_binding_genes, 
                                                                                     Rybp_binding_genes)), 
                              "111", ifelse(!PRC1_type_DE$X1 %in% c(Cbx7_binding_genes, Rybp_binding_genes),
                                            "100", 
                                            ifelse(!PRC1_type_DE$X1 %in% c(Cbx7_binding_genes, Ring1b_binding_genes),
                                                   "001",
                                                   ifelse(!PRC1_type_DE$X1 %in% c(Rybp_binding_genes, Ring1b_binding_genes),
                                                          "010",
                                                          ifelse(!PRC1_type_DE$X1 %in% Rybp_binding_genes,
                                                                 "110",
                                                                 ifelse(!PRC1_type_DE$X1 %in% Cbx7_binding_genes,
                                                                        "101",
                                                                        "011"
                                                                 )
                                                          )
                                                   )
                                            )
                              )
  )
  PRC1_type_DE$PRC1[PRC1_type_DE$PRC1 %in% c("001", "010", "011")] <- "Cbx7 or Rybp\nn=1037"
  PRC1_type_DE$PRC1[PRC1_type_DE$PRC1 == "100"] <- "Ring1b\nn=1666" 
  PRC1_type_DE$PRC1[PRC1_type_DE$PRC1 == "110"] <- "cPRC1\nn=612" 
  PRC1_type_DE$PRC1[PRC1_type_DE$PRC1 == "101"] <- "vPRC1\nn=580" 
  PRC1_type_DE$PRC1[PRC1_type_DE$PRC1 == "111"] <- "Both\nn=542" 
  PRC1_type_DE$PRC1 <- factor(PRC1_type_DE$PRC1,
                              levels = c("Cbx7 or Rybp\nn=1037", 
                                         "Ring1b\nn=1666",
                                         "vPRC1\nn=580",
                                         "cPRC1\nn=612", 
                                         "Both\nn=542"))
  
  ggplot(PRC1_type_DE, aes(x = PRC1, y = value, fill = X2)) +
    geom_hline(yintercept = 0, lty = 1, lwd = 1, col = "red") +
    geom_boxplot(outlier.size = 0, notch = T) +
    scale_fill_viridis_d(begin = 0.2, end = 0.84) +
    xlab("") + ylab("RNA-seq log2FC") +
    labs(fill = "Time") +
    ggpubr::theme_pubclean() +
    theme(axis.text = element_text(face = "bold", size = 12), 
          axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1), 
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 13))
  
  ggsave(filename = paste0("FigS5_RNA_DE_PRC1_groups.pdf"), 
         path = "../figS5/figs/",
         device = "pdf", width = 5, height = 5)
}

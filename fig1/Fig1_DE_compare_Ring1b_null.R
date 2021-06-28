# Rui Shao 2021 May
# Figure 1
# DE genes comparison

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")


# load previos conditional Ring1b knock-out DE results
source("F1_src_Ring1b_CKO_DE_res.R")

# ----------------------------------------------------------------------------------- #
Hox_genes <- mcols(gene.gr[grep("Hox", gene.gr$gene_name)])
Hox_genes$Type <- gsub("Hox(.).*", '\\1', Hox_genes$gene_name) %>% toupper()
Hox_genes$No <- gsub("Hox.", '', Hox_genes$gene_name) %>% as.numeric()
Hox_genes <- Hox_genes[order(Hox_genes$Type, as.numeric(Hox_genes$No)), ] #%>% as.matrix()

# for (i in LETTERS[1:4]) {
#   idx1 <- which(Hox_genes[, "Type"] == i)
#   idx2 <- Hox_genes[idx1, "No"] %>% as.numeric() %>% diff() #%>% ">"(1) %>% which()
#   idx2 <- rep(which(idx2 > 1), (idx2 - 1)[idx2 > 1])
#   Hox_genes <- rbind(matrix(NA,
#                             nrow = min(as.numeric(Hox_genes[idx1, "No"])) - 1, 
#                             ncol = ncol(Hox_genes)),
#                      .add_row(mat = Hox_genes[idx1, ], 
#                               .after_rows = idx2, 
#                               vals = matrix(NA, 
#                                             nrow = length(idx2), 
#                                             ncol = ncol(Hox_genes))),
#                      Hox_genes[-idx1, ]
#                      )
# }

Hox_genes <- as.data.frame(Hox_genes)
Hox_genes$Type <- factor(Hox_genes$Type, LETTERS[4:1])
Hox_genes$No <- factor(Hox_genes$No, 1:13)

# Hox_genes$log2FC_BAP1 <- rowMeans(log2FC_mat[Hox_genes$gene_id, c("P12_C24", "P12_C36")])
Hox_genes$log2FC_BAP1 <- log2FC_mat[Hox_genes$gene_id, "P12_C24"]
Hox_genes$pval_BAP1 <- res.list$P12_C36[Hox_genes$gene_id, "padj"] %>% log2()
Hox_genes$ME_BAP1 <- res.list$P12_C36[Hox_genes$gene_id, "baseMean"] %>% log2()


Hox_genes$log2FC_Ring1b_Tam <- GSE134053_DE_res[Hox_genes$gene_id, "log2FoldChange"]
Hox_genes$log2FC_Ring1b_Bla <- GSE132753_DESeq2_res[Hox_genes$gene_id, "LFC_PRC1CKO"]

# plot tile
g4.1 <- ggplot(Hox_genes, aes(x = No, y = Type)) +
  geom_point(aes(color = log2FC_BAP1), size = 7) +
  scale_color_viridis() +
  xlab("") + ylab("Hox clusters\n") + ggtitle("BAP1 P12C24") +
  labs(color = "log2FC") +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 12), 
        panel.grid.major.y = element_line(size = 0.4) )


g4.2 <- ggplot(Hox_genes, aes(x = No, y = Type)) +
  geom_point(aes(color = log2FC_Ring1b_Bla), size = 7) +
  scale_color_viridis(limits = c(0, max(Hox_genes$log2FC_Ring1b_Bla) )) +
  xlab("") + ylab("Hox clusters\n") + ggtitle("Ring1b CKO Blackledge et al.") +
  labs(color = "log2FC") +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 12), 
        panel.grid.major.y = element_line(size = 0.4) )


g4.3 <- ggplot(Hox_genes, aes(x = No, y = Type)) +
  geom_point(aes(color = log2FC_Ring1b_Tam), size = 7) +
  scale_color_viridis(limits = c(0, max(Hox_genes$log2FC_Ring1b_Tam) )) +
  xlab("") + ylab("Hox clusters\n") + ggtitle("Ring1b CKO Tamburri et al.") +
  labs(color = "log2FC") +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 12), 
        panel.grid.major.y = element_line(size = 0.4) )

ggsave(plot = grid.arrange(g4.1, g4.2, g4.3, ncol = 3),
       filename = paste0("Fig1_Hox_clusters_DE_log2FC.png"), 
       path = "figs/",
       device = "png", width = 18, height = 2)

# ----------------------------------------------------------------------------------- #
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



GSE132753_DE_genes <- rownames(GSE132753_DESeq2_res)[GSE132753_DESeq2_res$LFC_PRC1CKO > 1.5 &
                                                       GSE132753_DESeq2_res$LFC_PRC1CPM > 1.5 & 
                                                       GSE132753_DESeq2_res$padj_PRC1CKO < 0.05 & 
                                                       GSE132753_DESeq2_res$padj_PRC1CPM < 0.05]

GSE134053_DE_genes <- rownames(GSE134053_DE_res.list$I53S_OHT)[GSE134053_DE_res.list$Ring1b_OHT$log2FoldChange > 1.5 &
                                                                 GSE134053_DE_res.list$Ring1b_OHT$padj < 0.05 &
                                                                 GSE134053_DE_res.list$I53S_OHT$log2FoldChange > 1.5 &
                                                                 GSE134053_DE_res.list$I53S_OHT$padj < 0.05 &
                                                                 GSE134053_DE_res.list$I53S_WT_OHT$log2FoldChange > 1.5 &
                                                                 GSE134053_DE_res.list$I53S_WT_OHT$padj < 0.05 &
                                                                 GSE134053_DE_res.list$Ring1b_WT_OHT$log2FoldChange > 1.5 &
                                                                 GSE134053_DE_res.list$Ring1b_WT_OHT$padj < 0.05
                                                               ]

length(intersect.Vector(GSE132753_DE_genes, GSE134053_DE_genes))

length(intersect.Vector(genes_derepressed, GSE132753_DE_genes))
length(intersect.Vector(genes_derepressed, GSE134053_DE_genes))

length(intersect.Vector(genes_W_derepressed, GSE132753_DE_genes))
length(intersect.Vector(genes_W_derepressed, GSE134053_DE_genes))

# ----------------------------------------------------------------------------------- #

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

# ------------------------------Overlapped genes GO---------------------------------- #
source("../util/go.R")
enrichGeneSets(intersect.Vector(genes_derepressed, GSE132753_DE_genes),
               ontology = "BP")
ggsave(filename = paste0("FigS1_GO_BP_C1_Klose.png"), 
       path = "figs/",
       device = "png", width = 7, height = 6)

enrichGeneSets(intersect.Vector(genes_derepressed, GSE134053_DE_genes),
               ontology = "BP")
ggsave(filename = paste0("FigS1_GO_BP_C1_Pasini.png"), 
       path = "figs/",
       device = "png", width = 7, height = 6)

enrichGeneSets(intersect.Vector(GSE132753_DE_genes, GSE134053_DE_genes),
               ontology = "BP")
ggsave(filename = paste0("FigS1_GO_BP_Klose_Pasini.png"), 
       path = "figs/",
       device = "png", width = 7, height = 6)


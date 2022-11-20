
# --------------------------------------------------------------------------------------------------- #
# explore gene sets by derepression pattern along pulsing time


# --------------------------------------------------------------------------------------------------------- #
# add boxplots (b)
dat_lfc <- as.data.frame(log2FC_mat[H2Aub_repressed_genes_fig4, ])
dat_lfc$Class <- derepressed_gene_tab[H2Aub_repressed_genes_fig4, "derepression_rate"]

dat_lfc %>%
  dplyr::filter(complete.cases(.)) %>%
  reshape::melt() %>%
  ggplot(aes(x = Class, y = value, fill = variable)) +
  geom_hline(yintercept = 0, color = add.alpha("red", 0.7)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name = "Times", values = brewer.pal(4, "Blues")) +
  ylim(c(-3, 7)) + xlab("De-repression groups") + ylab("RNA-seq log2FC") + ggtitle("H2Aub repressed genes (n=1565)") +
  theme_setting

ggsave(filename = paste0("figS5_boxplot_log2FC_time_group_H2Aub_repressed.png"), 
       path = "../figS5/figs/",
       device = "png", width = 5, height = 4)


dat_lfc <- as.data.frame(log2FC_mat[BAP1_specific_genes, ])
dat_lfc$Class <- derepressed_gene_tab[BAP1_specific_genes, "derepression_rate"]

dat_lfc %>%
  dplyr::filter(complete.cases(.)) %>%
  reshape::melt() %>%
  ggplot(aes(x = Class, y = value, fill = variable)) +
  geom_hline(yintercept = 0, color = add.alpha("red", 0.7)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(name = "Times", values = c("grey80", "grey60", "grey40", "grey20")) +
  ylim(c(-3, 7)) + xlab("De-repression groups") + ylab("RNA-seq log2FC") + ggtitle("BAP1 specific genes (n=1327)") +
  theme_setting

ggsave(filename = paste0("figS5_boxplot_log2FC_time_group_BAP1_specific.png"), 
       path = "../figS5/figs/",
       device = "png", width = 5, height = 4)

# derepressed_gene_tab[, c(4:8, 14)] %>% 
#   as.data.frame() %>%
#   reshape::melt(id.vars = "Derepression") %>%
#   dplyr::mutate(value = as.numeric(as.character(value)),
#                 group = factor(Derepression, levels = c("Ring1b CKO", "BAP1 Specific"))) %>%
#   ggplot(aes(x = Derepression, y = log2(value + 1), fill = Derepression)) +
#   geom_boxplot(outlier.color = NA, width = 0.7) +
#   facet_grid(.~variable) +
#   ggpubr::stat_compare_means(comparisons = list(c("Ring1b CKO", "BAP1 Specific")),
#                              aes(label = ..p.signif..),
#                              method = "t.test", na.rm = T, vjust = 0.5) +
#   scale_fill_manual(values = colors_9[c(2,3,1)]) +
#   xlab("") + ylab("TSS ±1 kb Density (log2)") + ggtitle("") +
#   theme_setting +
#   theme(legend.position = "none",
#         strip.text = element_text(size = 13), 
#         plot.margin = margin(0,0,0,0))



# -------------------------------------- test Gene Ontology -------------------------------------- #
source("../util/go.R")
enrichGeneSets(derepressed_gene_tab$gene_id[derepressed_gene_tab$derepression_rate == "Early"], 
               method = "GO", ontology = "BP", is.GeneRatio = T, title = "Early", col.limits = c(0, 0.5)) +
  theme(legend.position = "none")
ggsave(filename = paste0("figS5_test_GO_Early_genes_BP_T.png"), 
       path = "../figS5/figs/",
       device = "png", width = 6.5, height = 5)

enrichGeneSets(derepressed_gene_tab$gene_id[derepressed_gene_tab$derepression_rate == "Middle"], 
               method = "GO", ontology = "BP", is.GeneRatio = T, title = "Middle", col.limits = c(0, 0.5)) +
  theme(legend.position = "none") + ylim(c(0, 0.07))
ggsave(filename = paste0("figS5_test_GO_Middle_genes_BP_T.png"), 
       path = "../figS5/figs/",
       device = "png", width = 6, height = 5)

enrichGeneSets(derepressed_gene_tab$gene_id[derepressed_gene_tab$derepression_rate == "Late"], 
               method = "GO", ontology = "BP", is.GeneRatio = T, title = "Late", col.limits = c(0, 0.5))
ggsave(filename = paste0("figS5_test_GO_Late_genes_BP_T.png"), 
       path = "../figS5/figs/",
       device = "png", width = 7, height = 5)


enriched_GO_list_heatmap(split(names(gene_cls_H2Aub), gene_cls_H2Aub),
                         title = "H2Aub repressed genes (n=1565)",
                         colorset = c("grey50", "cyan3"))
ggsave(filename = paste0("FigS5_GO_H2Aub_repressed_gene.png"), 
       path = "../figS5/figs",
       device = "png", width = 6.4, height = 4.5)


enriched_GO_list_heatmap(split(names(class_genes_rates[BAP1_specific_genes]), class_genes_rates[BAP1_specific_genes]),
                         title = "BAP1 specific genes (n=1327)",
                         colorset = c("grey50", "grey10"))
ggsave(filename = paste0("FigS5_GO_BAP1_specific_gene.png"), 
       path = "../figS5/figs",
       device = "png", width = 8, height = 4.5)

# enrichGeneSets(derepressed_gene_tab$gene_id[derepressed_gene_tab$derepression_rate == "Late" & derepressed_gene_tab$BAP1_specific], 
#                method = "GO", ontology = "BP", is.GeneRatio = T, title = "Late", col.limits = c(0, 0.5))



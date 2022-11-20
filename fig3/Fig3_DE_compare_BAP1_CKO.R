# filter DE log2FC results
log2FC_mat <- sapply(res.gene.list, function(x) x$log2FoldChange)
rownames(log2FC_mat) <- rownames(res.gene.list[[1]])
log2FC_mat <- log2FC_mat[!is.na(rowSums(log2FC_mat)), ] #%>% trim_quantile()


# --------------------------------- BAP1 KO DE ---------------------------------- #
gene_intersect_Conway <- intersect.Vector(rownames(log2FC_mat), Conway_res_table$ENSEMBL)
gene_intersect_Fursova <- intersect.Vector(rownames(log2FC_mat), Fursova_res_table$gene_id)
gene_intersect_Kolovos <- intersect.Vector(rownames(log2FC_mat), Kolovos_res_table$gene_id)


if (F) {
  smoothScatter(Conway_res_table$log2FoldChange[match(gene_intersect_Conway, Conway_res_table$ENSEMBL)],
                log2FC_mat[gene_intersect_Conway, 3])
  
  smoothScatter(Fursova_res_table$log2FoldChange[match(gene_intersect_Fursova, Fursova_res_table$gene_id)],
                log2FC_mat[gene_intersect_Fursova, 3])
  
  smoothScatter(Kolovos_res_table$log2FoldChange[match(gene_intersect_Kolovos, Kolovos_res_table$gene_id)],
                log2FC_mat[gene_intersect_Kolovos, 3])
  
  
  table(derepressed_gene_ids %in% Conway_res_list$up, derepressed_gene_ids %in% Kolovos_res_list$up)
  
  
  # BAP1 depression
  table(Conway_res_list$up %in% derepressed_gene_ids) #    324   165 
  table(Conway_res_list$down %in% derepressed_gene_ids) #   465   105 
  
  table(Fursova_res_list$up %in% derepressed_gene_ids) #   166   127 
  table(Fursova_res_list$down %in% derepressed_gene_ids) # 1185   440 
  
  table(Kolovos_res_list$up %in% derepressed_gene_ids) #   498   473 
  table(Kolovos_res_list$down %in% derepressed_gene_ids) # 1157   378 
  
  table(Conway_res_list$up %in% unique(Reduce(c, DE_gene.list[grep("down", names(DE_gene.list))])))
  table(Conway_res_list$down %in% unique(Reduce(c, DE_gene.list[grep("down", names(DE_gene.list))])))
  
  table(Fursova_res_list$up %in% unique(Reduce(c, DE_gene.list[grep("down", names(DE_gene.list))])))
  table(Fursova_res_list$down %in% unique(Reduce(c, DE_gene.list[grep("down", names(DE_gene.list))])))
  
  table(Kolovos_res_list$up %in% unique(Reduce(c, DE_gene.list[grep("down", names(DE_gene.list))])))
  table(Kolovos_res_list$down %in% unique(Reduce(c, DE_gene.list[grep("down", names(DE_gene.list))])))
  
  
  # Ring1b CKO
  table(Conway_res_list$up %in% Ring1b_CKO_genes) #     321   168 
  table(Conway_res_list$down %in% Ring1b_CKO_genes) #   512    58 
  
  table(Fursova_res_list$up %in% Ring1b_CKO_genes) #    102   191 
  table(Fursova_res_list$down %in% Ring1b_CKO_genes) #  1278   347 
  
  table(Kolovos_res_list$up %in% Ring1b_CKO_genes) #    604   367 
  table(Kolovos_res_list$down %in% Ring1b_CKO_genes) #  1280   255 
  
  # Ring1b CPM
  table(Conway_res_list$up %in% Ring1b_CPM_genes) #     160   110 
  table(Conway_res_list$down %in% Ring1b_CPM_genes) #   323    21 
  
  table(Fursova_res_list$up %in% Ring1b_CPM_genes) #    300   334 
  table(Fursova_res_list$down %in% Ring1b_CPM_genes) #  2553   413
  
  table(Kolovos_res_list$up %in% Ring1b_CPM_genes) #    478   169 
  table(Kolovos_res_list$down %in% Ring1b_CPM_genes) #  900   128 
}


Ring1b_CKO_genes <- unique(c(GSE132753_DE_genes_CKO, GSE134053_DE_genes_CKO)) # 2668
Ring1b_CPM_genes <- unique(c(GSE132753_DE_genes_CPM, GSE134053_DE_genes_CPM)) # 1932

repressed_gene_ids <- unique(Reduce(c, DE_gene.list[grep("down", names(DE_gene.list))])) # 1032

# survey
survey_BAP1 <- data.frame("up" = as.numeric(table(derepressed_gene_ids %in% Conway_res_list$up +
                                                    derepressed_gene_ids %in% Fursova_res_list$up +
                                                    derepressed_gene_ids %in% Kolovos_res_list$up)[-1]),
                          "down" = as.numeric(table(derepressed_gene_ids %in% Conway_res_list$down +
                                                      derepressed_gene_ids %in% Fursova_res_list$down +
                                                      derepressed_gene_ids %in% Kolovos_res_list$down)[-1]))

survey_Ring1b <- data.frame("up" = as.numeric(table(Ring1b_CKO_genes %in% Conway_res_list$up +
                                                      Ring1b_CKO_genes %in% Fursova_res_list$up +
                                                      Ring1b_CKO_genes %in% Kolovos_res_list$up)[-1]),
                            "down" = as.numeric(table(Ring1b_CKO_genes %in% Conway_res_list$down +
                                                        Ring1b_CKO_genes %in% Fursova_res_list$down +
                                                        Ring1b_CKO_genes %in% Kolovos_res_list$down)[-1]))

survey_BAP1_down <- data.frame("up" = c(as.numeric(table(repressed_gene_ids %in% Conway_res_list$up +
                                                           repressed_gene_ids %in% Fursova_res_list$up +
                                                           repressed_gene_ids %in% Kolovos_res_list$up)[-1]), 0),
                               "down" = c(as.numeric(table(repressed_gene_ids %in% Conway_res_list$down +
                                                             repressed_gene_ids %in% Fursova_res_list$down +
                                                             repressed_gene_ids %in% Kolovos_res_list$down)[-1]))
)

survey_BAP1_KO <- data.frame("up" = as.numeric(table(table(c(Conway_res_list$up, Fursova_res_list$up, Kolovos_res_list$up))))[1:3],
                             "down" = as.numeric(table(table(c(Conway_res_list$down, Fursova_res_list$down, Kolovos_res_list$down))))[1:3])


derepressed_gene_tab <- data.frame(gene_id = derepressed_gene_ids,
                                   BAP1_KO_up = derepressed_gene_ids %in% c(Conway_res_list$up, 
                                                                            Fursova_res_list$up,
                                                                            Kolovos_res_list$up),
                                   BAP1_KO_down = derepressed_gene_ids %in% c(Conway_res_list$down, 
                                                                              Fursova_res_list$down,
                                                                              Kolovos_res_list$down))

repressed_gene_tab <- data.frame(gene_id = repressed_gene_ids,
                                 BAP1_KO_up = repressed_gene_ids %in% c(Conway_res_list$up, 
                                                                        Fursova_res_list$up,
                                                                        Kolovos_res_list$up),
                                 BAP1_KO_down = repressed_gene_ids %in% c(Conway_res_list$down, 
                                                                          Fursova_res_list$down,
                                                                          Kolovos_res_list$down))

rbind(cbind(survey_Ring1b %>%
              as.matrix() %>% 
              reshape::melt(), group = "Ring1b CKO Up\nn = 2668"),
      cbind(survey_BAP1 %>%
              as.matrix() %>% 
              reshape::melt(), group = "BAP1 Pulse Up\nn = 2802"),
      cbind(survey_BAP1_down %>%
              as.matrix() %>% 
              reshape::melt() , group = "BAP1 Pulse Down\nn = 1391"),
      cbind(survey_BAP1_KO %>%
              as.matrix() %>% 
              reshape::melt(), group = "BAP1 CKO\nup = 1505\ndown = 2770")) %>%
  
  dplyr::mutate("BAP1 CKO" = factor(X2, c("up", "down"))) %>%
  ggplot(aes(x = X1, y = value, fill = `BAP1 CKO`)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 0.7) +
  facet_grid(group ~ ., scales = "free") +
  ggtitle("Differential Expression Survey") +
  xlab("Number of public datasets") + ylab("Gene overlap") +
  # ylim(c(0, 600)) +
  scale_x_continuous(breaks = 1:3, labels = 1:3) +
  scale_fill_manual(values = c("#A21821", "#796926")) +
  theme_setting +
  theme(strip.text = element_text(size = 13))

ggsave(filename = "FigS4_barplot_BAP1_KO_differential_expression_survey.pdf", 
       path = "../figS4/figs", device = "pdf", width = 7, height = 6)


# 
intersect_two_sets(unique(c(Conway_res_list$up, Fursova_res_list$up, Kolovos_res_list$up)), derepressed_gene_ids)
#  870  635 2544

intersect_two_sets(unique(c(Conway_res_list$down, Fursova_res_list$down, Kolovos_res_list$down)), derepressed_gene_ids)
# 2047  723 2456

intersect_two_sets(Ring1b_CKO_genes, derepressed_gene_ids)
intersect_two_sets(Ring1b_CKO_genes, repressed_gene_ids) # BAP1 induced down-regulation is irrelevant to PRC1





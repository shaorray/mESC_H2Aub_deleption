# --------------------------------------------------------------------------------------------------- #
# get Hox gene family
# gene.gr <- readRDS("../data/gene.gr.RData")
# gene.gr$gene_id <- gsub("\\..*", "", gene.gr$gene_id)
# names(gene.gr) <- gene.gr$gene_id

Hox_genes <- mcols(gene.gr[grepl("Hox", gene.gr$gene_name) & !grepl("Hox.*as$", gene.gr$gene_name)])
Hox_genes$Type <- gsub("Hox(.).*", '\\1', Hox_genes$gene_name) %>% toupper()
Hox_genes$No <- gsub("Hox.", '', Hox_genes$gene_name) %>% as.numeric()
Hox_genes <- Hox_genes[order(Hox_genes$Type, as.numeric(Hox_genes$No)), ] %>% as.data.frame()
Hox_genes$Type <- factor(Hox_genes$Type, LETTERS[4:1])
Hox_genes$No <- factor(Hox_genes$No, 1:13)

Hox_genes$H2Aub_repressed <- ifelse(Hox_genes$gene_id %in% derepressed_gene_ids, "De-repressed", "Unchanged")

Hox_genes$Ring1b <- ifelse(Hox_genes$gene_id %in% rownames(ChIP_TSS_B2_mat_log2FC)[which(ChIP_TSS_B2_mat_log2FC[, "Ring1b_P12"] < -1)],
                           "2x loss", "Stable")


g1 <- ggplot(Hox_genes, aes(x = No, y = Type)) +
  geom_point(aes(color = H2Aub_repressed), size = 5) +
  scale_color_manual(values = c("red3", "grey75")) +
  xlab("") + ylab("Hox clusters\n") + #ggtitle("BAP1 P12C24") +
  labs(color = "Response") +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 12), 
        panel.grid.major = element_blank(), 
        axis.ticks = element_line() )

g2 <- ggplot(Hox_genes, aes(x = No, y = Type)) +
  geom_point(aes(color = Ring1b), size = 5) +
  # scale_color_manual(values = c("red3", "grey75")) +
  xlab("") + ylab("Hox clusters\n") +
  labs(color = "Ring1b") +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 12), 
        panel.grid.major = element_blank(), 
        axis.ticks = element_line() )


ggsave(grid.arrange(g1, g2, nrow = 1), filename = "Fig_S3_Hox_clusters_derepression.png", 
       path = "../figS3/figs", 
       device = "png", width = 12, height = 3)


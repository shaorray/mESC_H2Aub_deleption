# arrange BAP1 log2FC by derepression time
log2FC_mat <- readRDS("data/log2FC_mat.RData")

log2FC_mat_std <- t(log2FC_mat[, 1:3]) %>% scale() %>% t()

# overview  ---------------------------------------------------------------------------------
p <- umap::umap(log2FC_mat)
png("../figS1/figs/FigS1_gene_log2FC_cluster_Umap.png", 
    width = 400, height = 400)
plot(p$layout, 
     xlim=c(-20, 30), 
     ylim=c(-30, 30), 
     col = log2FC_cls,
     pch = 19, cex = 0.1,
     xlab = "UMAP1",
     ylab = "UMAP2", 
     main = "Gene log2FC clustering (n = 16734)")

legend("topright", 
       c("Derepressed", "Weakly derepressed", "Unchanged", "Repressed"), 
       col = 1:4,
       pch = 19, cex = 0.8)
dev.off()

# early -> middle -> late
log2FC_mat2 <-
  log2FC_mat[intersect.Vector(c(genes_derepressed, genes_W_derepressed), use_gene_ids), ]
log2FC_mat2 <- log2FC_mat2[rowMins(log2FC_mat2) > 0, ]

max_up_time <- apply(log2FC_mat2, 1, which.max) # c("Early", "Middle", "Late")

# derepression time
derepression_time <- colSums(t(log2FC_mat2[, 1:3] / rowSums(log2FC_mat2[, 1:3])) * c(12, 24, 36))
derepression_cv <- rowSds(log2FC_mat2) / rowMeans(log2FC_mat2)


# plot heatmap ---------------------------------------------------------------------------------
cbind(log2FC_mat2, cls = max_up_time) %>% 
  as.data.frame() %>%
  arrange(cls, rowMeans(log2FC_mat2[, 1:3])) %>%
  cbind(id = seq_len(nrow(log2FC_mat2))) %>%
  # log2FC_mat[order(log2FC_cls, rowMeans(log2FC_mat[, 3:4])), ] %>%
  `rownames<-`(NULL) %>% 
  reshape::melt(id.vars = c("cls", "id")) %>% 
  ggplot(aes(x = variable, y = id, fill = value)) +
  geom_tile(width = 0.99) +
  facet_grid(cls ~ ., scales = "free", 
             space = "free", switch = "y", 
             labeller = labeller(cls = c("1"="Early", "2"="Middle", "3"="Mid-late", "4" = "Late"))) +
  labs(title = 'Derepressed genes (n = 2902)',  x = '', y = '') +
  scale_fill_gradientn(name = "Log2FC", 
                       colours = c("white", 'red'), 
                       values = (c(0, max(log2FC_mat2)) - min(log2FC_mat2)) /
                         diff(range(log2FC_mat2)) ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.x = element_text(size = 11, hjust = 1, angle = 45),
        axis.text.y = element_blank(), 
        legend.key.size = unit(13, "pt"), 
        plot.margin = unit(c(1, 1, 0, 1), "lines")) 


if (F) {
  library(Ternary)
  TernaryPlot(alab = 'P12', blab = 'P12C12', clab = 'P12C24')
  TernaryPoints(log2FC_mat2,
                col = colors_20[max_up_time],
                pch = 19, cex = 0.5
  )
  TernaryPoints(log2FC_mat2,
                col = colors_20[max_up_time],
                pch = 1, cex = 0.5
  )
  # genes that explain early up-regulation
  TernaryPoints(log2FC_mat2[intersect.Vector(rownames(log2FC_mat2),
                                             use_gene.gr$gene_id[use_gene.gr$PcG == "Target"]), ],
                pch = 19, cex = 0.5
  )
}






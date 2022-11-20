

# combine DE genes
derepressed_gene_ids_fig3 <- Reduce(c, DE_gene.list[grep("up", names(DE_gene.list))]) %>% unique()
derepressed_gene_ids_fig3 <- derepressed_gene_ids_fig3[derepressed_gene_ids_fig3 %ni% CpK_DE.gr$gene_id] # 355 out of 3179 genes are due to amber codon read-through

# extract first derepression time
derepressed_gene_times <- rep(1, length(derepressed_gene_ids_fig3))
for (i in 4:1) {
  idx <- derepressed_gene_ids_fig3 %in% DE_gene.list[grep("up", names(DE_gene.list))][[i]]
  derepressed_gene_times[idx] <- i
}

derepressed_gene_tab <- cbind(gene_id = derepressed_gene_ids_fig3,
                              time = derepressed_gene_times)

# derepressed_gene.gr <- gene.gr[derepressed_gene_ids_fig3]






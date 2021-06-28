library(DESeq2)
# kallisto est read count DE
dds_txRC <- DESeqDataSetFromMatrix(round(as.matrix(txRC_RNA)),
                                   colData = data.frame(condition = gsub("_R.", "", .design)),
                                   design = ~ condition)
sizeFactors(dds_txRC) <- SizeFactorCal(flyRPK_RNA[rowMeans(flyRPK_RNA) > 1, ])
dds_txRC <- DESeq(dds_txRC)

res_P12 <- results(dds_txRC, contrast = c("condition", "P12", "P0"))
res_P12_C12 <- results(dds_txRC, contrast = c("condition", "P12_C12", "P0"))
res_P12_C24 <- results(dds_txRC, contrast = c("condition", "P12_C24", "P0"))
res_P12_C36 <- results(dds_txRC, contrast = c("condition", "P12_C36", "P0"))

res.list <- list(P12 = res_P12,
                 P12_C12 = res_P12_C12,
                 P12_C24 = res_P12_C24,
                 P12_C36 = res_P12_C36)
saveRDS(res.list, "data/RNAseq_res.list.RData")

# list of DE genes -------------------------------------------------------------------------
# kallisto transcripts
DE_gene.list <- list("P12_up" = rownames(res_P12)[res_P12$log2FoldChange > 1.5 & res_P12$padj < 0.05],
                     "P12_down" = rownames(res_P12)[res_P12$log2FoldChange < (-1.5) & res_P12$padj < 0.05],
                     "P12_C12_up" = rownames(res_P12_C12)[res_P12_C12$log2FoldChange > 1.5 & res_P12_C12$padj < 0.05],
                     "P12_C12_down" = rownames(res_P12_C12)[res_P12_C12$log2FoldChange < (-1.5) & res_P12_C12$padj < 0.05],
                     "P12_C24_up" = rownames(res_P12_C24)[res_P12_C24$log2FoldChange > 1.5 & res_P12_C24$padj < 0.05],
                     "P12_C24_down" = rownames(res_P12_C24)[res_P12_C24$log2FoldChange < (-1.5) & res_P12_C24$padj < 0.05],
                     "P12_C36_up" = rownames(res_P12_C36)[res_P12_C36$log2FoldChange > 1.5 & res_P12_C36$padj < 0.05],
                     "P12_C36_down" = rownames(res_P12_C36)[res_P12_C36$log2FoldChange < (-1.5) & res_P12_C36$padj < 0.05])

DE_gene.list <- lapply(DE_gene.list, function(x) x[!is.na(x)])

saveRDS(DE_gene.list, "data/DE_gene.list.RData")


#  -------------------------------------------------------------------------
if (F) {
        # bam -> featureCounts gencode read count DE
        dds_gencode <- DESeqDataSetFromMatrix(round(as.matrix(count_table_rnaseq.norm)),
                                              colData = data.frame(condition = gsub("_R.", "", .design)),
                                              design = ~ condition)
        sizeFactors(dds_gencode) <- 1
        dds_gencode <- DESeq(dds_gencode)
        
        res_P12 <- results(dds_gencode, contrast = c("condition", "P12", "P0"))
        res_P12_C12 <- results(dds_gencode, contrast = c("condition", "P12_C12", "P0"))
        res_P12_C24 <- results(dds_gencode, contrast = c("condition", "P12_C24", "P0"))
        res_P12_C36 <- results(dds_gencode, contrast = c("condition", "P12_C36", "P0"))
        
        res_gencode.list <- list(P12 = res_P12,
                                 P12_C12 = res_P12_C12,
                                 P12_C24 = res_P12_C24,
                                 P12_C36 = res_P12_C36)
        saveRDS(res_gencode.list, "data/RNAseq.res_gencode.list.RData")
        
        # list of DE genes -------------------------------------------------------------------------
        # gencode exons
        DE_gencode.list <- with(res_gencode.list,
                                list("P12_up" = rownames(res_P12)[res_P12$log2FoldChange > 1.5 & res_P12$padj < 0.05],
                                     "P12_down" = rownames(res_P12)[res_P12$log2FoldChange < (-1.5) & res_P12$padj < 0.05],
                                     "P12_C12_up" = rownames(res_P12_C12)[res_P12_C12$log2FoldChange > 1.5 & res_P12_C12$padj < 0.05],
                                     "P12_C12_down" = rownames(res_P12_C12)[res_P12_C12$log2FoldChange < (-1.5) & res_P12_C12$padj < 0.05],
                                     "P12_C24_up" = rownames(res_P12_C24)[res_P12_C24$log2FoldChange > 1.5 & res_P12_C24$padj < 0.05],
                                     "P12_C24_down" = rownames(res_P12_C24)[res_P12_C24$log2FoldChange < (-1.5) & res_P12_C24$padj < 0.05],
                                     "P12_C36_up" = rownames(res_P12_C36)[res_P12_C36$log2FoldChange > 1.5 & res_P12_C36$padj < 0.05],
                                     "P12_C36_down" = rownames(res_P12_C36)[res_P12_C36$log2FoldChange < (-1.5) & res_P12_C36$padj < 0.05]))
        
        DE_gencode.list <- lapply(DE_gencode.list, function(x) x[!is.na(x)])
        
        saveRDS(DE_gencode.list, "data/DE_gencode.list.RData")
}

library(ggpubr)
res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                       keys = as.character(rownames(res_P12)),
                       keytype = "GENEID",
                       columns = "GENENAME")

g1 <- ggpubr::ggmaplot(res_P12[res$GENEID, ], genenames = res$GENENAME, 
                       top = 50, select.top.method = "fc", ggtheme = theme_pubclean(), 
                       fc = 1.5, fdr =  0.05, size = 0.5, ylim = c(-15, 15), 
                       main = "P12 vs P0")

g2 <- ggpubr::ggmaplot(res_P12_C12[res$GENEID, ], genenames = res$GENENAME, 
                       top = 50, select.top.method = "fc", ggtheme = theme_pubclean(),
                       fc = 1.5, fdr =  0.05, size = 0.5, ylim = c(-15, 15), 
                       main = "P12_C12 vs P0")

g3 <- ggpubr::ggmaplot(res_P12_C24[res$GENEID, ], genenames = res$GENENAME, 
                       top = 50, select.top.method = "fc", ggtheme = theme_pubclean(),
                       fc = 1.5, fdr =  0.05, size = 0.5, ylim = c(-15, 15), 
                       main = "P12_C24 vs P0")

g4 <- ggpubr::ggmaplot(res_P12_C36[res$GENEID, ], genenames = res$GENENAME, 
                       top = 50, select.top.method = "fc", ggtheme = theme_pubclean(),
                       fc = 1.5, fdr =  0.05, size = 0.5, ylim = c(-15, 15), 
                       main = "P12_C36 vs P0")

ggsave(grid.arrange(g1, g2, g3, g4, nrow = 1),
       filename = "RNAseq_gencode_MAplot_BAP1_pulse_chase.png",
       width = 16, height = 5, path = "figs")


# RNA-seq vs TT-seq ------------------------------------------------------------------------

plot(x = rowDiffs(log2(count_table_tx.norm[, 1:2])), 
     y = res.list$P12[rownames(count_table_tx.norm), "log2FoldChange"], 
     cex = 0.1, ylim = c(-15, 15), 
     ylab = "RNA-seq log2FC P12 vs P0",
     xlab = "TT-seq log2FC P12 vs P0")

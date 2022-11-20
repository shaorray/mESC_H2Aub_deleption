plot_scatter_n <- function(X, Y, gene_ids, target_ids = NA,
                           .xlab = "", .ylab = "", 
                           xlim = NULL, ylim = NULL, 
                           n.cluster = NA, rm.mid.cluster = T) {
  dat <- data.frame(x = X, y = Y) %>% 
    `rownames<-`(., gene_ids) %>%
    dplyr::filter(complete.cases(.) & is.finite(rowSums(.)) & !is.na(rowSums(.)))
  if (!is.null(xlim) & !is.null(ylim)) 
    dat <- dplyr::filter(dat, x > min(xlim) & x < max(xlim) & y > min(ylim) & y < max(ylim))
  # r <- single_variance_explained(dat[, 1], dat[, 2], is.cor = T) %>% round(3)
  
  g <- ggplot(dat, aes(x = x, y = y, color = get_dens(x, y))) +
    geom_point(cex = 0.5) +
    annotate("text", x = -Inf, y = Inf,
             hjust = -0.5, vjust = 1.2,
             label = paste0("n = ", nrow(dat))) +
    scale_x_continuous(name = .xlab, limits = xlim) +
    scale_y_continuous(name = .ylab, limits = ylim) +
    scale_color_gradientn(colours = c("grey80", "grey50", "grey30", "grey15", "grey10", "black")) +
    theme_setting +
    theme(legend.position = "none")
  
  if (!is.na(n.cluster)) {
    require(mclust)
    .Group <<- mclust::Mclust(dat[, 1:2], G = 4)$classification
    dat$Group <- .Group
    if (rm.mid.cluster) {
      group_coord <- data.frame(x = aggregate(dat$x, list(dat$Group), "median")$x,
                                y = aggregate(dat$y, list(dat$Group), "median")$x)
      group_border <- c(which.min(group_coord$x), which.min(group_coord$y), 
                        which.max(group_coord$x), which.max(group_coord$x))
      dat2 <<- dat[dat$Group %in% unique(group_border), ]
    } else {
      dat2 <- dat
    }
    g <- g + stat_ellipse(data = dat2, aes(x = x, y = y, group = Group))
  }
  
  if (!is.na(target_ids[1])) {
    dat_highlight <- dat[rownames(dat) %in% target_ids, ] %>% 
      dplyr::filter(complete.cases(.) &
                      x > (cut_baseMean) & 
                      !is.infinite(x) & 
                      !is.infinite(y))
    
    # add fraction distribution
    Group_precent <- 
      (as.matrix(table(rownames(dat2) %in% target_ids, dat2$Group))[2, ] /
         length(target_ids) * 100) %>% round(2) %>% paste0(., "%")
    
    Group_pos <- aggregate(dat2[, 1:2], list(dat2$Group), median) %>% data.frame()
    Group_pos$y <- aggregate(dat2[, 1:2], list(dat2$Group), min)[, "y"]
    
    g + geom_point(data = dat_highlight, aes(x = x, y = y),
                   size = 0.2, pch = 19, color = add.alpha("red", 0.5)) +
      geom_text(data = Group_pos, 
                aes(x = x, y = y, label = Group_precent), 
                color = "black", size = 5)
    
  } else {
    g
  }
}


# ------------------------------- load annotations -----------------------------------#
ChromHMM <- importRanges("/mnt/0E471D453D8EE463/genomeDir/chromHMM/mESC_E14_12_dense.annotated.mm10.bed")
ChromHMM$name <- gsub(".*_", "\\2", ChromHMM$name)

gene.mm10.gr <- importRanges("/mnt/0E471D453D8EE463/genomeDir/GENCODE/gencode.vM25.annotation.gtf")
gene.mm10.gr <- gene.mm10.gr[gene.mm10.gr$gene_type == "protein_coding" & gene.mm10.gr$type == "gene"]
gene.mm10.gr$gene_id <- gsub("\\..*", "", gene.mm10.gr$gene_id)
names(gene.mm10.gr) <- gene.mm10.gr$gene_id 
gene.mm10.gr <- gene.mm10.gr[!is.na(names(gene.mm10.gr))]
gene.mm10.gr$chromHMM <- ChromHMM$name[findOverlaps(promoters(gene.mm10.gr, upstream = 1, downstream = 0),
                                                    ChromHMM) %>% subjectHits()]

Bivalent_genes_chromHMM <- gene.mm10.gr$gene_id[gene.mm10.gr$chromHMM == "BivalentChromatin"]
Bivalent_genes_chromHMM_dr <- intersect.Vector(Bivalent_genes_chromHMM, derepressed_gene_ids)
Bivalent_genes_chromHMM_nc <- Bivalent_genes_chromHMM[Bivalent_genes_chromHMM %ni% derepressed_gene_ids]

Bivalent_genes_chromHMM_tss <- gene.mm10.gr$gene_id[findOverlaps(promoters(gene.mm10.gr),
                                                                 ChromHMM[ChromHMM$name == "BivalentChromatin"]) %>% 
                                                      countQueryHits() > 0]

# gene group intersections: H2Aub repressed, bivalent and PcG binding
intersect_three_sets(H2Aub_repressed_genes, Bivalent_genes_chromHMM_tss, PcG_binding_genes)
# 100  110  101  111  010  011  001 
# 253  278  223 1324 2065 1570 2478

if (F) {
  # DNA replication timing has limited connection with derepression
  Rep_phase.gr <- importRanges("/mnt/0E471D453D8EE463/genomeDir/Replication/CC_ChIP_EdU_Replication_continuous_impute.bed")
  
  mtch <- findOverlaps(promoters(gene.gr, upstream = 1, downstream = 0),
                       Rep_phase.gr)
  gene.gr$Rep_phase <- NA
  gene.gr$Rep_phase[queryHits(mtch)] <- Rep_phase.gr$name[subjectHits(mtch)]
}

# ------------------------------- chromatin state clustering ------------------------------------- #
derepressed_gene_ids <- Reduce(c, DE_gene.list[grep("up", names(DE_gene.list))]) %>% unique() %>%
  intersect.Vector(rownames(ChIP_TSS_B2_mat_cmb))


H3K4me3_tss <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2016_Chronis_ESC_H3K4me3_ChIP-Seq.mm9.bw", 
                        intervals = promoters(gene.gr[derepressed_gene_ids],
                                              upstream = 1000, downstream = 1000), fast = F)
H3K27me3_tss <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2016_Chronis_ESC_H3K27me3_ChIP-Seq.mm9.bw", 
                         intervals = promoters(gene.gr[derepressed_gene_ids],
                                               upstream = 1000, downstream = 1000), fast = F)
H2AZ_tss <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/2016_Surface_H2AZ_WT_ChIP_Seq.mm9.bw", 
                     intervals = promoters(gene.gr,
                                           upstream = 1000, downstream = 1000), fast = F)
H33_tss <- .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/GSM2582412_ESC_H3.3WT_YFP.bigwig", 
                    intervals = promoters(gene.gr,
                                          upstream = 1000, downstream = 1000), fast = F)
histone_dat <- cbind(H3K4me3 = log1p(H3K4me3_tss),
                     H3K27me3 = log1p(H3K27me3_tss),
                     H2AZ = log1p(H2AZ_tss),
                     H33 = log1p(H33_tss))


# [s] make a scatter plot of H2AZ and H3K27me3
plot_scatter_n(X = log2(H3K27me3_tss[, 1] + 1), Y = log2(H2AZ_tss[, 1] + 1), 
               gene_ids = gene.gr$gene_id, 
               target_ids = H2Aub_repressed_genes,
               .xlab = "log2 H3K27me3 (Chronis et al.)", 
               .ylab = "log2 H2A.Z  (Surface et al.)",
               xlim = c(0, 5), ylim = c(0, 7), 
               n.cluster = 4) + ggtitle("Derepressed Genes (n=1460) ~ TSS ±1 kb Density")

H3K27me3_H2AZ_tss_group <- .Group

ggsave(filename = "FigS4_scatter_Chronis_H2AZ_H3K27me3_derepressed_gene.png", 
       path = "../figS4/figs", device = "png", width = 5.5, height = 5)

# [s] make a scatter plot of H3K4me3 and H3K27me3
plot_scatter_n(X = log2(H3K27me3_tss[, 1] + 1), Y = log2(H3K4me3_tss[, 1] + 1), 
               gene_ids = gene.gr$gene_id,
               target_ids = H2Aub_repressed_genes,
               .xlab = "log2 H3K27me3 (Chronis et al.)", 
               .ylab = "log2 H3K4me3 (Chronis et al.)",
               xlim = c(0, 5), ylim = c(0, 8), 
               n.cluster = 4) + ggtitle("Derepressed Genes (n=1460) ~ TSS ±1 kb Density")

H3K27me3_H3K4me3_tss_group <- .Group

ggsave(filename = "FigS4_scatter_Chronis_H3K4me3_H3K27me3_derepressed_gene.png", 
       path = "../figS4/figs", device = "png", width = 5.5, height = 5)

# # H3K27me3 ~ H33 scatter plot 
# plot_scatter_n(X = log2(H3K27me3_tss[, 1] + 1), Y = log2(H33_tss[, 1] + 1), 
#              gene_ids = gene.gr$gene_id, target_ids = derepressed_gene_ids,
#              .xlab = "log2 H3K27me3 (Chronis et al.)", 
#              .ylab = "log2 H3.3  (Chen et al.)",
#              xlim = c(0, 5), ylim = c(0, 3), 
#              n.cluster = 4) + ggtitle("Derepressed Genes (n=1460) ~ TSS ±1kb Density")
# ggsave(filename = "FigS1_scatter_Chronis_H3.3_H3K27me3_derepressed_gene.pdf", 
#        path = "../figS1/figs", device = "pdf", width = 5.5, height = 5)

# # H2AZ ~ H33 scatter plot 
plot_scatter_n(X = log2(H2AZ_tss[, 1] + 1), Y = log2(H33_tss[, 1] + 1),
             gene_ids = gene.gr$gene_id, target_ids = derepressed_gene_ids,
             .xlab = "log2 H2A.Z  (Surface et al.)",
             .ylab = "log2 H3.3  (Chen et al.)",
             xlim = c(0, 7), ylim = c(0, 4),
             n.cluster = 3) + ggtitle("Derepressed Genes (n=2802) ~ TSS ±1kb density")
# 
# ggsave(filename = "FigS1_scatter_Chronis_H3.3_H2A.Z_derepressed_gene.pdf", 
#        path = "../figS1/figs", device = "pdf", width = 5.5, height = 5)

# H3K4me3 ~ H33 scatter plot
plot_scatter_n(X = log2(H3K4me3_tss[, 1] + 1), Y = log2(H33_tss[, 1] + 1),
               gene_ids = gene.gr$gene_id, target_ids = derepressed_gene_ids,
               .xlab = "log2 H3K4me3  (Chronis et al.)",
               .ylab = "log2 H3.3  (Chen et al.)",
               xlim = c(0, 7), ylim = c(0, 4),
               n.cluster = 3) + ggtitle("Derepressed Genes (n=2802) ~ TSS ±1kb density")

intersect_three_sets(H2Aub_repressed_genes,
                     names(H3K27me3_H3K4me3_tss_group[H3K27me3_H3K4me3_tss_group == 2]),
                     names(H3K27me3_H2AZ_tss_group[H3K27me3_H2AZ_tss_group == 2]))
# 100  110  101  111  010  011  001 
# 398  421   31 1228 3044 1899  189 

# group specification
# 1: not H2Aub repressed + bivalent
# 2: not H2Aub repressed + H3K4me3/H3K27me3

gene_group_1 <- Bivalent_genes_chromHMM_tss[!Bivalent_genes_chromHMM_tss %in% H2Aub_repressed_genes]

H3K27me3_H3K4me3_cluster <- names(H3K27me3_H3K4me3_tss_group)[H3K27me3_H3K4me3_tss_group == 2]
gene_group_2 <- H3K27me3_H3K4me3_cluster[!H3K27me3_H3K4me3_cluster %in% H2Aub_repressed_genes]

test_features <- gene_epi_features2[H3K27me3_H3K4me3_cluster, ]

test_features$idx <- 0
test_features$idx[rownames(test_features) %in% gene_group_2] <- 1
test_features <- test_features[complete.cases(test_features) & is.finite(rowSums(test_features)), ]
test_features$idx <- as.factor(test_features$idx)

library(caret)
library(randomForest)
library(e1071)
library(pROC)
library(gbm)

set.seed(1)
control <- trainControl(method = 'repeatedcv', 
                        number = 10, 
                        repeats = 1)

tunegrid <- expand.grid(.mtry = sqrt(ncol(test_features)))

inTraining <- createDataPartition(test_features$idx, p = .75, list = FALSE)
training <- test_features[ inTraining, ]
testing  <- test_features[-inTraining, ]

fit_rf_cls <- train(idx ~ ., 
                    data = training, 
                    method = 'rf', 
                    metric = 'Accuracy',
                    tuneGrid = tunegrid, 
                    trControl = control)
print(fit_rf_cls)

rf_cls_Imp <- varImp(fit_rf_cls, scale = FALSE)
plot(rf_cls_Imp)


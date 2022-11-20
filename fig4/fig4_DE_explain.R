# Rui Shao 2021 May
# Figure 4
# BAP1 pulse derepression explanation

source("F4_scr_get_features_density.r")

library(caret)
library(randomForest)
library(e1071)
library(pROC)
library(gbm)

# -------------------------------------- explain DE with MINUTE-ChIP data -------------------------------------- # 

# keep up-regulated genes
gene_intersect <- intersect.Vector(rownames(log2FC_mat)[rowMeans(log2FC_mat) > 0],
                                   intersect.Vector(rownames(ChIP_TSS_B1_mat),
                                                    rownames(ChIP_TSS_B2_mat)))
# remove BAP1 specific genes
gene_intersect <- gene_intersect[gene_intersect %ni% derepressed_gene_ids[derepressed_gene_ids %ni% H2Aub_repressed_genes]]

# use gene responses as proxies to represent H2Aub repression intensity 
repression_idx <- log2(pmax(0, rowMeans(geneRC_RNA[, 9:12]) - rowMeans(geneRC_RNA[, 1:3])) / rowMeans(geneRC_RNA[, 1:3]))

ChIP_features <- log2(ChIP_TSS_B1_mat[gene_intersect, grepl("_NT", colnames(ChIP_TSS_B1_mat)) & !grepl("IN_", colnames(ChIP_TSS_B1_mat))] + 1)
colnames(ChIP_features) <- gsub("_NT", "", colnames(ChIP_features))

ChIP_features2 <- log2(ChIP_TSS_B2_mat[gene_intersect, grep("(Cbx7|Rybp)_P0", colnames(ChIP_TSS_B2_mat))] + 1)
colnames(ChIP_features2) <- gsub("_P0", "", colnames(ChIP_features2))
ChIP_features <- cbind(ChIP_features, ChIP_features2, 
                       "RNA_RPK" = log2(geneRPK_RNA[gene_intersect, c("P32", "P33", "P41")] %>% rowMeans() + 1))
rm(ChIP_features2)

# linear explaination
R2_linear <- data.frame(Target = colnames(ChIP_features),
                           R2 = multi_variance_explained(ChIP_features[gene_intersect, ],
                                                         repression_idx[gene_intersect])
                        )
R2_linear$Target <- factor(R2_linear$Target, levels = R2_linear$Target[order(R2_linear$R2)])

# combine the response variable
ChIP_features$Y <- repression_idx[gene_intersect]
ChIP_features <- ChIP_features[complete.cases(ChIP_features) & is.finite(rowSums(ChIP_features)), ]

# train the same RF model
set.seed(1)
control <- trainControl(method = 'repeatedcv', 
                        number = 10, 
                        repeats = 1)

tunegrid <- expand.grid(.mtry = sqrt(ncol(ChIP_features)))

inTraining <- createDataPartition(ChIP_features$Y, p = .75, list = FALSE)
training <- ChIP_features[ inTraining, ]
testing  <- ChIP_features[-inTraining, ]

# linear model
fit <- lm(Y~., training)
summary(fit) # Adjusted R-squared:  0.364 
smoothScatter(testing$Y, predict(fit, testing))
cor(testing$Y, predict(fit, testing)) # 0.6353233

g1 <- ggplot(R2_linear, aes(x = Target, y = R2)) +
  geom_bar(stat = 'identity', width = 0.7, fill = "darkmagenta") +
  xlab("") + ylab(expression(R^2)) +
  coord_flip() +
  ggpubr::theme_pubr() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12))


g2 <- data.frame(log2FC = testing$Y, pred = predict(fit, testing)) %>%
  ggplot(aes(x = log2FC, y = pred, color = get_dens(log2FC, pred))) +
  geom_point() +
  annotate(geom = 'text', x = -Inf, y = Inf, 
           hjust = -0.5, vjust = 1, 
           label = paste0("r = ", round(cor(testing$Y, predict(fit, testing)), 3),
                          "\nn = ", nrow(testing))) +
  scale_color_viridis_c(end = 0.8) +
  xlab("Repression index test set") + ylab("Linear Model Predicted log2FC") + 
  ggpubr::theme_pubr() +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

# random forest
fit_rf_cls <- train(Y ~ ., 
                    data = training, 
                    method = 'rf', 
                    # metric = 'Accuracy',
                    tuneGrid = tunegrid, 
                    trControl = control)
print(fit_rf_cls) # Rsquared 0.3799298

rf_cls_Imp <- varImp(fit_rf_cls, scale = FALSE)
plot(rf_cls_Imp)

smoothScatter(testing$Y, predict(fit_rf_cls, testing))

# train with gbm 
fitControl <- trainControl(method = "cv", number = 10)

gbmGrid <-  expand.grid(interaction.depth = 9, 
                         n.trees = 50, 
                         shrinkage = 0.1,
                         n.minobsinnode = 80)

gbmFit <- train(Y ~ .,
                data = training, 
                method = "gbm", 
                trControl = fitControl, 
                tuneGrid = gbmGrid,
                verbose = F)
gbmFit; summary(gbmFit, las = 2) # Rsquared 0.3804834


gbmImp <- varImp(gbmFit, scale = FALSE)

gbm.Test <- predict(gbmFit, 
                    newdata = testing,
                    type="raw")

g3 <- data.frame(Imp = gbmImp$importance$Overall,
                 Target = gsub("`", "", rownames(gbmImp$importance))) %>%
  dplyr::mutate(Target = factor(Target, levels = Target[order(Imp)])) %>%
  ggplot(aes(x = Target, y = Imp)) +
  geom_bar(stat = 'identity', width = 0.7, fill = "gold4") +
  xlab("") + ylab("Importance") +
  coord_flip() +
  ggpubr::theme_pubr() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12))


g4 <- data.frame(log2FC = testing$Y, pred = gbm.Test) %>%
  ggplot(aes(x = log2FC, y = pred, color = get_dens(log2FC, pred))) +
  geom_point() +
  annotate(geom = 'text', x = -Inf, y = Inf, 
           hjust = -0.5, vjust = 1, 
           label = paste0("r = ", round(cor(testing$Y, gbm.Test), 3),
                          "\nn = ", nrow(testing))) +
  scale_color_viridis_c(end = 0.8) +
  xlab("Repression index test set") + ylab("GBM Predicted log2FC") + 
  ggpubr::theme_pubr() +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12))

ggsave(grid.arrange(g2, g1, g4, g3, nrow = 1),
       filename = "FigS6_ChIP_predict_RNA_log2FC_R2_Imp.png",
       path = "../figS6/figs", 
       height = 4, width = 16)

# ------------------------------------ Evaluate on certain gene groups ------------------------------------ #
# bivalent
R2_linear_bi <- data.frame(Target = colnames(ChIP_features)[-13], 
                           Gene_set = "Bivalent",
                           R2 = multi_variance_explained(ChIP_features[intersect.Vector(gene_intersect, bivalent_genes_chromHMM), -13], 
                                                         repression_idx[intersect.Vector(gene_intersect, bivalent_genes_chromHMM)]),
                           Cor_2 = multi_variance_explained(ChIP_features[intersect.Vector(gene_intersect, bivalent_genes_chromHMM), -13], 
                                                          rowMeans(log2FC_mat[intersect.Vector(gene_intersect, bivalent_genes_chromHMM), ]), is.cor = T),
                           Cor = multi_variance_explained(ChIP_features[intersect.Vector(gene_intersect, bivalent_genes_chromHMM), -13], 
                                                          repression_idx[intersect.Vector(gene_intersect, bivalent_genes_chromHMM)], is.cor = T)
                           )
R2_linear_bi$Target <- factor(R2_linear_bi$Target, levels = R2_linear_bi$Target[order(R2_linear$R2)])


R2_linear_PcG <- data.frame(Target = colnames(ChIP_features)[-13],
                            Gene_set = "PcG",
                            R2 = multi_variance_explained(ChIP_features[intersect.Vector(gene_intersect, PcG_binding_genes), -13], 
                                                          repression_idx[intersect.Vector(gene_intersect, PcG_binding_genes)]),
                            Cor_2 = multi_variance_explained(ChIP_features[intersect.Vector(gene_intersect, PcG_binding_genes), -13], 
                                                             rowMeans(log2FC_mat[intersect.Vector(gene_intersect, PcG_binding_genes), ]), is.cor = T),
                            Cor = multi_variance_explained(ChIP_features[intersect.Vector(gene_intersect, PcG_binding_genes), -13], 
                                                           repression_idx[intersect.Vector(gene_intersect, PcG_binding_genes)], is.cor = T)
                            )
R2_linear_PcG$Target <- factor(R2_linear_PcG$Target, levels = R2_linear_bi$Target[order(R2_linear$R2)])

R2_linear_cmb <- rbind(R2_linear_bi, R2_linear_PcG)


g1 <- ggplot(R2_linear_cmb, aes(x = Target, y = R2, group = Gene_set, fill = Gene_set)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  coord_flip() +
  xlab("") + ylab(expression(paste(R^2, " explaining repression index"))) +
  coord_flip() +
  scale_fill_manual(values = colors_9[6:7]) +
  ggpubr::theme_pubr() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12))

g2 <- ggplot(R2_linear_cmb, aes(x = Target, y = Cor, group = Gene_set, fill = Gene_set)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  coord_flip() +
  xlab("") + ylab("Correlation with repression index") +
  coord_flip() +
  scale_fill_manual(values = colors_9[6:7]) +
  ggpubr::theme_pubr() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12))

 
ggsave(grid.arrange(g1, g2, nrow = 1),
       filename = "Fig4_Bivalent_PcG_ChIP_RNA_log2FC_R2_Cor.png",
       path = "../fig4/figs", 
       height = 4, width = 8)




# bivalent gene derepressed with more H2Aub?

ChIP_features_bivalent <- ChIP_features[rownames(ChIP_features) %in% bivalent_genes_chromHMM, 
                                        c("Ezh2", "Ring1b", "H2Aub", "H3K27me3", "H3K4me3", "Pol II-NTD")]
ChIP_features_bivalent$Derepression <- ifelse(rownames(ChIP_features_bivalent) %in% H2Aub_repressed_genes,
                                              "True", "False")

reshape::melt(ChIP_features_bivalent, id.vars = "Derepression") %>%
  dplyr::mutate(variable = factor(variable, 
                                  levels = c("Ring1b", "Ezh2", "H2Aub", "H3K27me3", "H3K4me3", "Pol II-NTD")),
                Derepression = factor(Derepression, levels = c("True", "False"))) %>%
  ggplot(aes(x = Derepression, y = value, fill = Derepression)) +
  geom_boxplot(outlier.size = 0) +
  ggpubr::stat_compare_means(method = "t.test") +
  facet_grid(.~variable) +
  scale_fill_viridis_d(begin = 0.4, end = 0.9, direction = -1) +
  xlab("BAP1_PC de-repression") + ylab(expression(paste("log"[2], "TSS RPGC"))) + ggtitle("Bivalent genes (n=2685)") +
  theme_setting +
  theme(legend.position = "none",
        strip.text = element_text(size = 13))

ggsave(filename = paste0("FigS6_boxplot_bivalent_derepression_TSS_ChIP.png"), 
       path = "../figS6/figs/",
       device = "png", width = 12, height = 4.1) 

# --------------------------------------- predict gene response groups ------------------------------------- #
if (FALSE) {
  # with public data
  
  # gene_epi_features <- fetch_TU_feature(gene.gr[use_gene_ids], T)
  gene_epi_features2 <- fetch_TU_feature(gene.gr, is.fast = FALSE)
  gene_epi_features2[is.na(gene_epi_features2)] <- 0
  gene_epi_features2[gene_epi_features2 < 0] <- 0
  rownames(gene_epi_features2) <- gene.gr$gene_id
  
  # calculate feature importance with random forest
  gene_epi_features_cp <- gene_epi_features2
  gene_epi_features_cp$idx <- factor(log2FC_cls[rownames(gene_epi_features_cp)])
  gene_epi_features_cp[, -1] <- gene_epi_features_cp[, -1] %>% trim_quantile() %>% log1p()
  
  gene_epi_features_cp <- gene_epi_features_cp[complete.cases(gene_epi_features_cp), ]
  
  
  set.seed(1)
  control <- trainControl(method = 'repeatedcv', 
                          number = 10, 
                          repeats = 1)
  
  tunegrid <- expand.grid(.mtry = sqrt(ncol(gene_epi_features_cp)))
  
  inTraining <- createDataPartition(gene_epi_features_cp$idx, p = .75, list = FALSE)
  training <- gene_epi_features_cp[ inTraining, ]
  testing  <- gene_epi_features_cp[-inTraining, ]
  
  fit_rf_cls <- train(idx ~ ., 
                      data = training, 
                      method = 'rf', 
                      metric = 'Accuracy',
                      tuneGrid = tunegrid, 
                      trControl = control)
  print(fit_rf_cls)
  
  rf_cls_Imp <- varImp(fit_rf_cls, scale = FALSE)
  plot(rf_cls_Imp)
  
  
  # plot ROC
  gbm.probs <- predict(fit_rf_cls,
                       newdata = testing,
                       type="prob")
  plot(roc(testing[, "idx"],
           gbm.probs[, "1"]), col = 1, main = "Accuracy 0.554599")
  plot(roc(testing[, "idx"],
           gbm.probs[, "2"]), col = 2, add = T)
  plot(roc(testing[, "idx"],
           gbm.probs[, "3"]), col = 3, add = T)
  plot(roc(testing[, "idx"],
           gbm.probs[, "4"]), col = 4, add = T)
}
# ------------------------------------ De-repression time class ------------------------------------- #
# from 'Fig4_cluster_RNA_seq_log2FC.R'
gene_cls_H2Aub


gene_time_epi_features <- gene_epi_features2[names(gene_cls_H2Aub), ]
colnames(gene_time_epi_features)
gene_time_epi_features$derepression_time <- factor(gene_cls_H2Aub)

gbmFit_time <- train(derepression_time ~ .,
                     data = gene_time_epi_features, 
                     method = "gbm", 
                     trControl = fitControl, 
                     tuneGrid = gbmGrid,
                     verbose = F)
# Accuracy 0.7131181
gbmImp_time <- varImp(gbmFit_time, scale = FALSE)
plot(gbmImp_time)

c("m6A", "H33_YFP", "H3K36me3", "Med1", "CTCF", "cMyc", "Chd2")


# ----------------------------------- chromatin state explain --------------------------------------- #
TU.DE.gr <- readRDS("../data/TT_seq_TU_DE_mm10_gr.RData")

TU.DE.intergenic <- TU.DE.gr[TU.DE.gr$location == "intergenic"]

# add enhancer annotation
F5_enhancer <- importRanges("/mnt/0E471D453D8EE463/genomeDir/FANTOM5/GRCm38/F5.mm10.enhancers.bed")

TU.intergenic_as <- promoters(TU.DE.intergenic, upstream = 1000, downstream = 0)
levels(strand(TU.intergenic_as)) <- c("-", "+", "*")

mtch_dir <- findOverlaps(TU.DE.intergenic, TU.intergenic_as)
TU.DE.intergenic$direction <- ifelse(countQueryHits(mtch_dir) > 0,
                                     "Bidirectional", "Unidirectional" )
TU.DE.intergenic$enhancer <- ifelse(findOverlaps(F5_enhancer, TU.DE.intergenic) %>% countSubjectHits > 0,
                                    "Enhancer", "TX")
TU.DE.intergenic$enhancer_direction <- paste(TU.DE.intergenic$direction, TU.DE.intergenic$enhancer, sep = "_")

# chromHMM
TU.DE.intergenic$chromHMM <- ChromHMM$name[findOverlaps(promoters(TU.DE.intergenic, upstream = 1, downstream = 0),
                                                        ChromHMM) %>% subjectHits()]

TU.DE.coding <- TU.DE.gr[TU.DE.gr$location == "protein_coding", ]
TU.DE.coding$chromHMM <- ChromHMM$name[findOverlaps(promoters(TU.DE.coding, upstream = 1, downstream = 0),
                                                    ChromHMM) %>% subjectHits()]
# find coding TUs' neighbors

# gene.mm10.gr <- GenomicFeatures::genes(TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene)
# res_txdb <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
#                             keys = as.character(gene.mm10.gr$gene_id),
#                             keytype = "ENTREZID",
#                             columns = "GENEID")
# gene.mm10.gr$gene_id <- res_txdb$GENEID[match(gene.mm10.gr$gene_id, res_txdb$ENTREZID)]
# names(gene.mm10.gr) <- gene.mm10.gr$gene_id 
# gene.mm10.gr <- gene.mm10.gr[!is.na(names(gene.mm10.gr))]

mtch <- findOverlaps(gene.mm10.gr, TU.DE.coding)
TU.DE.coding$gene_id <- NA
TU.DE.coding$gene_id[subjectHits(mtch)] <- gene.mm10.gr$gene_id[queryHits(mtch)]
TU.DE.coding$log2FC_cls <- log2FC_cls[match(TU.DE.coding$gene_id, names(log2FC_cls))]
TU.DE.coding <- TU.DE.coding[!is.na(TU.DE.coding$log2FC_cls)]
TU.DE.coding_ext <- TU.DE.coding + 5e4


idx.DE.intergenic <- TU.DE.intergenic$log2FoldChange_rel_12h > 1 & TU.DE.intergenic$padj_rel_12h < 0.05
idx.DE.coding <- TU.DE.coding$log2FoldChange_rel_12h > 1 & TU.DE.coding$padj_rel_12h < 0.05

idx.DE.gene <- gene.mm10.gr$gene_id %in% derepressed_gene_ids
idx.DE.Ring1b_CKO <- gene.mm10.gr$gene_id %in% unique(c(GSE132753_DE_genes_CKO, GSE134053_DE_genes_CKO))


# test neighboring intergenic TU occurrence
mtch <- findOverlaps(TU.DE.intergenic, TU.DE.coding_ext[TU.DE.coding$log2FC_cls == 1])
plot(table(countSubjectHits(mtch)) / sum(countSubjectHits(mtch)), xlim = c(0, 20), type = 'l')
TU.DE.intergenic$chromHMM[countQueryHits(mtch) > 0 & idx.DE.intergenic] %>% table() %>% sort() %>%"/"(sum(countQueryHits(mtch) > 0 & idx.DE.intergenic, na.rm = T))

mtch <- findOverlaps(TU.DE.intergenic, TU.DE.coding_ext[TU.DE.coding$log2FC_cls == 2])
points(table(countSubjectHits(mtch)) / sum(countSubjectHits(mtch)), xlim = c(0, 20), type = 'l')
TU.DE.intergenic$chromHMM[countQueryHits(mtch) > 0 & idx.DE.intergenic] %>% table() %>% sort() %>%"/"(sum(countQueryHits(mtch) > 0 & idx.DE.intergenic, na.rm = T))

# plot intergenic TU chromHMM type
g1 <- data.frame(
  Type = rep(rownames(table(gene.mm10.gr$chromHMM)), 2),
  Freq = c(
    table(gene.mm10.gr$chromHMM, idx.DE.gene & !idx.DE.Ring1b_CKO)[, 2] /
      c(table(gene.mm10.gr$chromHMM)),
    table(gene.mm10.gr$chromHMM, idx.DE.gene & idx.DE.Ring1b_CKO)[, 2] /
      c(table(gene.mm10.gr$chromHMM))
  ),
  Group = rep(c("BAP1_specific", "Ring1b_CKO"), each = 11)
) %>%
  dplyr::mutate(Type = factor(Type, 
                              levels = rownames(table(gene.mm10.gr$chromHMM))[order(table(gene.mm10.gr$chromHMM, idx.DE.gene)[, 2] /
                                                    rowSums(table(gene.mm10.gr$chromHMM, idx.DE.gene)))])) %>%
  dplyr::mutate(Group = factor(Group, levels = c("Ring1b_CKO", "BAP1_specific"))) %>%
  ggplot(aes(x = Type, y = Freq, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(name = "", values = c("grey70", "lightblue")) +
  ylab("Frequency (Derepressed / All) (2780 / 21846)") + xlab("") +
  ggtitle("Coding TSS (BAP1_specific | Ring1b_CKO)") +
  ggpubr::theme_pubclean() +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12))

g2 <- table(gene.mm10.gr$chromHMM[idx.DE.gene], c("BAP1_specific", "Ring1b_CKO")[idx.DE.Ring1b_CKO[idx.DE.gene] + 1]) %>% 
  as.data.frame() %>%
  dplyr::mutate(Var1 = factor(Var1, 
                              levels = rownames(table(gene.mm10.gr$chromHMM))[order(table(gene.mm10.gr$chromHMM, idx.DE.gene)[, 2] /
                                                                                      rowSums(table(gene.mm10.gr$chromHMM, idx.DE.gene)))])) %>%
  dplyr::mutate(Var2 = factor(Var2, levels = c("Ring1b_CKO", "BAP1_specific"))) %>%
  ggplot(aes(x = Var1, y = Freq, group = Var2, fill = Var2)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(name = "", values = c("grey70", "lightblue")) +
  ylab("Number of Genes") + xlab("") +
  ggtitle("Coding TSS (BAP1_specific | Ring1b_CKO)") +
  ggpubr::theme_pubclean() +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12))

g3 <- data.frame(Type = rownames(table(TU.DE.intergenic$chromHMM)),
           Freq = table(TU.DE.intergenic$chromHMM, idx.DE.intergenic)[, 2] / 
             rowSums(table(TU.DE.intergenic$chromHMM, idx.DE.intergenic))) %>%
  dplyr::mutate(Type = factor(Type, levels = rownames(table(gene.mm10.gr$chromHMM))[order(table(gene.mm10.gr$chromHMM, idx.DE.gene)[, 2] /
                                                                                            rowSums(table(gene.mm10.gr$chromHMM, idx.DE.gene)))])) %>%
  ggplot(aes(x = Type, y = Freq)) +
  geom_bar(stat = "identity", width = 0.7, fill = "lightcyan3") +
  coord_flip() +
  ylab("Frequency (Derepressed / All) (1023 / 70667)") + xlab("") +
  ggtitle("Intergenic TSS") +
  ggpubr::theme_pubclean() +
  theme(legend.position = 'none',
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12))

ggsave(grid.arrange(g1, g2, g3, nrow = 1),
       filename = "FigS5_ChromHMM_Freq_TU_Derepression.pdf",
       path = "../figS5/figs", 
       height = 4, width = 18)


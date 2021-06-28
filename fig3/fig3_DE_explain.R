# Rui Shao 2021 May
# Figure 3
# BAP1 pulse derepression explanation

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
source("../util/getCoverage.R")
source("F3_scr_get_features_density.r")

# gene_epi_features <- fetch_TU_feature(gene.gr[use_gene_ids], T)
gene_epi_features2 <- fetch_TU_feature(gene.gr[use_gene_ids], F)
gene_epi_features2[is.na(gene_epi_features2)] <- 0
gene_epi_features2[gene_epi_features2 < 0] <- 0

gene_epi_features_cp <- gene_epi_features2
gene_epi_features_cp$idx <- factor(log2FC_cls[use_gene_ids])
gene_epi_features_cp[, -1] <- gene_epi_features_cp[, -1] %>% trim_quantile() %>% log1p()

gene_epi_features_cp <- gene_epi_features_cp[complete.cases(gene_epi_features_cp), ]

# ---------------calculate feature importance with random forest------------------ #
library(caret)
library(randomForest)
library(e1071)
library(pROC)

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



# --------------------------------------------------------------------------
hit_mat_90_2 <- hit_mat_90[, 1:119] %>% as.matrix() 
hit_mat_90_2$idx <- factor(log2FC_cls[use_gene_ids])

inTraining <- createDataPartition(hit_mat_90_2$idx, p = .75, list = FALSE)
training <- hit_mat_90_2[ inTraining, ]
testing  <- hit_mat_90_2[-inTraining, ]

fit_rf_cls <- train(idx ~ ., 
                    data = training, 
                    method = 'rf', 
                    metric = 'Accuracy',
                    tuneGrid = tunegrid, 
                    trControl = control)
print(fit_rf_cls)

rf_cls_Imp <- varImp(fit_rf_cls, scale = FALSE)
plot(rf_cls_Imp)



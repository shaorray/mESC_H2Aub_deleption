
source("../util/getCoverage.R")

# ------------------------------------- Processed bw files ----------------------------------------- #

gene.gr$gene_id <- gsub("\\..*", "", gene.gr$gene_id)
blacklist.gr <- importRanges("../data/peaks/blacklist_peak_mm9.bed")

# bigwig files
if (TRUE) {
  # TSS occupancy values 
  # batch 1, after input-bin normalization
  bw_files_1 <-
    list.files(path = "../data/bw/batch1/input_norm",
               pattern = ".*bw$",
               full.names = TRUE)
  
  # bw_files_1 <-
  #   list.files(path = "/mnt/E0767589767560E8/UPPMAX/PHILIP_LOF_20201122/bw",
  #              pattern = ".*bw$",
  #              full.names = TRUE)
  
  
  ChIP_TSS_B1_mat <- .countBW(bw_files  = bw_files_1,
                              intervals = promoters(gene.gr, upstream = 1000, downstream = 1000),
                              blacklist = blacklist.gr,
                              fast = FALSE ) 
  
  rownames(ChIP_TSS_B1_mat) <- gene.gr$gene_id
  colnames(ChIP_TSS_B1_mat) <- gsub(".*(input_norm/)(.*).bw", "\\2", bw_files_1)
  colnames(ChIP_TSS_B1_mat) <- gsub("Pol2", "Pol II", gsub("m3", "me3", colnames(ChIP_TSS_B1_mat)))
  
  
  # batch2, after input-background normalization
  bw_files_2 <- list.files(path = "../data/bw/batch2/input_norm",
                           pattern = "pool.bw$", full.names = TRUE)
  
  ChIP_TSS_B2_mat <- .countBW(bw_files = bw_files_2,
                              intervals = promoters(gene.gr, upstream = 1000, downstream = 1000),
                              blacklist = blacklist.gr, 
                              fast = FALSE)
  
  colnames(ChIP_TSS_B2_mat) <- gsub(".*input_norm\\/(.*)_pool.bw", "\\1", bw_files_2)
  rownames(ChIP_TSS_B2_mat) <- gene.gr$gene_id
  
  saveRDS(list("ChIP_TSS_B1_mat" = ChIP_TSS_B1_mat, 
               "ChIP_TSS_B2_mat" = ChIP_TSS_B2_mat),
          "../data/ChIP_TSS_input_bin_norm_batch1_batch2_list.RData")
  
  
  # ChIP occupancy fold-changes 
  # Batch 1
  ChIP_TSS_B1_mat_log2FC <- NULL
  target_names <- unique(gsub("(_.*)", "\\2", colnames(ChIP_TSS_B1_mat)))
  for (i in target_names) {
    tmp <- ChIP_TSS_B1_mat[, grep(paste0("^", i), colnames(ChIP_TSS_B1_mat))]
    idx <- grep(paste0(i, "_NT"), colnames(tmp))
    tmp <- as.matrix(tmp)
    ChIP_TSS_B1_mat_log2FC <- cbind(ChIP_TSS_B1_mat_log2FC, log2(tmp[, -idx] / tmp[, idx]))
  }
  rownames(ChIP_TSS_B1_mat_log2FC) <- rownames(ChIP_TSS_B1_mat)
  colnames(ChIP_TSS_B1_mat_log2FC) <- gsub("m3", "me3", colnames(ChIP_TSS_B1_mat_log2FC))
  colnames(ChIP_TSS_B1_mat_log2FC) <- gsub("Pol2", "Pol II", colnames(ChIP_TSS_B1_mat_log2FC))
  
  # Batch 2
  ChIP_TSS_B2_mat_log2FC <- NULL
  target_names <- unique(gsub("(_.*)", "\\2", colnames(ChIP_TSS_B2_mat)))
  for (i in target_names) {
    tmp <- ChIP_TSS_B2_mat[, grep(paste0("^", i), colnames(ChIP_TSS_B2_mat))]
    idx <- grep(paste0(i, "_P0"), colnames(tmp))
    tmp <- as.matrix(tmp)
    ChIP_TSS_B2_mat_log2FC <- cbind(ChIP_TSS_B2_mat_log2FC, log2(tmp[, -idx] / tmp[, idx]))
  }
  rownames(ChIP_TSS_B2_mat_log2FC) <- rownames(ChIP_TSS_B2_mat)
}

# ------------------------------------------- RPGC scaling factors ---------------------------------------------- #
# To allow comparison between different studies / batches / runs, 
# RPGC can be a universal scale beyond a particular MINUTE-ChIP. 
# Since bw files above were from input-bin normalized reads coverage, 
# only an average value of inputs is needed for RPGC conversion.

mm9_seqlengths <- seqlengths(getBSgenome("mm9"))
mm9_seqlengths <- mm9_seqlengths[names(mm9_seqlengths) %in% paste0("chr", c(1:19, c("X", "Y")))]

# batch 1
input_bam_files <- list.files("/mnt/0E471D453D8EE463/BAP1_MINUTE_ChIP/batch1/bam", "IN.*RC.*fltd.bam$", full.names = TRUE)
b1_RPGC_scale <- mean(unlist(sapply(input_bam_files, function(x) Rsamtools::countBam(x)[, "nucleotides"])) / sum(mm9_seqlengths))

ChIP_TSS_B1_mat <- ChIP_TSS_B1_mat / b1_RPGC_scale


# batch 2
input_bam_files <- list.files("/mnt/0E471D453D8EE463/BAP1_MINUTE_ChIP/batch2/bam", "sfINPUT_.*pooled.UCSC.bam$", full.names = TRUE)
b2_RPGC_scale <- mean(unlist(sapply(input_bam_files, function(x) Rsamtools::countBam(x)[, "nucleotides"])) / sum(mm9_seqlengths))

ChIP_TSS_B2_mat <- ChIP_TSS_B2_mat / b2_RPGC_scale


# ------------------------------------ Processed bam files ------------------------------------- #
# First batch
# bam input normalization
if (FALSE) {
  # bam input normalization
  bam_files_1 <-
    list.files(path = "/mnt/E0767589767560E8/UPPMAX/PHILIP_LOF_20201122/bam",
               pattern = ".*RC.mm9.fltd.bam$",
               full.names = T)
  
  ChIP_TSS_B1_BAM_mat <- .countBam(bam_files = bam_files_1,
                                   intervals = promoters(gene.gr, upstream = 2000, downstream = 2000),
                                   stranded = F, 
                                   paired.end = 'ignore')
  col_names <- gsub(".*bam/(.*)_RC.mm9.fltd.bam", "\\1", bam_files_1)
  col_names[grep("H3K4m3", col_names)] <- gsub("H3K4m3", "H3K4me3", col_names[grep("H3K4m3", col_names)])
  col_names[grep("Pol2-NTD", col_names)] <- gsub("Pol2-NTD", "Pol II", col_names[grep("Pol2-NTD", col_names)])
  col_names[grep("Pol2-S5p", col_names)] <- gsub("Pol2-S5p", "Pol II-S5p", col_names[grep("Pol2-S5p", col_names)])
  col_names[grep("Pol2-S2p", col_names)] <- gsub("Pol2-S2p", "Pol II-S2p", col_names[grep("Pol2-S2p", col_names)])
  col_names[grep("H3K27m3", col_names)] <- gsub("H3K27m3", "H3K27me3", col_names[grep("H3K27m3", col_names)])
  
  colnames(ChIP_TSS_B1_BAM_mat) <- col_names
  
  ChIP_TSS_B1_BAM_IN <- ChIP_TSS_B1_BAM_mat[, grep("^IN_", colnames(ChIP_TSS_B1_BAM_mat))[c(8,3,1,5,6,7,4,2,9,10)]] # reorder samples
  ChIP_TSS_B1_BAM_IN <- apply(ChIP_TSS_B1_BAM_IN, 2, function(x) ifelse(x == 0, median(x), x))  # fill zeros with medians
  
  ChIP_TSS_B1_BAM_mat_norm <- NA
  for (i in unique(gsub("(_.*)", "\\2", colnames(ChIP_TSS_B1_BAM_mat))) ) {
    ChIP_TSS_B1_BAM_mat_norm <- cbind(ChIP_TSS_B1_BAM_mat_norm,
                                      sweep(ChIP_TSS_B1_BAM_mat[, grep(paste0("^", i, "_"), 
                                                                       colnames(ChIP_TSS_B1_BAM_mat))[c(8,3,1,5,6,7,4,2,9,10)]],
                                            2, size_factor_cal(ChIP_TSS_B1_BAM_IN), "/")
    )
  }
  ChIP_TSS_B1_BAM_mat_norm[, c(1, grep("^IN_", colnames(ChIP_TSS_B1_BAM_mat_norm)))] <- NULL
  
  sf_P0_B1 <- size_factor_cal(ChIP_TSS_B1_BAM_mat[, grepl("NT", colnames(ChIP_TSS_B1_BAM_mat)) &
                                                    !grepl("IN", colnames(ChIP_TSS_B1_BAM_mat))])
  for (i in names(sf_P0_B1)) {
    idx <- grep(gsub("_NT", "", i), colnames(ChIP_TSS_B1_BAM_mat_norm))
    ChIP_TSS_B1_BAM_mat_norm[, idx] <- ChIP_TSS_B1_BAM_mat_norm[, idx] / sf_P0_B1[i]
  }
  
  # transform to log2FC
  ChIP_TSS_B1_BAM_mat_log2FC <- NULL
  for (i in unique(gsub("(_.*)", "\\2", colnames(ChIP_TSS_B1_BAM_mat_norm))) ) {
    idx <- grep(paste0("^", i, "_"), colnames(ChIP_TSS_B1_BAM_mat_norm))
    tmp <- trim_quantile(ChIP_TSS_B1_BAM_mat_norm[, idx])
    colnames(tmp) <- colnames(ChIP_TSS_B1_BAM_mat_norm[, idx])
    tmp <- log2(tmp / tmp[, 1])
    ChIP_TSS_B1_BAM_mat_log2FC <- cbind(ChIP_TSS_B1_BAM_mat_log2FC, tmp)
  }
  rownames(ChIP_TSS_B1_BAM_mat_log2FC) <- rownames(ChIP_TSS_B1_BAM_mat_norm)
  
  # additional comparison with E1 as reference
  ChIP_TSS_B1_BAM_mat_log2FC_E1 <- NULL
  for (i in unique(gsub("(_.*)", "\\2", colnames(ChIP_TSS_B1_BAM_mat_norm))) ) {
    idx <- grep(paste0("^", i, "_.*Ezh2i"), colnames(ChIP_TSS_B1_BAM_mat_norm))
    tmp <- trim_quantile(ChIP_TSS_B1_BAM_mat_norm[, idx])
    colnames(tmp) <- colnames(ChIP_TSS_B1_BAM_mat_norm[, idx])
    tmp <- log2(tmp / tmp[, 1])
    ChIP_TSS_B1_BAM_mat_log2FC_E1 <- cbind(ChIP_TSS_B1_BAM_mat_log2FC_E1, tmp)
  }
  rownames(ChIP_TSS_B1_BAM_mat_log2FC_E1) <- rownames(ChIP_TSS_B1_BAM_mat_norm)
  
  # additional comparison with B12 as reference
  ChIP_TSS_B1_BAM_mat_log2FC_B12 <- NULL
  for (i in unique(gsub("(_.*)", "\\2", colnames(ChIP_TSS_B1_BAM_mat_norm))) ) {
    idx <- grep(paste0("^", i, "_BAP1-12"), colnames(ChIP_TSS_B1_BAM_mat_norm))
    tmp <- trim_quantile(ChIP_TSS_B1_BAM_mat_norm[, idx])
    colnames(tmp) <- colnames(ChIP_TSS_B1_BAM_mat_norm[, idx])
    tmp <- log2(tmp / tmp[, 1])
    ChIP_TSS_B1_BAM_mat_log2FC_B12 <- cbind(ChIP_TSS_B1_BAM_mat_log2FC_B12, tmp)
  }
  rownames(ChIP_TSS_B1_BAM_mat_log2FC_B12) <- rownames(ChIP_TSS_B1_BAM_mat_norm)
  
}


# entire gene occupancy
if (FALSE) {
  ChIP_GENE_B1_BAM_mat <- .countBam(bam_files = bam_files_1,
                                    intervals = gene.gr,
                                    stranded = F, 
                                    paired.end = 'ignore')

  colnames(ChIP_GENE_B1_BAM_mat) <- col_names
  
  ChIP_GENE_B1_BAM_IN <- ChIP_GENE_B1_BAM_mat[, grep("^IN_", colnames(ChIP_GENE_B1_BAM_mat))[c(8,3,1,4,2,5,6,7,9,10)]] 
  ChIP_GENE_B1_BAM_IN <- apply(ChIP_GENE_B1_BAM_IN, 2, function(x) ifelse(x == 0, median(x), x)) # fill zeros with medians
  
  ChIP_GENE_B1_BAM_mat_norm <- NA
  for (i in unique(gsub("(_.*)", "\\2", colnames(ChIP_GENE_B1_BAM_mat))) ) {
    ChIP_GENE_B1_BAM_mat_norm <- cbind(ChIP_GENE_B1_BAM_mat_norm,
                                      sweep(ChIP_GENE_B1_BAM_mat[, grep(paste0("^", i, "_"), colnames(ChIP_GENE_B1_BAM_mat))[c(8,3,1,4,2,5,6,7,9,10)]],
                                            2, size_factor_cal(ChIP_GENE_B1_BAM_IN), "/")
    )
  }
  ChIP_GENE_B1_BAM_mat_norm[, c(1, grep("^IN_", colnames(ChIP_GENE_B1_BAM_mat_norm)))] <- NULL
  # ChIP_GENE_B1_BAM_mat_norm <- ChIP_GENE_B1_BAM_mat_norm[!apply(ChIP_GENE_B1_BAM_mat_norm, 1, function(x) any(x == 0)), ]
  
  sf_P0_B1_2 <- size_factor_cal(ChIP_GENE_B1_BAM_mat[, grepl("NT", colnames(ChIP_GENE_B1_BAM_mat)) &
                                                    !grepl("IN", colnames(ChIP_GENE_B1_BAM_mat))])
  for (i in names(sf_P0_B1_2)) {
    idx <- grep(gsub("NT", "", i), colnames(ChIP_GENE_B1_BAM_mat_norm))
    ChIP_GENE_B1_BAM_mat_norm[, idx] <- ChIP_GENE_B1_BAM_mat_norm[, idx] / sf_P0_B1_2[i]
  }
  
  # transform to log2FC
  ChIP_GENE_B1_BAM_mat_log2FC <- NULL
  for (i in unique(gsub("(_.*)", "\\2", colnames(ChIP_GENE_B1_BAM_mat_norm))) ) {
    idx <- grep(paste0("^", i, "_"), colnames(ChIP_GENE_B1_BAM_mat_norm))
    tmp <- trim_quantile(ChIP_GENE_B1_BAM_mat_norm[, idx])
    colnames(tmp) <- colnames(ChIP_GENE_B1_BAM_mat_norm[, idx])
    tmp <- log2(tmp / tmp[, 1])
    ChIP_GENE_B1_BAM_mat_log2FC <- cbind(ChIP_GENE_B1_BAM_mat_log2FC, tmp)
  }
  rownames(ChIP_GENE_B1_BAM_mat_log2FC) <- rownames(ChIP_GENE_B1_BAM_mat_norm)
}

# --------------------------------------------- Processed bam files --------------------------------------------- #
# Second batch
# ChIP manual normalization
# bam input normalization
# TSS
if (FALSE) {
  bam_files_2 <-
    list.files(path = "/mnt/E0767589767560E8/UPPMAX/PHILIP_LOF3_20210205R/bam",
               pattern = "pooled.mm9.bam$",
               full.names = T)
  
  ChIP_TSS_B2_BAM_mat <- .countBam(bam_files = bam_files_2,
                                   intervals = promoters(gene.gr, upstream = 2000, downstream = 2000),
                                   stranded = F, 
                                   paired.end = 'ignore')
  col_names <- gsub(".*bam/(.*)_pooled.mm9.bam", "\\1", bam_files_2)
  col_names[grep("Pol2", col_names)] <- gsub("Pol2", "Pol II", col_names[grep("Pol2", col_names)])
  colnames(ChIP_TSS_B2_BAM_mat) <- col_names
  
  ChIP_TSS_B2_BAM_IN <- ChIP_TSS_B2_BAM_mat[, grep("INPUT_", colnames(ChIP_TSS_B2_BAM_mat))] # fill zeros with medians
  ChIP_TSS_B2_BAM_IN <- apply(ChIP_TSS_B2_BAM_IN, 2, function(x) ifelse(x == 0, median(x), x))
  
  ChIP_TSS_B2_BAM_mat_norm <- NA
  for (i in unique(gsub("(_.*)", "\\2", colnames(ChIP_TSS_B2_BAM_mat))) ) {
    ChIP_TSS_B2_BAM_mat_norm <- cbind(ChIP_TSS_B2_BAM_mat_norm,
                                      sweep(ChIP_TSS_B2_BAM_mat[, grep(paste0("^", i, "_"), colnames(ChIP_TSS_B2_BAM_mat))],
                                            2, size_factor_cal(ChIP_TSS_B2_BAM_IN), "/")
    )
  }
  ChIP_TSS_B2_BAM_mat_norm[, c(1, grep("INPUT_", colnames(ChIP_TSS_B2_BAM_mat_norm)))] <- NULL
  # ChIP_TSS_B2_BAM_mat_norm <- ChIP_TSS_B2_BAM_mat_norm[!apply(ChIP_TSS_B2_BAM_mat_norm, 1, function(x) any(x == 0)), ]
  
  sf_P0_B2 <- size_factor_cal(ChIP_TSS_B2_BAM_mat[, grepl("P0", colnames(ChIP_TSS_B2_BAM_mat)) &
                                                    !grepl("INPUT", colnames(ChIP_TSS_B2_BAM_mat))])
  for (i in names(sf_P0_B2)) {
    idx <- grep(gsub("_P0", "", i), colnames(ChIP_TSS_B2_BAM_mat_norm))
    ChIP_TSS_B2_BAM_mat_norm[, idx] <- ChIP_TSS_B2_BAM_mat_norm[, idx] / sf_P0_B2[i]
  }
  
  # transform to log2FC
  ChIP_TSS_B2_BAM_mat_log2FC <- NULL
  for (i in unique(gsub("(_.*)", "\\2", colnames(ChIP_TSS_B2_BAM_mat_norm))) ) {
    idx <- grep(paste0("^", i, "_"), colnames(ChIP_TSS_B2_BAM_mat_norm))
    tmp <- trim_quantile(ChIP_TSS_B2_BAM_mat_norm[, idx])
    colnames(tmp) <- colnames(ChIP_TSS_B2_BAM_mat_norm[, idx])
    tmp <- log2(tmp / tmp[, 1])
    ChIP_TSS_B2_BAM_mat_log2FC <- cbind(ChIP_TSS_B2_BAM_mat_log2FC, tmp)
  }
  rownames(ChIP_TSS_B2_BAM_mat_log2FC) <- rownames(ChIP_TSS_B2_BAM_mat_norm)
}


# entire gene occupancy
if (FALSE) {
  ChIP_GENE_B2_BAM_mat <- .countBam(bam_files = bam_files_2,
                                    intervals = gene.gr,
                                    stranded = F, 
                                    paired.end = 'ignore')
  
  colnames(ChIP_GENE_B2_BAM_mat) <- col_names
  
  ChIP_GENE_B2_BAM_IN <- ChIP_GENE_B2_BAM_mat[, grep("INPUT_", colnames(ChIP_GENE_B2_BAM_mat))] 
  ChIP_GENE_B2_BAM_IN <- apply(ChIP_GENE_B2_BAM_IN, 2, function(x) ifelse(x == 0, median(x), x)) # fill zeros with medians

  ChIP_GENE_B2_BAM_mat_norm <- NA
  for (i in unique(gsub("(_.*)", "\\2", colnames(ChIP_GENE_B2_BAM_mat))) ) {
    ChIP_GENE_B2_BAM_mat_norm <- cbind(ChIP_GENE_B2_BAM_mat_norm,
                                      sweep(ChIP_GENE_B2_BAM_mat[, grep(paste0("^", i, "_"), colnames(ChIP_GENE_B2_BAM_mat))],
                                            2, size_factor_cal(ChIP_GENE_B2_BAM_IN), "/")
    )
  }
  ChIP_GENE_B2_BAM_mat_norm[, c(1, grep("INPUT_", colnames(ChIP_GENE_B2_BAM_mat_norm)))] <- NULL
  # ChIP_GENE_B2_BAM_mat_norm <- ChIP_GENE_B2_BAM_mat_norm[!apply(ChIP_GENE_B2_BAM_mat_norm, 1, function(x) any(x == 0)), ]
  
  sf_P0_B2_2 <- size_factor_cal(ChIP_GENE_B2_BAM_mat[, grepl("P0", colnames(ChIP_GENE_B2_BAM_mat)) &
                                                       !grepl("INPUT", colnames(ChIP_GENE_B2_BAM_mat))])
  for (i in names(sf_P0_B2_2)) {
    idx <- grep(gsub("_P0", "", i), colnames(ChIP_GENE_B2_BAM_mat_norm))
    ChIP_GENE_B2_BAM_mat_norm[, idx] <- ChIP_GENE_B2_BAM_mat_norm[, idx] / sf_P0_B2_2[i]
  }
  
  # transform to log2FC
  ChIP_GENE_B2_BAM_mat_log2FC <- NULL
  for (i in unique(gsub("(_.*)", "\\2", colnames(ChIP_GENE_B2_BAM_mat_norm))) ) {
    idx <- grep(paste0("^", i, "_"), colnames(ChIP_GENE_B2_BAM_mat_norm))
    tmp <- trim_quantile(ChIP_GENE_B2_BAM_mat_norm[, idx])
    colnames(tmp) <- colnames(ChIP_GENE_B2_BAM_mat_norm[, idx])
    tmp <- log2(tmp / tmp[, 1])
    ChIP_GENE_B2_BAM_mat_log2FC <- cbind(ChIP_GENE_B2_BAM_mat_log2FC, tmp)
  }
  rownames(ChIP_GENE_B2_BAM_mat_log2FC) <- rownames(ChIP_GENE_B2_BAM_mat_norm)
}

# additional reads
if (FALSE) {
  bam_files_2_1 <-
    list.files(path = "/mnt/0E471D453D8EE463/baseSpace/batch2/final/bam",
               pattern = "pooled.UCSC.bam$",
               full.names = T)
  
  ChIP_TSS_B2_1_BAM_mat <- .countBam(bam_files = bam_files_2_1,
                                   intervals = promoters(gene.gr, upstream = 2000, downstream = 2000),
                                   stranded = F, 
                                   paired.end = 'ignore')
  col_names <- gsub(".*bam/(.*)_pooled.UCSC.bam", "\\1", bam_files_2_1)
  col_names[grep("Pol2", col_names)] <- gsub("Pol2", "Pol II", col_names[grep("Pol2", col_names)])
  colnames(ChIP_TSS_B2_1_BAM_mat) <- col_names
  
  ChIP_TSS_B2_1_BAM_IN <- ChIP_TSS_B2_1_BAM_mat[, grep("INPUT_", colnames(ChIP_TSS_B2_1_BAM_mat))] # fill zeros with medians
  ChIP_TSS_B2_1_BAM_IN <- apply(ChIP_TSS_B2_1_BAM_IN, 2, function(x) ifelse(x == 0, median(x), x))
  
  ChIP_TSS_B2_1_BAM_mat_norm <- NA
  for (i in unique(gsub("(_.*)", "\\2", colnames(ChIP_TSS_B2_1_BAM_mat))) ) {
    ChIP_TSS_B2_1_BAM_mat_norm <- cbind(ChIP_TSS_B2_1_BAM_mat_norm,
                                      sweep(ChIP_TSS_B2_1_BAM_mat[, grep(paste0("^", i, "_"), colnames(ChIP_TSS_B2_1_BAM_mat))],
                                            2, size_factor_cal(ChIP_TSS_B2_1_BAM_IN), "/")
    )
  }
  ChIP_TSS_B2_1_BAM_mat_norm[, c(1, grep("INPUT_", colnames(ChIP_TSS_B2_1_BAM_mat_norm)))] <- NULL
  # ChIP_TSS_B2_1_BAM_mat_norm <- ChIP_TSS_B2_1_BAM_mat_norm[!apply(ChIP_TSS_B2_1_BAM_mat_norm, 1, function(x) any(x == 0)), ]
  
  sf_P0_B2 <- size_factor_cal(ChIP_TSS_B2_1_BAM_mat[, grepl("P0", colnames(ChIP_TSS_B2_1_BAM_mat)) &
                                                    !grepl("INPUT", colnames(ChIP_TSS_B2_1_BAM_mat))])
  for (i in names(sf_P0_B2)) {
    idx <- grep(gsub("_P0", "", i), colnames(ChIP_TSS_B2_1_BAM_mat_norm))
    ChIP_TSS_B2_1_BAM_mat_norm[, idx] <- ChIP_TSS_B2_1_BAM_mat_norm[, idx] / sf_P0_B2[i]
  }
  
  # transform to log2FC
  ChIP_TSS_B2_1_BAM_mat_log2FC <- NULL
  for (i in unique(gsub("(_.*)", "\\2", colnames(ChIP_TSS_B2_1_BAM_mat_norm))) ) {
    idx <- grep(paste0("^", i, "_"), colnames(ChIP_TSS_B2_1_BAM_mat_norm))
    tmp <- trim_quantile(ChIP_TSS_B2_1_BAM_mat_norm[, idx])
    colnames(tmp) <- colnames(ChIP_TSS_B2_1_BAM_mat_norm[, idx])
    tmp <- log2(tmp / tmp[, 1])
    ChIP_TSS_B2_1_BAM_mat_log2FC <- cbind(ChIP_TSS_B2_1_BAM_mat_log2FC, tmp)
  }
  rownames(ChIP_TSS_B2_1_BAM_mat_log2FC) <- rownames(ChIP_TSS_B2_1_BAM_mat_norm)
}



# ------------------------------------------- Public ChIP-seq data ---------------------------------------------- #

if (FALSE) {
  
  # ------------------------------- Blackledge, Ring1b CPM ------------------------------------- #
  bw_files <- list.files("/mnt/0E471D453D8EE463/GEO_bw/Ring1b_CKO", "Blackledge_PRC1CPM", full.names = T)
  sample_names <- gsub(".*PRC1CPM_(.*)_rep.*", "\\1", bw_files)
  
  ChIP_Blackledge_TSS <- .countBW(bw_files = bw_files,
                                  intervals = promoters(gene.gr, upstream = 1000, downstream = 1000),
                                  fast = FALSE) 
  
  ChIP_Blackledge_TSS_cmb <- sapply(unique(sample_names),
                                    function(x)
                                      rowMeans(all_ChIP_rerun_TSS_mat[, grep(x, sample_names)]))
  # TAM: tamoxifen-treated (72hr OHT) vs UNT; cell Ring1b(WT->I53A/D56K)fl/fl
  ChIP_Blackledge_TSS_cmb_log2FC <-
    sapply(unique(gsub( "(TAM|UNT)_", "", colnames(ChIP_Blackledge_TSS_cmb) )),
           function(x)
             ChIP_Blackledge_TSS_cmb[, grep(x, colnames(ChIP_Blackledge_TSS_cmb))] %>%
             log2() %>% rowDiffs() #%>% "*"(-1)
    ) %>% `rownames<-`(., rownames(ChIP_Blackledge_TSS_cmb)) 
  
  
  # --------------------------------------- Tamburri ------------------------------------------- #
  
  bw_files <- list.files("/mnt/0E471D453D8EE463/GEO_bw/Ring1b_CKO", "Tamburri", full.names = T)
  
  ChIP_Tamburri_TSS <- .countBW(bw_files = bw_files,
                                intervals = promoters(gene.gr, upstream = 1000, downstream = 1000),
                                fast = FALSE) 
  
  # I53S48: Ring1b catalytical dead + OHT (CPM)
  # WTOHT: Ring1b wild type + OHT (rescue)
  # PARETA: parental cell (control)
  # PAROHT: parental cell + OHT (CKO)
  
  # so use CPM (I53S48 vs WTOHT) or (I53S48 vs PARETA)
  ChIP_Tamburri_TSS <- ChIP_Tamburri_TSS[, -grep("I53SOHT-MTF2", colnames(ChIP_Tamburri_TSS))]
  targets <- unique(gsub(".*-(.*).mm9.bw", "\\1", colnames(ChIP_Tamburri_TSS)))
  
  ChIP_Tamburri_TSS_log2FC <- sapply(targets, function(x) {
    tmp <- ChIP_Tamburri_TSS[, grep(x, colnames(ChIP_Tamburri_TSS))]
    if (x == "MTF2") {
      log2(tmp[, grep("I53S48", colnames(tmp))] / tmp[, grep("WTOHT", colnames(tmp))])
    } else {
      log2(tmp[, grep("I53SOHT|I53S48", colnames(tmp))] / tmp[, grep("WTOHT", colnames(tmp))])
    }
  }) %>% `rownames<-`(., rownames(ChIP_Tamburri_TSS)) 
  
  # --------------------------------------- Kolovos ------------------------------------------- #
  bw_files <- list.files("/mnt/0E471D453D8EE463/GEO_BAP1/2020_Kolovos", 
                         pattern = "in_wt_mESCs.*bw$", recursive = T, full.names = T)
  
  ChIP_Kolovos_TSS <- .countBW(bw_files = bw_files,
                               intervals = promoters(gene.gr, upstream = 1000, downstream = 1000),
                               fast = FALSE)
  sample_names <- unique(gsub("_ChIP.*", "", colnames(ChIP_Kolovos_TSS)))
  
  ChIP_Kolovos_TSS <- sapply(sample_names, 
                             function(x) {
                               idx <- grep(x, colnames(ChIP_Kolovos_TSS))
                               if (length(idx) == 1) {
                                 out <- ChIP_Kolovos_TSS[, idx]
                               } else {
                                 out <- rowMeans(as.matrix(ChIP_Kolovos_TSS[, idx]))
                               }
                               out
                             }) 
  
  # --------------------------------------- Kweon ------------------------------------------- #
  bw_files <- list.files("/mnt/0E471D453D8EE463/GEO_BAP1/2019_Kweon", 
                         pattern = ".*bw$", recursive = T, full.names = T)
  
  ChIP_Kweon_TSS <- .countBW(bw_files = bw_files,
                             intervals = promoters(gene.gr, upstream = 1000, downstream = 1000),
                             fast = FALSE)
  
  # --------------------------------------- Conway ------------------------------------------- #
  bw_files <- list.files("/mnt/0E471D453D8EE463/GEO_BAP1/2021_Conway/", 
                         pattern = ".*bw$", recursive = T, full.names = T)
  
  ChIP_Conway_TSS <- .countBW(bw_files = bw_files,
                              intervals = promoters(gene.gr, upstream = 1000, downstream = 1000),
                              fast = FALSE)
  
  sample_names <- unique(gsub("_rep.", "", 
                              gsub("2021_Conway_(WT|WT_EV|pCAG)_(.*).mm9.bw", "\\2",
                                   colnames(ChIP_Conway_TSS))))
  
  ChIP_Conway_TSS <- sapply(sample_names, 
                            function(x) {
                              idx <- grep(x, colnames(ChIP_Conway_TSS))
                              if (length(idx) == 1) {
                                out <- ChIP_Conway_TSS[, idx]
                              } else {
                                out <- rowMeans(as.matrix(ChIP_Conway_TSS[, idx]))
                              }
                              out
                            }) 
  
  # --------------------------------------- Fursova ------------------------------------------- #
  bw_files <- list.files("/mnt/0E471D453D8EE463/GEO_BAP1/2021_Fursova", 
                         pattern = ".*bw$", recursive = T, full.names = T)
  
  ChIP_Fursova_TSS <- .countBW(bw_files = bw_files,
                               intervals = promoters(gene.gr, upstream = 1000, downstream = 1000),
                               fast = FALSE)
  
  sample_names <- unique(gsub("GSE161993_mESC_BAP1ff_(.*)_mm10.*", "\\1", colnames(ChIP_Fursova_TSS)))
  
  colnames(ChIP_Fursova_TSS) <- sample_names
}

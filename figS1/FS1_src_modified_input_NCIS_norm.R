# modified input NCIS normalization
# Rui Shao, Jun 2022

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
source("../util/getCoverage.R")
source("../util/norm_ChIP_utils.R")

# ---------------- ----------------- ----------------- ----------------- ----------------- #
# batch1 files
file_paths_rep <- c(list.files("/mnt/0E471D453D8EE463/BAP1_MINUTE_ChIP/batch1/bam", "RC.mm9.bam$", full.names = TRUE))
file_paths_rep_INPUT <- file_paths_rep[grep("IN_", file_paths_rep)]
file_paths_rep <- file_paths_rep[!grepl("INPUT", file_paths_rep)]

chip_targets <- unique(gsub(".*bam/(.*)_.*_.*", "\\1", file_paths_rep))
sample_conditions <- gsub(".*bam/IN_(.*)_RC.mm9.bam", "\\1", file_paths_rep_INPUT)

INPUT_Rle_list <- sapply(file_paths_rep_INPUT, 
                         function(file_path) {
                           bin_smooth_bam_to_Rle(file_path, bin_width = 1e5, stride_len = 1e3) # make sure no empty input region on genes
                         })
names(INPUT_Rle_list) <- sample_conditions

INPUT_Rle_mean_size <- mean(sapply(INPUT_Rle_list, function(x) sum(sum(x)) / sum(lengths(x)))) # to 1x genome coverage


# normalize to input coverage
dir.create("../bw_rep/batch1/input_norm/", recursive = TRUE)

for (chip in chip_targets) {
  tmp_files <- file_paths_rep[grep(paste0("/", chip, "_"), file_paths_rep)]
  tmp_conditions <- gsub(paste0(".*", chip, "_(.*)_RC.mm9.bam"), "\\1", tmp_files)
  tmp_Rle_list <- sapply(tmp_files,
                         function(x) 
                           coverage(readGAlignments(x, param = ScanBamParam(mapqFilter = 10)))[names(INPUT_Rle_list[[1]])]
  )
  names(tmp_Rle_list) <- tmp_conditions
  
  for (x in unique(tmp_conditions)) {
    export.bw(tmp_Rle_list[[x]] / INPUT_Rle_list[[x]] * INPUT_Rle_mean_size, 
              paste0("../bw_rep/batch1/input_norm/", chip, "_", x, ".bw")) 
  }
}


# ---------------- ----------------- ----------------- ----------------- ----------------- #
# batch2 files
file_paths_rep <- c(list.files("/mnt/0E471D453D8EE463/BAP1_MINUTE_ChIP/batch2/bam", "rep.*bam$", full.names = TRUE))
file_paths_rep_INPUT <- file_paths_rep[grep("INPUT", file_paths_rep)]
file_paths_rep <- file_paths_rep[!grepl("INPUT", file_paths_rep)]

sample_conditions <- unique(gsub(".*(P.*rep.).*", "\\1", file_paths_rep_INPUT))
chip_targets <- unique(gsub(".*bam/(.*)_P.*", "\\1", file_paths_rep))

# (combine two rounds of sequencing)
# load input Rle

INPUT_Rle_list <- sapply(sample_conditions, 
                         function(sample_condition) {
                           file_paths <- file_paths_rep_INPUT[grep(sample_condition, file_paths_rep_INPUT)]
                           bin_smooth_bam_to_Rle(file_paths, bin_width = 1e5, stride_len = 1e2) # make sure no empty input region on genes
                         })
INPUT_Rle_mean_size <- mean(sapply(INPUT_Rle_list, function(x) sum(sum(x)) / sum(lengths(x)))) # to 1x genome coverage

# estimate background ratio
bg_rep.size <- sapply(chip_targets, 
                      function(target) {
                        get_bg_est_res(file_paths_rep[grep(target, file_paths_rep)],
                                       file_paths_rep_INPUT,
                                       sample_conditions) }
                      )

# normalize to input coverage
dir.create("../bw_rep/batch2/input_norm/", recursive = TRUE)

for (chip in chip_targets) {
  tmp_files <- file_paths_rep[grep(paste0("/", chip, "_"), file_paths_rep)]
  tmp_conditions <- gsub(".*(P.*rep.).*", "\\1", tmp_files)
  tmp_Rle_list <- sapply(tmp_files,
                         function(x) 
                           coverage(readGAlignments(x, param = ScanBamParam(mapqFilter = 10)))[names(INPUT_Rle_list[[1]])]
                         )
  names(tmp_Rle_list) <- tmp_conditions
  
  tmp_bg_rep.size <- bg_rep.size[, chip] / mean(bg_rep.size[-(1:2), chip]) # stabilize within the same treatment
  tmp_bg_rep.size[1:2] <- 1 # skip bg normalization for control 
  
  for (x in unique(tmp_conditions)) {
    idx <- which(names(tmp_Rle_list) == x)
    export.bw(Reduce("+", tmp_Rle_list[idx]) / INPUT_Rle_list[[x]] * INPUT_Rle_mean_size / tmp_bg_rep.size[x], 
              paste0("../bw_rep/batch2/input_norm/", chip, "_", x, ".bw")) 
  }
}

# merge replicates
bw_file_list <- list.files("../bw_rep/batch2/input_norm", "rep..bw", full.names = TRUE)

for (chip in unique(gsub(".*input_norm/(.*)_.*_.*.bw", "\\1", bw_file_list)) ) { 
  tmp_files <- bw_file_list[grep(paste0("/", chip, "_"), bw_file_list)]
  tmp_names <- gsub(".*input_norm/(.*)_rep..bw", "\\1", tmp_files)
  sapply(unique(tmp_names), 
         function(target) { 
           tmp <- tmp_files[grep(paste0(target, "_"), tmp_files)]
           print(tmp)
           export.bw(Reduce("+", sapply(tmp, function(x) import.bw(x, as = "RleList"))),
                     paste0("../bw_rep/batch2/input_norm/", target, "_pool.bw"))
           message(paste(target, "done."))
         }
  )
}


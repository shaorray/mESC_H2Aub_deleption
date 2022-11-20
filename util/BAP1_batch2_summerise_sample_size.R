# sample size summary
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("../util/getCoverage.R")
source("../util/calibration/background_calibration.R")
library(Rsamtools)

# input reads
# input reads mapq30
# chip reads
# chip reads, mapq30
# chip bg reads
# bg reads
# NCIS size
# deepTools size

# ---------------- ----------------- ----------------- ----------------- -----------------
# files
file_paths_rep <- c(list.files("/mnt/0E471D453D8EE463/baseSpace/batch2/final/bam", "rep.*bam$", full.names = T))
file_paths_rep_INPUT <- file_paths_rep[grep("IN", file_paths_rep)]
file_paths_rep <- file_paths_rep[!grepl("IN", file_paths_rep)]

# 1. NCIS background -----------------------------------------------------------
get_bg_est_res <- function(file_paths_rep) {
  bg_rep.res <- list()
  for (chip in seq_along(file_paths_rep)) {
    target <- gsub(".*bam/(.*)_.*_.*.UCSC.bam", "\\1", chip)
    .rep <- gsub(".*_(.*_.*).UCSC.bam", "\\1", chip)
    
    input.path <- file_paths_rep_INPUT[grep(.rep, file_paths_rep_INPUT)]
    bg_rep.res <- c(bg_rep.res, 
                    list(NCIS.run(chip.path = chip, 
                                  input.path = input.path, 
                                  min.binsize = 1000,
                                  max.binsize = 50000,
                                  quant = 0.75)))
  }
  bg_rep.res
}

bg_rep.res <- get_bg_est_res(file_paths_rep)
names(bg_rep.res) <- gsub(".*/bam/(.*).UCSC.bam", "\\1", file_paths_rep)
saveRDS(bg_rep.res, "bg_norm/bg_rep.res.rds")

bg_rep.size <- unlist(lapply(bg_rep.res, function(x) x$est))
saveRDS(bg_rep.size, "bg_norm/bg_rep.size.rds")


# 2. Diaz et al. norm ----------------------------------------------------------------

source("../util/calibration/deeptools_scalefactor.R")

diy.res <- NULL
for (chip in file_paths_rep) {
  input.path <- file_paths_rep_INPUT[grep(gsub(".*_(.*_.*).UCSC.bam", "\\1", chip), file_paths_rep_INPUT)]
  
  diy.res <- rbind(diy.res, 
                   estimateScaleFacto_R(chip.path = chip, 
                                        input.path = input.path)
                   )
}

diy.res <- `names<-`(ifelse(diy.res[, 1] == 1, diy.res[, 2], 1 / diy.res[, 1]),
                     gsub(".*/bam/(.*).UCSC.bam", "\\1", file_paths_rep))


# 3. get read counts by different attributes -----------------------------------------


sbp1 <- ScanBamParam(mapqFilter = 30)
sbp2 <- ScanBamParam(flag = scanBamFlag(isPaired = TRUE))

input_rep_table <- NULL
for (i in seq_along(file_paths_rep_INPUT)) { 
  
  bg_tab <- data.frame("sample" = gsub(".*bam/(.*).UCSC.bam", "\\1", file_paths_rep_INPUT[i]),
                       
                       "input_reads_all" = Rsamtools::countBam(file_paths_rep_INPUT[i])$records,
                       
                       "input_reads_mapq30" = Rsamtools::countBam(file_paths_rep_INPUT[i], param = sbp1)$records,
                       
                       "input_reads_paired" = Rsamtools::countBam(file_paths_rep_INPUT[i], param = sbp2)$records 
  )
  
  input_rep_table <- rbind(input_rep_table, bg_tab)
}


rep_rc_tab <- NULL
for (i in file_paths_rep) { 
  sample_name <- gsub(".*bam/(.*).UCSC.bam$", "\\1", i)
  .rep <- gsub(".*_(.*_rep.).UCSC.bam", "\\1", i)
  
  bg.gr <- bg_rep.res[[sample_name]]$bg.bin
  
  bg_tab <- data.frame("sample" = gsub(".*bam/(.*).UCSC.bam", "\\1", i),
                       
                       "chip_reads_all" = Rsamtools::countBam(i)$records,
                       
                       "chip_reads_mapq30" = Rsamtools::countBam(i, param = sbp1)$records,
                       
                       "chip_reads_paired" = Rsamtools::countBam(i, param = sbp2)$records,
                       
                       "chip_bg_reads" = sum(.countBam(bam_files = i, intervals = bg.gr)),
                       
                       "chip_bg_length" = sum(width(bg.gr)),
                       
                       "bg_NCIS_size" = bg_rep.size[sample_name], 
                       
                       "bg_Diaz_et_al_size" = diy.res[sample_name])
  
  bg_tab <- cbind(bg_tab, input_rep_table[grep(.rep, input_rep_table[, 1]), ])
  rep_rc_tab <- rbind(rep_rc_tab, bg_tab)
}

write.table(rep_rc_tab, "bg_norm/BAP1_batch2_RC_summary.txt", quote = F, row.names = F)

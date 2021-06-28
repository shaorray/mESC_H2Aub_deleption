# Rui Shao
# 2021 May

# Filter artificial peaks in ChIP assays
# Input bw coverage files with false peaks in squared shape
# Ouput common overlaps of the squared outlier peaks

library(rtracklayer)
library(GenomicRanges)

# Files that contain artificial peaks
bw_files <- list.files("/mnt/E0767589767560E8/UPPMAX/BAP1_PC",
                       ".bw$", full.names = T) 

# Keep peaks with extreme values
get_outlier_peak <- function(bw_file, .quantile = 0.999){
  bw <- import(bw_file)
  reduce(bw[bw$score > quantile(bw$score, .quantile)])
}

peak_list <- sapply(bw_files, get_outlier_peak)

# View. Small sample size could decrease the number of peaks
peak_nums <- unname(Reduce("c", lapply(a, length)))
hist(peak_nums, main = "Sample outlier peaks number")

# Get the common peaks
intersect_gr <- function(gr_list) {
  if (length(gr_list) == 1) {
    return(gr_list[[1]])
  } else {
    GenomicRanges::intersect(gr_list[[1]], intersect_gr(gr_list[-1]))
  }
}

# Output result
peak_list_large <- peak_list[peak_nums > 1500]
peak_out <- intersect_gr(peak_list_large)

peak_out$score <- 0
export.bed(peak_out, "blacklist_peak_mm9.bed")

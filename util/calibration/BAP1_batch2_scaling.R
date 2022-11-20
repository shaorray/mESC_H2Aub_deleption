# background normalization with input sample


# load paths
chip.paths <- list.files("/mnt/0E471D453D8EE463/baseSpace/batch2/final/bam",
                         pattern = ".*bam$", full.names = T)
chip.paths <- chip.paths[-grep("pooled", chip.paths)]
input.paths <- chip.paths[grep("sfINPUT", chip.paths)]


sample_names <- unique(gsub(".*bam\\/(.*)_P.*", "\\1", chip.paths[!grepl("sfINPUT", chip.paths)]))

# ChIP/Input ratio list
res.list <- list()
for (chip in sample_names) {
  
  tmp.paths <- chip.paths[grep(chip, chip.paths)]
  
  tmp.list <- list()
  for (i in seq_along(tmp.paths)) {
    tmp.list <- c(tmp.list, 
                  list(NCIS.test(chip.path = tmp.paths[i], 
                                 input.path = input.paths[i], 
                                 max.binsize = 5000,
                                 quant = 0.75))
                  )
  }
  names(tmp.list) <- gsub(".*bam\\/(.*).UCSC.bam", "\\1", tmp.paths)
  
  res.list <- c(res.list, tmp.list)
}


bg.sf <- unlist(lapply(res.list2, function(x) x$est))

rc.ratio <- unlist(lapply(res.list2, function(x) x$r.seq.depth))

rc.sf <- SizeFactorCal(sapply(chip.paths[grepl("sfINPUT", chip.paths)], function(x) Rsamtools::idxstatsBam(x)[, 3]))

norm.factor <- bg.sf * rc.sf

write.table(cbind(names(norm.factor), unname(norm.factor)), 
            file = "/mnt/0E471D453D8EE463/baseSpace/batch2/bg_scaling.txt",
            row.names = F, col.names = F,
            quote = F)

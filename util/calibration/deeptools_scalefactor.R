# rebuild background estimation of scaling factor that has been applied in:
#   - Diaz et al (2012)
#   - deepTools: estimateScaleFactor()


estimateScaleFacto_R <- function(chip.path, 
                                 input.path,
                                 binsize = 1e4,
                                 chr.vec = NULL) {
  # process bam files
  sbp <- Rsamtools::ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = F,
                                                    isDuplicate = F), 
                                 what = c("rname", "pos"))
  
  chip.pos <- Rsamtools::scanBam(chip.path, param = sbp)
  idx <- !is.na(chip.pos[[1]]$pos) & !grepl("random", chip.pos[[1]]$rname)
  chip.pos <- lapply(chip.pos[[1]], function(x) x[idx])
  chip.pos <- split(chip.pos$pos, chip.pos$rname, drop = TRUE)
  
  input.pos <- Rsamtools::scanBam(input.path, param = sbp)
  idx <- !is.na(input.pos[[1]]$pos) & !grepl("random", input.pos[[1]]$rname)
  input.pos <- lapply(input.pos[[1]], function(x) x[idx])
  input.pos <- split(input.pos$pos, input.pos$rname, drop = TRUE)
  
  if (is.null(chr.vec)) {
    chip.name <- names(chip.pos)
    input.name <- names(input.pos)
    chr.vec <- intersect(chip.name, input.name)
    if (length(chr.vec) < length(chip.name)){
      cat("Control sample doesn't have chromosome", 
          setdiff(chip.name, chr.vec), 
          "which are in ChIP sample. These chromosomes are ignored.\n")
    }
    if (length(chr.vec) < length(input.name)){
      cat("ChIP sample doesn't have chromosome", 
          setdiff(chip.name, chr.vec), 
          "which are in control sample. These chromosomes are ignored.\n")
    }
  }
  
  chr.end.max <- sapply(chr.vec, 
                        function(chr) 
                          max( c(max(chip.pos[[chr]]),
                                 max(input.pos[[chr]])) )
  )
  
  bindata <- bin.data(chip.pos, 
                      input.pos, 
                      binsize, 
                      chr.vec = chr.vec, 
                      chr.end.max = chr.end.max)
  
  ind <- bindata$chip + bindata$input > 0
  
  res <- est.sizeFactorsSES(chip = bindata$chip[ind], 
                            input = bindata$input[ind])
}


est.sizeFactorsSES <- function(chip, input) {
  chip <- cumsum(sort(chip))
  input <- cumsum(sort(input))
  
  diff <- abs(chip / tail(chip, 1) - input / tail(input, 1))
  
  maxIndex <- which.max(diff)
  
  maxIndex <- as.integer(maxIndex * 0.8)
  while (maxIndex < length(chip)) {
    # in rare cases the maxIndex maps to a zero value.
    # In such cases, the next index is used until
    # a non zero value appears.
    cumSum <- c(chip[maxIndex], input[maxIndex])
    if (min(cumSum) > 0)
      break
    maxIndex <- maxIndex + 1
  }
  
  sizeFactorsSES <- min(cumSum) / cumSum
  
  return(sizeFactorsSES)
}

bin.data <- function(chip.pos, 
                     input.pos, 
                     binsize,
                     chr.vec = NULL, 
                     chr.end.max = NULL) { 
  chip <- list()
  input <- list()
  
  for(chr in chr.vec){
    bk <- seq(from = 0, by = binsize, 
              to = ceiling(chr.end.max[[chr]] / binsize) * binsize)
    chip[[chr]] <- hist(chip.pos[[chr]], breaks = bk, plot = FALSE)$counts
    input[[chr]] <- hist(input.pos[[chr]], breaks = bk, plot = FALSE)$counts
    
  }
  list(chip = unlist(chip, use.names = FALSE), 
       input = unlist(input, use.names = FALSE))
}
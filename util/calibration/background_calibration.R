# A method from "Normalization of ChIP-seq data with control"

if (!require("pacman")) install.packages("pacman")
pacman::p_load(Rsamtools, GenomicAlignments)


#' A function normalise sample size with background reads in the ChIPed
#' and input samples, by self-adapting bin size and low-cut of reads.
#' 
#' The figure is faceted according to scaling group.
#'
#' @param chip.path a file path of ChIP sample, BAM or BW.
#' @param input.path a file path of Input sample, BAM or BW
#' 
#' @param min.binsize start of bin scan size 
#' @param max.binsize end of bin scan size 
#' 
#' @param chr.vec a vector of chromasome names, e.g. c("chr1", "chr2")  
#' @param quant a quantile of the ceiling background bin counts   
#' 
#' @return A list of estmated ratio, estimated bg bin size, read cout ratio

NCIS.run <- function(chip.path, 
                     input.path, 
                     which.bg = NULL,
                     min.binsize = 100, 
                     max.binsize = 20000,  
                     min.stop.binsize = 100, 
                     chr.vec = NULL, 
                     quant = 0.75) {
  
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
  
  nchip <- sum(sapply(chr.vec, function(x) length(chip.pos[[x]])))
  ninput <- sum(sapply(chr.vec, function(x) length(input.pos[[x]])))
  r.seq.depth <- nchip / ninput
  
  sizevec <- rep(c(1,2,5), times = 3) * 10^rep(2:4, each=3)
  sizevec <- sizevec[sizevec <= max.binsize & sizevec >= min.binsize]
  norm.est <- rep(1000, length(sizevec))
  
  chr.end.max <- sapply(chr.vec, 
                        function(chr) 
                          max( c(max(chip.pos[[chr]]),
                                 max(input.pos[[chr]])) )
                        )
  
  # cache bin intervals
  empty.bin <- NULL
  bg.bin <- NULL
  
  binsize.est <- -1
  
  for(si in 1:length(sizevec)) {
    binsize <- sizevec[si]
    
    bindata <- bin.data(chip.pos, 
                        input.pos, 
                        binsize, 
                        chr.vec = chr.vec, 
                        chr.end.max = chr.end.max)
    ind <- bindata$chip + bindata$input > 0
    
    if (!is.null(which.bg)) {
      bin.tile <- tileGenome(seqlengths = chr.end.max, 
                             tilewidth = binsize, 
                             cut.last.tile.in.chrom = T)
      which.bg <- sortSeqlevels(which.bg)
      bin.tile <- sortSeqlevels(bin.tile)
      seqlengths(bin.tile) <- seqlengths(which.bg) <- pmax(seqlengths(bin.tile), seqlengths(which.bg))
      ind <- ind & countQueryHits(findOverlaps(bin.tile, which.bg)) > 0
    }
    
    res <- est.norm.quant.search(chip = bindata$chip[ind], 
                                 input = bindata$input[ind],
                                 quant = quant)
    
    norm.est[si] <- res$ratio
    
    # stopping criteria
    if (si > 1 & binsize.est < 0) {
      if (norm.est[si] >= norm.est[si-1]) {
        est <- norm.est[si-1]
        binsize.est <- sizevec[si-1]
      }
      if (si == length(sizevec) & binsize.est < 0) { # the end of binsize, no converge yet
        est <- norm.est[si]
        binsize.est <- sizevec[si]
      }
    } # end if(si>1 & binsize.est<0)
    
    if (binsize.est > 0 & binsize >= min.stop.binsize) {
      # extract bg bin
      if (!exists("bin.tile", inherits = F)) {
        bin.tile <- tileGenome(seqlengths = chr.end.max, 
                               tilewidth = binsize, 
                               cut.last.tile.in.chrom = T)
      }
      empty.bin <- IRanges::reduce(bin.tile[ind])
      bg.bin <- IRanges::reduce(bin.tile[which(ind)[res$bg.ind]])
      break
    }
  } 
  
  return(list(est = est, 
              binsize.est = binsize.est,
              r.seq.depth = r.seq.depth,
              empty.bin = empty.bin,
              bg.bin = bg.bin)
         )
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


est.norm.quant.search <- function(chip, input, quant = 0.75){
  
  total <- chip + input
  tbl <- table(total)
  total.count <- as.integer(names(tbl))
  cum.bc <- cumsum(tbl) # increment of bin number by #reads
  cum.prop <- cum.bc / length(total)
  
  if(cum.prop[1] >= quant) {
    threshold <- 1
  } else {
    # largest total before quant
    threshold <- max(total.count[cum.prop < quant])
  }
  ind <- total <= threshold
  chip.sum.low <- sum(chip[ind])
  input.sum.low <- sum(input[ind])
  
  .chip <- chip[!ind]
  .input <- input[!ind]
  .od <- order(.chip + .input)
  .chip <- .chip[.od]
  .input <- .input[.od]
  
  # start after threshold
  cum.bc.high <- 
    cum.bc[total.count > threshold] - 
    cum.bc[which(total.count == threshold)]
  
  all.cum.chip <- cumsum(.chip)
  all.cum.input <- cumsum(.input)
  
  # search from an arbitrary grade
  cum.chip <- c(0, all.cum.chip[cum.bc.high]) # number of read sums at each increment above the threshold
  cum.input <- c(0, all.cum.input[cum.bc.high])
  cum.ratio <- (chip.sum.low + cum.chip) / (input.sum.low + cum.input)
  
  if(sum(diff(cum.ratio) >= 0) == 0) 
    stop("No increase in cum.ratio, binding signal missing!")
  
  index <- min( (2:length(cum.ratio))[diff(cum.ratio) >= 0] )
  
  return(list(ratio = cum.ratio[index],
              index = index,
              bg.ind = which(total <= (index + threshold - 1) )
              )
         )
}


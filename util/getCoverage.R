# Rui Shao
# 2022 Jun

options(scipen=999)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(BiocManager, GenomicRanges, IRanges, dplyr, rtracklayer,
               GenomicAlignments, BSgenome, Rsamtools, foreach, doParallel)

registerDoParallel(cores = 4)

process_coverage <- function(file_name)
{ # this function converts input bam/bw file to coverage rle list
  # if input a paired-end bam file, coverage combines the reads strandness
  suffix <- toupper(sapply(strsplit(file_name,'\\.'), tail, 1))
  if (suffix %in% c('BW', 'BIGWIG') )
  {
    input_cov <- rtracklayer::import.bw(file_name, as = "RleList")
    if ("X" %in% names(input_cov)) 
    {
      names(input_cov) <- paste0("chr", names(input_cov))
    }
  } else if (suffix == 'BAM') 
  {
    input_cov <- coverage(GenomicAlignments::readGAlignments(file_name))
  } else {
    return(NULL)
  }
  input_cov <- input_cov[as.character(names(input_cov)) %in% paste0("chr", c(1:100, "X", "Y"))]
  return(input_cov)
}


bin_sums <- function(x, bin_width = NULL, stride = NULL, is_mean = FALSE)
{
  # Args:
  #   x:          a vector
  #   bin_width:  the window for sum
  #   stride:     bin step length
  #   add_name:   bin names
  
  if (is.null(bin_width)) bin_width = 1000
  if (is.null(stride) || bin_width < stride) stride = bin_width
  
  x = as.numeric(x)
  len = length(x)
  
  # bin ranges
  starts = c(1, seq_len((len - 1) %/% stride) * stride + 1)

  # extend bin window
  starts = starts - ((bin_width - stride) %/% 2) 
  ends = starts + bin_width - 1
  
  # trim
  starts[starts < 1] = 1
  ends[ends > len] = len
  
  bin_effective_len = ends - starts + 1
  
  # bin sums by stride
  if (any(is.na(x))) { # avoid NAs 
    na_c_sum = cumsum(is.na(x))
    bin_effective_len = bin_effective_len - (na_c_sum[ends] - na_c_sum[starts])
    x[is.na(x)] = 0
  }
  
  c_sums = cumsum(x)
  
  # output
  if (is_mean) {
    
    return((c_sums[ends] - c_sums[starts] + x[starts]) / bin_effective_len)
    
  } else {
    
    return(c_sums[ends] - c_sums[starts] + x[starts])
    
  }
}


get_cov_matrix <- function(intervals, input_file, bin_width, new_len, .df = 20)
{ # Args:
  # intervals: genomic intervals to extract coverage matrix
  # input_file: a file for process_coverage(), bw or bam
  # extract coverage from input interval, save as a matrix
  
  input_cov <- process_coverage(input_file)
  
  # intervals_order <- order(intervals)
  
  .seqnames <- names(input_cov)[names(input_cov) %in% seqnames(intervals)]
  out_matrix <- foreach(chr = .seqnames, .combine = rbind, .inorder = TRUE) %dopar%
    {
      tmp_gr <- intervals[seqnames(intervals) == chr]
      tmp_cov <- input_cov[[chr]]
      chr.end <- length(tmp_cov)
      gr.start <- min(start(tmp_gr))
      
      # padding ranges outside chromosome boundary
      if (gr.start < 1) {
        tmp_cov <- c(Rle(values = 0, lengths = abs(gr.start) + 1), tmp_cov)
        end(tmp_gr) <- end(tmp_gr) + abs(gr.start) + 1
        start(tmp_gr) <- start(tmp_gr) + abs(gr.start) + 1
      }
      
      if (max(end(tmp_gr)) > length(tmp_cov)) 
        tmp_cov <- c(tmp_cov, Rle(values = 0, lengths = max(end(tmp_gr)) - chr.end) )
      
      # Set a output new length if input intervals are different lengths, or resize with a bin width
      get_new_len <- function(n, tmp_cov, tmp_gr, .df, new_len) {
        
        cov_n <- as.numeric(tmp_cov[ranges(tmp_gr[n])])
        
        # remove invalid values 
        cov_n[is.na(cov_n) | is.infinite(cov_n)] <- 0
        
        if (all(cov_n == 0)) {
          return(rep(0, new_len))
        } else {
          x = bin_sums(cov_n, bin_width)
          if (.df > 1) x = smooth.spline(x = seq_along(x), 
                                         y = as.numeric(x), 
                                         df = .df)$y
          x = spline(x = seq_along(x),
                     y = x,
                     n = new_len, 
                     method = 'natural')$y
          x[x < 0] = 0
          return(x)
        }
      }   
      
      if (!(bin_width > 0) | is.na(bin_width) | is.null(bin_width)) bin_width = 50 # default
      if (!(new_len > 0) | is.na(new_len) | is.null(new_len)) new_len = 100
      
      out_mat <- matrix(0, nrow = length(tmp_gr), ncol = new_len)
      for (n in seq_along(tmp_gr)) {
        out_mat[n, ] <- get_new_len(n, tmp_cov, tmp_gr, .df, new_len)
      }
      out_mat / bin_width
    }
  
  if (as.character(strand(intervals[1])) != '*')
  {
    neg_strand <- which( as.character(strand(intervals)) == '-' )
    neg_mat <- out_matrix[neg_strand, ]
    out_matrix[neg_strand, ] <- neg_mat[, rev(seq_len(ncol(neg_mat)))]
  }
  out_matrix
}


## read file list
convert_coverage <- function(file_names,
                             intervals, 
                             bin_width = 50,
                             save_genome_cov = FALSE,
                             new_len = NA, 
                             sample_names = NA)
{ # Args:
  # file_names: BW or BAM files without specifying strandness
  # intervals: GRanges intervals, random chromosomes will be removed
  # save_genome_cov: if TRUE the full genome coverage Rle object will be in outputs
  # new_len: if NA this function will bin the interval with a fixed bin width,
  #          otherwise resize to the indicated new length
  # convert bw or bam files to a coverage array on input intervals
  
  intervals <- intervals[seqnames(intervals) %in% paste0("chr", c(1:100, "X", "Y"))]
  intervals <- sort(intervals)
  
  if (is.na(new_len))
  {
    out_array <- array(dim = c(length(intervals),
                               ceiling(width(intervals[1]) / bin_width),
                               length(file_names)))
  } else {
    message("Resize to the new length.")
    out_array <- array(dim = c(length(intervals),
                               new_len,
                               length(file_names)))
  }
  
  gcov_list <- list()
  for ( i in seq_along(file_names))
  {
    out_array[, , i] <- get_cov_matrix(intervals, input_file = file_names[i], bin_width, new_len)
    if (save_genome_cov) gcov_list <- c(gcov_list, list(input_cov))
  }
  
  if (any(is.na(sample_names))) sample_names <- seq_len(dim(out_array)[3])
  if (!is.null(names(intervals))) {
      row_names <- names(intervals)
    } else { 
      row_names <- seq_len(dim(out_array)[1])
    }
  dimnames(out_array) <- list(row_names, 
                              seq_len(dim(out_array)[2]), 
                              sample_names)
  
  if (save_genome_cov)
  {
    return( list(out_array = out_array, gcov_list = gcov_list) )
  } else
  {
    return( out_array )
  }
}


binBwBam <- function(file_names, 
                     which_genome = 'mm10', 
                     bin_width = 10000L,
                     stride = NULL, 
                     as_rle = FALSE)
{
  # this function convert bw/bam files to a binned matrix (or a Rle list)
  # returns mean coverage for the given bin width, by each stride.
  
  if (is.null(stride)) stride <- bin_width
  
  .get_coverage_vector <- function(file_name, bin_width, stride) {
    input_cov <- process_coverage(file_name)
    foreach(chr = names(input_cov), .combine = c) %dopar% {
      bin_sums(x = as.numeric(input_cov[[chr]]),
               bin_width = bin_width, 
               stride = stride, 
               is_mean = TRUE)
    }
  }
  
  if (as_rle) {
    bin_dat <- vector(mode = "list", length = length(file_names))
    
    for (i in seq_along(file_names)) {
      tmp_dat <- .get_coverage_vector(file_names[i], 
                                      bin_width = bin_width, 
                                      stride = stride)
      bin_dat[[i]] <- Rle(tmp_dat)
    }
    
    names(bin_dat) <- gsub(".*\\/(.*)\\.*", "\\1", file_names)
    
  } else {
    # get number of bins
    chr_bin <- tileGenome(seqlengths = seqlengths(BSgenome::getBSgenome(which_genome)), 
                          tilewidth = stride, 
                          cut.last.tile.in.chrom = TRUE) # return a GRanges object, not a list
    chr_bin <- chr_bin[seqnames(chr_bin) %in% paste0("chr", c(1:100, "X", "Y"))]
    
    # initialize a matrix
    bin_dat <- matrix(0, nrow = length(chr_bin), ncol = length(file_names))
    
    for (i in seq_along(file_names))
    {
      tmp_dat <- .get_coverage_vector(file_names[i],
                                      bin_width = bin_width,
                                      stride = stride)
      bin_dat[, i] <- tmp_dat
    } # end of for loop
    
    colnames(bin_dat) <- gsub(".*\\/(.*)\\.*", "\\1", file_names)
  }
  gc()
  
  bin_dat
}


get_reads_coverage <- function(bam.files, object.gr, targets = "all", flank.size, is.frag = FALSE)
{
  # sum reads coverage from bam files
  # output a list of Rle strings of each range
  # Agrs:
  #     bam.file: paired end strand specific bam file path
  #     object.gr: multiple intervals of gene or transcript
  #     targets: "full", "exon", "intron" (depends on Gene_input: complete gene reference)
  #     flank.size: flank both ends to find precise boundary or for plotting
  #     Gene_input: complete gene reference, GENCODE/ensembl
  #     is.frag: if coverage piles from paired-end reads fragments
  # refer to NGS.plot.r
  
  object.gr = object.gr + flank.size
  object.start = start(object.gr)
  object.end = end(object.gr)
  object.strand = as.character(strand(object.gr)) == "+"

  if(targets != "all")
    exon.gr = Gene_input[Gene_input$type=="exon"]
    # exon.gr = Gene_input[Gene_input$type=="exon" & Gene_input$transcript_id%in%object.gr$transcript_id]

  for(i in seq_along(bam.files))
  {
    file.name = gsub('*.bam','\\',strsplit(bam.files[i],'\\/')%>%unlist%>%tail(1))
    file.dir = gsub(paste0('*', file.name,'.bam'), '\\', bam.files[i])
    bam.index = paste0(file.dir,file.name,'.bam.bai')
    if(!file.exists(bam.index)) bam.index = Rsamtools::indexBam(files = bam.files[i])
    invisible(capture.output(is.paired.end <- testPairedEndBam(bam.files[i], index = bam.index)))

    .is.frag = is.frag
    if(!is.paired.end & is.frag) .is.frag = F
    # scanBamWhat: the info that need to be extracted from a bam file.
    sbw <- c('pos', 'qwidth', 'mapq', 'strand', 'rname',
             'mrnm', 'mpos', 'isize')
    sbp <- ScanBamParam(what = sbw, which = object.gr,
                        flag = scanBamFlag(isUnmappedQuery = FALSE,
                                           isSecondaryAlignment = FALSE,
                                           isFirstMateRead = is.paired.end)) # use first in pair reads for coverage

    # Scan bam file to retrieve short reads.
    sr.in.ranges <- scanBam(bam.files[i], param = sbp, index = bam.index)

    combind.list = function(a, b)
    {
      new = lapply(names(a), function(x) c(a[[x]],b[[x]]))
      names(new) = names(a)
      new
    }

    scanBamRevOrder <- function(org.gr, sbp) {
      # ScanBamParam re-arranges the input genomic ranges. Use range info to
      # construct a string vector to find the order to reverse it.
      org.grnames <- with(org.gr, paste(seqnames, start, end, sep=':'))
      sbw.gr <- as.data.frame(bamWhich(sbp))  # scan-bam-ed
      if('space' %in% names(sbw.gr)) {
        sbw.grnames <- with(sbw.gr, paste(space, start, end, sep=':'))
      } else if('group_name' %in% names(sbw.gr)) {
        sbw.grnames <- with(sbw.gr, paste(group_name, start, end, sep=':'))
      } else {
        stop("Cannot locate chromosome names in extracted short reads. Report
             this problem using issue tracking or discussion forum.\n")
      }
      match(org.grnames, sbw.grnames)
    }

    # Restore the original order.
    sr.in.ranges <- sr.in.ranges[scanBamRevOrder(as.data.frame(object.gr), sbp)]

    srg.temp = list()
    for(n in seq_along(object.gr))
    {
      # Filter short reads by mapping quality and insertion size.
      all.index <- with(sr.in.ranges[[n]], (mapq>3 & strand == as.character(strand(object.gr[n])) ))
      if(is.paired.end)
      {
        all.index <- all.index & with(sr.in.ranges[[n]], !is.na(isize))
        sr.in.ranges[[n]]$isize = with(sr.in.ranges[[n]], mpos-pos )
      }
      srg.temp <- c(srg.temp, list(lapply(sr.in.ranges[[n]], function(x) x=x[all.index] )) )
    }
    names(srg.temp) = names(sr.in.ranges)
    
    # combine multiple bam files
    if (!exists('srg.filtered'))
    {
      srg.filtered = srg.temp
    } else {
      srg.filtered = lapply(names(srg.temp),
                            function(x)
                              combind.list(srg.temp[[x]], srg.filtered[[x]]))
    }
  }
  
  cov.list = list()
  for (n in seq_along(object.gr))
  {
    if (targets == "all")
    {
      cov = with( srg.filtered[[n]],
                  coverage(IRanges(start = ifelse(rep(object.strand[n], length(pos)), pos, mpos) - object.start[n],
                                   width = ifelse(rep(.is.frag, length(isize)), abs(isize), 0) + qwidth ),
                           width = width(object.gr[n]) )
                )

    } else if (targets == "exon") {
      ov.index = with( srg.filtered[[n]],
                       findOverlaps(IRanges(start = pos, width = qwidth),
                                                       ranges(exon.gr[exon.gr$gene_id == object.gr[n]$gene_id])) %>%
                         countQueryHits() != 0 )

      cov = with( lapply(srg.filtered[[n]], function(x) x[ov.index]),
                  coverage(IRanges(start = pos - object.start[n],
                                   width = ifelse(.is.frag, abs(isize),0) + qwidth ),
                           width = width(object.gr[n]) )
                )

    } else if (targets == "intron") {
      ov.index = with(srg.filtered[[n]],
                      findOverlaps(IRanges(start = pos, width = qwidth),
                                   ranges(exon.gr[exon.gr$gene_id == object.gr[n]$gene_id])) %>%
                        countQueryHits()==0 )

      cov = with(lapply(srg.filtered[[n]], function(x) x[ov.index]),
                  coverage(IRanges(start = pos - object.start[n],
                                   width = ifelse(.is.frag, abs(isize),0) + qwidth ),
                           width = width(object.gr[n]) )
                )

    }

    if(!object.strand[n]) cov = rev(cov)
    cov.list = c(cov.list, cov)
  }
  return(cov.list)
}


resize_cov = function(cov.list, df=10, len=100)
{
  foreach(i = seq_along(cov.list), .combine = rbind) %dopar%
  {
    if (class(cov.list[i]) == "list") cov.tmp = cov.list[[i]] else cov.tmp = cov.list[i]
    cov = smooth.spline(x = seq_along(cov.tmp), y = as.numeric(cov.tmp), df = df)$y
    cov = spline(x = seq_along(cov), y = cov, n = len, method = 'natural')$y
    cov[cov < 0] = 0
    cov
  }
}


# read counting and coverage counting ------------------------------------------------------

.countBW <- function(bw_files, intervals, blacklist = NULL, fast = FALSE)
{ # Args:
  # extract intervals coverages from input bw files
  # output density table (row: intervals, column: bw_files)
  
  if (fast)
  {
    count_dat <- sapply(bw_files, function (file)
      unlist(summary(
        rtracklayer::BigWigFile(file), intervals
      ))$score)
  } else {
    count_dat <-
      foreach(i = seq_along(bw_files), .combine = cbind) %dopar%
      {
        bw_i <- rtracklayer::import(bw_files[i], format = 'BW', as = "Rle")
        foreach (chr = sort(unique(seqnames(intervals))), .combine = c) %dopar%
        {
          tmp.chr <- sort(intervals[as.character(seqnames(intervals)) == chr])
          if (!is.null(blacklist))
          {
            blacklist.chr <- blacklist[as.character(seqnames(blacklist)) == chr]
            if (length(blacklist.chr) > 0)
            {
              bw_i[[chr]][ranges(blacklist.chr)] <- 0
            }
          }
          sum(Views(bw_i[[chr]], ranges(tmp.chr))) / width(tmp.chr)
         }
      }
    if (length(bw_files) == 1) count_dat <- as.matrix(count_dat, ncol = 1)
    count_dat <- count_dat[rank(intervals), ]
  }
  count_dat <- as.data.frame(count_dat)
  rownames(count_dat) <- names(intervals)
  colnames(count_dat) <- gsub(".*\\/(.*)\\.*", "\\1", bw_files)
  return(count_dat)
}

.countBam <- function(bam_files, intervals, stranded = FALSE, paired.end = "ignore")
{ # extract count reads from bam files
  count_dat <- foreach ( bam_path = bam_files, .combine = cbind) %dopar%
  {
    if(!file.exists(paste0(bam_path,'.bai'))) indexBam(bam_path)
    read_counts <- bamsignals::bamCount(bam_path,
                                        intervals,
                                        paired.end = paired.end,
                                        ss = stranded, # split strandness
                                        verbose = FALSE)
    if(stranded) read_counts <- read_counts[1, ]
    return(read_counts)
  } %>% as.data.frame
  colnames(count_dat) <- gsub(".*\\/(.*)\\.*", "\\1", bam_files)
  if(is.null(names(intervals))) {
    rownames(count_dat) <- seq_along(intervals)
  } else {
    rownames(count_dat) <- names(intervals)
  }
  return(count_dat)
}

.coverBam <- function(bam_files, intervals, paired.end = "ignore", tlenFilter = c(0, 1000),
                      is.resize = FALSE, df = 20, len = 100)
{ # extract base-wise coverage of bam files reads
  # paired.end != "ignore", only first read in pairs
  # paired.end == "extend", fragment is treated as a single read
  registerDoParallel(cores = 2)
  cov_list <- foreach (bam_path = bam_files) %dopar%
  {
    if(!file.exists(paste0(bam_path,'.bai'))) indexBam(bam_path)
    cov.list = bamsignals::bamCoverage(bam_path,
                                       gr = intervals,
                                       paired.end = paired.end,
                                       tlenFilter = tlenFilter,
                                       verbose = FALSE)
    if (is.resize) {
      resize_cov(cov.list, df = df, len = len)
    } else {
      cov.list
    }
  }
  names(cov_list) <- gsub(".*\\/(.*)\\.*", "\\1", bam_files)
  return(cov_list)
}

extend_range <- function(x, up_flank, down_flank)
{
  .strand = as.character(strand(x))
  start(x) = start(x) - ifelse(.strand %in% c("+", "*"), up_flank, down_flank)
  end(x) = end(x) + ifelse(.strand %in% c("-", "*"), up_flank, down_flank)
  x
}

sandwich_Bam <- function(bam_files, 
                         intervals, 
                         paired_end = FALSE, 
                         stranded = FALSE,
                         flanks = NULL, # e.g. flanks = c(up=2000, down=2000),
                         new_lens = c(up=20, mid=100, down=40))
{ # extract interval read coverage with flanking regions from bam files
  # save to a list object for each bam file
  if (is.null(flanks) & length(unique(width(intervals))) > 1) stop("Interval without resizing must be the same width.")
  registerDoParallel(cores = 5)
  
  intervals <- keepSeqlevels(intervals, as.character(unique(seqnames(intervals))))
  
  if (!is.null(flanks)) 
  {
    intervals <- extend_range(intervals, flanks[1], flanks[2])
  } else {
    flanks <- c(0, 0)
  }
  
  cov_dat <- foreach(bam_path = bam_files) %dopar%
  {
      bam.index=paste0(bam_path,'.bai')
      if(!file.exists(bam.index)) bam.index = Rsamtools::indexBam(files = bam_path)
      if(paired_end)
      {
        sbw = c('pos', 'qwidth','strand','rname', 'mrnm', 'mpos', 'isize')
        flag = scanBamFlag(isFirstMateRead = paired_end, isSecondaryAlignment = FALSE)
      } else {
        sbw = c('pos', 'qwidth','strand','rname')
        flag=scanBamFlag(isSecondaryAlignment=F)
      }
      param = ScanBamParam(what = sbw, flag = flag, which = intervals)
      srg = scanBam(bam_path, param=param, index = bam.index)
      
      .cov_mat <- matrix(0, nrow = length(intervals), ncol = sum(new_lens))
      
      for (n in seq_along(intervals))
      {
        .strand = as.character(strand(intervals[n]))
        if (paired_end) {
          if (stranded)
          {
            all.index <- with(srg[[n]],
                              !is.na(mpos) & mrnm == rname & abs(isize) < 2000 & strand == .strand)
          } else {
            all.index <- with(srg[[n]], !is.na(mpos) & mrnm == rname & abs(isize) < 2000 )
          }
        } else {
          if (stranded)
          {
            all.index <- with(srg[[n]], !is.na(rname) & strand == .strand)
          } else {
            all.index <- with(srg[[n]], !is.na(rname) )
          }
        }
        

        srg.tmp <- lapply(srg[[n]], function(x) x = x[all.index])
        .cov = with(srg.tmp, coverage(IRanges(start = pos - start(intervals[n]),
                                              width = qwidth ),
                                      width = width(intervals[n]) ))
        if ( .strand == '-') .cov = rev(.cov)
        
        if (sum(flanks) > 0)
        {
          .cov_resize = NULL
          if (flanks[1] > 0) .cov_resize = spline(x = seq_len(flanks[1]),
                                                  y = .cov[seq_len(flanks[1])],
                                                  n = new_lens[1],
                                                  method = 'natural')$y
          
          .cov_resize = c(.cov_resize, spline(x = seq_len(width(intervals[n]) - sum(flanks)),
                                              y = .cov[(sum(flanks[1]) + 1) : (width(intervals[n]) - sum(flanks[2]))],
                                              n = new_lens[2],
                                              method = 'natural')$y)
          
          if (flanks[2] > 0) .cov_resize = c(.cov_resize, spline(x = seq_len(sum(flanks[2])),
                                                                 y = .cov[(width(intervals[n]) - sum(flanks[2]) + 1) : width(intervals[n])],
                                                                 n = new_lens[3],
                                                                 method = 'natural')$y)
        }
        .cov_resize[.cov_resize < 0] = 0
        .cov_mat[n, ] = .cov_resize
      }
      rownames(.cov_mat) <- names(intervals)
      .cov_mat
    } # end of processing all bam files
  names(cov_dat) <- gsub(".*\\/(.*)\\.*", "\\1", bam_files)
  return(cov_dat)
}


sandwich_BW <- function(bw_files, 
                        intervals,
                        flanks = NULL, # e.g. flanks = c(up=2000, down=2000),
                        new_lens = c(up=100, mid=50, down=200))
{ # extract interval read coverage with flanking regions from bw files
  # save to a list object for each bw file
  
  if (is.null(flanks) & length(unique(width(intervals))) > 1) 
    stop("Interval without resizing must be the same width.")
  
  registerDoParallel(cores = 2)
  
  intervals <- keepSeqlevels(intervals, as.character(unique(seqnames(intervals))))
  
  if (!is.null(flanks)) 
  {
    intervals <- extend_range(intervals, flanks[1], flanks[2])
  } else {
    flanks <- c(0, 0)
  }
  
  resize_flank <- function(.cov, flanks, new_lens = c(10, 50, 10)) {
    
    if (any(is.na(.cov))) .cov[is.na(.cov)] <- 0 # remove NAs
    
    if (sum(.cov) == 0) return(rep(0, sum(new_lens)))
    
    interval_width <- length(.cov)
    
     c(spline(x = seq_len(flanks[1]),
              y = .cov[seq_len(flanks[1])],
              n = new_lens[1],
              method = 'natural')$y,
       
       spline(x = seq_len(interval_width - sum(flanks)),
              y = .cov[(sum(flanks[1]) + 1) : (interval_width - sum(flanks[2]))],
              n = new_lens[2],
              method = 'natural')$y,
       
       spline(x = seq_len(sum(flanks[2])),
              y = .cov[(interval_width - sum(flanks[2]) + 1) : interval_width],
              n = new_lens[3],
              method = 'natural')$y)
    
  }
  
  cov_mat_list <- foreach(bw_path = bw_files) %dopar% {
    .cov <- rtracklayer::import.bw(bw_path, as = "RleList")
    
    chrs <- intersect(unique(seqnames(intervals)), names(.cov))
    
    .cov_mat = NULL
    for (chr in chrs) {
      tmp_intervals <- intervals[seqnames(intervals) == chr]
      tmp_cov <- .cov[[chr]]
      chr.end <- length(tmp_cov)
      min_start <- min(start(tmp_intervals))
      
      # padding ranges outside chromosome boundary
      if (min_start < 1) {
        tmp_cov <- c(Rle(values = 0, lengths = abs(min_start) + 1), tmp_cov)
        end(tmp_intervals) <- end(tmp_intervals) + abs(min_start) + 1
        start(tmp_intervals) <- start(tmp_intervals) + abs(min_start) + 1
      }
      if (max(end(tmp_intervals)) > chr.end) 
        tmp_cov <- c(tmp_cov, Rle(values = 0, lengths = max(end(tmp_intervals)) - chr.end) )
      
      # get interval coverage
      tmp_cov <- as.vector(Views(tmp_cov, ranges(tmp_intervals)))
      tmp_strands <- as.character(strand(tmp_intervals))
      tmp_cov[tmp_strands == "-"] <- sapply(tmp_cov[tmp_strands == "-"], rev)
      
      .cov_resize <- matrix(0, nrow = length(tmp_intervals), ncol = sum(new_lens))
      for (i in seq_along(tmp_cov)) {
        .cov_resize[i, ] <- resize_flank(.cov = tmp_cov[[i]], 
                                         flanks = flanks, 
                                         new_lens = new_lens)
      }
      .cov_resize[.cov_resize < 0 | is.na(.cov_resize)] <- 0
      .cov_mat <- rbind(.cov_mat, `rownames<-`(.cov_resize, names(tmp_cov)))
    }
    
    message(bw_path)
    return(.cov_mat)
  } # end of processing all bw files
  
  names(cov_mat_list) <- gsub(".*\\/(.*)\\.*", "\\1", bw_files)
  return(cov_mat_list)
}


pile_interval_coverage <- function(query.list, intervals, out_width, is.strand = TRUE)
{ # this function converts GRange objects in query.list to coverage matrix upon subject intervals
  # e.g. find motif coverage on intervals
  cov.list <- list()
  for (i in seq_along(query.list))
  {
    query <- query.list[[i]]
    if (is.strand) {
      cov.mat.p <- matrix(0, nrow = length(intervals), ncol = out_width)
      cov.mat.m <- matrix(0, nrow = length(intervals), ncol = out_width)
      
      mtch.p <- findOverlaps(query, intervals, ignore.strand = FALSE) # plus overlap
      levels(strand(query)) <- c('-', '+', '*') %>% as.factor
      mtch.m <- findOverlaps(query, intervals, ignore.strand = FALSE) # minus overlap
      
      s.p <- split(queryHits(mtch.p), list(subjectHits(mtch.p)))
      s.m <- split(queryHits(mtch.m), list(subjectHits(mtch.m)))
      
      for (i in names(s.p))
      {
        i.n <- as.numeric(i)
        start.p <- start(query[ s.p[[i]] ]) - start(intervals[i.n])
        width.p <- width(query[ s.p[[i]] ])
        
        cov.mat.p[i.n, ] <- coverage( IRanges( start = start.p, width = width.p ),
                                      width = width(intervals[i.n]) ) %>% spline(n=out_width) %>% '$'(y)
      }
      
      for ( i in names(s.m))
      {
        i.n <- as.numeric(i)
        start.m <- start(query[ s.m[[i]] ]) - start(intervals[i.n])
        width.m <- width(query[ s.m[[i]] ])
        
        cov.mat.m[i.n, ] <- coverage( IRanges( start = start.m, width = width.m ),
                                      width = width(intervals[i.n]) ) %>% spline(n=out_width) %>% '$'(y)
      }
      
      minus.intevals <- as.character(strand(intervals)) == '-'
      cov.mat.p[minus.intevals, ] <- rev(cov.mat.p[minus.intevals, ])
      cov.mat.m[minus.intevals, ] <- rev(cov.mat.m[minus.intevals, ])
      
      cov.list <- c(cov.list, list(list('plus' = cov.mat.p, 'minus' = cov.mat.m)))
    } else { # if not stranded
      message("Calculate without strand.")
      cov.mat <- matrix(0, nrow = length(intervals), ncol = out_width)
      mtch <- findOverlaps(query, intervals, ignore.strand = TRUE)
      
      m.split <- split(queryHits(mtch), list(subjectHits(mtch)))
      for (i in names(m.split))
      {
        i.n <- as.numeric(i)
        .start <- start(query[ m.split[[i]] ]) - start(intervals[i.n])
        .width <- width(query[ m.split[[i]] ])
        
        cov.mat[i.n, ] <- coverage( IRanges( start = .start, width = .width ),
                                    width = width(intervals[i.n]) ) %>% spline(n=out_width) %>% '$'(y)
      }
      minus.intevals <- as.character(strand(intervals)) == '-'
      cov.mat[minus.intevals, ] <- rev(cov.mat[minus.intevals, ])
      cov.list <- c(cov.list, list(list('query_unstranded' = cov.mat)))
    }
  }
  names(cov.list) <- names(query.list)
  return(cov.list)
}

view_coverage <- function(bam.file.iist, 
                          bam.sample.names = NULL,
                          interval, 
                          stranded = TRUE,
                          is.fragment = FALSE,
                          size.factor.list = NULL, 
                          log_scale = TRUE, 
                          low_cut = 2,
                          smoothen = FALSE, 
                          df = 300,
                          bin_width = 1,
                          cov_color_list = NULL,
                          annotation_list = NULL, 
                          anno_box_height = 2, 
                          anno_text_cex = 1.5, 
                          anno_text_offset = -1) {
  
  # Initialize parameters
  if (is.null(names(bam.file.iist))) 
    names(bam.file.iist) <- seq_along(bam.file.iist)
  if (is.null(bam.sample.names)) 
    bam.sample.names <- names(bam.file.iist)
  if (is.null(size.factor.list)) 
    size.factor.list <- lapply(bam.file.iist, function(x) rep(1, length(x)))
  if (any(names(size.factor.list) != names(bam.file.iist))) 
    names(size.factor.list) <- names(bam.file.iist)
  if (is.null(cov_color_list)) 
    cov_color_list <- lapply(bam.file.iist, function(x) c("plus"='#bd5734', "minus"='#feb236'))
  
  combine_rep <- function(cov.list, size.factors, log_scale, 
                          smoothen = FALSE, df = 300) { 
    # get the mean normalised coverage
    out <- Reduce("+", sapply(seq_along(cov.list), 
                              function(i) cov.list[[i]] / size.factors[i]))
    out <- as.numeric(out / length(cov.list))
    if (smoothen) out <- smooth.spline(seq_along(out), y = out, df = df)$y
    out[out < 0] <- 0
    if (log_scale) {
      out[out < low_cut] <- 1
      return(log10(out))
    }
    out
  }
  
  # Initialize an empty list
  coverage_list <- vector(mode = "list", length = length(bam.file.iist))
  names(coverage_list) <- names(bam.file.iist)
  
  # Get coverage and combine replicates
  for (sample in names(bam.file.iist)) 
  {
    bam_files = bam.file.iist[[sample]]
    plus.cov = 
      combine_rep(
        cov.list = unlist(lapply(bam_files, get_reads_coverage, object.gr = interval, flank.size=0, is.frag = is.fragment)),
        size.factors = size.factor.list[[sample]], 
        log_scale = log_scale, 
        smoothen = smoothen, 
        df = df)
    if (stranded) 
    {
      levels(strand(interval)) = c("-", "+", "*")
      minus.cov = 
        combine_rep(
          cov.list = unlist(lapply(bam_files, get_reads_coverage, object.gr = interval, flank.size=0, is.frag = is.fragment)),
          size.factors = size.factor.list[[sample]], 
          log_scale = log_scale, 
          smoothen = smoothen, 
          df = df)
      minus.cov = rev(minus.cov)
      levels(strand(interval)) = c("-", "+", "*")
      # save interval coverage into list for plotting
      coverage_list[[sample]] = list("plus" = bin_sums(plus.cov, bin_width = bin_width) / bin_width, 
                                     "minus" = bin_sums(minus.cov, bin_width = bin_width) / bin_width)
    } else {
      coverage_list[[sample]] = list("plus" = bin_sums(plus.cov, bin_width = bin_width) / bin_width)
    }
  }
  
  # Plot start
  par(mfrow = c(length(coverage_list) + length(annotation_list), 1))
  par(mar = c(1,5,1,1))
  
  x_lims = c(0, length(coverage_list[[1]][[1]]))
  if (stranded) {
    max(unlist(lapply(coverage_list, function(x) max(x[["minus"]]))))
    y_lims = c(-max(unlist(lapply(coverage_list, function(x) max(x[["minus"]])))),
               max(unlist(lapply(coverage_list, function(x) max(x[["plus"]])))))
  } else {
    y_lims = c(0, max(unlist(lapply(coverage_list, function(x) max(x[["plus"]])))))
  }
  y_off_set = y_lims[2] / 30
  
  for (sample in names(bam.file.iist)) 
  {
    plot(x_lims, y_lims, type = "n",bty="n", xaxt='n', xlab="", yaxt = "n", ylab = sample, cex.lab = 1.2)
    axis(side = 2, lwd=1, cex.axis=1)
    barplot(height = coverage_list[[sample]][["plus"]], col=cov_color_list[[sample]][1],
            add = TRUE, axes = FALSE, space = 0, border = NA, offset = y_off_set)
    if (stranded)
      barplot(height = -coverage_list[[sample]][["minus"]], col=cov_color_list[[sample]][2], 
              add = TRUE, axes = FALSE, space = 0, border = NA, offset = -y_off_set)
  }
  
  for (anno_name in names(annotation_list)) {
    # Extract annotations in the plotting window
    if (length(annotation_list[[anno_name]] == 0)) next()
    anno = annotation_list[[anno_name]]
    anno = anno[countQueryHits(findOverlaps(anno, interval, ignore.strand = TRUE)) > 0]
    
    if (anno_name == "Gene") {
      exons = reduce(anno[anno$type == "exon"]) 
      genes = anno[anno$type == "gene"]
      genes.plus.anno = genes[strand(genes) == "+"]
      genes.minus.anno = genes[strand(genes) == "-"]
      
      exon.plus = intersect(IRanges(start = 1, width = width(interval)), 
                            shift(ranges(exons[as.character(strand(exons)) == "+"]),
                                  shift = -start(interval)))
      exon.minus = intersect(IRanges(start = 1, width = width(interval)), 
                             shift(ranges(exons[as.character(strand(exons)) == "-"]),
                                   shift = -start(interval)))
      genes.plus = intersect(IRanges(start = 1, width = width(interval)), 
                             shift(ranges(genes[as.character(strand(genes)) == "+"]),
                                   shift = -start(interval)))
      genes.minus = intersect(IRanges(start = 1, width = width(interval)), 
                              shift(ranges(genes[as.character(strand(genes)) == "-"]),
                                    shift = -start(interval)))
      
      # gene ranges
      start.plus = start(genes.plus)
      start.minus = end(genes.minus)
      mid.plus = start(resize(genes.plus, fix = "center", 1))
      mid.minus = start(resize(genes.minus, fix = "center", 1))
      
      plot(NULL, xlim = c(1, width(interval)), ylim = c(-14, 14), 
           xaxt="n", yaxt="n", bty="n",pch="",ylab="",xlab="", main="", sub="")
      
      if (length(genes.plus) > 0) {
        rect(xleft = start(exon.plus), ybottom = 4, 
             xright = end(exon.plus), ytop = 7,
             col = '#0000B2', border = NA)
        rect(xleft = start(genes.plus), ybottom = 5.25, 
             xright = end(genes.plus), ytop = 5.75,
             col = '#0000B2', border = NA)
      }
      if (length(genes.minus) > 0) {
        rect(xleft = start(exon.minus), ybottom = -7, 
             xright = end(exon.minus), ytop = -4,
             col = '#0000B2', border = NA)
        rect(xleft = start(genes.minus), ybottom = -5.75, 
             xright = end(genes.minus), ytop = -5.25,
             col = '#0000B2', border = NA)
      }
      
      # gene names
      label.plus = genes.plus.anno$gene_name
      label.minus = genes.minus.anno$gene_name
      
      # remove embeded genes
      if (length(genes.plus) < sum(strand(genes) == "+")) 
        label.plus = genes.plus.anno$gene_name[countQueryHits(findOverlaps(genes.plus.anno, genes.plus.anno, type = "within")) < 2]
      if (length(genes.minus) < sum(strand(genes) == "-")) 
        label.minus = genes.plus.anno$gene_name[countQueryHits(findOverlaps(genes.minus.anno, genes.minus.anno, type = "within")) < 2]
      
    } else if (anno_name == "TU") {
      
      anno.plus = anno[strand(anno) == "+"]
      anno.minus = anno[strand(anno) == "-"]
      
      range.plus = intersect(IRanges(start = 1, width = width(interval)), 
                             shift(ranges(anno.plus), shift = -start(interval)))
      range.minus = intersect(IRanges(start = 1, width = width(interval)), 
                              shift(ranges(anno.minus), shift = -start(interval)))
      
      if (length(range.plus) < length(anno.plus)) anno.plus = anno.plus[anno.plus$location %ni% c("dsRNA", "usRNA")]
      if (length(range.minus) < length(anno.minus)) anno.minus = anno.minus[anno.minus$location %ni% c("dsRNA", "usRNA")]
      
      # TU ranges
      start.plus = start(range.plus)
      start.minus = end(range.minus)
      mid.plus = start(resize(range.plus, fix = "center", 1))
      mid.minus = start(resize(range.minus, fix = "center", 1))
      
      # TU rectangular
      plot(NULL, xlim = c(1, width(interval)), ylim = c(-14, 14), 
           xaxt="n", yaxt="n", bty="n",pch="",ylab="",xlab="", main="", sub="")
      
      if (length(range.plus) > 0) 
        rect(xleft = start(range.plus), ybottom = 4.5, 
             xright = end(range.plus), ytop = 6.5,
             col = '#800800', border = NA)
      
      if (length(range.minus) > 0) 
        rect(xleft = start(range.minus), ybottom = -6.5, 
             xright = end(range.minus), ytop = -4.5,
             col = '#800800', border = NA)
      
      # TU names
      label.plus = anno.plus$location
      label.minus = anno.minus$location
    }
    
    # TSS arrow and TU name
    if (length(label.plus) > 0) 
    {
      text(x = mid.plus, y = 2.5 - anno_text_offset, 
           labels = label.plus, cex = anno_text_cex)
      arrows(x0 = start.plus, y0 = 8, 
             x1 = start.plus + width(interval) / 75, y1 = 8, 
             code = 2, angle = 10, cex = 0.1, length = 0.05)
      arrows(x0 = start.plus, y0 = 5.5, 
             x1 = start.plus, y1 = 7.9, length = 0)
    }
    
    if (length(label.minus) > 0) 
    {
      text(x = mid.minus, y = -8.5 - anno_text_offset,
           labels = label.minus, cex = anno_text_cex)
      arrows(x0 = start.minus, y0 = -8, 
             x1 = start.minus - width(interval) / 75, y1 = -8, 
             code = 2, angle = 10, cex = 0.1, length = 0.05)
      arrows(x0 = start.minus, y0 = -5.5, 
             x1 = start.minus, y1 = -7.9, length = 0)
    }
  }
}



insertion_size_coverage <- function(bam.files, object.gr, is.matrix = FALSE, is_CI = FALSE)
{
  # bin pair-end reads insertion size coverage from bam files
  # all input bam files from a certain sample will be merged
  # output a list of Rle strings of each intervel
  # Agrs:
  #     bam.file: paired end strand specific bam file path
  #     object.gr: multiple intervals of gene or transcript
  
  combind.list = function(a, b)
  {
    out = lapply(names(a), function(x) c(a[[x]], b[[x]]))
    names(out) = names(a)
    out
  }
  
  for (i in seq_along(bam.files))
  {
    file.name = gsub('*.bam','\\', strsplit(bam.files[i],'\\/') %>% unlist %>% tail(1))
    file.dir = gsub(paste0('*', file.name, '.bam'), '\\', bam.files[i])
    bam.index = paste0(file.dir, file.name, '.bam.bai')
    if (!file.exists(bam.index)) bam.index = Rsamtools::indexBam(files = bam.files[i])
    
    # scanBamWhat: the info that need to be extracted from a bam file.
    sbw = c('pos', 'qwidth', 'mapq', 'strand', 'mpos', 'isize')
    sbp = ScanBamParam(what = sbw, 
                       which = object.gr,
                       flag = scanBamFlag(isProperPair = TRUE,
                                          isUnmappedQuery = FALSE,
                                          isSecondaryAlignment = FALSE,
                                          isFirstMateRead = TRUE)) # use first in pair reads for coverage
    
    # Scan bam file to retrieve short reads.
    sr.in.ranges = scanBam(bam.files[i], param = sbp, index = bam.index)
    
    srg.temp = foreach(n = seq_along(object.gr)) %dopar%
      {
        # Filter short reads by mapping quality and insertion size.
        all.index <- with(sr.in.ranges[[n]], 
                          mapq > 3 & !is.na(isize) & abs(isize) < 1000 & abs(isize) > 1) 
        lapply(sr.in.ranges[[n]], function(x) x = x[all.index] )
      }
    names(srg.temp) = names(sr.in.ranges)
    
    # combine multiple bam files
    if (!exists('srg.filtered'))
    {
      srg.filtered = srg.temp
    } else {
      srg.filtered = lapply(names(srg.temp),
                            function(x) 
                              combind.list(srg.temp[[x]], srg.filtered[[x]]) 
      )
    }
    
    scanBamRevOrder <- function(org.gr, sbp) {
      # ScanBamParam re-arranges the input genomic ranges. Use range info to
      # construct a string vector to find the order to reverse it.
      org.grnames <- with(org.gr, paste(seqnames, start, end, sep=':'))
      sbw.gr <- as.data.frame(bamWhich(sbp))  # scan-bam-ed
      if('space' %in% names(sbw.gr)) {
        sbw.grnames <- with(sbw.gr, paste(space, start, end, sep=':'))
      } else if('group_name' %in% names(sbw.gr)) {
        sbw.grnames <- with(sbw.gr, paste(group_name, start, end, sep=':'))
      } else {
        stop("Cannot locate chromosome names in extracted short reads. Report
             this problem using issue tracking or discussion forum.\n")
      }
      match(org.grnames, sbw.grnames)
    }
    
    # Restore the original order.
    srg.filtered <- srg.filtered[scanBamRevOrder(as.data.frame(object.gr), sbp)]
    
  }
  
  if (is_CI) {
    out_CI = NULL
    for (n in seq_along(object.gr)) {
      
      if (length(srg.filtered[[n]][[1]]) < 10) {
        # continue if there are reads on the current interval
        out_CI = c(out_CI, NA)
      } else {
        read_widths = abs(srg.filtered[[n]]$isize)
        out_CI = c(out_CI, sum(read_widths > 120 & read_widths < 220) / sum(read_widths > 20 & read_widths < 80))
      }
    }
    return(out_CI)
  }
  
  if (is.matrix) {
    
    n.pos = width(object.gr[1]) 
    read_mat = matrix(0, nrow = 200, ncol = n.pos)
    
    for (n in seq_along(object.gr)) {
      mat = matrix(0, nrow = 200, ncol = n.pos)
      
      if (length(srg.filtered[[n]][[1]]) > 0) {
        # continue if there are reads on the current interval
        .strand = as.character(strand(object.gr[n]))
        read_starts = with(srg.filtered[[n]], ifelse(strand == "+", pos, mpos)) - start(object.gr[n])
        if (.strand == "-") {
          read_starts = width(object.gr[n]) - abs(srg.filtered[[n]]$isize) - read_starts
        } 
        
        read_widths = abs(srg.filtered[[n]]$isize)
        
        for (i in seq_along(read_starts)) {
          .width = read_widths[i]
          if (.width <= 220 & .width > 20) {
            .range = read_starts[i] + seq_len(.width) - 1
            .range = .range[.range > 0 & .range <= n.pos]
            mat[.width - 20, .range] = 1
          }
        }
      } 
      read_mat = read_mat + mat
    }
    
    return(read_mat)
    
  } else {
    out = foreach(n = seq_along(object.gr)) %dopar% {
      n.pos = width(object.gr[n]) 
      cov = rep(NA, n.pos)
      if (length(srg.filtered[[n]][[1]]) > 0) {
        # continue if there are reads on the current interval
        .strand = as.character(strand(object.gr[n]))
        
        read_starts = with(srg.filtered[[n]], ifelse(strand == "+", pos, mpos)) - start(object.gr[n])
        read_gr = IRanges(start = read_starts, width = abs(srg.filtered[[n]]$isize))
        read_cov <- coverage(read_gr, weight = 1, width = width(object.gr[n]))
        intsertion_cov <- coverage(read_gr, weight = abs(srg.filtered[[n]]$isize), width = width(object.gr[n]))
        cov <- intsertion_cov / read_cov
        if (.strand == "-") cov = rev(cov) 
      }
      cov
    }
    return(out)
  }
}

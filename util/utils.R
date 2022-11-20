# utilities
pacman::p_load(BiocManager, rtracklayer, Rsamtools,
               GenomicRanges, IRanges,
               org.Mm.eg.db, biomaRt, AnnotationFilter,
               foreach, doParallel, 
               dplyr, tibble, matrixStats,
               ggplot2, grid, gridExtra, viridis, RColorBrewer)

registerDoParallel(cores = 4)

# ------------------------------------- ChIP target colors ------------------------------------- #
sample_colors <- c("H2Aub" = colors_20[16], 
                   "H3K27me3" = colors_20[6],
                   "H3K27ac" = colors_20[3],
                   "H3K4me3" =  colors_20[2],
                   "Ring1b" =  colors_20[15],
                   "Ezh2" =  colors_20[5],
                   "Pol II" =  colors_20[10],
                   "Pol II-S5p" =  colors_20[8], 
                   "Pol II-S2p" =  colors_20[11],
                   "Cbx7" = colors_9[1],
                   "Rybp" = colors_9[3],
                   "TT-seq" = colors_9[9])

# Input data --------------------------------------------------------------------------------
importRanges <- function(data, seq_lvls = paste0('chr',c(1:21,'X','Y')))
{
  format <- toupper(sapply(strsplit(data,'\\.'), tail, 1))
  import_input <- rtracklayer::import(data, format = format)
  seqlevelsStyle(import_input) <- "UCSC"
  import_input <- import_input[ seqnames(import_input) %in% seq_lvls ]
  seqlevels(import_input) <- as.character(unique(seqnames(import_input)))
  return(import_input)
}

readKallistoResults <- function( file_paths, sample_names )
{
  txCounts <- NULL
  for(i in file_paths)
  {
    temp.table <- read.table(file = i, header=T, colClasses=c("character", "numeric", "numeric", "numeric", "numeric"))
    txCounts <- cbind(txCounts, temp.table[, 4])
  }
  rownames(txCounts) <- unname(temp.table$target_id)
  colnames(txCounts) <- sample_names[seq_along(file_paths)]
  txCounts <- keepOneTx(txCounts) # the max in tx variants
  return(txCounts)
}

keepOneTx <- function(count_table, rowname_gene_id = F, is_gene_sum = F)
{ # process kallisto counts on genecode fasta reference
  # save only one count for each protein_coding gene and lincRNA 

  count_table <- count_table[rowSums(count_table) > 0, ]
  trim_type <- function(x) gsub("(.*\\|)(.*)(\\|$)", "\\2", x)
  gene_types <- rownames(count_table) %>% trim_type
  count_table <- count_table[gene_types %in% c('protein_coding', 'lincRNA'), ]
  out_tx_id <- T
  trim_gene <- function(x) lapply(strsplit(x, '\\|'),
                                  function(y) gsub("*(\\..*)",
                                                   "\\2",
                                                   y[ifelse(out_tx_id, 1, 2)])
                                  ) %>% unlist
  tx_ids <- rownames(count_table) %>% trim_gene
  out_tx_id <- F
  gene_ids <- rownames(count_table) %>% trim_gene
  idx <- split(tx_ids, list(gene_ids))

  rownames(count_table) <- tx_ids

  row_sums <- rowSums(count_table)
  doParallel::registerDoParallel(cores = 8)

  new_table <- foreach(i = idx, .combine = rbind) %dopar%
    {
      if(length(i) == 1)
      {
        matrix(count_table[i, ],
               nrow = 1,
               dimnames = list(i, colnames(count_table)))
      } else {
        tmp <- data.frame(count_table[i, ])
        if (is_gene_sum) {
          matrix(colSums(tmp),
                 nrow = 1,
                 dimnames = list(i[1], colnames(count_table)))
        } else {
          as.matrix(tmp[which.max(row_sums[i]), ])
        }
      }
    }

  colnames(new_table) <- colnames(count_table)
  if (rowname_gene_id) rownames(new_table) <- names(idx)

  return(new_table)
}

# stats tools --------------------------------------------------------------------------

size_factor_cal <- function(gene.counts, ref_expr = NULL)
{
  # estimate relative size factor comparing to reference TU gene counts
  # Args:
  # gene.counts: genes by row and samples by column
  # ref_expr: relative size factors to a reference count table
  gene.counts = gene.counts[!apply(gene.counts, 1, function(x) any(is.na(x)) ) & rowSums(gene.counts) > 0,]
  if (!is.null(ref_expr)){
    gene.counts = cbind(ref_expr, gene.counts)
    # geometric mean
    geoMeans = apply(gene.counts, 1, function(x) exp(sum(log(x[x > 0]), na.rm = TRUE) / length(x)) )
    # quotients list to median (size factor)
    quoMedian = apply(sweep(gene.counts, 1, geoMeans,'/'), 2, median)
    return(quoMedian[-1] / quoMedian[1])
  } else{
    geoMeans=apply(gene.counts, 1, function(x) exp(sum(log(x[x > 0]), na.rm=T) / length(x)) )
    return(apply(sweep(gene.counts, 1, geoMeans,'/'), 2, median))
  }
}

rowMax <- function(mat) apply(mat, 1, max)

geoMeans <- function(mat)
{
  if (is.null( dim(mat) )) {
    mat
  } else {
    exp( log(rowProds(mat)) / ncol(mat))
  }
}

max_min_norm <- function(x) {
  # perform max-min normalization, transform variable to the range (0, 1)
  if (!is.null(dim(x))) {
    cbind(max_min_norm(x[, 1]), max_min_norm(x[, -1]))
  } else {
    x[is.infinite(x)] = NA
    .min = min(x, na.rm = T)
    return((x - .min) / (max(x, na.rm = T) - .min) )
  }
}

groupConsec <- function( seq ) # -> order
{ # input a sequence of classes
  # sum up consecutive numbers
  vals <- seq[1]
  val_group <- 1
  `%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))
  for ( i in seq_along(seq)[-1] )
  {
    if (seq[i-1] == seq[i])
    {
      val_group[length(val_group)] %+=% 1
    } else
    {
      vals <- c( vals, seq[i] )
      val_group <- c(val_group, 1)
    }
  }
  names(val_group) <- vals
  return(val_group)
}

kink_index <- function(x, method = "slope", .p = 0.05)
{
  # this function returns outlier indexes of the input vector
  # to retrieve highly enriched loci, e.g. HOMER uses slope = 1 to get super enhancers
  
  x_ord <- order(x)
  
  if (method == "slope")
  {
    x[is.na(x)] <- 0
    slopes_above_step <- diff(x[x_ord]) > 1
    stopifnot(sum(slopes_above_step) > 0)
    return(x_ord[seq(min(which(slopes_above_step)) + 1, length(x))])
  }
  
  if (method == "weight")
  {
    x[is.na(x)] <- 0
    x[is.infinite(x)] <- max(x[!is.infinite(x)])
    weight_above_mean <- which.min(cumsum(diff(x[x_ord]) - mean(diff(x[x_ord]))))
    return(x_ord[seq(weight_above_mean, length(x))])
  }
  
  if (method == "poisson")
  {
    lambda <- mean(x[is.finite(x)], na.rm = TRUE)
    scaling <- 10 / lambda # remove sample normalization
    x <- x * scaling
    lambda <- lambda * scaling
    p <- ppois(x, lambda, lower.tail = FALSE)
    return(which(p < .p))
  }
}

single_variance_explained <- function(X, Y, is.cor = FALSE)
{ 
  stopifnot(length(X) == length(Y))
  rm.idx = is.na(X) | is.infinite(X) | is.na(Y) | is.infinite(Y)
  X = X[!rm.idx]
  Y = Y[!rm.idx]
  # output correlation
  if (is.cor) return(cor(X, Y))
  # or reduced variance: 1 - covariance / variance 
  1 - var(Y - X) / var(Y) 
}

multi_variance_explained <- function(X, Y, is.cor = F)
{ # r_squared decomposition of regression through origin
  # X: n-row matrix with k features, Y: n responses
  keep.idx = !apply(cbind(X, Y), 1, function(x) any(is.na(x) | is.infinite(x)))
  X = X[keep.idx, ]; Y = Y[keep.idx]
  stopifnot(var(Y) > 0 & length(Y) == nrow(X))
  if (is.cor) return(c(cor(Y, X)))
  X_tilde = cbind(1, scale(X)) # add an intercept error term ε
  Y_tilde = scale(Y)
  beta_hat = solve(t(X_tilde) %*% X_tilde) %*% t(X_tilde) %*% Y_tilde
  # let the error term ε=0
  # then R2 = sum(Y_hat^2) / sum(Y_tilde^2), where sum(Y_tilde^2) = length(Y) - 1
  # R2 decomposition is:
  c(beta_hat * t(cov(Y_tilde, X_tilde)))[-1]
}


trim_quantile <- function(x, q = 0.995, lower = TRUE, upper = TRUE) {
  # trim outliers with the indicated quantile
  .trim <- function(x, q = 0.995) {
    if (!is.null(dim(x))) {
      cbind(.trim(x[, 1], q), .trim(x[, -1], q))
    } else {
      if (upper)
        x[x >= quantile(x, q, na.rm = TRUE)] = quantile(x, q, na.rm = TRUE)
      if (lower)
        x[x <= quantile(x, 1-q, na.rm = TRUE)] = quantile(x, 1-q, na.rm = TRUE)
      return(x)
    }
  }
  `colnames<-`(.trim(x, q), colnames(x))
}


quantile_mean <- function(x, q = 0.01, lower = TRUE, upper = TRUE) {
  
  if (!is.null(dim(x))) {
    c(quantile_mean(x[, 1], q, lower, upper),
      quantile_mean(x[, -1], q, lower, upper))
  } else {
    min_x = quantile(x, q)
    max_x = quantile(x, 1 - q)
    
    if (lower)  x = x[x >= min_x]
    if (upper) x = x[x <= max_x]
    
    mean(x, na.rm = TRUE)
  }
}

.add_row <- function(mat, .after_rows, vals) {
  # this function adds new cases to indicated rows
  # mat:          matrix object
  # .after_rows:  multiple row indexes
  # vals:         matrix object in the same dimension as mat for the new cases
  if (length(.after_rows) > 1 && .after_rows[1] != .after_rows[2]) {
    # if consecutive cases
    consec <- .after_rows == .after_rows[1]
    mat_out <- rbind(mat[seq_len(.after_rows[1]), ], 
                     vals[consec, ],
                     .add_row(mat = mat[(.after_rows[1] + 1) : nrow(mat), ],
                              .after_rows = .after_rows[!consec] - .after_rows[1], 
                              vals = vals[!consec, ])
    )
  } else {
    mat_out <- rbind(mat[seq_len(.after_rows[1]), ], 
                     vals, 
                     mat[(.after_rows[1] + 1) : nrow(mat), ])
  }
  mat_out
}


exclude.Vector <- function(v1, v2) v1[v1 %ni% v2]

intersect_list <- function(.list) {
  if (length(.list) > 2) {
    tmp <- c(unlist(.list[1]), intersect_list(.list[-1]))
    tmp[duplicated(tmp)]
  } else if (length(.list) == 2) {
    unlist(.list)[duplicated(unlist(.list))]
  } else {
    unlist(.list)
  }
}

intersect_two_sets <- function(s1, s2) {
  # return number of s1_unique, common, s2_unique
  c(sum(!s1 %in% s2), sum(s1 %in% s2), sum(!s2 %in% s1))
}

intersect_three_sets <- function(s1, s2, s3, show_list = FALSE) {
  # return in the order
  # 1 0 0
  # 1 1 0
  # 1 0 1
  # 1 1 1
  # 0 1 0
  # 0 1 1
  # 0 0 1
  
  if (!show_list) {
    `names<-`(c(sum(!s1 %in% c(s2, s3)),
                sum(!s1[s1 %in% s2] %in% s3), 
                sum(!s1[s1 %in% s3] %in% s2),
                sum(s1[s1 %in% s2] %in% s3),
                sum(!s2 %in% c(s1, s3)),
                sum(!s2[s2 %in% s3] %in% s1), 
                sum(!s3 %in% c(s1, s2))
    ),
    c("100", "110", "101", "111", "010", "011", "001"))
  } else {
    s_list <- list(s1, s2, s3)
    sapply(c("100", "110", "101", "111", "010", "011", "001"), 
           function(x) {
             is_x <- as.logical(as.numeric(unlist(strsplit(x, ""))))
             exclude.Vector(intersect_list(s_list[is_x]), unlist(s_list[!is_x]))
           })
  }
}


# sequence analysis ------------------------------------------------------------------------

interval_pattern_hits <- function(intervals,
                                  pattern, 
                                  which_genome = 'mm10', 
                                  weight_pos = NA,
                                  to_coverage = F, 
                                  out_width = 200, ...)
{
  # convert interval to DNA sequence to count the pattern hits
  if (which_genome == 'mm10')
    genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  if (which_genome == 'mm9')
    genome <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9
  
  # remove ranges beyond genome 
  seq_lengths <- seqlengths(genome)
  for (chr in seqlevels(intervals))
  {
    chr_idx <- seqnames(intervals) == chr
    gaps_to_end <- end(intervals[chr_idx]) - seq_lengths[chr] 
    idx <- gaps_to_end > 0
    if (sum(idx) > 0) 
      end(intervals[chr_idx][idx]) = start(intervals[chr_idx][idx]) = seq_lengths[chr] 
  }
  seqs <- getSeq(genome, intervals)
  
  # output one hot matrix
  if (pattern == "one_hot")
    return(sapply(seqs, oneHot_nt))
  
  # output the number of matches
  if (!to_coverage & is.na(weight_pos)) 
    return(Biostrings::vcountPattern(pattern = pattern, subject = seqs))
  
  query.list = Biostrings::vmatchPattern(pattern = pattern, subject = seqs)
  if (!is.na(weight_pos)) {
    dists = mclapply(seq_along(query.list), 
                     function(n) 
                       sum(1 / (distance(query.list[[n]], IRanges(start = weight_pos[n], width = 0))^2 + 1))
    )
    return(unlist(dists))
  }
  # otherwise get pattern coverage
  cov.list = sapply(seq_along(query.list), 
                    function(n) coverage(query.list[[n]], 
                                         width = width(intervals[n])))
  cov.mat = Reduce("rbind", 
                   lapply(cov.list, 
                          function(y) spline(x = seq_along(y), 
                                             y = y,
                                             n = out_width)$y
                   )
  )
  return(cov.mat)
}

oneHot_nt <- function(dna) 
{
  nt <- c("A", "C", "G", "T")
  dna %>% as.character() %>% strsplit(split='') %>% unlist() -> dna
  mat <- matrix(0, nrow = 4, ncol = length(dna), dimnames = list(nt, seq_along(dna)))
  mtch <- match(dna, nt)
  mat[matrix(c(mtch, seq_along(dna)), ncol=2)] <- 1
  return(mat)
}

# plot utilities ----------------------------------------------------------------------------
plot_Vennerable <- function (list_1, list_2, 
                             name_1, name_2,
                             color_set = c(1, 2, 3), 
                             color_Text = c(1, 2), ...)
{
  pacman::p_load(GeneOverlap, RColorBrewer, grid)
  
  Sets = list(list_1, list_2)
  names(Sets) = c(name_1, name_2)
  p = Vennerable::compute.Venn(Vennerable::Venn(Sets = Sets))
  gp = Vennerable::VennThemes(p)
  gp$Face$`10`$fill = color_set[1]
  gp$Face$`11`$fill = color_set[2]
  gp$Face$`01`$fill = color_set[3]
  gp$FaceText$`10`$fontsize = 30
  gp$FaceText$`11`$fontsize = 30
  gp$FaceText$`01`$fontsize = 30
  gp$Set$Set1$lwd = 0
  gp$Set$Set2$lwd = 0
  gp$SetText$Set1$fontsize = 32
  gp$SetText$Set2$fontsize = 32
  gp$SetText$Set1$col = color_Text[1]
  gp$SetText$Set2$col = color_Text[2]
  # plot and add pval
  plot(p, gp = gp, show = list(Universe=F), ...)
  
  grid.text(
    paste("Jaccard =", round(length(intersect.Vector(list_1, list_2)) / 
                               length(unique(c(list_1, list_2))), 2)),
    x = 0.5, y=0.95, gp=gpar(cex=2))
}

get_dens <- function(X, Y, n.grid = 100) {
  # from https://slowkow.com/notes/ggplot2-color-by-density/
  dens <- MASS::kde2d(X, Y, n = n.grid)
  ix <- findInterval(X, dens$x)
  iy <- findInterval(Y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
  # X_cut <- cut(X, seq(min(X), max(X), diff(range(X)) / n.grid ) )
  # Y_cut <- cut(Y, seq(min(Y), max(Y), diff(range(Y)) / n.grid ) )
  # den <- as.numeric(sqrt(table(X_cut)[X_cut] * table(Y_cut)[Y_cut]))
}

plain <- function(x, ...) {
  format(x, ..., scientific = FALSE, drop0trailing = TRUE)
}

plot_ma <- function(baseMean,
                    log2FC,
                    pval,
                    cut_baseMean = 5, 
                    DE_log2FC = 1.5,
                    up_num = 1,
                    down_num = 1,
                    xlab = expression("log"[10] ~ "baseMean"),
                    ylab = expression("log"[2] ~ "fold-change"),
                    title = "", 
                    p_col = "grey50") {
  
  dat <- data.frame("x" = log10(baseMean), 
                    "y" = log2FC, #- median(log2FC, na.rm = T), 
                    "z" = pval) %>% 
    dplyr::filter(complete.cases(.) & is.finite(rowSums(.))) %>%
    dplyr::mutate("density" = get_dens(x, y, n.grid = 200)) 
  
  dat_DE <- dat[which(dat$x > log10(cut_baseMean) & abs(dat$y) > DE_log2FC & abs(dat$z) < 0.05), ]
  
  g <- ggplot(dat, aes(x = x, y = y, color = density)) +
    geom_point(size = 1) + 
    geom_hline(yintercept = c(-DE_log2FC, 0, DE_log2FC),
               lwd = c(0.2, 0.5, 0.2), lty = 2) +
    scale_color_gradient(low = add.alpha("grey50", 0.5),
                         high = p_col) +
    ylim(c(-10, 10)) +
    xlab(xlab) +
    ylab(ylab) +
    annotate("text", x = Inf, y = Inf, size = 5,
             label = paste("Up =", up_num),
             vjust = 1.3, hjust = 1.1) + 
    annotate("text", x = Inf, y = -Inf, size = 5,
             label = paste("Down =", down_num),
             vjust = -1.3, hjust = 1.1) +
    ggtitle(title) +
    theme_setting +
    theme(legend.position = "none")
  g + geom_point(data = dat_DE,
                 mapping = aes(x = x, y = y), 
                 color = add.alpha(ifelse(dat_DE$y > 0, "red", "blue"), alpha = 0.5),
                 size = 0.5, pch = 2)
}

add.alpha <- function(col, alpha=1)
{
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb) / 255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}

plot_umap <- function(mat) {
  require(ggplot2)
  require(umap)
  
  dat = `colnames<-`(as.data.frame(umap(mat)$layout), c("UMAP1", "UMAP2"))
  
  ggplot(dat, aes(x = UMAP1, y = UMAP2, color = get_dens(UMAP1, UMAP2))) +
    geom_point() +
    scale_color_viridis_c() +
    xlim(range(dat$UMAP1) * 1.3) + ylim(range(dat$UMAP2) * 1.3) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_text(hjust = 0, vjust = 0, size = 13), 
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
}

bin_smoothen <- function(x, bin_size = 6) {
  
  # this function smooths the input with mean values in a sliding window 
  # Args:
  # data: a vector, or a matrix smooth by column
  # bin_size: averaging window size
  
  if (!is.null(dim(x))) {
    
    cbind(bin_smoothen(x[, 1], bin_size), bin_smoothen(x[, -1], bin_size))
    
  } else {
    
    # padding the start and end
    x <- c(rep(NA, bin_size %/% 2), x, rep(NA, bin_size %/% 2 + bin_size %% 2))
    x_na <- is.na(x) # skip NAs for averaging
    x[is.na(x)] <- 0
    
    c_sums = cumsum(x)
    b_sums = c_sums[seq(bin_size, length(x) - 1)] - c(0, c_sums[seq(1, length(x) - bin_size - 1)])
    
    c_nas = cumsum(x_na)
    b_nas = c_nas[seq(bin_size, length(x) - 1)] - c(0, c_nas[seq(1, length(x) - bin_size - 1)])
    
    # return
    b_sums / c(rep(bin_size, length(x) - bin_size) - b_nas)
  }
}


plot_freq_intersect <- function(dat, .by = "Group", .levels = NA, .split = "Category", .color = 1:10, top_n = 10) {
  # Args:
  #   dat: gene_id (unique)
  #   .by: row feature
  #   .levels: dot plot row names 
  #   .split: bar split feature
  
  require(ggplot2)
  require(cowplot)
  
  if (is.na(.levels[1])) {
    .levels <- seq_len(as.character(nchar(dat[1, .by])))
  }
  
  # limit less frequent groups
  top_groups <- names(tail(sort(table(dat[, .by])), top_n))
  dat <- dat[dat[[.by]] %in% top_groups, ]
  
  dat$Group <- factor(dat[, .by], levels = names(sort(table(dat[, .by]), decreasing = T)))
  dat$Type <- dat[, .split]
  
  # barplot
  dat_g1 <- dplyr::count(dat, Group, Type, sort = TRUE)
  
  frac_tbl <- table(dat[, "Group"], dat[, "Type"])
  frac_tbl <- frac_tbl / rowSums(frac_tbl)
  
  dat_g1_text <- data.frame(Group = rownames(frac_tbl),
                            Fraction = paste0(round(frac_tbl[, 1], 2) * 100, "%"),
                            Type = colnames(frac_tbl)[1])
  
  Group_n <- dplyr::count(dat, Group, sort = TRUE)
  dat_g1_text$n  <- Group_n[match(dat_g1_text$Group, Group_n$Group), "n"]
  
  g1 <- ggplot(dat_g1, aes(x = Group, y = n, fill = Type)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(data = dat_g1_text, 
              aes(x = Group, y = n, label = Fraction),
              vjust = -0.25, hjust = 0.5,
              size = 3) +
    xlab("") + ylab("Number of genes") +
    ggpubr::theme_pubclean() +
    scale_fill_manual(name = .split, values = c("grey70", .color)) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.margin = margin(c(0,1,1,1)))
  
  # dot plot
  dat_tile <- NULL
  for (m in seq_along(levels(dat$Group))) {
    tmp <- cbind(x = m, y = .levels, 
                 color = ifelse(strsplit(levels(dat$Group)[m], "")[[1]] == 1, 1, 0))
    dat_tile <- rbind(dat_tile, tmp)
  }
  
  dat_tile <- as.data.frame(dat_tile) 
  dat_tile$x <- factor(dat_tile$x, levels = order(as.numeric(dat_tile$x)))
  dat_tile$y <- factor(dat_tile$y, levels = rev(.levels))
  
  g2 <- ggplot(dat_tile, aes(x = x, y = y, color = factor(color), group = x)) +
    geom_point(size = 4) +
    scale_color_manual(values = c("grey90", .color)) +
    xlab("") + ylab("") +
    ggpubr::theme_pubclean() +
    theme(panel.grid = element_blank(), 
          legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.margin = margin(c(1,1,0,1)))
  
  # add lines
  dat_seg <- NULL
  for (i in unique(dat_tile$x)) {
    tmp <- dat_tile[dat_tile$x == i, ]
    tmp <- tmp[tmp$color == 1, ]
    if (nrow(tmp) > 1) {
      dat_seg <- rbind(dat_seg,
                       data.frame(x = tmp$x[-nrow(tmp)], y = tmp$y[-nrow(tmp)],
                                  xend = tmp$x[-1], yend = tmp$y[-1]))
    }
  }
  
  if (!is.null(dat_seg)) {
    levels(dat_seg$x) <- levels(dat_seg$xend) <- levels(dat_tile$x)
    levels(dat_seg$y) <- levels(dat_seg$yend) <- levels(dat_tile$y)
    dat_seg$color <- 1
    
    g2 <- g2 + geom_segment(data = dat_seg, aes(x = x, y = y, xend = xend, yend = yend), lty = 2)
  }
  
  cowplot::plot_grid(g1, g2, ncol = 1, align = "v", rel_heights = c(1, 0.4))
}

# graphic settings ---------------------------------------------------------------------------
theme_setting <- list(theme_minimal(),
                      theme(legend.position = 'right',
                            legend.justification = c(0, 1),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            # panel.border = element_rect(colour = "black", fill=NA, size=1),
                            plot.title = element_text(size = 14, face="bold", vjust = 0.5),
                            axis.text=element_text(size=12, face = "bold"),
                            axis.title=element_text(size=15,face="plain"),
                            legend.title = element_text(size=14, face="bold"),
                            legend.text = element_text(size=12, face="plain"),
                            axis.ticks.x = element_line(size = 1), 
                            axis.ticks.y = element_line(size = 1),
                            # axis.ticks.length = unit(20, "pt"),
                            panel.border = element_rect(colour = "black", fill=NA, size=1.2))
                      )

colors_9 <- c('#59C7EB', '#FEA090', '#9AA0A7', '#077187', '#0A9086', '#3E5496', '#E0607E', '#8E2043', '#EFDC60')

colors_20 <- c('#bf4040', '#cc8800', '#808000', '#408000', '#006622', '#2d8659', '#5c8a8a',
               '#0073e6', '#4d4dff', '#5353c6', '#4d0099', '#660080', '#602060', '#c6538c',
               '#99003d', '#cc5200', '#b32400', '#663300', '#b3b300', '#4d9900')

colors_n <- c('#c90000', '#c94600', '#c99700', '#6f9c00', '#009c56', '#00838f',
              '#0550b3', '#3e0080', '#560659', '#adadad')


# extensions -------------------------------------------------------------------------------
`%+=%` = function(e1, e2) eval.parent(substitute(e1 <- e1 + e2))
`%-=%` = function(e1, e2) eval.parent(substitute(e1 <- e1 - e2))
`%*=%` = function(e1, e2) eval.parent(substitute(e1 <- e1 * e2))
`%/=%` = function(e1, e2) eval.parent(substitute(e1 <- e1 / e2))
`%c=%` = function(e1, e2) eval.parent(substitute(e1 <- c(e1, e2)))
`%ni%` = Negate(`%in%`)

# chromosome length ------------------------------------------------------------------------
seqlengths_mm9 <- c("chr1" =	197195432,
                    "chr2" =	181748087,
                    "chr3" =	159599783,
                    "chr4" =	155630120,
                    "chr5" =	152537259,
                    "chr6" =	149517037,
                    "chr7" =	152524553,
                    "chr8" =	131738871,
                    "chr9" =	124076172,
                    "chr10" =	129993255,
                    "chr11" =	121843856,
                    "chr12" =	121257530,
                    "chr13" =	120284312,
                    "chr14" =	125194864,
                    "chr15" =	103494974,
                    "chr16" =	98319150,
                    "chr17" =	95272651,
                    "chr18" =	90772031,
                    "chr19" =	61342430,
                    "chrX" =	166650296,
                    "chrY" =	15902555)

chr.mm9.gr <- GRanges(seqnames = paste0("chr", c(1:19, "X", "Y")),
                      IRanges(start = 1, width = c(197195432,
                                                   181748087,
                                                   159599783,
                                                   155630120,
                                                   152537259,
                                                   149517037,
                                                   152524553,
                                                   131738871,
                                                   124076172,
                                                   129993255,
                                                   121843856,
                                                   121257530,
                                                   120284312,
                                                   125194864,
                                                   103494974,
                                                   98319150,
                                                   95272651,
                                                   90772031,
                                                   61342430,
                                                   166650296,
                                                   15902555)),
                      strand = "*")
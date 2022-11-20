
# ------------------------------------------- Gene DE ------------------------------------------- #

if (FALSE) {
  require(DESeq2)
  
  filenames <- sort(list.files('../../mESC_Ezh2i_BAP1/data/kallisto_output', full.names = T)) # count with combined reference GENCODE vM20 and ncRNA annotation
  sampleNewName <- gsub(".*/", "\\2", filenames)
  
  # read count
  count_table <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                    as = 'matrix', what = "est_counts")
  
  colnames(count_table) <- sampleNewName
  
  # save sample sizes
  LRNA.sizefactor <- size_factor_cal(count_table[grepl("chrS(2|4|5|8)", rownames(count_table)), ])
  saveRDS(LRNA.sizefactor, "../fig1/data/LRNA.sizefactor.rds")
  
  # get DE results
  txCounts <- count_table[grepl("^ENS", rownames(count_table)), ] %>%
    keepOneTx(rowname_gene_id = T, is_gene_sum = T)
  
  dds <- DESeqDataSetFromMatrix(round(as.matrix(txCounts)),
                                colData = data.frame(condition = gsub("_(1|2)$", "\\2", sampleNewName) ),
                                design = ~ condition)
  dds <- DESeq(dds)
  
  get_res <- function(dds) {
    res_E0B6 <- results(dds, contrast = c("condition", "mESC_Ezh2i_0d_BAP1_6h_2", "mESC_Ezh2i_0d_BAP1_0h_8"))
    res_E0B12 <- results(dds, contrast = c("condition", "mESC_Ezh2i_0d_BAP1_12h_4", "mESC_Ezh2i_0d_BAP1_0h_8"))
    
    res_E1B6 <- results(dds, contrast = c("condition", "mESC_Ezh2i_1d_BAP1_6h_1", "mESC_Ezh2i_0d_BAP1_0h_8"))
    res_E1B12 <- results(dds, contrast = c("condition", "mESC_Ezh2i_1d_BAP1_12h_3", "mESC_Ezh2i_0d_BAP1_0h_8"))
    
    res_E1B0 <- results(dds, contrast = c("condition", "mESC_Ezh2i_1d_BAP1_0h_5", "mESC_Ezh2i_0d_BAP1_0h_8"))
    res_E2B0 <- results(dds, contrast = c("condition", "mESC_Ezh2i_2d_BAP1_0h_6", "mESC_Ezh2i_0d_BAP1_0h_8"))
    res_E7B0 <- results(dds, contrast = c("condition", "mESC_Ezh2i_7d_BAP1_0h_7", "mESC_Ezh2i_0d_BAP1_0h_8"))
    
    res.list <- list(res_E0B6, res_E0B12, res_E1B6, res_E1B12, res_E1B0, res_E2B0, res_E7B0)
    names(res.list) <- c("res_E0B6", "res_E0B12", "res_E1B6", "res_E1B12",
                         "res_E1B0", "res_E2B0", "res_E7B0")
    return(res.list)
  }
  
  TTseq_res.list <- get_res(dds)
  saveRDS(TTseq_res.list, "../fig1/data/TTseq_res.list.RData")
  
  # differential expression above the thresholds
  get_DE_genes <- function(res, direction = "up", lFC = 1.5, padj = 0.05) {
    idx <- is.finite(res$log2FoldChange) & !is.na(res$padj)
    if (direction == "up") {
      out <- rownames(res)[res$log2FoldChange > lFC & res$padj < padj]
    } else {
      out <- rownames(res)[res$log2FoldChange < -lFC & res$padj < padj]
    }
    out[!is.na(out)]
  }
  
  DE_tt.gene.list <- with(TTseq_res.list,
                          list("E0B6_up" = get_DE_genes(res_E0B6, "up"),
                               "E0B6_down" = get_DE_genes(res_E0B6, "down"),
                               "E0B12_up" = get_DE_genes(res_E0B12, "up"), 
                               "E0B12_down" = get_DE_genes(res_E0B12, "down"), 
                               "E1B0_up" = get_DE_genes(res_E1B0, "up"), 
                               "E1B0_down" = get_DE_genes(res_E1B0, "down"), 
                               "E2B0_up" = get_DE_genes(res_E2B0, "up"), 
                               "E2B0_down" = get_DE_genes(res_E2B0, "down"), 
                               "E7B0_up" = get_DE_genes(res_E7B0, "up"), 
                               "E7B0_down" = get_DE_genes(res_E7B0, "down"), 
                               "E1B12_up" = get_DE_genes(res_E1B12, "up"), 
                               "E1B12_down" = get_DE_genes(res_E1B12, "down") ))
  
  lengths(DE_tt.gene.list)
  
  # up
  intersect_two_sets(DE_tt.gene.list$E0B6_up, DE_tt.gene.list$E0B12_up) #  8 189 471
  intersect_two_sets(DE_tt.gene.list$E0B12_up, DE_tt.gene.list$E1B12_up) # 86 574 189
  
  # down
  intersect_two_sets(DE_tt.gene.list$E0B6_down, DE_tt.gene.list$E0B12_down) # 20  21 202
  intersect_two_sets(DE_tt.gene.list$E0B12_down, DE_tt.gene.list$E1B12_down) # 117 106  91
}


TTseq_res.list <- readRDS("../fig1/data/TTseq_res.list.RData")


if (TRUE) {
  message("TT-seq DE MAplots")
  g_list <- with(TTseq_res.list,
                 list(plot_ma(baseMean = res_E0B6$baseMean,
                              log2FC = res_E0B6$log2FoldChange, 
                              pval = res_E0B6$padj, 
                              DE_log2FC = 1.5,
                              title = "\nTT-seq: P6"),
                      
                      plot_ma(baseMean = res_E0B12$baseMean,
                              log2FC = res_E0B12$log2FoldChange, 
                              pval = res_E0B12$padj, 
                              DE_log2FC = 1.5, 
                              title = "\nTT-seq: P12"),
                      
                      plot_ma(baseMean = res_E1B12$baseMean,
                              log2FC = res_E1B12$log2FoldChange, 
                              pval = res_E1B12$padj, 
                              DE_log2FC = 1.5,
                              title = "\nTT-seq: Ezh2i 1d + P12"), 
                      
                      plot_ma(baseMean = res_E1B0$baseMean,
                              log2FC = res_E1B0$log2FoldChange, 
                              pval = res_E1B0$padj, 
                              DE_log2FC = 1.5,
                              title = "\nTT-seq: Ezh2i 1d"),
                      
                      plot_ma(baseMean = res_E2B0$baseMean,
                              log2FC = res_E2B0$log2FoldChange, 
                              pval = res_E2B0$padj, 
                              DE_log2FC = 1.5,
                              title = "\nTT-seq: Ezh2i 2d"),
                      
                      plot_ma(baseMean = res_E7B0$baseMean,
                              log2FC = res_E7B0$log2FoldChange, 
                              pval = res_E7B0$padj, 
                              DE_log2FC = 1.5,
                              title = "\nTT-seq: Ezh2i 7d")
                 ))
  
  ggsave(plot = do.call(grid.arrange,
                        c(g_list, ncol = 3)),
         filename = paste0("Fig4_MAplot_TT-seq_gene_est_count_DESeq2.png"),
         path = "../fig4/figs/",
         device = "png", width = 10, height = 6)
}


# ------------------------------------------- TU DE ------------------------------------------- #
if (FALSE) {
  Count_RNA_reads <- function(TU.list, bam_files)
  {
    require(Rsubread)
    locations <- unique(TU.list[[1]]$location)
    TU.out <- GRanges()
    for (loc in locations)
    {
      TU.tmp <- reduce(
        Reduce(c, lapply(TU.list, function(TU) TU[TU$location == loc]) )
      )
      TU.tmp$location <- loc
      TU.out %c=% TU.tmp
    }
    # append spike-in normalised read counts to TU GRanges output
    mcols(TU.out) %c=% featureCounts(files = bam_files,
                                     annot.ext = data.frame(GeneID = seq_along(TU.out),
                                                            Chr = seqnames(TU.out),
                                                            Start = start(TU.out),
                                                            End = end(TU.out),
                                                            Strand = strand(TU.out)),
                                     isPairedEnd=TRUE, nthreads = 8)$counts
    return(TU.out)
  }
  
  LRNA.sizefactor <- readRDS("../fig1/data/LRNA.sizefactor.rds")
  
  if (FALSE) {
    # mm9 annotations
    TU.list <- lapply(list.files("../data/TU_anno/mm9", pattern = "TU_filter.*", full.names = TRUE), 
                      importRanges)
    
    # count reads on annotated TUs
    bam_files <- list.files("/mnt/0E471D453D8EE463/TT_seq_data/RS20201112/bam_mm9",
                            pattern = ".bam$", full.names = TRUE)
  }
  
  if (TRUE) {
    # mm10 annotations
    TU.list <- lapply(list.files("../data/TU_anno/mm10", pattern = ".*", full.names = TRUE), 
                      importRanges)
    
    # count reads on annotated TUs
    bam_files <- list.files("/mnt/0E471D453D8EE463/TT_seq_data/RS20201112/bam_mm10",
                            pattern = "out.bam$", full.names = TRUE)
  }
  
  # featureCounts
  TU.counts <- Count_RNA_reads(TU.list = TU.list, bam_files = bam_files)
  
  TU.count_table <- sweep(as.matrix(mcols(TU.counts)[, -1]), 2, LRNA.sizefactor, "/")
  
  TU.count_table <- sapply(1:8, function(x) {
    rowMeans(TU.count_table[, c(2*x-1, 2*x)])
  })
  
  colnames(TU.count_table) <- gsub("mESC_(.*h)_.*", "\\1", names(LRNA.sizefactor[(1:8)*2-1]))
  
  require(DESeq2)
  dds <- DESeqDataSetFromMatrix(as.matrix(mcols(TU.counts)[, -1]),
                                colData = data.frame(condition = gsub("mESC_(.*h)_.*", "\\1", names(LRNA.sizefactor)) ),
                                design = ~ condition)
  dds_relative <- DESeq(dds) # relative DE uses internal normalization
  
  sizeFactors(dds) <- LRNA.sizefactor
  dds_spikein <- DESeq(dds) # absolute DE uses spike-in normalization
  
  res_rel_6h <- results(dds_relative, contrast = c("condition", "Ezh2i_0d_BAP1_6h", "Ezh2i_0d_BAP1_0h"))
  res_sp_6h <- results(dds_spikein, contrast = c("condition", "Ezh2i_0d_BAP1_6h", "Ezh2i_0d_BAP1_0h"))
  
  res_rel_12h <- results(dds_relative, contrast = c("condition", "Ezh2i_0d_BAP1_12h", "Ezh2i_0d_BAP1_0h"))
  res_sp_12h <- results(dds_spikein, contrast = c("condition", "Ezh2i_0d_BAP1_12h", "Ezh2i_0d_BAP1_0h"))
  
  res_dat <- data.frame("location" = mcols(TU.counts)[, 1],
                        "baseMean_rel" = res_rel_6h[, "baseMean"],
                        "log2FoldChange_rel_6h" = res_rel_6h[, "log2FoldChange"],
                        "padj_rel_6h" = res_rel_6h[, "padj"],
                        "log2FoldChange_sp_6h" = res_sp_6h[, "log2FoldChange"],
                        "padj_sp_6h" = res_sp_6h[, "padj"],
                        
                        "log2FoldChange_rel_12h" = res_rel_12h[, "log2FoldChange"],
                        "padj_rel_12h" = res_rel_12h[, "padj"],
                        "log2FoldChange_sp_12h" = res_sp_12h[, "log2FoldChange"],
                        "padj_sp_12h" = res_sp_12h[, "padj"]
  )
  
  TU.DE.gr <- TU.counts
  mcols(TU.DE.gr) <- res_dat
  
  saveRDS(TU.DE.gr, "../data/TT_seq_TU_DE_mm10_gr.RData")
}


# ------------------------------------- RNA-seq vs TT-seq ------------------------------------ #

if (FALSE) {
  message("RNA-seq vs TT-seq")
  
  # get TT-seq read counts
  bam_files <- list.files("/mnt/0E471D453D8EE463/TT_seq_data/RS20201112/bam_mm10/",
                          ".bam$", full.names = T)
  
  gencode.mm9 <-
    importRanges("/mnt/0E471D453D8EE463/genomeDir/GENCODE/gencode.mouse.v1.annotation.gtf")
  gencode.mm9 <- gencode.mm9[gencode.mm9$type == "exon" &
                               gencode.mm9$gene_type %in% c("protein_coding", "lincRNA")]
  
  count_table_ttseq <- Rsubread::featureCounts(
    bam_files,
    annot.ext = data.frame(
      GeneID = gencode.mm9$gene_id,
      Chr = as.character(seqnames(gencode.mm9)),
      Start = start(gencode.mm9),
      End = end(gencode.mm9),
      Strand = strand(gencode.mm9)
    ),
    isPairedEnd = TRUE,
    strandSpecific = 1,
    nthreads = 16
  )$counts
  
  
  # count_table_ttseq.norm <- sweep(count_table_ttseq, 2, LRNA.sizefactor, "/")
  count_table_ttseq.norm <- sweep(count_table_ttseq, 2, size_factor_cal(count_table_ttseq), "/")
  
  ttRC.log2FC <- log2(rowMeans(count_table_ttseq.norm[, 3:4]) / rowMeans(count_table_ttseq.norm[, 1:2]))


  # get RNA-seq read counts
  RNAseq_RC.log2FC <- log2(rowMeans(exon_RC_RNA[, 4:6]) / rowMeans(exon_RC_RNA[, 1:3])) # "exon_RC_RNA" is from "F1_src_load_RNA_read_count.R"
  
  
  cmp.log2FC <- data.frame(TT_seq = ttRC.log2FC,
                           RNA_seq = RNAseq_RC.log2FC) %>%
    dplyr::filter(complete.cases(.) & RNA_seq > -10 & RNA_seq < 10 & is.finite(rowSums(.)))
  
  ggplot(cmp.log2FC, aes(x = RNA_seq, y = TT_seq, color = get_dens(RNA_seq, TT_seq))) +
    geom_point(size = 1, alpha = 0.2) +
    scale_color_gradientn(name = "Density", colors = rev(colors_n)) +
    xlab("RNA-seq log2FC P12 vs P0") + 
    ylab("TT-seq log2FC P12 vs P0") + 
    annotate(geom = "text", 
             x = -6, y = 11,
             hjust = "left", vjust = "top",
             label = paste0("r = ", round(cor(cmp.log2FC)[1, 2], 3),
                            "\nn = ", nrow(cmp.log2FC))) +
    theme_setting +
    theme(legend.position = "none")
  
  ggsave(filename = "FigS4_scatter_RNA_TT_seq_log2FC_correlation.png", 
         path = "../figS4/figs",
         width = 5, height = 4.5, device = "png")
}

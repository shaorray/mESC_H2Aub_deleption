# Rui Shao 2022 Aug
# Figure 4
# TU DE neighboring effect

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
source("../util/getCoverage.R")

# functions ----------------------------------------------------------------------------------------------
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
  mcols(TU.out) %c=% featureCounts(bam_files,
                                   annot.ext = data.frame(GeneID = seq_along(TU.out),
                                                          Chr = seqnames(TU.out),
                                                          Start = start(TU.out),
                                                          End = end(TU.out),
                                                          Strand = strand(TU.out)),
                                   isPairedEnd=TRUE, nthreads = 8)$counts
  return(TU.out)
}

# mm9 read counts ----------------------------------------------------------------------------------------
if (T) {
  
  TU.list <- list("P0" = importRanges("../data/TU_anno/mm9/TU_filter+mESC_Ezh2i_0d_BAP1_0h_8_mm9.gff3"),
                  "P6" = importRanges("../data/TU_anno/mm10/TU_filter+mESC_Ezh2i_0d_BAP1_6h_2.gff3"),
                  "P12" = importRanges("../data/TU_anno/mm10/TU_filter+mESC_Ezh2i_0d_BAP1_12h_4.gff3"),
                  "E1" = importRanges("../data/TU_anno/mm10/TU_filter+mESC_Ezh2i_1d_BAP1_0h_5.gff3"),
                  "E1P12" = importRanges("../data/TU_anno/mm10/TU_filter+mESC_Ezh2i_1d_BAP1_12h_3.gff3"),
                  "E2" = importRanges("../data/TU_anno/mm10/TU_filter+mESC_Ezh2i_2d_BAP1_0h_6.gff3"),
                  "E7" = importRanges("../data/TU_anno/mm10/TU_filter+mESC_Ezh2i_7d_BAP1_0h_7.gff3") )
  
  bam_files <- list.files("/mnt/0E471D453D8EE463/TT_seq_data/RS20201112/bam_mm9/",
                          pattern = ".bam$", full.names = TRUE)
  
  # featureCounts
  TU.counts <- Count_RNA_reads(TU.list = TU.list, bam_files = bam_files)
  
  LRNA.sizefactor <- readRDS("../fig1/data/LRNA.sizefactor.rds")
  readCounts <- mcols(TU.counts)[, -1]
  colnames(readCounts) <- gsub("\\.Aligned.*", "\\1", colnames(readCounts))
  colnames(readCounts) <- gsub("\\.", "\\_", colnames(readCounts))
  readCounts <- t(t(as.matrix(readCounts[, names(LRNA.sizefactor)])) / LRNA.sizefactor)
  mcols(TU.counts)[, -1] <- readCounts
  
  # use gene interval instead of TU annotations
  TU.coding.gr <- gene.gr[countSubjectHits(findOverlaps(TU.counts, gene.gr, minoverlap = 500)) > 0]
  names(TU.coding.gr) <- TU.coding.gr$gene_id
  TU.coding.gr <- TU.coding.gr[TU.coding.gr$gene_type == "protein_coding"]
  coding.readCounts <- .countBam(bam_files, TU.coding.gr, stranded = TRUE)
  colnames(coding.readCounts) <- gsub("mESC_(.*).Aligned.sortedByCoord.out.bam", "\\1", colnames(coding.readCounts))
  
  library(DESeq2)
  # gene
  dds_gene <- DESeqDataSetFromMatrix(coding.readCounts,
                                   colData = data.frame(condition = gsub("(.*)_._.$", "\\1", colnames(coding.readCounts)) ),
                                   design = ~ condition)
  dds_gene <- DESeq(dds_gene)
  
  res_gene_P12 <- results(dds_gene, contrast = c("condition", "Ezh2i_0d_BAP1_12h", "Ezh2i_0d_BAP1_0h"))
  res_gene_E1P12 <- results(dds_gene, contrast = c("condition", "Ezh2i_1d_BAP1_12h", "Ezh2i_0d_BAP1_0h"))
  
  gene_DE_res <- data.frame("location" = "protein_coding",
                            "gene_id" = TU.coding.gr$gene_id,
                            "baseMean" = res_gene_P12[, "baseMean"],
                            "log2FoldChange_P12" = res_gene_P12[, "log2FoldChange"],
                            "padj_P12" = res_gene_P12[, "padj"],
                            "log2FoldChange_E1P12" = res_gene_E1P12[, "log2FoldChange"],
                            "padj_E1P12" = res_gene_E1P12[, "padj"]) 
  mcols(TU.coding.gr) <- gene_DE_res
  
  # TU
  dds_TU <- DESeqDataSetFromMatrix(round(as.matrix(readCounts)),
                                     colData = data.frame(condition = gsub("mESC_(.*)_._.$", "\\1", colnames(readCounts)) ),
                                     design = ~ condition)
  dds_TU <- DESeq(dds_TU)
  
  res_P12 <- results(dds_TU, contrast = c("condition", "Ezh2i_0d_BAP1_12h", "Ezh2i_0d_BAP1_0h"))
  res_E1P12 <- results(dds_TU, contrast = c("condition", "Ezh2i_1d_BAP1_12h", "Ezh2i_0d_BAP1_0h"))
 
  TU_DE_res <- data.frame("location" = mcols(TU.counts)[, 1],
                          "gene_id" = TU.counts$gene_id,
                          "baseMean" = res_P12[, "baseMean"],
                          "log2FoldChange_P12" = res_P12[, "log2FoldChange"],
                          "padj_P12" = res_P12[, "padj"],
                          "log2FoldChange_E1P12" = res_E1P12[, "log2FoldChange"],
                          "padj_E1P12" = res_E1P12[, "padj"]) 
  TU_DE_res.gr <- TU.counts
  mcols(TU_DE_res.gr) <- TU_DE_res
}


# -----------------------------------------------------------------------------------------------------
# gene-intergenic TU neighboring
TU.intergenic.gr <- TU_DE_res.gr[TU_DE_res.gr$location == "intergenic"]
TU.intergenic.gr$id <- seq_along(TU.intergenic.gr)

gene.intergenic.pairs <- foreach (i = seqlevels(TU.coding.gr), .combine = rbind) %dopar% {
  # extract intergenic neighbors for each chromosome
  genes.tmp = TU.coding.gr[seqnames(TU.coding.gr) == i]
  genes.tmp = genes.tmp[order(start(genes.tmp))]
  TU.tmp = TU.intergenic.gr[seqnames(TU.intergenic.gr) == i]
  TU.strands = as.character(strand(TU.tmp)) == "+"
  
  out_table <- NULL
  for (j in seq_along(genes.tmp)) {
    gene.strand = as.character(strand(genes.tmp[j])) == "+"
    gaps = start(TU.tmp) - start(genes.tmp[j])
    
    same.strand = (TU.strands == gene.strand)
    
    if (length(genes.tmp) == 1) {
      up.tu = .ifelse(gene.strand, gaps < 0, gaps > 0)
      down.tu = .ifelse(gene.strand, gaps > 0, gaps < 0)
    } else if (j == 1) {
      up.tu =  .ifelse(gene.strand, 
                       gaps < 0, 
                       gaps > 0 & start(TU.tmp) < start(genes.tmp[2]) )
      down.tu = .ifelse(!gene.strand, 
                        gaps < 0,
                        gaps > 0 & start(TU.tmp) < start(genes.tmp[2]))
    } else if (j == length(genes.tmp)) {
      up.tu =  .ifelse(gene.strand, 
                       start(TU.tmp) > start(genes.tmp[j-1]) & gaps < 0, 
                       gaps > 0)
      down.tu = .ifelse(!gene.strand,
                        start(TU.tmp) > start(genes.tmp[j-1]) & gaps < 0, 
                        gaps > 0)
    } else {
      up.tu =  .ifelse(gene.strand, 
                       gaps < 0 & start(TU.tmp) > start(genes.tmp[j-1]), 
                       gaps > 0 & start(TU.tmp) < start(genes.tmp[j+1]))
      down.tu = .ifelse(!gene.strand, 
                        gaps < 0 & start(TU.tmp) > start(genes.tmp[j-1]),
                        gaps > 0 & start(TU.tmp) < start(genes.tmp[j+1]))
    }
    up.s = same.strand & up.tu
    up.as = !same.strand & up.tu
    down.s = same.strand & down.tu
    down.as = !same.strand & down.tu
    
    out_table <- rbind(out_table,
                       data.frame("gene_id" = rep(genes.tmp[j]$gene_id, sum(up.s + up.as + down.s + down.as)),
                                  "strand" = rep(strand(genes.tmp[j]), sum(up.s + up.as + down.s + down.as)),
                                  "type" = c(rep("up.sense", sum(up.s)), 
                                             rep("up.antisense", sum(up.as)),
                                             rep("down.sense", sum(down.s)),
                                             rep("down.antisense", sum(down.as)) ),
                                  "gap" = c(gaps[up.s], gaps[up.as], gaps[down.s], gaps[down.as]),
                                  "id" = c(TU.tmp$id[up.s], TU.tmp$id[up.as], TU.tmp$id[down.s], TU.tmp$id[down.as])))
  }
  out_table
}

# gap of TU TSS to gene boundary
gaps <- start(promoters(TU.intergenic.gr[gene.intergenic.pairs$id], upstream = 0, downstream = 0)) - 
  cbind(start(TU.coding.gr[gene.intergenic.pairs$gene_id]),
        end(TU.coding.gr[gene.intergenic.pairs$gene_id]))

idx <- ifelse(gene.intergenic.pairs$strand == "+", 
              ifelse(grepl("up", gene.intergenic.pairs$type), 1, 2), 
              ifelse(grepl("up", gene.intergenic.pairs$type), 2, 1))

gene.intergenic.pairs$gap_boundary <- 
  sapply(seq_along(idx),
         function(x) gaps[x, idx[x]]
  ) * ifelse(gene.intergenic.pairs$strand == "+", 1, -1) / 1e3

tmp_table <- cbind(mcols(TU.coding.gr[gene.intergenic.pairs$gene_id])[, 3:7], 
                   mcols(TU.intergenic.gr[gene.intergenic.pairs$id])[, 3:7])
colnames(tmp_table) <- paste0(c(rep("Gene_", 5), rep("TU_", 5)), colnames(tmp_table))
gene.intergenic.pairs <- cbind(gene.intergenic.pairs, tmp_table)
rm(tmp_table)

# show co-expression ------------------------------------------------------------------
gene.intergenic.pairs <-
  gene.intergenic.pairs[!apply(gene.intergenic.pairs[, grep("log2FoldChange", colnames(gene.intergenic.pairs))], 1, 
                               function(x) any(is.na(x))), ]

gap_breaks <- c(min(gene.intergenic.pairs$gap_boundary), 
                c(-400, -200, -100, -75, -50, -30, -20, -10, 0, 10, 20, 30, 50, 75, 100, 200, 400),# * 1000,
                max(gene.intergenic.pairs$gap_boundary))

gene.intergenic.pairs$gap_class <- cut(gene.intergenic.pairs$gap_boundary, 
                                       breaks = gap_breaks, 
                                       labels = 1:18) %>% as.numeric()
gene.intergenic.pairs <- gene.intergenic.pairs[!is.na(gene.intergenic.pairs$gap_class), ]

# plot by TU log2FC correlation
if (TRUE) {
  mat <- NULL
  for( i in 1:18 ) {
    tmp.pairs <- gene.intergenic.pairs[as.numeric(gene.intergenic.pairs$gap_class) == i, ]
    
    PcG.idx <- tmp.pairs$gene_id %in% Ring1b_enriched_genes
    
    tmp.pairs <- tmp.pairs[, grep("log2", colnames(gene.intergenic.pairs))] %>% as.matrix()
    mat <- rbind(mat, c("Pos" = i, 
                        "cor_P12_PcG" = cor(tmp.pairs[PcG.idx, 1], tmp.pairs[PcG.idx, 3]),
                        "cor_E1P12_PcG" = cor(tmp.pairs[PcG.idx, 2], tmp.pairs[PcG.idx, 4]),
                        "cor_2i_non" = cor(tmp.pairs[!PcG.idx, 1], tmp.pairs[!PcG.idx, 3]),
                        "cor_mTORi_non" = cor(tmp.pairs[!PcG.idx, 2], tmp.pairs[!PcG.idx, 4])
    ))
  }
  mat[10:18, 1] <- (10:18) + 3 # add gene box
  mat <- rbind(mat[1:9, ], matrix(c(10:12, rep(NA, 12)), nrow = 3), mat[10:18, ])
  
  dat_mat <- reshape2::melt(data = as.data.frame(mat[, 1:5]), id = "Pos")
  dat_mat$Sample <- rep(c(rep("P12", 21), rep("E1P12", 21)), 2)
  dat_mat$PcG <- factor(c(rep("PcG", 21*2), rep("Other", 21*2)),levels = c("PcG", "Other"))
  dat_mat$line_group <- factor(cumsum(is.na(dat_mat$value)))
  
  ggplot(dat_mat[dat_mat$Sample == "P12", ], aes(x = Pos, y = value, color = PcG, group = line_group)) +
    geom_hline(yintercept = 0, lty = 2, col = "grey85") +
    # geom_vline(xintercept = 10, col = "grey85") +
    geom_rect(aes(xmin = 10, xmax = 12, ymin = -Inf, ymax = 0),
              fill = "#0000B2", linetype = 0) +
    geom_rect(aes(xmin = 9, xmax = 13, ymin = 0, ymax = Inf),
              fill = "grey95", linetype = 0) +
    geom_rect(aes(xmin = 10, xmax = 12, ymin = 0, ymax = Inf),
              fill = "grey75", linetype = 0) +
    # geom_point() + 
    geom_line(lwd = 2) +
    
    scale_x_continuous(name = "\nDistance to gene (kb)", 
                       breaks = c(1, 5, 8, 10, 12, 14, 17, 21), 
                       labels = c(gap_breaks[c(2, 5, 8)], c("5\'", "3\'"), gap_breaks[c(12, 15, 18)]) ) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4),
                       labels = c(0, 0.2, 0.4), 
                       limits = c(-0.1, 0.4)) +
    scale_color_manual(name = "Intergenic", values = colors_9[7:6]) +
    ylab("Transcription log2FC correlation") +
    theme_setting + 
    theme(strip.text = element_blank(),
          axis.ticks = element_line(), 
          panel.grid.minor = element_blank(), 
          panel.spacing = unit(1, "lines"),
          legend.position = "top")
}

table(unique(gene.intergenic.pairs$gene_id) %in% Ring1b_enriched_genes) #  Other: 9788   Ring1b_enriched: 759 


ggsave(filename = "FigS8_TU_gene_log2FC_position_correlation_direction.png",
       path = "../figS8/figs", device = "png", width = 5, height = 4 )

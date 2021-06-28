setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")

# mm9 genes
gene.gr <- GenomicFeatures::genes(TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene)
res <- biomaRt::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                       keys = as.character(gene.gr$gene_id),
                       keytype = "ENTREZID",
                       columns = c("GENEID", "GENENAME", "GENEBIOTYPE"))
gene.gr$gene_id <- res$GENEID[match(names(gene.gr), res$ENTREZID)]
gene.gr$gene_type <- res$GENEBIOTYPE[match(names(gene.gr), res$ENTREZID)]
gene.gr$gene_name <- res$GENENAME[match(names(gene.gr), res$ENTREZID)]
names(gene.gr) <- gene.gr$gene_id
gene.gr <- gene.gr[!is.na(names(gene.gr))]
saveRDS(gene.gr, "../data/gene.gr.RData")

# TT-seq TX counts -------------------------------------------------------------------------------------
# load kallisto counts, mm10
filenames <- sort(list.files('../data/kallisto_output', full.names = T)) # count with combined reference GENCODE vM20 and ncRNA annotation
sampleNewName <- gsub(".*/", "\\2", filenames)

count_table <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                  as = 'matrix', what = "est_counts")
colnames(count_table) <- sampleNewName

txRC <- count_table[grepl("^ENS", rownames(count_table)), ] %>%
  keepOneTx(rowname_gene_id = T)

tuRC <- count_table[!grepl("^chrS|^ENS", rownames(count_table)), ]
tuRC <- tuRC[rowSums(tuRC) > 0, ]

spRC <- count_table[grep("^chrS", rownames(count_table)), ]

merge_replicates <- function(dat) {
  idx = seq_len(ncol(dat) %/% 2) * 2
  (dat[, idx - 1] + dat[, idx]) / 2
}

LRNA.sizefactor <- spRC[1:4, ] %>% merge_replicates() %>% SizeFactorCal()
# saveRDS(LRNA.sizefactor, "data/LRNA.sizefactor.RData")

txRC.norm <- sweep(merge_replicates(txRC), 2, LRNA.sizefactor, "/")

if (T) {
  sample.sizefactor <- SizeFactorCal(txRC.norm[apply(txRC.norm, 1, function(x) mean(x) > 1), ])
  dat <- data.frame(sample_size = sample.sizefactor,
                    sample_num = paste0("#", (c(8, 4, 2, 5, 3, 1, 6, 7))),
                    treatment = gsub("mESC_(.*)_1", "\\1", names(sample.sizefactor)))
  
  dat$treatment <- factor(dat$treatment, levels = dat$treatment[order(dat$sample_num)])
  ggplot(dat, aes(x = sample_num, y = sample_size, fill = treatment)) +
    geom_bar(stat = "identity") +
    ggtitle("mRNA (n = 17401)") +
    theme_minimal()
  ggsave("sample_sizes_mRNA.png", path = "figs", width = 5, height = 4)
  
  
  tuRC.norm <- sweep(merge_replicates(tuRC), 2, LRNA.sizefactor, "/")
  nc.sample.sizefactor <- SizeFactorCal(tuRC.norm[apply(tuRC.norm, 1, function(x) mean(x) > 0.1), ])
  
  dat <- data.frame(sample_size = nc.sample.sizefactor,
                    sample_num = paste0("#", (c(8, 4, 2, 5, 3, 1, 6, 7))),
                    treatment = gsub("mESC_(.*)_1", "\\1", names(nc.sample.sizefactor)))
  
  dat$treatment <- factor(dat$treatment, levels = dat$treatment[order(dat$sample_num)])
  ggplot(dat, aes(x = sample_num, y = sample_size, fill = treatment)) +
    geom_bar(stat = "identity") +
    ggtitle("all ncRNA annotations (n = 39497)") +
    theme_minimal()
  # ggsave("sample_sizes_ncRNA.png", path = "figs", width = 5, height = 4)
}

if (F) {
  # get TT-seq read counts
  bam_files <- list.files("/mnt/0E471D453D8EE463/TT_seq_data/RS20201112/bam_mm9/",
                          ".bam$", full.names = T)
  
  gencode.mm9 <- importRanges("/mnt/0E471D453D8EE463/genomeDir/GENCODE/gencode.mouse.v1.annotation.gtf")
  gencode.mm9 <- gencode.mm9[gencode.mm9$type == "exon" & 
                               gencode.mm9$gene_type %in% c("protein_coding", "lincRNA")]
  
  
  count_table_ttseq <- Rsubread::featureCounts(bam_files,
                                               annot.ext = data.frame(GeneID = gencode.mm9$gene_id,
                                                                      Chr = as.character(seqnames(gencode.mm9)),
                                                                      Start = start(gencode.mm9),
                                                                      End = end(gencode.mm9),
                                                                      Strand = strand(gencode.mm9)),
                                               isPairedEnd=TRUE, strandSpecific = 1, nthreads = 16)$counts
  
  count_table_ttseq.norm <- sweep(merge_replicates(count_table_ttseq), 2, LRNA.sizefactor, "/")
}

# load ChIP-seq data
bw_files <- list.files("/mnt/E0767589767560E8/UPPMAX/Ezh2i_BAP1_ChIP",
                        ".bw$", full.names = T)
sample_names <- gsub(".*ChIP\\/(.*)", "\\1", bw_files)

sample_table <- data.frame(Target = gsub("(_BAP1|_Ezh2i|-NTD|_NT|_Trp).*", "", sample_names),
                           Treatment = gsub("_RC.fltd.bw", "", sample_names))

source("../util/getCoverage.R")
count_table_ChIP <- .countBW(bw_files = bw_files, 
                             intervals = promoters(gene.gr, upstream = 2000, downstream = 2000), 
                             fast = T)
count_table_ChIP <- count_table_ChIP[match(rownames(count_table_ttseq),
                                           res$GENEID[match(rownames(count_table_ChIP), 
                                                            res$ENTREZID)]), ]

sample_ChIP_sizes <- t(matrix(size_factor_cal(count_table_ChIP), nrow = 9))
colnames(sample_ChIP_sizes) <- sample_table$Treatment[1:9]
rownames(sample_ChIP_sizes) <- unique(sample_table$Target)

pheatmap::pheatmap(sample_ChIP_sizes[, -8] / sample_ChIP_sizes[, 8],
                   cluster_rows = F, cluster_cols = F,
                   main = "Relative changes to control")

if (F) {
  # correlate ChIP changes
  log2FC_list <- list(log2(count_table_ttseq.norm[, c(3,2,6,5,4,7,8)] / count_table_ttseq.norm[, 1]))
  for (i in unique(sample_table$Target)) {
    tmp_counts <- count_table_ChIP[, sample_table$Target == i]
    tmp_FC <- (log2(tmp_counts[, c(3,1,4,2,5,6,7)] / tmp_counts[, 8]))
    log2FC_list <- c(log2FC_list, list(tmp_FC))
  }
  
  Treatment = c("BAP1_6h", "BAP1_12h", "BAP1_6h_Ezh2_1d", "BAP1_12h_Ezh2_1d", 
                "Ezh2_1d", "Ezh2_2d", "Ezh2_7d")
  
  g_list <- list()
  for (i in 1:8) {
    for (j in 1:7) {
      dat <- data.frame(Tx_log2FC = log2FC_list[[1]][, j],
                        ChIP_log2FC = log2FC_list[[i+1]][, j])
      dat <- dat[apply(dat, 1, function(x) all(!is.infinite(x) & !is.na(x))), ]
      
      g_list <- c(g_list, 
                  list(ggplot(dat, aes(x = ChIP_log2FC, y = Tx_log2FC)) +
                         # geom_point(aes(color = with(dat, get_dens(ChIP_log2FC, Tx_log2FC)))) +
                         geom_hex(bins = 50) +
                         scale_fill_gradient(low="grey75", high=colors_9[i]) +
                         ggtitle(paste(unique(sample_table$Target)[i], Treatment[j])) +
                         theme_minimal() +
                         theme(legend.position = "none")))
    }
  }
  
  ggsave(do.call(grid.arrange, c(g_list, nrow = 4)),
         filename = "cor_all_Tx_ChIP.png",
         width = 40, height = 10, path = "figs")
  
  # estimate kinetics by bins
  sample_bw_bins <- binBwBam(file_names = bw_files)
  colnames(sample_bw_bins) <- paste(sample_table$Target, sample_table$Treatment, sep = "_")
  
  estimate_kinetics <- function(dat, point) {
    point = cbind(1, point)
    dat = log(dat)
    (solve(t(point) %*% point) %*% t(point) %*% t(dat))[2, ]
  }
  
  bin_H2Aub_kinetics <- estimate_kinetics(sample_bw_bins[, c(8, 3, 1)], c(0, 6, 12))
  bin_H3K27m3_kinetics <- estimate_kinetics(sample_bw_bins[, c(17, 14:15)], c(0, 1, 2))
  
  # estimate kinetics by genes
  gene_H2Aub_kinetics <- estimate_kinetics(count_table_ChIP[, c(8, 3, 1)], c(0, 6, 12))
  gene_H3K27m3_kinetics <- estimate_kinetics(count_table_ChIP[, c(17, 14:15)], c(0, 1, 2))
}


# RNA-seq BGI ------------------------------------------------------------------------------------
filenames <- sort(list.files('../data/kallisto_output_BGI/', full.names = T)) # count with combined reference GENCODE vM20 and ncRNA annotation
sampleNewName <- gsub(".*/", "\\2", filenames)

.design <- c("Ni_R1", "Ni_R2", "Ni_R3", "Pr_R1", "Pr_R2", "Pr_R3",
             "Nip_R1", "Nip_R2", "Nip_R3", "Prp_R1", "Prp_R2", "Prp_R3",
             "P0_R1", "P0_R2", "P0_R3",
             "P12_R1", "P12_R2", "P12_R3",
             "P12_C12_R1", "P12_C12_R2", "P12_C12_R3",
             "P12_C24_R1", "P12_C24_R2", "P12_C24_R3",
             "P12_C36_R1", "P12_C36_R2", "P12_C36_R3")
names(.design) <- 
  sampleNewName[order(gsub("(^.).*", "\\1", sampleNewName), as.numeric(gsub("^.", "", sampleNewName)))]
  
count_table_RNA <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                  as = 'matrix', what = "tpm")
colnames(count_table_RNA) <- sampleNewName
count_table_RNA <- count_table_RNA[, names(.design)]

txRPK_RNA <- count_table_RNA[grepl("^ENS", rownames(count_table_RNA)), ] %>%
  keepOneTx(rowname_gene_id = T)

flyRPK_RNA <- count_table_RNA[grep("^FBtr", rownames(count_table_RNA)), ]

SizeFactorCal(flyRPK_RNA[rowMeans(flyRPK_RNA) > 1, ]) # spike-in drosophila genome 

# get read counts for DE test
count_table_RNA_rc <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                      as = 'matrix', what = "est_counts")
colnames(count_table_RNA_rc) <- sampleNewName

txRC_RNA <- count_table_RNA_rc[grepl("^ENS", rownames(count_table_RNA_rc)), names(.design)] %>%
            keepOneTx(rowname_gene_id = T)
saveRDS(txRC_RNA, "data/txRC_RNA_seq.RData")

if (F) {
  # get RNA-seq read counts
  bam_files <- list.files("/mnt/E0767589767560E8/UPPMAX/BGI/mESC/bam_mm9",
                          ".bam$", full.names = T)
  
  count_table_rnaseq <- Rsubread::featureCounts(bam_files,
                                                annot.ext = data.frame(GeneID = gencode.mm9$gene_id,
                                                                       Chr = as.character(seqnames(gencode.mm9)),
                                                                       Start = start(gencode.mm9),
                                                                       End = end(gencode.mm9),
                                                                       Strand = strand(gencode.mm9)),
                                                isPairedEnd=TRUE, 
                                                strandSpecific = 1,
                                                nthreads = 16)$counts
  colnames(count_table_rnaseq) <- sampleNewName
  count_table_rnaseq.norm <- sweep(count_table_rnaseq, 
                                   MARGIN = 2, 
                                   SizeFactorCal(flyRPK_RNA[rowMeans(flyRPK_RNA) > 1, ]), 
                                   "/")
  rownames(count_table_rnaseq.norm) <- gsub("\\..*", "", rownames(count_table_rnaseq.norm))
}

source("run_DESeq.R")


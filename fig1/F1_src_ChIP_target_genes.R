
message("get PcG target genes")

# signal enrichment ---------------------------------------------------------------------------
message("Ring1b target genes")
Ring1b_enriched_genes <-
  kink_index(ChIP_TSS_B1_mat[, "Ring1b_NT"],
             method = "poisson", .p = 0.01) %>%
  "["(rownames(ChIP_TSS_B1_mat), .) # 2402, at 0.470823

message("Ezh2 target genes")
Ezh2_enriched_genes <-
  kink_index(ChIP_TSS_B1_mat[, "Ezh2_NT"],
             method = "poisson", .p = 0.01) %>%
  "["(rownames(ChIP_TSS_B1_mat), .) # 2328, at 0.2260815

PcG_enriched_genes <- intersect.Vector(Ring1b_enriched_genes, Ezh2_enriched_genes) # 1480
saveRDS(PcG_enriched_genes, "data/PcG_enriched_genes.rds")


# show stat cutoff illustrative plots

g1 <- data.frame(y = ChIP_TSS_B1_mat[, "Ring1b_NT"]) %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::mutate(y = sort(y),
                x = seq_along(y),
                color = y > 0.470823) %>%
  ggplot(aes(x = x, y = y)) +
  geom_rect(mapping = aes(xmin = -Inf,
                          xmax = Inf,
                          ymin = -Inf,
                          ymax = 0.470823), 
            fill = add.alpha("grey85", 0.5)) +
  geom_point(aes(color = color, pch = color, size = color)) +
  scale_color_manual(values = c("black", "blue3")) +
  scale_shape_manual(values = c(5, 1)) +
  scale_size_manual(values = c(0.5, 2)) +
  xlab("Order of genes") + ylab("TSS ±1 kb density") +
  ggtitle("Ring1b enriched genes (n=2402)") +
  theme_bw() +
  theme_setting +
  theme(panel.grid = element_blank(), 
        legend.position = "none") 

g2 <- data.frame(y = ChIP_TSS_B1_mat[, "Ezh2_NT"]) %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::mutate(y = sort(y),
                x = seq_along(y),
                color = y > 0.2260815) %>%
  ggplot(aes(x = x, y = y)) +
  geom_rect(mapping = aes(xmin = -Inf,
                          xmax = Inf,
                          ymin = -Inf,
                          ymax = 0.2260815), 
            fill = add.alpha("grey85", 0.5)) +
  geom_point(aes(color = color, pch = color, size = color)) +
  scale_color_manual(values = c("black", "blue3")) +
  scale_shape_manual(values = c(5, 1)) +
  scale_size_manual(values = c(0.5, 2)) +
  xlab("Order of genes") + ylab("TSS ±1 kb density") +
  ggtitle("Ezh2 enriched genes (n=2328)") +
  theme_bw() +
  theme_setting +
  theme(panel.grid = element_blank(), 
        legend.position = "none") 

ggsave(filename = "../figS1/figs/FigS1_Ring1b_Ezh2_enriched_genes_ppois_001.png",
       grid.arrange(g1, g2, nrow = 1), device = "png",
       width = 7, height = 3.5)


# H2Aub binding genes
H2Aub_enriched_genes <-
  rownames(ChIP_TSS_B1_mat)[kink_index(ChIP_TSS_B1_mat[, "H2Aub_NT"] + ChIP_TSS_B2_mat[, "H2Aub_P0"], 
                                       method = "poisson", 0.01)] # 1699


# GSM3100852_Ring1b_EZH2_Peaks  ---------------------------------------------------------------

if (FALSE) {
  GSM3100852_Ring1b_EZH2_Peaks <- importRanges("../data/peaks/GSM3100852_Ring1b_EZH2_Peaks.mm9.bed")
  GSM3100852_Ring1b_EZH2_Peaks <- GSM3100852_Ring1b_EZH2_Peaks[GSM3100852_Ring1b_EZH2_Peaks$name == "both"]
  PcG_GSM3100852_genes <- gene.gr$gene_id[countQueryHits(findOverlaps(promoters(gene.gr), 
                                                                      GSM3100852_Ring1b_EZH2_Peaks)) > 0]
}

# Bivalent target genes -----------------------------------------------------------------------

if (FALSE) {
  
  H3K4me3_enriched_genes <- rownames(ChIP_TSS_B1_mat)[kink_index(ChIP_TSS_B1_mat[, "H3K4m3_NT"], 
                                                                 method = "poisson", .p = 0.01)] # 2053
  
  H3K27me3_enriched_genes <-
    rownames(ChIP_TSS_B1_mat)[kink_index(ChIP_TSS_B1_mat[, "H3K27m3_NT"], 
                                         method = "poisson", 0.01)] # 2846
  
  # H2AZ binding genes
  H2AZ_enriched_genes <-
    .countBW(bw_files = "/mnt/0E471D453D8EE463/GEO_bw/Histone/2016_Surface_H2AZ_WT_ChIP_Seq.mm9.bw", 
             intervals = promoters(gene.gr), fast = T) %>% 
    unlist() %>% 
    kink_index(method = "poisson", 0.05) %>% 
    "["(names(gene.gr), .) # 4665
  
  # H33 binding genes
  H33_enriched_genes <-
    .countBW("/mnt/0E471D453D8EE463/GEO_bw/Histone/GSM2582412_ESC_H3.3WT_YFP.bigwig", 
             intervals = promoters(gene.gr), fast = T) %>%
    unlist() %>% 
    kink_index(method = "poisson", 0.05) %>% 
    "["(names(gene.gr), .) # 4481
}



# PRC1 types ----------------------------------------------------------------------------------

if (FALSE) {
  Ring1b_binding_genes <-
    rownames(ChIP_TSS_B1_mat)[kink_index(ChIP_TSS_B1_mat[, "Ring1b_NT"] + ChIP_TSS_B2_mat[, "Ring1b_P0"], 
                                         method = "poisson", 0.05)] # 3173
  Ezh2_binding_genes <-
    rownames(ChIP_TSS_B1_mat)[kink_index(ChIP_TSS_B1_mat[, "Ezh2_NT"] + ChIP_TSS_B2_mat[, "Ezh2_P0"], 
                                         method = "poisson", 0.05)] # 3306
  Cbx7_binding_genes <-
    rownames(ChIP_TSS_B2_mat)[kink_index(ChIP_TSS_B2_mat[, "Cbx7_P0"], 
                                             method = "poisson", 0.05)] # 3700
  Rybp_binding_genes <-
    rownames(ChIP_TSS_B2_mat)[kink_index(ChIP_TSS_B2_mat[, "Rybp_P0"], 
                                             method = "poisson", 0.05)] # 4155
  
  Rsf1_peaks <- importRanges("../data/peaks/Rsf1_peaks.narrowPeak")
  Rsf1_peaks <- Rsf1_peaks[Rsf1_peaks$qValue > 7]
  Rsf1_binding_genes <- names(gene.gr)[findOverlaps(gene.gr, Rsf1_peaks) %>% countQueryHits() > 0] # 6500
  
  PcG_binding_genes <- intersect.Vector(Ring1b_binding_genes, Ezh2_binding_genes) # 2093
  saveRDS(PcG_binding_genes, "data/PcG_binding_genes.rds")
}


# H2Aub repressed genes ------------------------------------------------------------------------------
message("get derepressed genes")
# Ring1b CKO overlapped genes + PcG bound genes
derepressed_gene_ids <- unique(unlist(DE_gene.list[grep("up", names(DE_gene.list))]))

H2Aub_repressed_genes <- derepressed_gene_ids[derepressed_gene_ids %in% c(Ring1b_CKO_genes, PcG_binding_genes)]
saveRDS(H2Aub_repressed_genes, "../fig1/data/H2Aub_repressed_genes.rds")

message("BAP1 specific; overlapped; PRC1 CKO")
intersect_two_sets(derepressed_gene_ids, Ring1b_CKO_genes)


if (FALSE) {
  derepressed_gene_tab <- data.frame(BAP1_specific = ifelse(derepressed_gene_ids %in% Ring1b_CKO_genes, "Ring1b_CKO", "BAP1_specific"),
                                     PcG_binding = ifelse(derepressed_gene_ids %in% PcG_binding_genes, "High", "Low")) %>%
    dplyr::count(BAP1_specific, PcG_binding,sort = TRUE)
  
  derepressed_gene_tab$BAP1_specific <- factor(derepressed_gene_tab$BAP1_specific, 
                                               levels = c("Ring1b_CKO", "BAP1_specific"))
  derepressed_gene_tab$PcG_binding <- factor(derepressed_gene_tab$PcG_binding, 
                                             levels = c("High", "Low"))
  
  pacman::p_load(ggalluvial)
  
  ggplot(derepressed_gene_tab,
         aes(y = n, axis1 = BAP1_specific, axis2 = PcG_binding)) +
    geom_alluvium(aes(fill = BAP1_specific), width = 1/14) +
    geom_stratum(width = 1/10, fill = "grey45", color = "white", lwd = 1.5) +
    scale_x_discrete(limits = c("Up-regulation\noverlap", "PcG binding"), 
                     expand = c(.05, .1)) +
    geom_label(stat = "stratum", infer.label = TRUE, label.size = 0, 
               size = 4, label.padding = unit(0.2, "lines")) +
    scale_fill_manual(values = c("cyan", "grey50")) +
    ggtitle("Derepressed genes (n=2802)") +
    theme_setting +
    theme(panel.border = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  
  ggsave(filename = "FigS2_intergenic_nascent_TUs_overlaps.pdf", 
         path = "../figS2/figs", width = 5, height = 5)
}

# highly expressed genes ------------------------------------------------------------------------------
message("get highly expressed genes")
filenames <- sort(list.files('../data/kallisto_output_BGI/', full.names = T)) # count with combined reference GENCODE vM20 and ncRNA annotation
sampleNewName <- gsub(".*/", "\\2", filenames)

.design <- c("P0_R1", "P0_R2", "P0_R3",
             "P12_R1", "P12_R2", "P12_R3",
             "P12_C12_R1", "P12_C12_R2", "P12_C12_R3",
             "P12_C24_R1", "P12_C24_R2", "P12_C24_R3",
             "P12_C36_R1", "P12_C36_R2", "P12_C36_R3")
names(.design) <- 
  sampleNewName[order(gsub("(^.).*", "\\1", sampleNewName), as.numeric(gsub("^.", "", sampleNewName)))]

# get read counts for DESeq2
count_table_RNA_tpm <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                          as = 'matrix', what = "tpm")
gene_tpm_RNA <- count_table_RNA_tpm[grepl("^ENS", rownames(count_table_RNA_tpm)), names(.design)] %>%
  keepOneTx(rowname_gene_id = TRUE, is_gene_sum = TRUE)


highly_expressed_genes <- rownames(gene_tpm_RNA)[rowMeans(gene_tpm_RNA[, 1:3]) > quantile(rowMeans(gene_tpm_RNA[, 1:3]), 0.95)]
highly_expressed_genes <- highly_expressed_genes[!highly_expressed_genes %in% derepressed_gene_ids]

saveRDS(highly_expressed_genes, "data/highly_expressed_genes.rds")



# ---------------------------------- ChIP log2FC cross chromHMM states# ---------------------------------- #

ChromHMM <- importRanges("../data/mESC_E14_12_dense.annotated.mm10.bed")
ChromHMM$name <- gsub(".*_", "\\2", ChromHMM$name)

ChromHMM_mm9 <- liftOver(ChromHMM, import.chain("../data/liftOver_chains/mm10ToMm9.over.chain")) %>% unlist()

biv_idx <- countSubjectHits(findOverlaps(ChromHMM_mm9[ChromHMM_mm9$name == "BivalentChromatin"], 
                                         promoters(gene.gr, upstream = 200, downstream = 200))) > 0
bivalent_genes_chromHMM <- gene.gr$gene_id[biv_idx] # 4737






setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../util/utils.R")
source("../util/getCoverage.R")

gene.gr <- readRDS("../data/gene.gr.RData")
gene.gr$gene_id <- gsub("\\..*", "", gene.gr$gene_id)
names(gene.gr) <- gene.gr$gene_id

# --------------------------------------------------------------------------------------------------- #
# test chromatin TAD
library(readxl)
PAD <- readxl::read_xlsx("../data/mmc3_Polycomb group proteins regulate chromatin architecture in mouse oocytes and early embryos.xlsx")

PAD.gr <- 
  PAD[-(1:2), ] %>% 
  `colnames<-`(., c("chr", "start", "end")) %>% 
  as.data.frame() %>%
  mutate(strand = "*", score = 0) %>%
  GenomicRanges::makeGRangesFromDataFrame()

TAD_mm9.gr <- import("../data/TAD.Ren.total.HindIII.combined.domain_mm9.bed")
TAD_mm9.gr$H3K27me3_b1 <- .countBW("../bw_rep/batch1/input_norm/H3K27m3_NT.bw", TAD_mm9.gr)[, 1]
TAD_mm9.gr$H3K27me3_b2 <- .countBW("../bw_rep/batch2/input_norm/H3K27me3_P0_pool.bw", TAD_mm9.gr)[, 1]
TAD_mm9.gr$H3K27me3_enriched <- FALSE
TAD_mm9.gr$H3K27me3_enriched[intersect(kink_index(TAD_mm9.gr$H3K27me3_b1, method = "log-normal", .p = 0.1),
                                       kink_index(TAD_mm9.gr$H3K27me3_b2, method = "log-normal", .p = 0.1))] <- TRUE

# annotate TAD: 1. with genes 2. overlap with PAD 3. # of de-repressed genes
TAD_mm9.gr$num_gene <- findOverlaps(gene.gr, TAD_mm9.gr) %>% countSubjectHits()
TAD_mm9.gr$num_H2Aub_gene <- findOverlaps(gene.gr[H2Aub_repressed_genes], TAD_mm9.gr) %>% countSubjectHits()

mtch <- findOverlaps(TAD_mm9.gr, PAD.gr)
mtch_list <- split(subjectHits(mtch), queryHits(mtch))
TAD_mm9.gr$PAD_ov <- 0
PAD_ov <- sapply(names(mtch_list), 
                 function(i) {
                   sum(width(PAD.gr[mtch_list[[i]]])) / width(TAD_mm9.gr[as.numeric(i)])})
TAD_mm9.gr$PAD_ov[as.numeric(names(mtch_list))] <- PAD_ov
  

.get_combination <- function(len) {
  out <- cbind(rep(seq_len(len), len),
               rep(seq_len(len), each = len))
  out[out[, 1] < out[, 2], ]
}

H2Aub_repressed_genes_2 <- intersect(H2Aub_repressed_genes, rownames(log2FC_mat))
mtch <- findOverlaps(TAD_mm9.gr, gene.gr[H2Aub_repressed_genes_2])
tad_pair_log2FC_tab <- NULL
for (i in unique(queryHits(mtch))) {
  tmp <- mtch[queryHits(mtch) == i]
  if (length(tmp) > 1) {
    mat <- matrix(.get_combination(length(tmp)), ncol = 2)
    tad_pair_log2FC_tab <- rbind(tad_pair_log2FC_tab,
                                 cbind(mean(log2FC_mat[H2Aub_repressed_genes_2[subjectHits(tmp)[mat[, 1]]], ]),
                                       mean(log2FC_mat[H2Aub_repressed_genes_2[subjectHits(tmp)[mat[, 2]]], ]),
                                       ifelse(TAD_mm9.gr$PAD_ov[i] > 0.5, "PAD", "Other") )
                                 )
  } 
}


tad_pair_log2FC_tab <- tad_pair_log2FC_tab %>%
  as.data.frame() %>%
  `colnames<-`(., c("Gene", "Pair", "PAD")) %>%
  mutate(Gene = as.numeric(as.character(Gene)),
         Pair = as.numeric(as.character(Pair)))

ggplot(tad_pair_log2FC_tab, aes(x = Gene, y = Pair)) +
  geom_point(color = "blue3") +
  geom_point(data = tad_pair_log2FC_tab[tad_pair_log2FC_tab$PAD == "PAD", ], 
             mapping = aes(x = Gene, y = Pair), color = "red") +
  annotate("text", x = -Inf, y = Inf,
           hjust = -0.5, vjust = 1.2, color = "blue3",
           label = paste0("All pairs\n",
                          "r = ", round(with(tad_pair_log2FC_tab, cor(Gene, Pair)), 3),
                          "\nn = ", nrow(tad_pair_log2FC_tab))) +
  annotate("text", x = Inf, y = Inf,
           hjust = 1.2, vjust = 1.2, color = 'red',
           label = paste0("PAD pairs\n",
                          "r = ", round(with(tad_pair_log2FC_tab[tad_pair_log2FC_tab$PAD == "PAD", ], cor(Gene, Pair)), 3),
                          "\nn = ", sum(tad_pair_log2FC_tab$PAD == "PAD"))) +
  ggtitle("Gene de-repression within TAD") + xlab("Gene_1 RNA-seq log2FC") + ylab("Gene_2 RNA-seq log2FC") +
  theme_setting +
  theme(legend.position = "none")

ggsave("Fig_S8_TAD_derepression_association.png", 
       path = "../figS8/figs", 
       device = "png", width = 5, height = 4)


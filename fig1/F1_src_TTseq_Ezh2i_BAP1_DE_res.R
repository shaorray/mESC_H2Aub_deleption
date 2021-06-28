library(DESeq2)

load_TTseq_DE_res <- function() {
  
  filenames <- sort(list.files('../data/kallisto_output', full.names = T)) # count with combined reference GENCODE vM20 and ncRNA annotation
  sampleNewName <- gsub(".*/", "\\2", filenames)
  
  # TPM
  TPM_table <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                    as = 'matrix', what = "tpm")
  
  colnames(TPM_table) <- sampleNewName
  
  txTPM <- TPM_table[grepl("^ENS", rownames(TPM_table)), ] %>%
    keepOneTx(rowname_gene_id = T)
  spTPM <- TPM_table[grepl("chrS(2|4|5|8)", rownames(TPM_table)), ]
  txTPM_norm <- sweep(txTPM, 2, size_factor_cal(spTPM), "/")
  TT_seq_txTPM_norm_cmb <<- sapply(unique(gsub("_(1|2)$", "", colnames(txTPM_norm))),
                                      function(x) 
                                        txTPM_norm[, grep(x, colnames(txTPM_norm))] %>%
                                        rowMeans())
  # read count
  count_table <- SummarizedExperiment::readKallisto(paste0(filenames, "/abundance.tsv"),
                                                    as = 'matrix', what = "est_counts")
  
  colnames(count_table) <- sampleNewName
  
  txCounts <- count_table[grepl("^ENS", rownames(count_table)), ] %>%
    keepOneTx(rowname_gene_id = T)
  
  dds <- DESeqDataSetFromMatrix(round(as.matrix(txCounts)),
                                colData = data.frame(condition = gsub("_(1|2)$", "\\2", sampleNewName) ),
                                design = ~ condition)
  dds <- DESeq(dds)
  
  get_res <- function(dds) {
    res_E0B6 <<- results(dds, contrast = c("condition", "mESC_Ezh2i_0d_BAP1_6h_2", "mESC_Ezh2i_0d_BAP1_0h_8"))
    res_E0B12 <<- results(dds, contrast = c("condition", "mESC_Ezh2i_0d_BAP1_12h_4", "mESC_Ezh2i_0d_BAP1_0h_8"))
    
    res_E1B6 <<- results(dds, contrast = c("condition", "mESC_Ezh2i_1d_BAP1_6h_1", "mESC_Ezh2i_0d_BAP1_0h_8"))
    res_E1B12 <<- results(dds, contrast = c("condition", "mESC_Ezh2i_1d_BAP1_12h_3", "mESC_Ezh2i_0d_BAP1_0h_8"))
    
    res_E1B0 <<- results(dds, contrast = c("condition", "mESC_Ezh2i_1d_BAP1_0h_5", "mESC_Ezh2i_0d_BAP1_0h_8"))
    res_E2B0 <<- results(dds, contrast = c("condition", "mESC_Ezh2i_2d_BAP1_0h_6", "mESC_Ezh2i_0d_BAP1_0h_8"))
    res_E7B0 <<- results(dds, contrast = c("condition", "mESC_Ezh2i_7d_BAP1_0h_7", "mESC_Ezh2i_0d_BAP1_0h_8"))
    
    res.list <- list(res_E0B6, res_E0B12, res_E1B6, res_E1B12, res_E1B0, res_E2B0, res_E7B0)
    names(res.list) <- c("res_E0B6", "res_E0B12", "res_E1B6", "res_E1B12",
                         "res_E1B0", "res_E2B0", "res_E7B0")
    return(res.list)
  }
  TTseq_res.list <<- get_res(dds)
  
  # TU DE
  tuCounts <- count_table[!grepl("^ENS|^chrS", rownames(count_table)), ]
  
  dds <- DESeqDataSetFromMatrix(round(as.matrix(tuCounts)),
                                colData = data.frame(condition = gsub("_(1|2)$", "\\2", sampleNewName) ),
                                design = ~ condition)
  dds <- DESeq(dds)
  
  TTseq_tu_res.list <<- get_res(dds)
}

load_TTseq_DE_res()


if (F) {
  g_list <- list(plot_ma(baseMean = res_E0B6$baseMean,
                         log2FC = res_E0B6$log2FoldChange, 
                         pval = res_E0B6$padj, 
                         title = "TU: P6"),
                 
                 plot_ma(baseMean = res_E0B12$baseMean,
                         log2FC = res_E0B12$log2FoldChange, 
                         pval = res_E0B12$padj, 
                         title = "TU: P12"),
                 
                 plot_ma(baseMean = res_E1B0$baseMean,
                         log2FC = res_E1B0$log2FoldChange, 
                         pval = res_E1B0$padj, 
                         title = "TU: Ezh2i 1d"),
                 
                 plot_ma(baseMean = res_E2B0$baseMean,
                         log2FC = res_E2B0$log2FoldChange, 
                         pval = res_E2B0$padj, 
                         title = "TU: Ezh2i 2d"),
                 
                 plot_ma(baseMean = res_E7B0$baseMean,
                         log2FC = res_E7B0$log2FoldChange, 
                         pval = res_E7B0$padj, 
                         title = "TU: Ezh2i 7d"),
                 
                 plot_ma(baseMean = res_E1B6$baseMean,
                         log2FC = res_E1B6$log2FoldChange, 
                         pval = res_E1B6$padj, 
                         title = "TU: Ezh2i 1d P6"),
                 
                 plot_ma(baseMean = res_E1B12$baseMean,
                         log2FC = res_E1B12$log2FoldChange, 
                         pval = res_E1B12$padj, 
                         title = "TU: Ezh2i 1d P12"))
  
  
  ggsave(plot = do.call(grid.arrange,
                        c(g_list,
                          ncol = 4)),
         filename = paste0("FigS1_TU_DE_TT-seq.png"),
         path = "../figS1/figs/",
         device = "png", width = 10, height = 6)
  
  
}

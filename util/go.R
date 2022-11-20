library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ggpubr)

enrichGeneSets <- function(gene_id, 
                           method = "GO", 
                           title = NULL, 
                           ontology = "BP", 
                           top_n_term = 10, 
                           colorset = viridis::viridis(10),
                           is.GeneRatio = T,
                           col.limits = c(0, 1))
{ # Agrs:
  # gene_id: ENSEMBL gene id of interest
  # ontology: BP = Biological Process; CC = Cellular Component; MF = Molecular Function
  
  gene_id <- as.character(gene_id)
  all_gene_ids <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)$gene_id
  
  # convert to ENTREZ_GENE_ID
  if (grepl("^EN", gene_id[1])) {
    res <- biomaRt::select(org.Mm.eg.db,
                           keys = gene_id,
                           keytype = "ENSEMBL", 
                           columns = "ENTREZID")
    gene_id <- res$ENTREZID[match(gene_id, res$ENSEMBL)]
  } else {
    res <- biomaRt::select(org.Mm.eg.db,
                           keys = gene_id,
                           keytype = "SYMBOL", 
                           columns = "ENTREZID")
    gene_id <- res$ENTREZID[match(gene_id, res$SYMBOL)]
  }
  
  stopifnot(length(gene_id) > 0) 
  
  
  if (method == 'GO')
  { 
    ego <- enrichGO(gene          = gene_id,
                    universe      = all_gene_ids,
                    OrgDb         = org.Mm.eg.db,
                    ont           = ontology,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
  }
  
  if (method == 'DAVID')
  {
    # if(!require(RDAVIDWebService)) BiocManager::install("RDAVIDWebService")
    # if(!require(rJava)) install.packages("rJava")
    # library(rJava)
    # library(RDAVIDWebService)
    # 
    # ego <- enrichDAVID(gene          = gene_id,
    #                    idType        = "ENTREZ_GENE_ID",
    #                    universe      = all_gene_ids,
    #                    minGSSize     = 10,
    #                    maxGSSize     = 500,
    #                    annotation    = "GOTERM_BP_FAT",
    #                    pvalueCutoff  = 0.05,
    #                    pAdjustMethod = "BH", 
    #                    qvalueCutoff  = 0.2)
  }
  
  if (method == 'enricher')
  { # Cell Marker
    cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Mouse_cell_markers.txt') %>%
      tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
      dplyr::select(cellMarker, geneID) %>%
      dplyr::mutate(geneID = strsplit(geneID, ', '))
    
    ego <- enricher(gene             = gene_id,
                    pvalueCutoff     = 0.05,
                    pAdjustMethod    = "BH",
                    universe         = all_gene_ids,
                    minGSSize        = 10,
                    maxGSSize        = 500,
                    qvalueCutoff     = 0.2,
                    TERM2GENE        = cell_markers)
  }
  
  if (method == 'KEGG')
  {
    ego <- enrichKEGG(gene              = gene_id, 
                      organism          = "mmu", 
                      keyType           = "kegg",
                      pvalueCutoff      = 0.05, 
                      pAdjustMethod     = "BH", 
                      universe          = all_gene_ids,
                      minGSSize         = 10, 
                      maxGSSize         = 500, 
                      qvalueCutoff      = 0.2,
                      use_internal_data = FALSE)
  }
  
  if (method == 'MKEGG')
  {
    ego <- enrichMKEGG(gene              = gene_id, 
                       organism          = "mmu", 
                       keyType           = "kegg",
                       pvalueCutoff      = 0.05, 
                       pAdjustMethod     = "BH", 
                       universe          = all_gene_ids,
                       minGSSize         = 10, 
                       maxGSSize         = 500, 
                       qvalueCutoff      = 0.2)
  }
  
  if (is.null(title)) title = method
  
  idx <- order(ego@result$Count, decreasing = T) %>% head(top_n_term)
  dat <- with(ego@result[idx, ], 
              data.frame(Description = Description, 
                         GO = ID,
                         p.adjust = p.adjust, 
                         GeneRatio = Count / length(ego@gene),
                         Count = Count
              )
  )
  dat$log.p.adjust <- dat$p.adjust %>% as.numeric() %>% log() %>% "*"(-1)
  dat$GeneRatio <- as.numeric(dat$GeneRatio)
  
  if (is.GeneRatio) {
    ggdotchart(dat, x = "Description", y = "GeneRatio",
               color = "p.adjust",                           # Color by groups
               sorting = "descending",                       # Sort value in descending order
               add = "segments",                             # Add segments from y = 0 to dots
               rotate = T,
               dot.size = 7,                                 # Large dot size
               label = round(dat$Count),                     # Add counts as dot labels
               font.label = list(color = "grey", size = 12, 
                                 vjust = 0.5)                # Adjust label parameters
    ) + 
      scale_color_gradientn(colours = colorset, limits = col.limits) + 
      theme_setting +
      theme(legend.position = "right", 
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 12),
            plot.margin = unit(c(1,0,1,0), "cm")) +
      xlab("") + 
      ggtitle(title)
  } else {
    ggdotchart(dat, 
               x = "Description", 
               y = "log.p.adjust",
               color = "GeneRatio",
               sorting = "descending",                       
               add = "segments", 
               rotate = T,
               dot.size = 7,
               label = round(dat$Count),
               font.label = list(color = "grey", 
                                 face = "bold",
                                 size = 12, 
                                 vjust = 0.5)
    ) + 
      scale_color_gradientn(colours = colorset, limits = col.limits) + 
      theme_setting +
      theme(legend.position = "right", 
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 12),
            plot.margin = unit(c(1,0,1,0), "cm")) +
      xlab("") +ylab("-log(p.adjust)") + 
      ggtitle(title)
  }
}


enriched_GO_list_heatmap <- function(gene_id_list, 
                                     title = NULL, 
                                     ontology = "BP", 
                                     top_n_term = 10, 
                                     colorset = viridis::viridis(10),
                                     is.GeneRatio = T) {
  
  all_gene_ids <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)$gene_id
  
  # convert to ENTREZ_GENE_ID
  res <- biomaRt::select(org.Mm.eg.db,
                         keys = all_gene_ids,
                         keytype = "ENTREZID", 
                         columns = c("ENSEMBL", "SYMBOL"))
  
  ego_list <- lapply(gene_id_list, function(gene_id) {
    
    if (grepl("^EN", gene_id)) {
      gene_res <- res$ENSEMBL
    } else {
      gene_res <- res$SYMBOL
    }
    
    enrichGO(gene          = res$ENTREZID[match(as.character(gene_id), gene_res)],
             universe      = all_gene_ids,
             OrgDb         = org.Mm.eg.db,
             ont           = ontology,
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.01,
             qvalueCutoff  = 0.05)
  })
  
  # keep go terms appear in all list gene groups
  go_ids <- lapply(ego_list, function(ego) {
    ego@result$ID
  }) %>% 
    unlist() %>% 
    table() %>% 
    as.data.frame() %>% 
    dplyr::filter(Freq == length(ego_list)) %>%
    dplyr::select('.') %>% 
    as.matrix() %>% 
    as.character()
  
  go_ratios <- lapply(ego_list, function(ego) {
    with(ego@result[go_ids, ], lapply(strsplit(GeneRatio, "\\/"), 
                                      function(x) {x = as.numeric(x)
                                      x[1] / x[2]}) %>% unlist())
  }) %>% Reduce(cbind, .) %>% 
    `rownames<-`(., go_ids)
  
  
  idx <- order(rowSums(go_ratios), decreasing = T) %>% head(top_n_term) %>% rev()
  go_ids <- unique(go_ids[idx])
  
  dat <- lapply(names(ego_list), function(sample_name) {
    ego <- ego_list[[sample_name]]
    with(ego@result[go_ids, ], 
         data.frame(Sample = sample_name,
                    Description = Description,
                    GO = ID, 
                    p.adjust = p.adjust,
                    Gene_Ratio = lapply(strsplit(GeneRatio, "\\/"), 
                                       function(x) {x = as.numeric(x)
                                       x[1] / x[2]}) %>% unlist()
                    )
         )
  }) %>% Reduce(rbind, .)
  
  dat$Description <- factor(dat$Description, levels = ego_list[[1]]@result[go_ids, "Description"])
  dat$`-log p.adj` <- dat$p.adjust %>% as.numeric() %>% log() %>% "*"(-1)
  dat$Gene_Ratio <- as.numeric(dat$Gene_Ratio)
  
  ggplot(dat, aes(x = Sample, y = Description)) +
    geom_tile(fill = NA) +
    geom_point(aes(size = Gene_Ratio, color = `-log p.adj`)) +
    scale_color_gradientn(colours = colorset) + 
    guides(size = guide_legend(title = "Gene Ratio")) +
    xlab("") + ylab("") + ggtitle(title) +
    theme_setting +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

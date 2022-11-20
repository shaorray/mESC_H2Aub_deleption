

# "log2FC_mat" is from "Fig2_gene_derepression_and_ChIP_targets.R"
# ---------------------------------- cluster log2FC from DESeq2 -------------------------------- #
order_cls <- function(cls, val) {
  factor(cls, levels = order(aggregate(val, list(cls), mean)[, 2], decreasing = T)) %>%
    as.numeric()
}

if (FALSE) {
  # kmeans
  set.seed(1)
  log2FC_cls <- kmeans(log2FC_mat, centers = 4)$cluster
  
  log2FC_cls <- order_cls(log2FC_cls, rowMeans(log2FC_mat))
  names(log2FC_cls) <- rownames(log2FC_mat)
  
  genes_derepressed <- rownames(log2FC_mat)[log2FC_cls == 1]
  genes_W_derepressed <- rownames(log2FC_mat)[log2FC_cls == 2]
  genes_unchanged <- rownames(log2FC_mat)[log2FC_cls == 3]
  genes_repressed <- rownames(log2FC_mat)[log2FC_cls == 4]
}


if (FALSE) {
  # overlap with PcG target peaks
  use_gene.gr <- gene.gr[use_gene_ids]
  gene_cls_H2Aub <- log2FC_cls[use_gene_ids]
  use_gene.gr$Type <- c("Derepressed", "Weakly\nDerepressed", "Unchanged", "Repressed")[gene_cls_H2Aub]
  
  use_gene.gr$PcG <- "None"
  use_gene.gr$PcG[use_gene.gr$gene_id %in% PcG_binding_genes] <- "Target" # "PcG_binding_genes" is from "F2_src_ChIP_target_genes.R"
  
  PcG_counts <- mcols(use_gene.gr) %>% as.data.frame() %>%
    dplyr::count(Type, PcG, sort = TRUE)
  
  PcG_counts$Type <- factor(PcG_counts$Type, 
                            levels = c("Derepressed", "Weakly\nDerepressed", "Unchanged", "Repressed"))
  PcG_counts$PcG <- factor(PcG_counts$PcG,
                           levels = c("Target", "None"))
  
  library("ggforce")
  dat_ggforce <- PcG_counts  %>%
    ggforce::gather_set_data(1:2) %>%
    arrange(x, Type, desc(PcG))
  
  ggplot(dat_ggforce, aes(x = x, id = id, split = y, value = n)) +
    geom_parallel_sets(aes(fill = Type), alpha = 0.8, axis.width = 0.2,
                       n=100, strength = 0.5) +
    geom_parallel_sets_axes(axis.width = 0.25, fill = "gray90",  size = 0) +
    geom_parallel_sets_labels(colour = 'gray35', size = 0, angle = 0, fontface="bold") +
    scale_fill_manual(values  = colors_n[c(2,3,6,10)]) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x  = element_blank()
    )
  
  ggsave(filename = "Fig2_DE_group_PcG_overlaps.pdf", 
         path = "../fig2/figs", width = 2.5, height = 5)
}



# -------------------------------------- cluster by derepression time --------------------------------------- #

derepressed_gene_ids_fig4 <- Reduce(c, DE_gene.list[grep("up", names(DE_gene.list))]) %>% 
  unique() %>% intersect.Vector(., names(gene.gr)) %>% intersect.Vector(., rownames(geneRC_RNA.norm.rlog)) 
length(derepressed_gene_ids_fig4) # 2892

BAP1_specific_genes <- derepressed_gene_ids_fig4[derepressed_gene_ids_fig4 %ni% c(GSE132753_DE_genes_CKO, GSE134053_DE_genes_CKO)]
length(BAP1_specific_genes) # 1327

H2Aub_repressed_genes_fig4 <- derepressed_gene_ids_fig4[derepressed_gene_ids_fig4 %ni% BAP1_specific_genes]


# derepression rate is defined by linear coefficient given the 5 PC RNA-seq read counts (P0, P12, P12C12, P12C24, P12C36)
get_derepression_rate <- function(mat) {
  mat = as.matrix(mat)
  X = cbind(1, seq_len(nrow(mat)))
  (solve(t(X) %*% X) %*% t(X) %*% mat)[2, 1]
}

get_derepression_rel_residual <- function(mat) {
  mat = as.matrix(mat)
  X = cbind(1, seq_len(nrow(mat)))
  sum(abs((X %*% (solve(t(X) %*% X) %*% t(X) %*% mat)) - mat)) / mean(mat)
}


# direct calculation from RPK
RNA_derepression_rate <- apply(geneRC_RNA.norm.rlog[derepressed_gene_ids_fig4, ], 1, get_derepression_rate)

class_genes_rates <- `names<-`(cut(RNA_derepression_rate, 
                                   breaks = c(-Inf, 0.07, 0.3, Inf), labels = c("Early", "Middle", "Late")),
                               names(RNA_derepression_rate))

# combine with depressed genes
derepressed_gene_times <- rep(1, length(derepressed_gene_ids_fig4))
for (i in 4:1) {
  idx <- derepressed_gene_ids_fig4 %in% DE_gene.list[grep("up", names(DE_gene.list))][[i]]
  derepressed_gene_times[idx] <- i
}

derepressed_gene_tab <- data.frame(gene_id = derepressed_gene_ids_fig4,
                                   BAP1_specific = derepressed_gene_ids_fig4 %in% BAP1_specific_genes,
                                   emergence_time = derepressed_gene_times,
                                   derepression_rate = class_genes_rates)

derepressed_gene_tab$derepression_rate <- factor(derepressed_gene_tab$derepression_rate,
                                                 levels = c("Early", "Middle", "Late"))

# prepare heatmap parameters, set order and group by derepression groups
row_order <- NULL
for (i in rev(c("Early", "Middle", "Late"))) {
  tmp_idx <- which(derepressed_gene_tab$derepression_rate == i)
  # row_order %c=% tmp_idx[cov_mat_H2Aub[[1]][tmp_idx, ] %>% 
  #                          trim_quantile() %>%
  #                          log1p() %>%
  #                          rowMeans(., na.rm = T) %>% 
  #                          order(.)]
  row_order %c=% tmp_idx[log2FC_mat[tmp_idx, ] %>%
                           rowMedians(., na.rm = T) %>%
                           order(.)]
}
.breaks = derepressed_gene_tab[, 3] %>% 
  table() %>% unname() %>% rev() %>% cumsum()

# --------------------------------------------------------------------------------------------- #
# # plot heatmap
gene_cls_H2Aub <- class_genes_rates[H2Aub_repressed_genes_fig4]
log2FC_mat_derepressed_fig4 <- log2FC_mat[H2Aub_repressed_genes_fig4, ]


log2FC_mat_de <- sapply(res.gene.list, function(x) x$log2FoldChange)
rownames(log2FC_mat_de) <- rownames(res.gene.list[[1]])
log2FC_mat_de <- log2FC_mat_de[derepressed_gene_ids, ]

log2FC_mat_de[row_order, ] %>%
  `rownames<-`(NULL) %>% 
  reshape::melt() %>% 
  ggplot(aes(x = X2, y = X1, fill = value)) +
  geom_tile() +
  labs(title = '',  x = '', y = '') +
  scale_fill_gradientn(name = expression(log[2]~FC), 
                       colours = c('blue', "white", 'red', 'red4'), 
                       values = (c(min(log2FC_mat_de), 0, max(log2FC_mat_de) / 3, max(log2FC_mat_de)) -
                                   min(log2FC_mat_de)) /
                         diff(range(log2FC_mat_de)),
                       guide = guide_colorbar(title.position = "top", 
                                              title.hjust = 0.5)
  ) +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), 
        strip.text.y = element_blank(), 
        legend.position = "bottom",
        legend.direction = "horizontal", 
        legend.key.size = unit(13, "pt"), 
        plot.margin = unit(c(0, 1, 0, 1), "lines")) +
  scale_y_break(c(.breaks[1], .breaks[1])) +
  scale_y_break(c(.breaks[2], .breaks[2])) 

ggsave(filename = paste0("../figs3/figs/FigS3_tmp_TSS_Heatmap_log2FC_derepressed_gene.png"), 
       device = "png", width = 1.5, height = 8)

derepressed_gene_RPK <- txRPK_RNA.norm[derepressed_gene_ids, c("P32", "P33", "P41")] %>% rowMeans()

(derepressed_gene_RPK[row_order] + 0.1) %>%
  log2() %>%
  as.matrix(ncol = 1) %>%
  `rownames<-`(NULL) %>% 
  reshape::melt() %>% 
  ggplot(aes(x = X2, y = X1, fill = value)) +
  geom_tile() +
  labs(title = '',  x = '', y = '') +
  scale_fill_gradientn(name = expression(log[2]~RPK), 
                       colours = c('brown', "white", 'deepskyblue', 'deepskyblue4'), 
                       values = (c(min(derepressed_gene_RPK), 0, max(derepressed_gene_RPK) / 1.5, max(derepressed_gene_RPK)) -
                                   min(derepressed_gene_RPK)) /
                         diff(range(derepressed_gene_RPK)),
                       guide = guide_colorbar(title.position = "top", 
                                              title.hjust = 0.5)) +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # plot.margin = margin(0,1,0, 1, "lines"),
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), 
        strip.text.y = element_blank(), 
        legend.position = "bottom",
        legend.direction = "horizontal", 
        legend.key.size = unit(13, "pt"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "white")) +
  scale_y_break(c(.breaks[1], .breaks[1])) +
  scale_y_break(c(.breaks[2], .breaks[2])) 


ggsave(filename = paste0("../fig4/figs/Fig4_tmp_TSS_Heatmap_log2_RPK_derepressed_gene.png"), 
       device = "png", width = 1.5, height = 8)
ggsave(filename = paste0("../fig4/figs/Fig4_tmp_TSS_Heatmap_log2_RPK_derepressed_gene2.png"), 
       device = "png", width = 0.7, height = 8)

# -------------------------------------------------------------------------------------------- #
# plot log2FC heatmap


cbind(log2FC_mat_derepressed_fig4, 
      cls = gene_cls_H2Aub) %>% 
  as.data.frame() %>%
  arrange(cls, rowMedians(log2FC_mat_derepressed_fig4)) %>%
  cbind(id = seq_len(nrow(log2FC_mat_derepressed_fig4))) %>%
  `rownames<-`(NULL) %>% 
  reshape::melt(id.vars = c("cls", "id")) %>% 
  ggplot(aes(x = variable, y = id, fill = value)) +
  geom_tile(width = 0.95) +
  facet_grid(cls ~ ., scales = "free", space = "free") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = paste0('H2Aub repressed genes (n=', nrow(log2FC_mat_derepressed_fig4), ')'), 
       x = '', y = '') +
  scale_fill_gradientn(name = "Log2FC", 
                       colours = c('blue', "white", 'red'), 
                       values = (c(min(log2FC_mat_derepressed_fig4), 0, max(log2FC_mat_derepressed_fig4)) - min(log2FC_mat_derepressed_fig4)) /
                         diff(range(log2FC_mat_derepressed_fig4)) ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, size = 1),
        axis.ticks.x = element_blank(), 
        axis.text.x = element_text(size = 10, hjust = 1, angle = 45),
        axis.text.y = element_blank(), 
        strip.text.y = element_blank(), 
        legend.position = "none",
        legend.direction = "horizontal", 
        legend.key.size = unit(13, "pt"), 
        plot.margin = unit(c(1, 1, 0, 1), "lines")) 

ggsave(filename = paste0("Fig4_RNAseq_log2FC_class_Heatmap_no_legend.png"), 
       path = "../fig4/figs/",
       device = "png", width = 3.2, height = 5)

# RPK color band
used_gene_RPK <- geneRPK_RNA[H2Aub_repressed_genes_fig4, c("P32", "P33", "P41")] %>% rowMeans()
used_gene_RPK[used_gene_RPK < 0.004] <- 0
used_gene_RPK <- used_gene_RPK %>% trim_quantile() %>% log2() 

used_gene_RPK %>% 
  cbind(cls = gene_cls_H2Aub) %>%
  as.data.frame() %>%
  arrange(cls, rowMedians(log2FC_mat_derepressed_fig4)) %>%
  cbind(id = seq_along(used_gene_RPK)) %>%
  `rownames<-`(NULL) %>% 
  reshape::melt(id.vars = c("cls", "id")) %>% 
  
  ggplot(aes(x = variable, y = id, fill = value)) +
  geom_tile(width = 0.95) +
  facet_grid(cls ~ ., scales = "free", space = "free") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  
  labs(title = '',  x = '', y = '') +
  scale_fill_gradientn(name = expression(log[2]~RPK),
                       colours = brewer.pal(8, "YlGnBu")) +
  # scale_fill_gradientn(name = expression(log[2]~RPK), 
  #                      colours = c('brown', "white", 'deepskyblue', 'deepskyblue4'),
  #                      values = (c(min(used_gene_RPK), 0, max(used_gene_RPK) / 1.5, max(used_gene_RPK)) -
  #                                  min(used_gene_RPK)) /
  #                        diff(range(used_gene_RPK)),
  #                      guide = guide_colorbar(title.position = "top",
  #                                             title.hjust = 0.5)) +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0,1,0, 1, "lines"),
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), 
        strip.text.y = element_blank(), 
        legend.position = "bottom",
        legend.direction = "horizontal", 
        legend.key.size = unit(13, "pt"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "white"))

ggsave(filename = paste0("Fig4_Heatmap_All_Genes_RPK.png"), 
       path = "../fig4/figs",
       device = "png", width = 1, height = 8)

ggsave(filename = paste0("Fig4_Heatmap_All_Genes_RPK_legend.png"), 
       path = "../fig4/figs",
       device = "png", width = 3, height = 8)


# --------------------------------------- add TSS coverage --------------------------------------- #
# "cov_mats" are from "F2_src_gene_tss_coverage_heatmap.R"

plot_TSS_heatmap_2 <- function(mat, .cls, ref_mat,
                               .title = "", 
                               .cut_outlier = 0.995,
                               .col = "black", 
                               bg_val = 0.1,
                               .show_legend = F) 
{
  
  mat <- trim_quantile(mat, q = .cut_outlier)
  mat <- mat - bg_val
  mat[mat < 0] <- 0
  .max <- max(mat)
  .min <- 0
  
  mat %>%
    cbind(cls = .cls) %>% 
    as.data.frame() %>%
    arrange(cls, rowMedians(ref_mat)) %>%
    cbind(id = seq_len(nrow(mat))) %>%
    `rownames<-`(NULL) %>% 
    reshape::melt(id.vars = c("cls", "id")) %>% 
    dplyr::mutate(variable = as.numeric(variable)) %>%
    ggplot(aes(x = variable, y = id, fill = value)) +
    geom_tile() +
    facet_grid(cls ~ ., scales = "free", space = "free") +
    
    scale_x_continuous(expand = c(0, 0),
                       breaks = c(1, 51, 100),
                       # labels = c('', '', '')
                       labels = c('-15 kb', 'TSS', '+15 kb')
    ) +
    scale_y_continuous(expand = c(0,0)) +
    
    scale_fill_gradientn(name = "",
                         colours = c(add.alpha("white", 0.9), .col, add.alpha("black", 0.9)),
                         limits = c(.min, .max),
                         breaks = c(.min, (.min + .max) / 2, .max),
                         labels = round(c(0, .max / 2, .max), 1),
                         guide = guide_colorbar(label = TRUE,
                                                draw.ulim = TRUE,
                                                draw.llim = TRUE,
                                                frame.colour = "black",
                                                ticks = TRUE,
                                                nbin = 10,
                                                label.position = "bottom",
                                                barwidth = 6.8,
                                                barheight = 0.8,
                                                direction = 'horizontal')
                         ) +
    
    labs(title = .title,  x = '', y = '') +
    
    theme_minimal() +
    theme(panel.border = element_rect(fill = NA, size = 1),
          panel.grid = element_blank(),
          axis.ticks.x = element_line(), 
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = unit(c(0, 1, 0, 1), "lines"),
          legend.position = "bottom", 
          strip.text.y = element_blank(), 
          legend.justification = "center")
}

ggsave(filename = paste0("Fig4_TSS_Heatmap_Genes_Group_H2Aub.png"), 
       plot = plot_TSS_heatmap_2(mat = cov_mat_H2Aub_1,  
                                 .cls = gene_cls_H2Aub, 
                                 ref_mat = cov_mat_H2Aub_2, 
                                 bg_val = 0.1,
                                 .col = sample_colors["H2Aub"]),
       path = "../fig4/figs/", dpi = 320,
       device = "png", width = 2, height = 5)


ggsave(filename = paste0("Fig4_TSS_Heatmap_Genes_Group_H3K27me3.png"), 
       plot = plot_TSS_heatmap_2(mat = cov_mat_H3K27me3_1,  
                                 .cls = gene_cls_H2Aub,
                                 ref_mat = cov_mat_H2Aub_2, 
                                 bg_val = 0.1,
                                 .col = sample_colors["H3K27me3"]),
       path = "../fig4/figs/", dpi = 320,
       device = "png", width = 2, height = 5)

ggsave(filename = paste0("Fig4_TSS_Heatmap_Genes_Group_Ring1b.png"), 
       plot = plot_TSS_heatmap_2(mat = cov_mat_Ring1b_1,  
                                 .cls = gene_cls_H2Aub, 
                                 .cut_outlier = 0.99,
                                 ref_mat = cov_mat_H2Aub_2, 
                                 bg_val = 0.1,
                                 .col = sample_colors["Ring1b"]),
       path = "../fig4/figs/", dpi = 320,
       device = "png", width = 2, height = 5)

ggsave(filename = paste0("Fig4_TSS_Heatmap_Genes_Group_Ezh2.png"), 
       plot = plot_TSS_heatmap_2(mat = cov_mat_Ezh2_1,  
                                 .cls = gene_cls_H2Aub,
                                 ref_mat = cov_mat_H2Aub_2, 
                                 bg_val = 0.1,
                                 .col = sample_colors["Ezh2"]),
       path = "../fig4/figs/", dpi = 320,
       device = "png", width = 2, height = 5)

ggsave(filename = paste0("Fig4_TSS_Heatmap_Genes_Group_Pol2.png"), 
       plot = plot_TSS_heatmap_2(mat = cov_mat_Pol2_1,  
                                 .cls = gene_cls_H2Aub, 
                                 .cut_outlier = 0.98,
                                 ref_mat = cov_mat_H2Aub_2, 
                                 bg_val = 0.08,
                                 .col = sample_colors["Pol II"]),
       path = "../fig4/figs/", dpi = 320,
       device = "png", width = 2, height = 5)

ggsave(filename = paste0("Fig4_TSS_Heatmap_Genes_Group_Pol2S5p.png"), 
       plot = plot_TSS_heatmap_2(mat = cov_mat_Pol2S5p_2,  
                                 .cls = gene_cls_H2Aub,
                                 ref_mat = cov_mat_H2Aub_2, 
                                 .cut_outlier = 0.99,
                                 bg_val = 0.05,
                                 .col = sample_colors["Pol II-S5p"]),
       path = "../fig4/figs/", dpi = 320,
       device = "png", width = 2, height = 5)

ggsave(filename = paste0("Fig4_TSS_Heatmap_Genes_Group_Cbx7.png"), 
       plot = plot_TSS_heatmap_2(mat = cov_mat_Cbx7_2,  
                                 .cls = gene_cls_H2Aub,
                                 ref_mat = cov_mat_H2Aub_2, 
                                 .cut_outlier = 0.98,
                                 bg_val = 0.1,
                                 .col = sample_colors["Cbx7"]),
       path = "../fig4/figs/", dpi = 320,
       device = "png", width = 2, height = 5)

ggsave(filename = paste0("Fig4_TSS_Heatmap_Genes_Group_Rybp.png"), 
       plot = plot_TSS_heatmap_2(mat = cov_mat_Rybp_2,  
                                 .cls = gene_cls_H2Aub,
                                 ref_mat = cov_mat_H2Aub_2, 
                                 .cut_outlier = 0.99,
                                 bg_val = 0.1,
                                 .col = sample_colors["Rybp"]),
       path = "../fig4/figs/", dpi = 320,
       device = "png", width = 2, height = 5)


# add top coverage profile
plot_TSS_coverage_2 <- function(mat, 
                                bg_val,
                                .cls,  
                                .max,
                                .col = 1:4) {
  
  mat <- trim_quantile(mat, q = 0.995)
  mat <- mat - bg_val
  mat[mat < 0] <- 0
  .min <- 0
  
  mat %>%
    cbind(cls = .cls) %>% 
    as.data.frame() %>%
    cbind(id = seq_len(nrow(mat))) %>%
    `rownames<-`(NULL) %>% 
    reshape::melt(id.vars = c("cls", "id")) %>% 
    group_by(cls, variable) %>%
    dplyr::summarize(Mean = mean(value, na.rm=TRUE)) %>%
    dplyr::mutate(cls = factor(cls)) %>%
    ggplot(aes(x = variable, y = Mean, group = cls, color = cls)) +
    geom_line(size = 0.8) + 
    scale_y_continuous(limits = c(.min, .max)) +
    scale_color_manual(values = .col) +
    xlab("") + ylab(paste("")) +
    theme_setting +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          legend.position = "none")
}

gene_cls_H2Aub <- class_genes_rates[H2Aub_repressed_genes_fig4]

ggsave(filename = paste0("Fig4_TSS_Coverage_Profile_Genes_Group_H2Aub.png"), 
       plot = plot_TSS_coverage_2(cov_mat_H2Aub_1,
                                  bg_val = 0.1, 
                                  .col = colors_n[c(2,3,6)],
                                  .cls = gene_cls_H2Aub,
                                  .max = 0.4),
       path = "../fig4/figs/",
       device = "png", width = 3, height = 1.5)

ggsave(filename = paste0("Fig4_TSS_Coverage_Profile_Genes_Group_H3K27me3.png"), 
       plot = plot_TSS_coverage_2(cov_mat_H3K27me3_1,
                                  bg_val = 0.1,
                                  .col = colors_n[c(2,3,6)],
                                  .cls = gene_cls_H2Aub, 
                                  .max = 0.4),
       path = "../fig4/figs/",
       device = "png", width = 3, height = 1.5)

ggsave(filename = paste0("Fig4_TSS_Coverage_Profile_Genes_Group_Ring1b.png"), 
       plot = plot_TSS_coverage_2(cov_mat_Ring1b_1,
                                  bg_val = 0.1,
                                  .col = colors_n[c(2,3,6)],
                                  .cls = gene_cls_H2Aub, 
                                  .max = .8),
       path = "../fig4/figs/",
       device = "png", width = 3, height = 1.5)

ggsave(filename = paste0("Fig4_TSS_Coverage_Profile_Genes_Group_Ezh2.png"), 
       plot = plot_TSS_coverage_2(cov_mat_Ezh2_1,
                                  bg_val = 0.1,
                                  .col = colors_n[c(2,3,6)],
                                  .cls = gene_cls_H2Aub, 
                                  .max = .3),
       path = "../fig4/figs/",
       device = "png", width = 3, height = 1.5)

ggsave(filename = paste0("Fig4_TSS_Coverage_Profile_Genes_Group_Cbx7.png"), 
       plot = plot_TSS_coverage_2(cov_mat_Cbx7_2,
                                  bg_val = 0.1,
                                  .col = colors_n[c(2,3,6)],
                                  .cls = gene_cls_H2Aub, 
                                  .max = 1.5),
       path = "../fig4/figs/",
       device = "png", width = 3, height = 1.5)

ggsave(filename = paste0("Fig4_TSS_Coverage_Profile_Genes_Group_Rybp.png"), 
       plot = plot_TSS_coverage_2(cov_mat_Rybp_2,
                                  bg_val = 0.1,
                                  .col = colors_n[c(2,3,6)],
                                  .cls = gene_cls_H2Aub, 
                                  .max = 1.2),
       path = "../fig4/figs/",
       device = "png", width = 3, height = 1.5)


ggsave(filename = paste0("Fig4_TSS_Coverage_Profile_Genes_Group_Pol2.png"), 
       plot = plot_TSS_coverage_2(cov_mat_Pol2_1[, 3:98],
                                  bg_val = 0.1,
                                  .col = colors_n[c(2,3,6)],
                                  .cls = gene_cls_H2Aub, 
                                  .max = 0.1),
       path = "../fig4/figs/",
       device = "png", width = 3, height = 1.5)


ggsave(filename = paste0("Fig4_TSS_Coverage_Profile_Genes_Group_Pol2S5p.png"), 
       plot = plot_TSS_coverage_2(cov_mat_Pol2S5p_1[, 3:98],
                                  bg_val = 0.1,
                                  .col = colors_n[c(2,3,6)],
                                  .cls = gene_cls_H2Aub, 
                                  .max = .1),
       path = "../fig4/figs/",
       device = "png", width = 3, height = 1.5)




# --------------------------------------- Plot ecdf by DE groups --------------------------------- #
if (FALSE) {
  df <- data.frame(ChIP_TSS_mat_P0[, grep("NT|P0_ALL", colnames(all_ChIP_TSS_mat))], 
                   groups = c("Derepressed", "Weakly Derepressed", "Unchanged", "Weakly Repressed")[log2FC_cls])
  df$groups <- factor(df$groups, c("Derepressed", "Weakly Derepressed", "Unchanged", "Weakly Repressed"))
  
  g3.1 <- ggplot(df, aes(log2(H2Aub_NT_RC), colour = groups)) +
    stat_ecdf(lwd = 1) +
    theme_classic() + 
    scale_color_viridis_d(option = "E", direction = -1) +
    xlim(c(-3, 3)) +
    xlab("log2 TSS density") + ylab("ECDF") +
    ggtitle("H2Aub  ") +
    theme_setting +
    theme(legend.position = "none")
  
  g3.2 <- ggplot(df, aes(log2(H3K27m3_NT_RC), colour = groups)) +
    stat_ecdf(lwd = 1) +
    theme_classic() + 
    scale_color_viridis_d(option = "E", direction = -1) +
    xlim(c(-3, 3)) +
    xlab("log2 TSS density") + ylab("ECDF") +
    ggtitle("H3K27me3  ") +
    theme_setting +
    theme(legend.position = "none")
  
  g3.3 <- ggplot(df, aes(log2(Cbx7_P0_ALL), colour = groups)) +
    stat_ecdf(lwd = 1) +
    theme_classic() + 
    scale_color_viridis_d(option = "E", direction = -1) +
    xlim(c(-3, 4)) +
    xlab("log2 TSS density") + ylab("ECDF") +
    ggtitle("Cbx7  ") +
    theme_setting +
    theme(legend.position = "none")
  
  g3.4 <- ggplot(df, aes(log2(H3K4m3_NT_RC), colour = groups)) +
    stat_ecdf(lwd = 1) +
    theme_classic() + 
    scale_color_viridis_d(option = "E", direction = -1) +
    xlim(c(-3, 3)) +
    xlab("log2 TSS density") + ylab("ECDF") +
    ggtitle("H3K4me3  ") +
    theme_setting +
    theme(legend.position = "none")
  
  g3.5 <- ggplot(df, aes(log2(Pol2.S2p_NT_RC), colour = groups)) +
    stat_ecdf(lwd = 1) +
    theme_classic() + 
    scale_color_viridis_d(option = "E", direction = -1) +
    xlim(c(-4, 4)) +
    xlab("log2 TSS density") + ylab("ECDF") +
    ggtitle("Pol2-S2p  ") +
    theme_setting +
    theme(legend.position = "none")
  
  g3.6 <- ggplot(df, aes(log2(Ring1b_P0_ALL), colour = groups)) +
    stat_ecdf(lwd = 1) +
    theme_classic() + 
    scale_color_viridis_d(option = "E", direction = -1) +
    xlim(c(-3, 6)) +
    xlab("log2 TSS density") + ylab("ECDF") +
    ggtitle("Ring1b  ") +
    theme_setting +
    theme(legend.position = "none")
  
  g3.7 <- ggplot(df, aes(log2(Ezh2_NT_RC), colour = groups)) +
    stat_ecdf(lwd = 1) +
    theme_classic() + 
    scale_color_viridis_d(option = "E", direction = -1) +
    xlim(c(-3, 4)) +
    xlab("log2 TSS density") + ylab("ECDF") +
    ggtitle("Ezh2  ") +
    theme_setting +
    theme(legend.position = "none")
  
  g3.8 <- ggplot(df, aes(log2(Rybp_P0_ALL), colour = groups)) +
    stat_ecdf(lwd = 1) +
    theme_classic() + 
    scale_color_viridis_d(option = "E", direction = -1) +
    xlim(c(-3, 4)) +
    xlab("log2 TSS density") + ylab("ECDF") +
    ggtitle("Rybp  ") +
    theme_setting +
    theme(legend.position = "none")
  
  g3.9 <- ggplot(df, aes(log2(Pol2.NTD_NT_RC), colour = groups)) +
    stat_ecdf(lwd = 1) +
    theme_classic() + 
    scale_color_viridis_d(option = "E", direction = -1) +
    xlim(c(-4, 4)) +
    xlab("log2 TSS density") + ylab("ECDF") +
    ggtitle("Pol2  ") +
    theme_setting +
    theme(legend.position = "none")
  
  g3.10 <- ggplot(df, aes(log2(Pol2.S5p_NT_RC), colour = groups)) +
    stat_ecdf(lwd = 1) +
    theme_classic() + 
    scale_color_viridis_d(option = "E", direction = -1) +
    xlim(c(-4, 4)) +
    xlab("log2 TSS density") + ylab("ECDF") +
    ggtitle("Pol2-S5p  ") +
    theme_setting +
    theme(legend.position = "none")
  
  ggsave(plot = grid.arrange(g3.1, g3.2, g3.3, g3.4, g3.5, 
                             g3.6, g3.7, g3.8, g3.9, g3.10, ncol = 5),
         filename = paste0("Fig1_ECDF_ChIP_DE_groups.png"), 
         path = "figs/",
         device = "png", width = 15, height = 6)
}
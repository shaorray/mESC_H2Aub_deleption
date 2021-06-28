# transcription factor binding sites
BiocManager::install(c("JASPAR2018", "TFBSTools"))
# remotes::install_github("da-bar/JASPAR2020")
require(JASPAR2018)
# require(JASPAR2020)
require(TFBSTools)

require(foreach)
require(doParallel)
registerDoParallel(cores = 3)
# -------------------------------------------------------------
# get mm9 gene intervals
require(dplyr)
require(GenomicRanges)
require(TxDb.Mmusculus.UCSC.mm9.knownGene)
require(GenomeInfoDb)
require(org.Mm.eg.db)
require(Biostrings)
require(BSgenome.Mmusculus.UCSC.mm9)

use_gene.gr <- promoters(gene.gr[use_gene_ids], upstream = 500, downstream = 100)

mm9_prom.seq <- getSeq(BSgenome.Mmusculus.UCSC.mm9, use_gene.gr)

# -------------------------------------------------------------

# using species human (9606) or mouse (10090)
pfm <- TFBSTools::getMatrixSet(JASPAR2018, list(species=10090))

# JASPAR2018 database doesn't have Nangog PFM
# download directly http://jaspar2020.genereg.net/matrix/UN0383.1/
pfm_Nanog <- readJASPARMatrix("../data/TFBS/UN0383.1.jaspar", matrixClass = "PFM")

hit_mat_75 <- foreach( i = seq_along(pfm), .combine = cbind ) %dopar%
  {
    lapply(mm9_prom.seq, 
           function(s)
             matchPWM(as.matrix(pfm[[1]]), s, min.score="75%") %>%
             length)
  } %>% as.data.frame

hit_mat_90 <- foreach( i = c(pfm, pfm_Nanog), .combine = cbind ) %dopar%
  {
    lapply(mm9_prom.seq, 
           function(s)
           {
             matches <- searchSeq(toPWM(i), s, min.score = "90%")
             if (length(matches) > 0)
             {
               paste(start(matches), end(matches)) %>% unique %>% length
             } else {
               return(0)
             }
           })
    # matchPWM(as.matrix(pfm[[i]]), s, min.score="75%") %>%
    # length)
  } %>% as.data.frame

colnames(hit_mat_75) <- colnames(hit_mat_90) <- 
  lapply(pfm, function(x) x@name) %>% unlist() %>% unname()
hit_mat_75$gene_id <- hit_mat_90$gene_id <- use_gene.gr$gene_id
hit_mat_75$gene_symbol <- hit_mat_90$gene_symbol <- use_gene.gr$gene_name

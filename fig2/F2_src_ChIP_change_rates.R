
solve_rate <- function(y, t = c(0, 6, 12)) {
  # first order kinetics solved by linear regression 
  # with function log(y) ~ lambda * t + log(y_0)
  # where t in {0, 6, 12} or {0, 12, 24, 36}
  X <- cbind(1, t)
  solve(t(X) %*% X) %*% t(X) %*% y
}


# get change rates
get_change_rates <- function(.name) {
  # calculate log-scale binding and binding change kinetics by 12 hours interval
  if (sum(grepl(paste0(.name, ".*_RC"), colnames(all_ChIP_TSS_mat)))) {
    ChIP_TSS_mat_b1 <- all_ChIP_TSS_mat[use_gene_ids, grep(paste0(.name, ".*_RC"),
                                                           colnames(all_ChIP_TSS_mat))] %>%
      as.matrix() %>% log()
    
    TSS_rates_b1 <- foreach(i = seq_len(nrow(ChIP_TSS_mat_b1)), 
                            .combine = cbind) %dopar% 
      {
        solve_rate(ChIP_TSS_mat_b1[i, ], t = c(0, 6, 12))
      } %>% t() %>% `colnames<-`(., c("NT", "Rate"))
  } else {
    TSS_rates_b1 <- `colnames<-`(matrix(NA, nrow = length(use_gene_ids), ncol = 2),
                                 c("NT", "Rate"))
  }
  
  ChIP_TSS_mat_b2 <- all_ChIP_rerun_TSS_mat[, grep(.name, colnames(all_ChIP_rerun_TSS_mat))] %>%
    as.matrix() %>% log()
  
  TSS_rates_b2 <- foreach(i = seq_len(nrow(ChIP_TSS_mat_b2)), 
                          .combine = rbind) %dopar% 
    {
      c(solve_rate(ChIP_TSS_mat_b2[i, 1:4], t = c(0, 0, 12, 12)),
        solve_rate(ChIP_TSS_mat_b2[i, 3:6], t = c(0, 0, 12, 12)),
        solve_rate(ChIP_TSS_mat_b2[i, 5:8], t = c(0, 0, 12, 12)),
        solve_rate(ChIP_TSS_mat_b2[i, 7:10], t = c(0, 0, 12, 12)))
    } %>% `colnames<-`(., c("P0", "P12_Rate", "P12", "P12C12_Rate", 
                            "P12C12", "P12C24_Rate", "P12C24", "P12C36_Rate"))
  # return
  out <- `rownames<-`(cbind(TSS_rates_b1, TSS_rates_b2), rownames(ChIP_TSS_mat_b2))
  out[is.infinite(out)] <- NA
  data.frame(out)
}




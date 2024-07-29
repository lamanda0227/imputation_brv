library(dplyr)

prepare_input <- function(fname, case_control_ref){
    if(!is.na(fname)) {
        aggregate_geno <- data.table::fread(fname, header = TRUE)
        aggregate_geno <- left_join(case_control_ref, aggregate_geno) %>% select(-FID, -IID, -assigned_id) %>% as.matrix()
        class(aggregate_geno) <- "numeric"

        num_var <- which(aggregate_geno != 0) %>% length()
        if(num_var >= 2){
            return(aggregate_geno)
        } else {
            return("LESS THAN TWO VAR PER GENE")
        }
    } else {
        return("FILE MISSING") 
    }
}

BRV <- function(fname, case_control_ref, y, sample, gene, es, causal_prop, prev, out_folder){
    if(sample == "hrc_topmed_exome") fname <- gsub(".raw", "_aggregate.raw", fname)
    X <- prepare_input(fname, case_control_ref)
    
    ## count number of rare variants for each individual
    if(!is.matrix(X)){
        cat("X is not matrix", file = sprintf("%s/%s_prev%.1f_es%.1f_causal%.1f_%s_maf001.csv", out_folder, gene, prev, es, causal_prop, sample))
        return(list(result_df = "X is not matrix", sample_name = sample))
    }
    
    if(is.na(X)[1]){
        result_df <- "X is NA"
        cat("X is NA", file = sprintf("%s/%s_prev%.1f_es%.1f_causal%.1f_%s_maf001.csv", out_folder, gene, prev, es, causal_prop,sample))
    } else {
        X.fit <- glm(y ~ X, family = "binomial")
                      
        if(dim(coef(summary(X.fit)))[1] != 1){
          zstat <- coef(summary(X.fit))[2, 3]
          pval <- coef(summary(X.fit))[2, 4]
        }else{
          zstat <- NA
          pval <- 1
        }
        
        result_df <-  data.frame(gene = gene, sample = sample, es = es, causal_prop = causal_prop, prev = prev, zstat = zstat, pval = pval, maf = 0.01)
        data.table::fwrite(result_df, file = sprintf("%s/%s_prev%.1f_es%.1f_causal%.1f_%s_maf001.csv", out_folder, gene, prev, es, causal_prop,sample))
    }
    
  return(list(result_df = result_df, sample_name = sample))
}

library(dscrutils)
library(genio)
library(dplyr)

### FUNCTIONS ###
return_exome_mat <- function(geno, case_control_ref){
    geno <- geno %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("IID")
    geno <- geno %>% mutate(idx = c(1:nrow(geno)))

    sampled_indv <- left_join(case_control_ref, geno %>% select(IID, idx), by = c("real_id" = "IID")) %>% pull(idx)
    geno <- geno[sampled_indv,] %>% select(-IID, -idx) %>% as.matrix()
    class(geno) <- "numeric"

    num_var <- apply(geno, 1, function(row) sum(row >= 1))
    num_var <- which(num_var >= 2) %>% length()
    if(num_var >= 2){
        return(geno)
    } else {
        return("ERROR")
    }
}

return_topmed_mat <- function(fname, case_control_ref){
    if(file.exists(fname)){
        topmed_geno <- data.table::fread(fname)

        topmed_geno <- topmed_geno %>%
            select(-c(1,3:6)) %>%
            t() %>%
            data.frame()
        topmed_geno <- topmed_geno %>%
            tibble::rownames_to_column("FID_IID") %>%
            tidyr::separate(FID_IID, c("Zero", "FID", "IID")) %>%
            select(-Zero) %>%
            mutate(idx = c(1:nrow(topmed_geno)))

        sampled_indv <- left_join(case_control_ref, topmed_geno %>% select(IID, idx),
                                  by = c("real_id"="IID")) %>% pull(idx)

        topmed_X <- topmed_geno[sampled_indv,] %>% select(-FID, -IID, -idx) %>% as.matrix()
        class(topmed_X) <- "numeric"
        topmed_X <- 2 - topmed_X

        num_var <- apply(topmed_X, 1, function(row) sum(row >= 1))
        num_var <- which(num_var >= 2) %>% length()
        if(num_var >= 2){
            return(topmed_X)
        } else {
            return("ERROR")
        }
     } else {
        return("FILE MISSING")
    }
}
                         
subset_X <- function(X, maf){
    pool_maf <- apply(X, 2, function(x) mean(x, na.rm = TRUE)/2)
    maf_select <- which(pool_maf < maf)
    X_select <- X[,maf_select]
    
    X_new  = rowSums(X_select, na.rm=TRUE) %>% as.matrix()
    X.fit <- glm(y ~ X_new, family = "binomial")
                      
    if(dim(coef(summary(X.fit)))[1] != 1){
      zstat <- coef(summary(X.fit))[2, 3]
      pval <- coef(summary(X.fit))[2, 4]
    }else{
      zstat <- NA
      pval <- 1
    }
                     
    return(data.frame(zstat = zstat, pval = pval, maf = maf)) 
}

BRV <- function(X, y, sample, es, causal_prop, prev){
    ## cout number of rare variants for each individual
    if(is.na(y)[1]){
      cat("y is NA", file = sprintf("%s/%s_%s.csv", out_folder, gene, sample))
      return(list(result_df = "y is NA", sample_name = sample))
    }
    if(!is.matrix(X)){
        cat("X is not matrix", file = sprintf("%s/%s_%s.csv", out_folder, gene, sample))
        return(list(result_df = "X is not matrix", sample_name = sample))
    }
    
    if(is.na(X)[1]){
        result_df <- "X is NA"
        cat("X is NA", file = sprintf("%s/%s_%s.csv", out_folder, gene, sample))
    } else {
        brv_result_001 <- subset_X(X, 0.01)
        brv_result_0005 <- subset_X(X, 0.005)
        brv_result_0001 <- subset_X(X, 0.001)
        
        result_df <- rbind(brv_result_001, brv_result_0005, brv_result_0001) %>% 
            mutate(sample = sample, es = es, causal_prop = causal_prop, prev = prev)
        # data.table::fwrite(result_df, sprintf("%s/%s_%s.csv", out_folder, gene, sample))
    }
    
  return(result_df)
}

### MAIN ###
# retrieve start and end values from command line
start <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
end <- as.numeric(commandArgs(trailingOnly = TRUE)[2])

# If end is greater than 3194, assign it to 3194
if (end > 3194) {
  end <- 3194
}

# Output the start and end values
print(paste("Start:", start))
print(paste("End:", end))

dsc_dir = '/home/tl3031/project/imputation-rvtest/workflows/imputation_aggregated_analysis/dsc_rsq03_maf001/bin_pheno_168206id_brv'
maf_dir <- "~/project/imputation-rvtest/workflows/imputation_aggregated_analysis/rsq_missing_variant_comp/maf001"
setwd('/home/tl3031/project/imputation-rvtest/workflows/imputation_aggregated_analysis/dsc_rsq03_maf001/bin_pheno_168206id_brv')

result <- dscrutils::dscquery(dsc_dir, targets = c("parse_input",
                                                   "parse_input.idx",
                                                   "fixed_effect_prop.deletrious_effect",
                                                   "fixed_effect_prop.causal_del",
                                                   "bin_phenotype",
                                                   "bin_phenotype.prevalence"),
                             module.output.files = c("parse_input", "bin_phenotype"))

result_filtered <- result %>% filter(fixed_effect_prop.deletrious_effect == 1.5, fixed_effect_prop.causal_del == "c(1,1)", bin_phenotype.prevalence == 0.1)

result_df_fname <- sprintf("%s/brv_result_df_%i_%i.csv", maf_dir, start, end)
if (file.exists(result_df_fname)) file.remove(result_df_fname)
fname <- file(result_df_fname, "w")
writeLines("dataset,zstat,pval", fname)
close(fname)

for(i in c(start:end)){
    print(i)
    gene <- readRDS(sprintf("./%s.rds", result_filtered[i,]$parse_input.output.file))$gene %>% strsplit(.,"_") %>% unlist()
    gene <- gene[[2]]
    
    pheno_rds <- readRDS(sprintf("./%s.rds", result_filtered[i,]$bin_phenotype.output.file))
    y <- pheno_rds$y
    case_control_ref <- pheno_rds$case_control_ref
    
    # read exome mat
    bim_fname <- sprintf("%s/%s_exome.bim", maf_dir, gene)
    fam_fname <- sprintf("%s/%s_exome.fam", maf_dir, gene)
    bed_fname <- sprintf("%s/%s_exome.bed", maf_dir, gene)
    if(file.exists(bim_fname)){ 
        bim <- read_bim(bim_fname)
        fam <- read_fam(fam_fname)
        exome <- read_bed(bed_fname, bim$id, fam$id)
        exome_X <- return_exome_mat(exome, case_control_ref)
        exome_result <- BRV(exome_X, y, "exome", 1.5, 1, 0.1)
    } else {
        exome_result <- list(zstat = NA, pval = NA)
    }
    
    # read topmed mat
    topmed_X <- return_topmed_mat(sprintf("%s/%s.traw", maf_dir, gene), case_control_ref)    
    topmed_result <- BRV(topmed_X, y, "topmed", 1.5, 1, 0.1)
    
    result_subdf <- rbind(topmed_result, exome_result)
    write.table(result_subdf, file = result_df_fname, append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
}
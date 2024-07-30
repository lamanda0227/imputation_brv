library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

gene_start <- as.integer(args[1])
gene_end <- as.integer(args[2])
bin_phenotype_idx <- as.integer(args[3])
rsq <- as.integer(args[4])
maf <- as.numeric(args[5])
out_folder <- args[6]

maf_c <- gsub("\\.", "", as.character(maf))
maf_c <- paste0("maf", maf_c)

format_imputed <- function(fname){
    df <- fread(fname, header = TRUE)
    df_names <- colnames(df)[7: dim(df)[2]] %>% stringr::str_split(pattern = "_", simplify = TRUE)
    colnames(df)[7:dim(df)[2]] <- df_names[,1]

    fid_iid <- df$IID %>% stringr::str_split(pattern = "_", simplify = TRUE)
    df$FID <- as.integer(fid_iid[,1])
    df$IID <- as.integer(fid_iid[,2])

    df_X <- df %>% arrange(FID, IID) %>% select(-c(3:6)) %>% as.matrix()
    df_X[, 3:ncol(df_X)] <- 2 - df_X[, 3:ncol(df_X)]
    
    return(as.data.frame(df_X))
}

format_exome <- function(fname){
    df <- fread(fname, header = TRUE)
    df_names <- colnames(df)[7: dim(df)[2]] %>% stringr::str_split(pattern = "_", simplify = TRUE)
    colnames(df)[7:dim(df)[2]] <- df_names[,1]

    df_X <- df %>% arrange(FID, IID) %>% select(-c(3:6))
    return(df_X)
}

prepare_input <- function(fname, case_control_ref, sample){
    if(!is.na(fname)) {
        if(sample == "exome"){
            geno <- format_exome(fname)
        } else if (sample %in% c("hrc", "topmed")) {
            geno <- format_imputed(fname)
        } else {
            geno <- fread(fname, header = TRUE)
        }
        
        geno <- left_join(case_control_ref, geno) %>% select(-FID, -IID, -assigned_id) %>% as.matrix()
        geno[is.na(geno)] <- 0
        aggregate_geno <- rowSums(geno)
        
        non_zero_var <- rowSums(geno >= 0.5)
        num_ids <- sum(non_zero_var >= 2)

        if(num_ids >= 2){
            return(as.matrix(aggregate_geno))
        } else {
            return("LESS THAN 2 VAR PER ID")
        }
    } else {
        return("FILE MISSING")
    }
}

BRV <- function(fname, case_control_ref, y, sample, gene, es, causal_prop, prev, out_folder, maf){
    maf_c <- gsub("\\.", "", as.character(maf))
    maf_c <- paste0("maf", maf_c)
    
    if(sample %in% c("exome", "hrc", "topmed")){
        fname <- gsub("_aggregate.txt", ".raw", fname)
    } else {
        fname <- gsub("_aggregate.raw", ".raw", fname)
    }
    
    print(fname)
    X <- prepare_input(fname, case_control_ref, sample)

    if(!is.matrix(X)){
        cat("X is not matrix", file = sprintf("%s/%s_prev%.1f_es%.1f_causal%.1f_%s_%s.csv", 
                                              out_folder, gene, prev, es, causal_prop, sample, maf_c))
        return(list(result_df = "X is not matrix", sample_name = sample))
    }

    if(is.na(X)[1]){
        result_df <- "X is NA"
        cat("X is NA", file = sprintf("%s/%s_prev%.1f_es%.1f_causal%.1f_%s_%s.csv", 
                                      out_folder, gene, prev, es, causal_prop, sample, maf_c))
        return(list(result_df = "X is NA", sample_name = sample))
    } else {
        X.fit <- glm(y ~ X, family = "binomial")

        if(dim(coef(summary(X.fit)))[1] != 1){
          zstat <- coef(summary(X.fit))[2, 3]
          pval <- coef(summary(X.fit))[2, 4]
        }else{
          zstat <- NA
          pval <- 1
        }

        result_df <-  data.frame(gene = gene, sample = sample, es = es, 
                                 causal_prop = causal_prop, prev = prev, 
                                 zstat = zstat, pval = pval, maf = maf)
        
        fwrite(result_df, file = sprintf("%s/%s_prev%.1f_es%.1f_causal%.1f_%s_%s.csv", 
                                         out_folder, gene, prev, es, causal_prop, sample, maf_c))
    }

  return(list(result_df = result_df, sample_name = sample))
}

for(gene_idx in c(gene_start:gene_end)){
    print(gene_idx)
    result_dir <- "~/project/git/imputation_brv/workflow/dsc_pipeline_rsq03/bin_pheno_168206id_brv_rsq03"
    
    parse_input_fname <- sprintf("%s/parse_input/parse_input_%i.rds", result_dir, gene_idx)
    parse_input <- readRDS(parse_input_fname)

    exome_fname <- gsub("maf001", maf_c, parse_input$exome_fname)
    if(rsq == 8) exome_fname <- gsub("rsq03", "rsq08", exome_fname)
    if(!file.exists(exome_fname)) exome_fname <- NA

    hrc_fname <- gsub("maf001", maf_c, parse_input$hrc_fname)
    if(rsq == 8) hrc_fname <- gsub("rsq03", "rsq08", hrc_fname)
    if(!file.exists(hrc_fname)) hrc_fname <- NA

    topmed_fname <- gsub("maf001", maf_c, parse_input$topmed_fname)
    if(rsq == 8) topmed_fname <- gsub("rsq03", "rsq08", topmed_fname)
    if(!file.exists(topmed_fname)) topmed_fname <- NA

    hrc_topmed_fname <- gsub("maf001", maf_c, parse_input$hrc_topmed_fname)
    if(rsq == 8) hrc_topmed_fname <- gsub("rsq03", "rsq08", hrc_topmed_fname)
    if(!file.exists(hrc_topmed_fname)) hrc_topmed_fname <- NA

    hrc_topmed_exome_fname <- gsub("maf001", maf_c, parse_input$hrc_topmed_exome_fname)
    if(rsq == 8) hrc_topmed_exome_fname <- gsub("rsq03", "rsq08", hrc_topmed_exome_fname)
    if(!file.exists(hrc_topmed_exome_fname)) hrc_topmed_exome_fname <- NA

    bin_phenotype_fname <- sprintf("%s/bin_phenotype/parse_input_%i_bin_phenotype_%i.rds", result_dir, gene_idx, bin_phenotype_idx)
    bin_phenotype <- readRDS(bin_phenotype_fname)

    y <- bin_phenotype$y
    prev <- bin_phenotype$prev
    es <- bin_phenotype$es
    causal_prop <- bin_phenotype$causal_prop
    case_control_ref <- bin_phenotype$case_control_ref

    BRV(exome_fname, case_control_ref, y, "exome", parse_input$gene, es, causal_prop, prev, out_folder, maf)
    BRV(hrc_fname, case_control_ref, y, "hrc", parse_input$gene, es, causal_prop, prev, out_folder, maf)
    BRV(topmed_fname, case_control_ref, y, "topmed", parse_input$gene, es, causal_prop, prev, out_folder, maf)
    BRV(hrc_topmed_fname, case_control_ref, y, "hrc_topmed", parse_input$gene, es, causal_prop, prev, out_folder, maf)
    BRV(hrc_topmed_exome_fname, case_control_ref, y, "hrc_topmed_exome", parse_input$gene, es, causal_prop, prev, out_folder, maf)
}

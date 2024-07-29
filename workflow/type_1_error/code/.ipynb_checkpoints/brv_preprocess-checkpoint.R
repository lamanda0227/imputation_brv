library(dplyr)
library(data.table)
library(genio)

return_exome_mat <- function(fname, case_control_ref){
    if(file.exists(fname)){
        bim_fname <- paste0(fname, '.bim')
        fam_fname <- paste0(fname, '.fam')
        bed_fname <- paste0(fname, '.bed')

        bim <- read_bim(bim_fname)
        fam <- read_fam(fam_fname)
        exome <- read_bed(bed_fname, bim$id, fam$id)

        geno <- exome %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("IID")
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
}

return_exome_mod_mat <- function(fname, case_control_ref){
    if (length(file.exists(fname)) != 1){
        geno <- fname %>% as.matrix()
        class(geno) <- "numeric"

        num_var <- apply(geno, 1, function(row) sum(row >= 1))
        num_var <- which(num_var >= 2) %>% length()
        if(num_var >= 2){
            return(geno)
        } else {
            return("ERROR")
        }
    }
    else if(length(file.exists(fname)) == 1 && file.exists(fname) == TRUE) {
        geno <- data.table::fread(fname, header = TRUE)
        geno <- geno %>% 
                    t() %>% 
                    as.data.frame() %>% 
                    tibble::rownames_to_column("IID")
        geno <- geno %>% mutate(idx = c(1:nrow(geno)))
        sampled_indv <- left_join(case_control_ref, geno %>% select(IID, idx), by = c("real_id"="IID")) %>% pull(idx)

        geno <- geno[sampled_indv,] %>% select(-IID, -idx) %>% as.matrix()
        class(geno) <- "numeric"

        num_var <- apply(geno, 1, function(row) sum(row >= 1))
        num_var <- which(num_var >= 2) %>% length()
        if(num_var >= 2){
            return(geno)
        } else {
            return("ERROR")
        }
    } else {
        return("FILE MISSING") 
    }
}

return_hrc_mat <- function(fname, case_control_ref){
    if(file.exists(fname)){
        hrc_geno <- data.table::fread(fname)

        hrc_geno <- hrc_geno %>% 
            select(-c(1,3:6)) %>% 
            t() %>% 
            data.frame()        
        hrc_geno <- hrc_geno %>% 
            tibble::rownames_to_column("FID_IID") %>% 
            tidyr::separate(FID_IID, c("FID", "IID")) %>%
            mutate(idx = c(1:nrow(hrc_geno)))

        sampled_indv <- left_join(case_control_ref, hrc_geno %>% select(IID, idx), 
                                  by = c("real_id"="IID")) %>% pull(idx)

        hrc_X <- hrc_geno[sampled_indv,] %>% select(-FID, -IID, -idx) %>% as.matrix()
        class(hrc_X) <- "numeric"
        hrc_X <- 2 - hrc_X

        num_var <- apply(hrc_X, 1, function(row) sum(row >= 1))
        num_var <- which(num_var >= 2) %>% length()
        if(num_var >= 2){
            return(hrc_X)
        } else {
            return("ERROR")
        }
    } else {
        return("FILE MISSING")
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
                     
brv_preprocess <- function(fname, case_control_ref, sample, maf){
    if(maf == 0.005){
        fname <- gsub("maf001", "maf0005", fname)
        if(sample == "hrc_topmed_exome"){
            X <- return_exome_mod_mat(fname, case_control_ref)
        } else if (sample == "exome") {
            X <- return_exome_mat(fname, case_control_ref)
        } else if(!is.na(fname)){
            if(grepl("hrc_topmed", sample)) {
                X <- return_exome_mod_mat(fname, case_control_ref)
            } else if (grepl("hrc", sample)){
                X <- return_hrc_mat(fname, case_control_ref)
            } else if (grepl("topmed", sample)) {
                X <- return_topmed_mat(fname, case_control_ref)
            }
        } else {
            X <- NA
        }
    } else if (maf == 0.001) {
        fname <- gsub("maf001", "maf0001", fname)
        if(sample == "hrc_topmed_exome"){
            X <- return_exome_mod_mat(fname, case_control_ref)
        } else if (sample == "exome") {
            X <- return_exome_mat(fname, case_control_ref)
        } else if(!is.na(fname)){
            if(grepl("hrc_topmed", sample)) {
                X <- return_exome_mod_mat(fname, case_control_ref)
            } else if (grepl("hrc", sample)){
                X <- return_hrc_mat(fname, case_control_ref)
            } else if (grepl("topmed", sample)) {
                X <- return_topmed_mat(fname, case_control_ref)
            }
        } else {
            X <- NA
        }
        
    } else {
        if(sample == "hrc_topmed_exome"){
            X <- fname
        } else if (sample == "exome") {
            X <- return_exome_mat(fname, case_control_ref)
        } else if(!is.na(fname)){
            if(grepl("hrc_topmed", sample)) {
                X <- return_exome_mod_mat(fname, case_control_ref)
            } else if (grepl("hrc", sample)){
                X <- return_hrc_mat(fname, case_control_ref)
            } else if (grepl("topmed", sample)) {
                X <- return_topmed_mat(fname, case_control_ref)
            }
        } else {
            X <- NA
        }
        
    }
        
    return(list(X=X, sample=sample, maf=maf))
} 
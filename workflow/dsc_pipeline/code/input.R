library(dplyr)
library(tidyverse)
library(genio)
library(data.table)

format_cadd_dir <- function(dir){
    pre <- strsplit(dir, "_(?!.*_)", perl=TRUE)[[1]][1]
    sub <- strsplit(dir, "_(?!.*_)", perl=TRUE)[[1]][2]

    cadd_dir <- sprintf("%s_cadd_%s", pre, sub)
    return(cadd_dir)    
}

return_cadd_name <- function(dir, fname){
    new_dir <- format_cadd_dir(dir)
    new_fname <- paste0(new_dir, "/", fname, ".traw")
    if(!file.exists(new_fname)){
        new_fname <- NA
    }
    
    return(new_fname)
}

parse_input <- function(exome_folder, 
                        hrc_topmed_exome_folder, hrc_topmed_folder, 
                        hrc_folder, topmed_folder, 
                        idx){
    
    flst <- list.files(hrc_topmed_exome_folder, pattern=".traw") %>% tools::file_path_sans_ext()
    fname <- flst[[idx]]
    
    hrc_topmed_exome_fname <- paste0(hrc_topmed_exome_folder, "/", fname, ".traw")
    
    exome_file <- list.files(exome_folder, pattern = paste0(fname,".log"))
    if(length(exome_file) != 0){
        bim_fname <- paste0(exome_folder, "/", fname, '.bim')
        fam_fname <- paste0(exome_folder, "/", fname, '.fam')
        bed_fname <- paste0(exome_folder, "/", fname, '.bed')
    
        bim <- read_bim(bim_fname)
        fam <- read_fam(fam_fname)
        exome <- read_bed(bed_fname, bim$id, fam$id)
    } else {
        exome <- NA
        bim_fname <- NA
        fam_fname <- NA
    }
    
    hrc_fname <- paste0(hrc_folder, "/", fname, ".traw")
    if(!file.exists(hrc_fname)){
        hrc_fname <- NA
    }
    
    topmed_fname <- paste0(topmed_folder, "/", fname, ".traw")
    if(!file.exists(topmed_fname)){
        hrc_fname <- NA
    }
    
    hrc_topmed_fname <- paste0(hrc_topmed_folder, "/", fname, ".traw")
    if(!file.exists(hrc_topmed_fname)){
        hrc_topmed_fname <- NA
    }
    
    cadd_fnames <- parse_input_cadd(exome_folder, 
                                          hrc_topmed_exome_folder, hrc_topmed_folder, 
                                          hrc_folder, topmed_folder, 
                                          fname)
    
    return(list(gene = fname,
                bim_fname = bim_fname,
                exome = exome,
                exome_cadd = cadd_fnames$exome_cadd,
                hrc_topmed_exome_fname = hrc_topmed_exome_fname,
                hrc_topmed_exome_cadd_fname = cadd_fnames$hrc_topmed_exome_cadd_fname, 
                hrc_fname = hrc_fname,
                hrc_cadd_fname = cadd_fnames$hrc_cadd_fname,
                topmed_fname = topmed_fname,
                topmed_cadd_fname = cadd_fnames$topmed_cadd_fname,
                hrc_topmed_fname = hrc_topmed_fname,
                hrc_topmed_cadd_fname = cadd_fnames$hrc_topmed_cadd_fname))
}

parse_input_cadd <- function(exome_folder, 
                             hrc_topmed_exome_folder, hrc_topmed_folder, 
                             hrc_folder, topmed_folder, 
                             fname){
    exome_cadd_dir <- paste0(exome_folder, "_cadd")
    exome_file <- list.files(exome_cadd_dir, pattern = paste0(fname,".log"))
    if(length(exome_file) != 0){
        bim_fname <- paste0(exome_cadd_dir, "/", fname, '.bim')
        fam_fname <- paste0(exome_cadd_dir, "/", fname, '.fam')
        bed_fname <- paste0(exome_cadd_dir, "/", fname, '.bed')
    
        bim <- read_bim(bim_fname)
        fam <- read_fam(fam_fname)
        exome <- read_bed(bed_fname, bim$id, fam$id)
    } else {
        exome <- NA
        bim_fname <- NA
        fam_fname <- NA
    }
    
    hrc_topmed_exome_cadd_fname <- return_cadd_name(hrc_topmed_exome_folder, fname)
    hrc_topmed_cadd_fname <- return_cadd_name(hrc_topmed_folder, fname)
    hrc_cadd_fname <- return_cadd_name(hrc_folder, fname)
    topmed_cadd_fname <- return_cadd_name(topmed_folder, fname)
    
    return(list(exome_cadd = exome,
                hrc_topmed_exome_cadd_fname = hrc_topmed_exome_cadd_fname,
                hrc_topmed_cadd_fname = hrc_topmed_cadd_fname,
                hrc_cadd_fname = hrc_cadd_fname,
                topmed_cadd_fname = topmed_cadd_fname))
}

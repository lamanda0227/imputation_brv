library(dplyr)
library(tidyverse)
library(genio)
library(data.table)

parse_input <- function(exome_folder, 
                        hrc_topmed_exome_folder, hrc_topmed_folder, 
                        hrc_folder, topmed_folder, 
                        idx){
    
    flst <- list.files(hrc_topmed_exome_folder, pattern="_aggregate.raw") %>% tools::file_path_sans_ext()
    fname <- gsub("_aggregate", "", flst[[idx]])
    
    hrc_topmed_exome_fname <- paste0(hrc_topmed_exome_folder, "/", fname, ".raw")
    
    exome_fname <- paste0(exome_folder, "/", fname, "_aggregate.txt")
    if(!file.exists(exome_fname)){
        exome_fname <- NA
    }
    
    hrc_fname <- paste0(hrc_folder, "/", fname, "_aggregate.txt")
    if(!file.exists(hrc_fname)){
        hrc_fname <- NA
    }
    
    topmed_fname <- paste0(topmed_folder, "/", fname, "_aggregate.txt")
    if(!file.exists(topmed_fname)){
        hrc_fname <- NA
    }
    
    hrc_topmed_fname <- paste0(hrc_topmed_folder, "/", fname, "_aggregate.raw")
    if(!file.exists(hrc_topmed_fname)){
        hrc_topmed_fname <- NA
    }
    
    return(list(gene = fname,
                exome_fname = exome_fname,
                hrc_fname = hrc_fname,
                topmed_fname = topmed_fname,
                hrc_topmed_fname = hrc_topmed_fname,
                hrc_topmed_exome_fname = hrc_topmed_exome_fname))
}

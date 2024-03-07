library(dplyr)
library(tidyverse)
library(genio)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
from=as.numeric(args[1])
to=as.numeric(args[2])
rsq=as.numeric(args[4])

if(args[3]=="05"){
    maf=0.005
} else if (args[3] == "01"){
    maf = 0.001
} else {
    maf = 0.01
}

maf_c <- gsub("\\.", "", as.character(maf))

exome_folder = sprintf("~/project/imputation-rvtest/workflows/imputation_aggregated_analysis/maf%s/exome", maf_c)
hrc_folder = sprintf("~/project/imputation-rvtest/workflows/imputation_aggregated_analysis/maf%s/hrc_rsq0%d", maf_c, rsq)
topmed_folder = sprintf("~/project/imputation-rvtest/workflows/imputation_aggregated_analysis/maf%s/topmed_rsq0%d", maf_c, rsq)

annot_fname = sprintf("~/project/imputation-rvtest/analysis/imputation_aggregated_analysis/hrc_topmed/hrc_topmed_168206ids_rsq0%d_maf%s_annot.csv.gz", rsq, maf_c)

annot <- fread(annot_fname)
genes <- annot$Gene.refGene %>% unique()

counter = from
for(g in genes[c(from:to)]){
    print(counter)
    counter = counter + 1
    
    chr <- annot %>% filter(Gene.refGene == g) %>% pull(Chr) %>% unique()
    fname <- paste0("chr", chr, "_", g)
    print(fname)
    
    hrc_snplist <- annot %>% filter(Gene.refGene == g) %>% filter(source == "hrc") %>% pull(ID_hg19)
    topmed_snplist <- annot %>% filter(Gene.refGene == g) %>% filter(source == "topmed") %>% pull(ID)

    hrc_fname <- paste0(hrc_folder, "/", fname, ".traw")
    if(file.exists(hrc_fname)){
        hrc <- fread(hrc_fname)
        hrc_names <- colnames(hrc)[7: dim(hrc)[2]] %>% stringr::str_split(pattern = "_", simplify = TRUE)
        colnames(hrc)[7:dim(hrc)[2]] <- hrc_names[,2]
        
        hrc_gene <- hrc %>% filter(SNP %in% hrc_snplist)
        hrc_gene <- hrc_gene %>% remove_rownames() %>% select(-c("CHR", "(C)M", "POS", "COUNTED", "ALT", "SNP"))
        hrc_gene <- 2-hrc_gene %>% select(sort(colnames(hrc_gene)))  
        bed_full <- hrc_gene
    }

    topmed_fname <- paste0(topmed_folder, "/", fname, ".traw")
    if(file.exists(topmed_fname)){
        topmed <- fread(topmed_fname)
        topmed_names <- colnames(topmed)[7: dim(topmed)[2]] %>% stringr::str_split(pattern = "_", simplify = TRUE)
        colnames(topmed)[7:dim(topmed)[2]] <- topmed_names[,3]
        
        topmed_gene <- topmed %>% filter(SNP %in% topmed_snplist)
        topmed_gene <- topmed_gene %>% remove_rownames() %>% select(-c("CHR", "(C)M", "POS", "COUNTED", "ALT", "SNP"))
        topmed_gene <- 2 - topmed_gene %>% select(sort(colnames(topmed_gene)))
        
        if(exists("bed_full")){
            bed_full <- rbind(bed_full, topmed_gene)
        } else {
            bed_full <- topmed_gene
        }
    }
    fwrite(bed_full, 
           sprintf("/home/tl3031/project/imputation-rvtest/workflows/imputation_aggregated_analysis/maf%s/hrc_topmed_rsq0%d/%s.traw", 
                   maf_c, rsq, fname))
}

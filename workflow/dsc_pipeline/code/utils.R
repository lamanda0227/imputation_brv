library(genio)
library(dplyr)

mk_dir <- function(directory){
  if (file.exists(directory) == FALSE){
    mk_dir_cmd <- paste('mkdir', directory)
    run_cmd(mk_dir_cmd)
  }
}

run_cmd <- function(cmd_str, shell_exec="/bin/bash", fout='', ferr='', quit_on_error=TRUE, ...) {
  if (ferr!=FALSE) {
    write("Running shell command:", stdout())
    write(cmd_str, stdout())
  }
  out <- system2(shell_exec, args = c("-c", shQuote(cmd_str)), stdout = fout, stderr=ferr, ...)
  if (out != 0 && quit_on_error && fout != TRUE && ferr != TRUE){
    stop(paste(strsplit(cmd_str, " +")[[1]][1], "command failed (returned a non-zero exit status)"))
  }
  return(out)
}

# write out R files as plink binary files
write_out_plink <- function(r_data, map, out_dir, out_name, num_chr = 1) {
  # mk_dir(out_dir)
  
  ## Write out fam file
  fam_data <- r_data %>% select(c(1:5), ncol(r_data))
  colnames(fam_data) <- c('fam', 'id', 'pat', 'mat', 'sex', 'pheno')
  
  if(fam_data$pheno[1] %in% c(0, 1)){
      fam_data$pheno <- fam_data$pheno + 1
  }  
  write.table(fam_data, paste0(out_dir, out_name, ".fam"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  ## Write out map file
  map <- map %>% rename('chr' = 'Chromosome',
                        'pos' = 'Position.bp',
                        'id' = 'ID',
                        'ref' = 'REF',
                        'alt' = 'ALT',
                        'posg' = 'morgan') %>%
                  select(chr, id, posg, pos, alt, ref)
  write.table(map, paste0(out_dir, out_name, ".bim"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  ## Write out bed file
  genotype <- r_data %>% select(-c(1:5, ncol(r_data))) %>%
    as.matrix %>%
    t
  genotype[is.na(genotype)] <- 0
  row.names(genotype) <- map$id
  colnames(genotype) <- fam_data$id
  write_bed(paste0(out_dir, out_name), genotype)
}

create_grm_pop <- function(sample){
  num_pop <- sample %>% pull(IID) %>% length()
  pop_grm <- matrix(0, num_pop, num_pop)
  diag(pop_grm) <- 1
  
  return(pop_grm)
}

create_grm <- function(sample){
  num_ind <- sample %>% pull(IID) %>% length()
  grm <- matrix(0, num_ind, num_ind)
  
  i = 1
  while (i < num_ind){
    fam_id <- sample$FID[i]
    fam_size <- sample$FID[sample$FID == fam_id] %>% length()
    
    grm[i:(i+fam_size-1), i:(i+fam_size-1)] = 0.25  # sibling kinship
    grm[i, i+1] = grm[i+1, i] = 0  # parent kinship
    
    i = i + fam_size
  }
  diag(grm) <- 0.5  # self kinship
  grm <- grm * 2
  
  return(grm)
}

mean_impute <- function(geno){
  f <- apply(geno, 2, function(x) mean(x,na.rm = TRUE))
  for (i in 1:length(f)) geno[,i][which(is.na(geno[,i]))] <- f[i]
  return(geno)
}

parse_x_y <- function(sample, sample_type, map, out_dir){
  # writing plink files susie_glm
  plink_name <- tempfile(pattern = "plink_", out_dir)
  write_out_plink(sample, map, out_dir, plink_name)

  X_name <- tempfile(pattern = paste0("fit_", sample_type, "_X_"), out_dir)
  y_name <- tempfile(pattern = paste0("fit_", sample_type, "_y_"), out_dir)
  grm_name <- tempfile(pattern = paste0("fit_", sample_type, "_K_"), out_dir)
  
  X <- sample[,-c(1:5, ncol(sample))] %>% as.matrix()
  X <- mean_impute(X)
  y <- sample['trait'] %>% as.matrix()
  
  grm <- matrix()
  if (sample_type=='pop'){
    K <- create_grm_pop(sample) %>% as.matrix()
  } else {
    K <- create_grm(sample) %>% as.matrix()
  }
  
  X_name <- tempfile(pattern = paste0("fit_", sample_type, "_X_"), out_dir)
  X_name <- paste0(X_name, '.txt')
  y_name <- tempfile(pattern = paste0("fit_", sample_type, "_y_"), out_dir)
  y_name <- paste0(y_name, '.txt')
  K_name <- tempfile(pattern = paste0("fit_", sample_type, "_K_"), out_dir)
  K_name <- paste0(K_name, '.txt')
  
  write.table(X, X_name, quote = FALSE, sep = '\t')
  write.table(y, y_name, quote = FALSE, sep = '\t')
  write.table(K, K_name, quote = FALSE, sep = '\t')
  
  return(list(X = X,
              y = y,
              X_name = X_name,
              y_name = y_name,
              K_name = K_name,
              plink_name = paste0(plink_name, '.bim'),
              susie_glm_name = paste0('susie_glm_', tail(unlist(strsplit(plink_name, '_')), 1), '.txt')))
}

create_group_def_file <- function(num_var_per_group, map){
    num_var <- nrow(map)
    
    group_idx <- c(1: floor(num_var / num_var_per_group)) %>% rep(num_var_per_group)
    
    diff <- num_var - length(group_idx)
    group_idx <- c(group_idx, rep(ceiling(num_var / num_var_per_group), diff)) %>% sample()
    
    group_def_file <-
        map %>% 
        mutate(SetID = paste0("Set", group_idx), Weight = 1) %>% 
        select(SetID, Chromosome, Position.bp, REF, ALT, Weight) %>%
        arrange(Chromosome, Position.bp)
    return(group_def_file)
}
 
# generate_causal_idx <- function(causal_set_idx, causal_prop, gt_one_odds_prop){
#     mat <- expand.grid(c(1:1000), c(1:length(causal_set_idx)))
#     mat <- mat %>% mutate(causal = 0, idx = c(1:nrow(mat))) %>% rename(rep_idx = Var1, var_idx = Var2)

#     causal_reps <- sample(c(1:1000), 1000 * causal_prop)
#     causal_reps_idx <- mat %>% filter(rep_idx %in% causal_reps) %>% pull(idx)
    
#     all_gt_one_idx <- sample(causal_reps_idx, length(causal_reps_idx) * gt_one_odds_prop)
#     all_lt_one_idx <- setdiff(causal_reps_idx, all_gt_one_idx)
#     mat$causal[mat$idx %in% all_gt_one_idx] <- 1
#     mat$causal[mat$idx %in% all_lt_one_idx] <- -1
    
#     causal_reps_idx_selected <- sample(c(1:1000), 1)
#     causal_mat <- mat %>% filter(rep_idx == causal_reps_idx_selected)
    
#     gt_one_idx <- causal_mat %>% filter(causal == 1) %>% pull(idx)
#     lt_one_idx <- causal_mat %>% filter(causal == -1) %>% pull(idx)
    
#     return(list(gt_one_idx = gt_one_idx, 
#                 lt_one_idx = lt_one_idx))
# }
             
generate_causal_idx <- function(causal_set_idx, causal_prop, gt_one_odds_prop){
    causal_set_idx_df <- data.frame()
    
    if(length(causal_set_idx) < 100){
        padding <- rep(0, 100 - length(causal_set_idx))
        causal_set_idx_df <- data.frame(causal_idx = c(padding, causal_set_idx) %>% sample(), idx = c(1:100))
    } else {
        causal_set_idx_df <- data.frame(causal_idx = causal_set_idx %>% sample(), idx = c(1:length(causal_set_idx)))
    }
    
    gt_one <- causal_set_idx_df %>% pull(idx) %>% sample(100 * causal_prop * gt_one_odds_prop)
    lt_one <- causal_set_idx_df %>% filter(! idx%in% gt_one) %>% pull(idx) %>% sample(100 * causal_prop * (1 - gt_one_odds_prop))
        
    gt_one_idx <- causal_set_idx_df %>% filter(idx %in% gt_one) %>% filter(causal_idx != 0) %>% pull(causal_idx)
    lt_one_idx <- causal_set_idx_df %>% filter(idx %in% lt_one) %>% filter(causal_idx != 0) %>% pull(causal_idx)

    return(list(gt_one_idx = gt_one_idx,
               lt_one_idx = lt_one_idx))
}

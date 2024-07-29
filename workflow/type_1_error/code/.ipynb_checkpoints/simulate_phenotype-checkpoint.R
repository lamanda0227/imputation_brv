library(dplyr)

logistic_generator_matrix <- function(geno, prev, beta_eff){
    print(dim(geno))
    print(length(beta_eff))
    geno[is.na(geno)] <- 0
    beta_tot <-  as.matrix(geno) %*% as.vector(beta_eff)
    alpha <- log(prev/(1-prev))
    epsilon <- rnorm(length(beta_tot))
    trait_l_bin <- 1 / (1 + exp(-(alpha + beta_tot + epsilon)))
    samp_func <- function(X) (sample(0:1, size=1, prob=c(1-X,X)))
    traits_bin<-apply(trait_l_bin %>% bind_cols(), MARGIN=c(1,2), samp_func) %>% 
        `colnames<-`(paste0('D', 1:1)) 
    return(as.vector(traits_bin))
}

generate_fixed_effect_prop <- function(num_snp, deletrious_effect, protective_effect, causal_prop, deleterious_prop){
    # simulating effect
    effect <- rep(0, num_snp)

    ### select the index of positive and negative snps from each group
    probs <- runif(length(effect))
    names(probs) <- c(1:length(probs))
    gt_one <- which(probs < (causal_prop * deleterious_prop)) %>% names() %>% as.numeric()

    if(length(gt_one)==0){
      lt_one <- which(probs < (causal_prop * (1-deleterious_prop))) %>% names() %>% as.numeric()
    } else if (deleterious_prop == 0.5){
        lt_one <- sample(gt_one, (length(gt_one)/2))
        gt_one <- setdiff(gt_one, lt_one)
    }else{
      lt_one <- which(probs[-gt_one] < (causal_prop * (1-deleterious_prop))) %>% names() %>% as.numeric()
    }
    
    if(length(gt_one) == 0){
        effect[lt_one] <- -log(deletrious_effect)
    } else {
        effect[gt_one] <- log(deletrious_effect)
        effect[lt_one] <- log(protective_effect)
    }
    
    # return(list(effect = effect, es = deletrious_effect, causal_prop = causal_prop))
    return(effect)
}

simulate_bin_phenotype <- function(hrc_topmed_exome_fname, deletrious_effect, protective_effect, causal_prop, deleterious_prop, prevalence, case_ss, control_ss){   
    hrc_topmed_exome <- data.table::fread(hrc_topmed_exome_fname, header = TRUE)
    
    effect <- generate_fixed_effect_prop(ncol(hrc_topmed_exome) - 2, deletrious_effect, protective_effect, causal_prop, deleterious_prop)
    sample_with_trait <- data.frame(hrc_topmed_exome,
                                    trait = logistic_generator_matrix(hrc_topmed_exome[,-c(1:2)] %>% as.matrix(), prevalence, effect))
    
    case <- subset(sample_with_trait, trait == 1)
    control <- subset(sample_with_trait, trait == 0)
    
    num_case <- nrow(case)
    num_control <- nrow(control)
    
    if(num_control == 0){
        return(list(y = NA,
                    exome_X = NA,
                    exon_X = NA,
                    exome_add_var_X = NA,
                    hrc_topmed_X = NA,
                    hrc_X = NA,
                    topmed_X = NA,
                    case_control_ref = NA))
    }
    
    if (num_case < case_ss & num_control < control_ss){
        case_a <- case
        case_b <- data.frame(control %>% select(-trait),
                             trait = logistic_generator_matrix(control[,-c(1:2)] %>% select(-trait), prevalence, effect)) %>%
                    filter(trait == 1) %>%
                    sample_n(case_ss - num_case)
        case <- rbind(case_a, case_b)
        
        control_a <- control
        control_b <- data.frame(case %>% select(-trait),
                                trait = logistic_generator_matrix(case[,-c(1:2)] %>% select(-trait), prevalence, effect)) %>%
                        filter(trait == 0) %>%
                        sample_n(control_ss - num_control, replace = TRUE)
        control <- rbind(control_a, control_b)
                
    } else if (num_case < case_ss){
        case_a <- case
        case_b <- data.frame(control %>% select(-trait),
                             trait = logistic_generator_matrix(control[,-c(1:2)] %>% select(-trait), prevalence, effect)) %>%
                    filter(trait == 1) %>%
                    sample_n(case_ss - num_case)
        case <- rbind(case_a, case_b)
        control <- control %>% sample_n(control_ss)
    } else if (num_control < control_ss) {
        control_a <- control
        control_b <- data.frame(case %>% select(-trait),
                                trait = logistic_generator_matrix(case[,-c(1:2)] %>% select(-trait), prevalence, effect)) %>%
                        filter(trait == 0) %>%
                        sample_n(control_ss - num_control, replace = TRUE)
        control <- rbind(control_a, control_b)
        case <- case %>% sample_n(case_ss)
    } else {
        case <- case %>% sample_n(case_ss)
        control <- control %>% sample_n(control_ss)
    }

    # case$real_id <- rownames(case)
    case$assigned_id <- c(1:case_ss)
    # control$real_id <- rownames(control)
    control$assigned_id <- c(1:control_ss) + case_ss
    
    sample <- rbind(case, control)
    rownames(sample) <- sample$assigned_id
    y <- sample %>% pull(trait)
    case_control_ref <- sample %>% select(FID, IID, assigned_id)

    return(list(prev = prevalence,
                es = deletrious_effect, 
                causal_prop = causal_prop,
                y = y,
                case_control_ref = case_control_ref))
}


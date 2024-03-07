library(dplyr)

logistic_generator_matrix <- function(geno, prev, beta_eff){
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

simulate_bin_phenotype <- function(effect, prevalence, case_ss, control_ss, exome_add_var_fname){    
    exome_add_var <- data.table::fread(exome_add_var_fname, header = TRUE) %>% t()
    sample_with_trait <- data.frame(exome_add_var,
                                    trait = logistic_generator_matrix(exome_add_var, prevalence, effect))
    
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
                             trait = logistic_generator_matrix(control %>% select(-trait), prevalence, effect)) %>%
                    filter(trait == 1) %>%
                    sample_n(case_ss - num_case)
        case <- rbind(case_a, case_b)
        
        control_a <- control
        control_b <- data.frame(case %>% select(-trait),
                                trait = logistic_generator_matrix(case %>% select(-trait), prevalence, effect)) %>%
                        filter(trait == 0) %>%
                        sample_n(control_ss - num_control, replace = TRUE)
        control <- rbind(control_a, control_b)
                
    } else if (num_case < case_ss){
        case_a <- case
        case_b <- data.frame(control %>% select(-trait),
                             trait = logistic_generator_matrix(control %>% select(-trait), prevalence, effect)) %>%
                    filter(trait == 1) %>%
                    sample_n(case_ss - num_case)
        case <- rbind(case_a, case_b)
        control <- control %>% sample_n(control_ss)
    } else if (num_control < control_ss) {
        control_a <- control
        control_b <- data.frame(case %>% select(-trait),
                                trait = logistic_generator_matrix(case %>% select(-trait), prevalence, effect)) %>%
                        filter(trait == 0) %>%
                        sample_n(control_ss - num_control, replace = TRUE)
        control <- rbind(control_a, control_b)
        case <- case %>% sample_n(case_ss)
    } else {
        case <- case %>% sample_n(case_ss)
        control <- control %>% sample_n(control_ss)
    }

    case$real_id <- rownames(case)
    case$assigned_id <- c(1:case_ss)
    control$real_id <- rownames(control)
    control$assigned_id <- c(1:control_ss) + case_ss
    
    sample <- rbind(case, control)
    rownames(sample) <- sample$assigned_id
    X <- sample %>% select(-trait, -real_id, -assigned_id) # hrc_topmed_exome
    y <- sample %>% pull(trait)
    case_control_ref <- sample %>% select(real_id, assigned_id)

    return(list(prev = prevalence,
                y = y,
                X = as.matrix(X),
                case_control_ref = case_control_ref))
}


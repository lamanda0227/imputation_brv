library(dplyr)
library(genio)

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
    
    return(list(effect = effect, es = deletrious_effect, causal_prop = causal_prop))
}
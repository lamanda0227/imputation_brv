library(dplyr)

my_check <- function(y, X, sample){
  ## Internal function to check argument:s y, X

  # y as numeric vector
  if (!is.vector(y) || mode(y) != "numeric")
    stop("Sorry, argument 'y' must be a numeric vector")
  # no misssing values in y
  if (any(is.na(y))) 
    stop("No missing data allowed in argument 'y' ")
  # binary values (0, 1) in y
  if (!all(y %in% c(0, 1)))
    stop("Argument 'y' must contain only 0 and 1")
  # X as matrix or data frame
  if(!is.matrix(X) & !is.data.frame(X))
    stop("Argument 'X' must be a matrix or data.frame")    
  # compatibility between X and y
  if (nrow(X) != length(y)) 
    stop("'X' and 'y' have different lengths")
  # force X as matrix
  if (!is.matrix(X)) X = as.matrix(X)

  # results
  list(y=y, X=X)
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

BRV <- function(X, y, sample, gene, out_folder, es, causal_prop, prev){
    ## cout number of rare variants for each individual
    if(is.na(y)[1]){
      cat("y is NA", file = sprintf("%s/%s_%s.csv", out_folder, gene, sample))
      return(result_df <- "y is NA")
    }
    if(!is.matrix(X)){
        cat("X is not matrix", file = sprintf("%s/%s_%s.csv", out_folder, gene, sample))
        return(result_df <- "X is not matrix")
    }
    
    if(is.na(X)[1]){
        result_df <- "X is NA"
        cat("X is NA", sprintf("%s/%s_%s.csv", out_folder, gene, sample))
    } else {
        brv_result_001 <- subset_X(X, 0.01)
        brv_result_0005 <- subset_X(X, 0.005)
        brv_result_0001 <- subset_X(X, 0.001)
        
        result_df <- rbind(brv_result_001, brv_result_0005, brv_result_0001) %>% 
            mutate(gene = gene, sample = sample, es = es, causal_prop = causal_prop, prev = prev)
        data.table::fwrite(result_df, sprintf("%s/%s_%s.csv", out_folder, gene, sample))
    }
    
  return(result_df)
}

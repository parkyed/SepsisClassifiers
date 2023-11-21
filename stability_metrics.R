# source: https://github.com/nogueirs/JMLR2018/blob/master/R/getStability.R

getStability <- function(X,alpha=0.05) {
  ## the input X is a binary matrix of size M*d where:
  ## M is the number of bootstrap replicates
  ## d is the total number of features
  ## alpha is the level of significance (e.g. if alpha=0.05, we will get 95% confidence intervals)
  ## it's an optional argument and is set to 5% by default
  ### first we compute the stability
  
  M<-nrow(X)
  d<-ncol(X)
  hatPF<-colMeans(X)
  kbar<-sum(hatPF)
  v_rand=(kbar/d)*(1-kbar/d)
  stability<-1-(M/(M-1))*mean(hatPF*(1-hatPF))/v_rand ## this is the stability estimate
  
  ## then we compute the variance of the estimate
  ki<-rowSums(X)
  phi_i<-rep(0,M)
  for(i in 1:M){ 
    phi_i[i]<-(1/v_rand)*((1/d)*sum(X[i,]*hatPF)-(ki[i]*kbar)/d^2-(stability/2)*((2*kbar*ki[i])/d^2-ki[i]/d-kbar/d+1))
  }
  phi_bar=mean(phi_i)
  var_stab=(4/M^2)*sum((phi_i-phi_bar)^2) ## this is the variance of the stability estimate
  
  ## then we calculate lower and upper limits of the confidence intervals
  z<-qnorm(1-alpha/2) # this is the standard normal cumulative inverse at a level 1-alpha/2
  upper<-stability+z*sqrt(var_stab) ## the upper bound of the (1-alpha) confidence interval
  lower<-stability-z*sqrt(var_stab) ## the lower bound of the (1-alpha) confidence interval
  
  return(list("stability"=stability,"variance"=var_stab,"lower"=lower,"upper"=upper))
  
}

# functions to calculate the correlation matrix of a matrix - two alternative packages

scaleCorrelationCoop <- function(mat){
  scaled.mat <- scale(mat, center = TRUE, scale = FALSE)
  pcorr.mat <- abs(pcor(scaled.mat, inplace = TRUE))
  return(pcorr.mat)
}


scaleCorrelationRfast <- function(mat){
  scaled.mat <- scale(mat, center = TRUE, scale = FALSE)
  pcorr.mat <- abs(cora(scaled.mat, large = FALSE))
  return(pcorr.mat)
}

# function to create a binary matrix from a list of feature subsets and a full list of the features they are taken from

createBinaryMatrix <- function(sel_f_list, full_f_list){
  
  #' creates a binary matrix from the list of selected features and the full feature list
  #'
  #' @param sel_f_list a list of lists of IDs of the selected features. sub lists are list of characters
  #' @param full_f_list the full list of feature names (characters) they were originally selected from, in original order - i.e. training set column names
  
  # get indices of selected features
  
  (f_idx <- map(sel_f_list, function(x) {match(x, full_f_list)}))
  
  # get the number of genes in the list
  
  (n_genes <-  length(full_f_list))
  
  # setNames names the list elements with their index, to capture the original indices of the list elements
  # Filter(Negate(is.null),...) then removes any list elements that have NULL values, but allows preservation of original indices as the names
  # then melt this to a dataframe, where the names, i.e. original list indices become the second column of d2
  
  (f_list_melt <- melt(Filter(Negate(is.null), setNames(f_idx, seq_along(f_idx)))))
  
  # create a matrix of zeros to become the binary matrix - row for each list of selected features, cols for the full feature set
  
  (bin_mat <- matrix(0, nrow=length(f_idx), ncol=max(n_genes)))
  
  # cbind(as.numeric(d2...)) recreates a matrix, swapping the columns, so the row index is the first column, and column index is second
  # then m1 matrix can be sliced on this, i.e. on the coordinates where the value is replaced by 1
  
  bin_mat[cbind(as.numeric(f_list_melt[,2]), f_list_melt[,1])] <- 1
  
  # add back the variable names as column names
  colnames(bin_mat) <- full_f_list
  
  return(bin_mat)
}


## alternative, maybe shorter, from stabm package
# get the list of names of all the features: this was elsewhere in the stabm package as reused by different functions
#(F.all = colnames(sim.mat))

# create a zeros matrix, with number of columns equal to the total number of features over the datasets, rows = number of subsamples
#(Z = matrix(0, nrow = n, ncol = p))

# code for function to create the selected features index, from the list of features selected each run and the full feature list
#    for (i in seq_along(features)) {
#      Z[i, ] = as.numeric(F.all %in% features[[i]])


# Zucknick stability measure

# params for manual testing
#featues <- noise.results.features[[1]]
#C <- noise.corr
#threshold <- 0.5
#(F1 <- noise.results.features[[1]][[1]])
#(F2 <- noise.results.features[[1]][[2]])

getStabilityZucknick = function(features, C, threshold) {
  
  #' Calculates zucknick stability measure
  #' 
  #' @param features the list of the sets of selected features
  #' @param C similarity matrix of all features based on all datasets (i.e before filtering)
  #' @param threshold threshold for indicating which features are similar or not, two features are similar if entry in sim.mat is >=threshold
  
  # get the union of all features selected over all samples, and remoe NA values
  (F.all <- unique(unlist(features)))
  F.all <- F.all[!is.na(F.all)]
  
  # subset C corrlation matrix on the union of all features selected
  ## NOTE - should add check to match order of F.all to C sim.mat rows and columns, as may not be the same
  (sim.mat <- C[F.all,F.all])
  
  # set values in C below threshold to zero
  sim.mat[sim.mat < threshold] <- 0.0  
  
  # calculate score for each pairwise comparison of feature sets
  scores <- zucknickScores(features = features, F.all = F.all, sim.mat = sim.mat)
  
  # calculate mean score over pairwise feature set comparisons
  score <- mean(scores)
  
  return(score)
}

# function to calculate the scores for the pairwise feature sets

zucknickScores <- function(features, F.all, sim.mat){
  
  measureScoreHelper(features = features,
                     
                     measureFun= function(F1, F2) {
                       
                       # calculation the union of F1, and F2, the denominator in the score, check > 0
                       (lu = length(union(F1, F2)))
                       if (lu == 0) {
                         return(NA_real_)
                       }
                       
                       # get the indices of the features selected in the feature sets being compared
                       (indices.1 = which(F.all %in% F1))
                       (indices.2 = which(F.all %in% F2))
                       
                       # calc the indices of the features in F1, but not in F2 and vice versa
                       (indices.1not2 = which(F.all %in% setdiff(F1, F2)))
                       (indices.2not1 = which(F.all %in% setdiff(F2, F1)))
                       
                       # slice similarity matrix to get thresholded correlation between genes in F2 but not F1, with all genes in F1; average
                       add.sim1 = 0
                       if (length(indices.2not1) > 0) {
                         sim.part = sim.mat[indices.2not1, indices.1, drop = FALSE]
                         add.sim1 = sum(sim.part) / length(indices.2)
                       }
                       
                       # slice similarity matrix to get thresholded correlation between genes in F1 but not F2, with all genes in F2; average
                       add.sim2 = 0
                       if (length(indices.1not2) > 0) {
                         sim.part = sim.mat[indices.1not2, indices.2, drop = FALSE]
                         add.sim2 = sum(sim.part) / length(indices.1)
                       }
                       
                       # calculate intersection (numerator of jaccard), and add the similarity values
                       (res = length(intersect(F1, F2)) + add.sim1 + add.sim2)
                       
                       # calculate score
                       (res = res / lu)
                       
                       return(res)
                     })
}


# function to enumerate out the pairwise comparisons

measureScoreHelper = function(features, measureFun) {
  
  # calculate the number of subsamples n
  (n = length(features))
  scores = unlist(lapply(1:(n - 1), function(i) {
    sapply((i + 1):n, function(j) {
      measureFun(features[[i]], features[[j]])
    })
  }))
  return(scores)
}


# function for regular jaccard
## UPDATE THIS AS ABOVE TO GET THE FULL FUNCTION FOR COMPARISON

measureJaccard <-  function(F1, F2) {
  
  # calculation the union of F1, and F2, the denominator in the score, check > 0
  (lu = length(union(F1, F2)))
  if (lu == 0) {
    return(NA_real_)
  }
  
  linsect = length(intersect(F1, F2))
  
  # calculate score
  (res = linsect / lu)
  
  return(res)
}


# Sechidis stability measure
# function to calculate the absolute value of the correlation
# with coop package - inplace=FALSE faster but memory intensive - see documentation

## updated function for sechidis stability score, include the univariate filter as part of the feature selection process

getStabilitySechidisOriginalSetSize = function(features, C, threshold) {
  
  #' Calculates sechidis stability measure, normalising with reference to the full original dataset size prior to any univariate filtering
  #' Has not impact if mapped to a series of datasets the same size
  #' 
  #' @param features the list of the sets of selected features
  #' @param C similarity matrix of all features based on all datasets (i.e before filtering) e.g. pearson corr, values btw 0 and 1, must be named with the elements that appear in the features vectors
  #' @param threshold threshold for indicating which features are similar or not, two features are similar if entry in sim.mat is >=threshold
  
  # get the total number of features
  (p = ncol(C))
  
  # calculate the number of subsamples n, and the number of features p over the datset
  (n = length(features))
  
  # get names of all the features from the column headings of C
  F.all <- colnames(C)
  
  # create a zeros matrix, with number of columns equal to the total number of features over the datasets, rows = number of subsamples
  (Z = matrix(0, nrow = n, ncol = p))
  
  # code for function to create the selected features index, from the list of features selected each run and the full feature list
  for (i in seq_along(features)) {
    Z[i, ] = as.numeric(F.all %in% features[[i]])
  }
  
  # check to ensure that the mean number of features selected is not zero (i.e. no features selected) or equal to the total number of features
  ns.mean = mean(rowSums(Z))
  if (ns.mean == 0 || ns.mean == p) {
    return(NA_real_)
  }
  
  # set values in C below threshold to zero
  C[C < threshold] <- 0.0
  
  # calculate S, the covariance matrix of Z
  S = covar(Z)
  
  # calculate the elements used in the denominator - k bar, k squared b, and construct the matrix big sigma
  k.bar <-  mean(rowSums(Z))
  k2.bar <- mean(rowSums(Z)^2)
  diag.element <-  k.bar/p * (1  - k.bar/p)
  off.diag.element <-  (k2.bar - k.bar) / (p^2 - p) - k.bar^2 / p^2
  Sigma0 <-  matrix(off.diag.element, nrow = p, ncol = p)
  diag(Sigma0) <-  diag.element
  
  # new section: transform matrices to sparse format to speed up computations
  C <- Matrix(C, sparse = T)
  S <- Matrix(S, sparse = T)
  
  # construct the numerator and denominator of the term in the score
  num = sum(diag(crossprod(t(C), S)))
  denom = sum(diag(crossprod(t(C), Sigma0)))
  
  # make sure the denominator is not zero, then calculate score
  if (denom == 0) {
    return(NA_real_)
  } else {
    score = 1 - num / denom
    return(score)
  }
}


# function for sechidis stability score

getStabilitySechidis = function(Z, C, threshold) {
  
  #' Calculates sechidis stability measure, taking the number of features from the dimensions of the input dataset
  #' 
  #' @param Z the matrix of selected features - n runs * p features, each row represents a feature selection
  #' @param C similarity matrix of all features based on all datasets (i.e. boostraps), e.g. pearson corr, values btw 0 and 1, must be named with the elements that appear in the features vectors
  #' @param threshold threshold for indicating which features are similar or not, two features are similar if entry in sim.mat is >=threshold
  
  # calculate the number of subsamples n, and the number of features p over the datset
  (n = nrow(Z))
  (p = ncol(Z))
  
  # check to ensure that the mean number of features selected is not zero (i.e. no features selected) or equal to the total number of features
  ns.mean = mean(rowSums(Z))
  if (ns.mean == 0 || ns.mean == p) {
    return(NA_real_)
  }
  
  # set values in C below threshold to zero
  C[C < threshold] <- 0.0
  
  # calculate S, the covariance matrix of Z
  S = covar(Z)
  
  # calculate the elements used in the denominator - k bar, k squared b, and construct the matrix big sigma
  k.bar <-  mean(rowSums(Z))
  k2.bar <- mean(rowSums(Z)^2)
  diag.element <-  k.bar/p * (1  - k.bar/p)
  off.diag.element <-  (k2.bar - k.bar) / (p^2 - p) - k.bar^2 / p^2
  Sigma0 <-  matrix(off.diag.element, nrow = p, ncol = p)
  diag(Sigma0) <-  diag.element
  
  # new section: transform matrices to sparse format to speed up computations
  C <- Matrix(C, sparse = T)
  S <- Matrix(S, sparse = T)
  
  # construct the numerator and denominator of the term in the score
  num = sum(diag(crossprod(t(C), S)))
  denom = sum(diag(crossprod(t(C), Sigma0)))
  
  # make sure the denominator is not zero, then calculate score
  if (denom == 0) {
    return(NA_real_)
  } else {
    score = 1 - num / denom
    return(score)
  }
}


# logistic regression functions using glmnet

# glmnet: useful sources
# https://glmnet.stanford.edu/articles/glmnet.html#logistic-regression-family-binomial-
# https://stackoverflow.com/questions/48709211/cv-glmnet-and-leave-one-out-cv
# https://stackoverflow.com/questions/70287268/obtain-variable-selection-order-from-glmnet
# https://slowkow.com/notes/sparse-matrix/ 


cv.lr.feature.selection <- function(x,y,k,a){
  
  #' Performs feature selection using penalised logistic regression, with a k-fold cross validation loop
  #'
  #' @param x the input matrix of counts, training set, either vst or rlog transformed. row = samples, columns = features
  #' @param y the class labels
  #' @param k the number of folds in the cross validation loop, can be set to loocv
  #' @param a alpha, balancing parameter between lasso (a=1) and ridge (a=0)
  
  # determine number of folds
  if(k=="loocv"){k <- length(y)}
  
  # fit cross validated model using loocv and mis classification error as loss
  cvfit <- cv.glmnet(x, y, family='binomial', alpha=a, type.measure = 'class', nfolds = k, grouped=FALSE, keep = TRUE, parallel = FALSE)
  
  # score the predictions on the validation set against true labels ($fit.preval are the validation set predictions as log odds)
  perf.scores.val <- assess.glmnet(cvfit$fit.preval, newy = y, family = "binomial")
  
  # create dataframe of lambda, misclassification error, auc and number of features
  (lambda.features.df <- data.frame('lambda' = cvfit$lambda, 'class.error' = cvfit$cvm, 'auc' = perf.scores.val$auc, 'num.features' = cvfit$nzero))  
  
  # get the best auc value
  (best.auc  <- max(lambda.features.df$auc))
  
  # filter the dataframe to only the lambda values that maximised auc
  (lambda.features.df.fil <- dplyr::filter(lambda.features.df, auc == best.auc))
  
  # select the maximum value of lambda within the filtered set => most regularised model, that maximises auc
  (best.lambda  <- max(lambda.features.df.fil$lambda))
  
  # create coefficients object from the model fitted on the full training dataset, at the optimal value of lambda
  (sparse.coefs.matrix   <- coef(cvfit$glmnet.fit, s = best.lambda)) 
  
  # check that the coefs matrix is correct size
  if(!identical(length(sparse.coefs.matrix[-1]), ncol(x))){stop()}
  
  # get the indices of the selected features, removing the intercept which is indexed at zero
  (fs.idx <- sparse.coefs.matrix@i[-1])
  
  # return the error and the features selected by the model
  return(list('results.df' = lambda.features.df, 'lambda.val'= best.lambda, 'auc.val' = best.auc, 'fs.idx'=fs.idx))
}


cv.lr.feature.selection.feat.range <- function(x,y,k,a,fmin,fmax){
  
  #' Performs feature selection using penalised logistic regression, with a k-fold cross validation loop, selecting number of features in a user specified range
  #'
  #' @param x the input matrix of counts, training set, either vst or rlog transformed. row = samples, columns = features
  #' @param y the class labels
  #' @param k the number of folds in the cross validation loop, can be set to loocv
  #' @param a alpha, balancing parameter between lasso (a=1) and ridge (a=0)
  #' @param fmin the minimum number of features to select
  #' @param fmax the maximum number of features to select
  
  # determine number of folds
  if(k=="loocv"){k <- length(y)}
  
  # fit cross validated model using loocv and mis classification error as loss
  cvfit <- cv.glmnet(x, y, family='binomial', alpha=a, type.measure = 'class', nfolds = k, grouped=FALSE, keep = TRUE, parallel = FALSE)
  
  # score the predictions on the validation set against true labels ($fit.preval are the validation set predictions as log odds)
  
  perf.scores.val <- assess.glmnet(cvfit$fit.preval, newy = y, family = "binomial")  # contains auc
  (cm.val <- confusion.glmnet(cvfit$fit.preval, newy = y, family = "binomial"))      # gives confusion matrices
  (perf.scores.val.f1 <- map(cm.val, scoreConfusionMatrix) %>% unlist())             # f1 score from confusion matrices
  
  # create dataframe of lambda, misclassification error, auc, f1 score, and number of features
  (lambda.features.df <- data.frame('lambda' = cvfit$lambda, 'class.error' = cvfit$cvm, 'auc' = perf.scores.val$auc, 'f1' = perf.scores.val.f1, 'num.features' = cvfit$nzero))  
  
  # if min and max number of features supplied, filter dataframe to feature numbers within that range
  
  if(!is.na(fmin) && !is.na(fmax) && fmin<=fmax){
    
    # filter the dataframe to only the lambda values that result in a number of features in the fmin fmax range
    (lambda.features.df <- dplyr::filter(lambda.features.df, num.features >= fmin & num.features <= fmax))  
  }
  # get the best f1 score
  (best.f1  <- max(lambda.features.df$f1))
  
  # filter the dataframe to only the lambda values that maximised f1
  (lambda.features.df.fil <- dplyr::filter(lambda.features.df, f1 == best.f1))
  
  # select the minimum value of lambda within the filtered set => model with number of features closest to the fmax
  (best.lambda  <- min(lambda.features.df.fil$lambda))
  
  # create coefficients object from the model fitted on the full training dataset, at the optimal value of lambda
  (sparse.coefs.matrix   <- coef(cvfit$glmnet.fit, s = best.lambda)) 
  
  # check that the coefs matrix is correct size
  if(!identical(length(sparse.coefs.matrix[-1]), ncol(x))){stop()}
  
  # get the indices of the selected features, removing the intercept which is indexed at zero
  (fs.idx <- sparse.coefs.matrix@i[-1])
  
  # return the error and the features selected by the model
  return(list('results.df' = lambda.features.df, 'lambda.val'= best.lambda, 'f1.val' = best.f1, 'fs.idx'=fs.idx))
}

scoreConfusionMatrix <- function(cm){
  #' Helper function to convert a confusion matrix to f1 score
  #'
  #' @param cm confusion matrix of the form output by confusion.glmnet
  
  fpn <- function(x){tryCatch(x, silent = TRUE,error = function(e){return(0)})}
  
  (fn = fpn(cm['-1','1']))  
  (tp = fpn(cm['1','1']))  
  (fp = fpn(cm['1','-1']))   
  (tn = fpn(cm['-1','-1']))
  
  (recall = tp / (tp + fn))
  (precision = tp / (tp + fp))
  (f1 = (2*precision*recall)/(precision + recall))
  
  return(f1)
}


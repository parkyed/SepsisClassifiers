
rfNFeatModel <- function(x, y, n.feature, ntree){
  
  # train a model
  (rf.model <- randomForest(x = x, y = y, importance=TRUE, ntree=ntree))
  
  # variable importance, just pull out the gini base measure
  (rf.feat.imp <- importance(rf.model)[, c('MeanDecreaseGini'), drop=FALSE] %>% as.data.frame(.) %>% arrange(., desc(MeanDecreaseGini)))
  
  # get the feature names of the best.n.features based on list ranked mean decrease in gini coefficient
  (best.n.features <- rownames(rf.feat.imp[1:n.feature,, drop=FALSE]))
  
  # get the index of the selected features based on the original feature set
  (fs.idx = which(colnames(x) %in% best.n.features))
  
  # subset the original dataset to just the selected features
  (x.selected.features = x[,fs.idx])                   # check: all(best.n.features %in% colnames(x.selected.features))
  
  # retrain a model on the selected features only
  (rf.model.retrain = randomForest(x = x.selected.features, y = y, importance=FALSE, ntree=ntree))
  
  # predict on the out of bag examples
  (rf.pred = predict(rf.model.retrain, type='response', drop=FALSE))
  
  # performance scores
  (rf.scores = caret::confusionMatrix(data = rf.pred, reference = y)$byClass)
  
  return(list('scores' = rf.scores, 'best.n.features'= best.n.features, 'f1.val' = rf.scores['F1'], 'fs.idx'=fs.idx))
}

# basic rf cv model to determine optimum number of features and extract the index of the selected feature set
# note: results are different depending on whether recursive is set to true or false. no way in this library to determine the best feature set if recursive is set to true, so have to retrain and just use this number to extract the most important features

rfRFEModel <- function(x, y, n.folds, ntree){
  
  #' Performs feature selection using Random Forest with recursive feature elimination (using the randomForest package built in simplified rfe)
  #' Applies k-fold cross validation to select the optimal number of parameters; parameters eliminated based on mean decrease in gini coeff over forest
  
  #' @param x the input matrix of counts, training set, either vst or rlog transformed. row = samples, columns = features
  #' @param y the class labels
  #' @param n.folds the number of folds in the cross validation loop
  
  # set correct number of folds for leave-one-out-cv
  if(n.folds=="loocv"){n.folds <- length(y)}
  
  # train cross validated model for the optimum number of features
  rfcv = rfcv(trainx = x, trainy = y, cv.fold=n.folds, recursive = TRUE, ntree=ntree)
  
  # get the performance scores for each number of features tested
  rf.cv.res <-map_dfr(rfcv$predicted, function(x){caret::confusionMatrix(data = x, reference = y)$byClass[c('Sensitivity', 'Specificity', 'Precision', 'F1')]})
  rf.cv.res$n.features <- as.integer(unlist(names(rfcv$predicted)))
  
  # determine best hyperparameters and corresponding accuracy score, using number of features and accuracy to decide
  rf.cv.res.fil <- rf.cv.res %>%
    
    dplyr::filter(F1 == max(F1, na.rm = TRUE)) %>%     # filter to best classification performance, remove NAs
    dplyr::filter(n.features == min(n.features))       # get minimum number of features for given performance => single line of results
  
  (best.n.features = rf.cv.res.fil$n.features)
  (best.f1.score = rf.cv.res.fil$F1)
  
  # retrain a model without cross validation, and extract the n.features with the top feature importance
  (rf.retrain <- randomForest(x = x, y = y, importance=TRUE, ntree=ntree))
  
  # variable importance, just pull out the gini base measure
  (rf.feat.imp <- importance(rf.retrain)[, c('MeanDecreaseGini'), drop=FALSE] %>% as.data.frame(.) %>% arrange(., desc(MeanDecreaseGini)))
  
  # get the feature names of the best.n.features based on list ranked mean decrease in gini coefficient
  (best.features <- rownames(rf.feat.imp[1:best.n.features,, drop=FALSE]))
  
  # get the index of the selected features based on the original feature set
  (fs.idx = which(colnames(x) %in% best.features))
  
  return(list('results.df' = rf.cv.res, 'best.n.features'= best.n.features, 'f1.val' = best.f1.score, 'fs.idx'=fs.idx))
}

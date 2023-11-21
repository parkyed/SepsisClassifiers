
normalize <- function(x) {
  
  #' normalise a matrix of count data to values between 0 and 1
  
  return ((x - min(x)) / (max(x) - min(x)))
}


createFoldDataFrame <- function(n.folds, x.mat, y.true){
  
  #' Create a data frame with the train, and validation sets for each fold in the cross validation scheme
  
  # Create the index for k equally size folds
  fold.idx <- cut(seq(1,nrow(x.mat)),breaks=n.folds,labels=FALSE)
  
  # create the required number of folds based on the index
  svm.cv.folds = lapply(c(1:n.folds), createFoldData, fold.idx=fold.idx, x.mat=x.mat, y.true=y.true)
  
  # format as a dataframe
  svm.cv.folds.df <- as.data.frame(do.call(rbind, svm.cv.folds))
  
  return(svm.cv.folds.df)
}


createFoldData <- function(fold, fold.idx, x.mat, y.true){
  
  #' Function creates the x.train, y.train, x.val and y.val sets for a single fold based on a fold index
  #' Returns a list of length 4, with x.train, y.train, x.val, and y.val as the four list items
  
  val.idx <- which(fold.idx==fold, arr.ind=TRUE)
  x.val  <- x.mat[val.idx, , drop=FALSE]
  y.val  <- y.true[val.idx]
  x.train <- x.mat[-val.idx, ,  drop=FALSE]
  y.train <- y.true[-val.idx]
  
  return(list(x.train=x.train, y.train=y.train, x.val=x.val, y.val=y.val))
}


svmTrainPredict <- function(x.train, y.train, x.val,lambda,...){
  
  #' Function trains an L1 penalised SVM model using the LiblineaR package and predicts on a validation set, for a given cost value
  #' Returns a vector of predictions
  #' 
  #' @param lambda the value of lambda to use in the model
  #' @param x.train count data set for this fold
  #' @param y.train class labels for this fold
  
  # train model
  model = LiblineaR(data = x.train, target = y.train, type = 5, cost = lambda)
  
  # predict on the validation set (in case of loocv, this is just one example)
  y.pred = predict(model, newx = x.val, proba = FALSE, decisionValues = FALSE)$predictions
  
  return(y.pred = y.pred)
}


svmPenalisedCV <- function(lambda1, svm.cv.folds.df, x.mat, y.true){
  
  #' Cross validation function that runs train and predict over vectorised set of folds, for a single value of lambda
  #' Returns the performance on the validation set in the cross validation loop - accuracy, f1 score
  #' Retrains a model on the full data set to identify the indices of the features selected, and the number of features
  #'
  #' @param lambda1 the value of lambda to use in this run of the cross validation over all the folds
  #' @param svm.cv.folds.df data frame with the train, and validation sets for each fold in the cross validation scheme
  #' @param x.mat the full training set count values used to retrain the model and select the features for this value of lambda
  #' @param y.true the full training set class labels used to retrain the model and select the features for this value of lambda
  
  # map all folds to the train / predict function to return the predictions on the validation set, unlist to create pred list over all samples
  y.pred <- pmap(svm.cv.folds.df, svmTrainPredict, lambda=lambda1) %>% unlist(.)
  
  # performance evaluation: score the predictions vs. the true labels
  cm = caret::confusionMatrix(factor(y.pred, levels=c(1, -1)), factor(y.true, levels=c(1, -1)))
  accuracy = cm$overall['Accuracy']
  f1 = cm$byClass['F1']
  sensitivity = cm$byClass['Sensitivity']
  specificity = cm$byClass['Specificity']
  precision = cm$byClass['Precision']
  
  # feature selection: refit the model on the full training dataset and extract the selected features
  model.retrain = LiblineaR(data = x.mat, target = y.true, type = 5, cost = lambda1)
  fs.idx = which(model.retrain$W!=0)
  n.features    = length(fs.idx)
  
  # return performance metrics and selected feature information
  return(list(lambda = lambda1, accuracy = accuracy, f1.score = f1, sensitivity = sensitivity, specificity = specificity, precision= precision, n.features = n.features, fs.idx = fs.idx))
}


svmCVModel <- function(x.matrix, y.labels, n.folds, param.grid){
  
  #' Performs feature selection using L1 regularised linear support vector machine
  #' Applies k-fold cross validation to select the best value of the regularisation parameter lambda, using classification accuracy
  #' A maximum number of features can be selected to constrain the amount of regularisation
  
  #' @param x.matrix the input matrix of counts, training set, either vst or rlog transformed. row = samples, columns = features
  #' @param y.labels the class labels
  #' @param n.folds the number of folds in the cross validation loop
  #' @param param.grid the vector of values of lambda to test in the cross validation loop
  #' @param max.features the maximum number of features allowed in the output feature set
  
  # set correct number of folds for leave-one-out-cv
  if(n.folds=="loocv"){n.folds <- length(y.labels)}  
  
  # create the folds dataframe
  cv.folds.df <- createFoldDataFrame(n.folds, x.matrix, y.labels)
  
  # map cross valuation function to all values of lambda - doing the grid search, convert to dataframe
  svm.cv.res <- map(param.grid, svmPenalisedCV, cv.folds.df, x.matrix, y.labels) %>% do.call(rbind, .) %>% as.data.frame(.)
  
  # convert numeric columns to type numeric
  svm.cv.res[ , c(1:3)] <- apply(svm.cv.res[ , c(1:3)], 2, function(x) as.numeric(x))
  
  # determine best hyperparameters and corresponding accuracy score, using number of features and accuracy to decide
  svm.cv.res.fil <- svm.cv.res %>%
    
    dplyr::filter(f1.score == max(f1.score, na.rm = TRUE)) %>%     # filter to best classification performance, remove NAs
    dplyr::filter(lambda == min(lambda))                           # get minimum value of lambda, most regularised, for given performance => single line of results
  
  best.lambda = svm.cv.res.fil$lambda
  best.f1.score = svm.cv.res.fil$f1.score
  fs.idx = unlist(svm.cv.res.fil$fs.idx)
  
  # return the error, features selected by the model as an index and by ID
  
  return(list('results.df' = svm.cv.res[-8], 'lambda.val'= best.lambda, 'f1.val' = best.f1.score, 'fs.idx'=fs.idx))
}


# code for testing component functions with single values of lambda
#lambda1 = 0.1
#lambda = 0.1
#pmap(svm.cv.folds.df, svmTrainPredict, lambda=lambda1) %>% unlist(.)
#svmPenalisedCV(lambda1, svm.cv.folds.df, x.mat, y.true)

# alternative apply / map functions for the dataframe to the svmFunction
#svm.ll.res.map2 <- map2(svm.cv.folds.df$x.train, svm.cv.folds.df$y.train, svmFunction)
#svm.ll.res.mapp <- mapply(svmFunction, svm.cv.folds.df$x.train, svm.cv.folds.df$y.train, USE.NAMES = TRUE)

# performance scoring code - not used, as doesn't deal with NA values in the confusion matrix

# performance evaluation: score the predictions vs. the true labels
#(cm = as.matrix(table(y.pred, y.true)))
#(tp = cm[2,2])
#(fp = cm[2,1])
#(fn = cm[1,2])
#(tn = cm[1,1])
#(accuracy = (tp + tn)/sum(cm))
#(sensitivity = tp / (tp + fn))
# #(specificity = tn / (tn + fp))
# (recall = sensitivity)
# (precision = tp / (tp + fp))
# (f1 = (2*precision*recall)/(precision + recall))



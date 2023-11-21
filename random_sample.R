# functions to create a % random stratified sample of a data matrix and labels, stratified by the label classes

sampleFunction <- function(xy.mat, pcent.samp){
  
  #' Helper function that creates the stratified sample of the dataframe
  
  samples.output <- xy.mat %>% dplyr::group_by(class) %>%
    dplyr::slice_sample(prop=pcent.samp) %>%
    ungroup() %>%
    tibble::column_to_rownames(var = "analysisID")
  return(samples.output)
}

createRandomSample <- function(x, y, n.samp, pcent.samp){
  
  #' Creates n.samp stratified random samples of the dataset, using pcent.samp % of samples
  #' 
  #' @param x the input matrix of counts, training set, either vst or rlog transformed. row = samples, columns = features
  #' @param y the class labels
  #' @param n.samp the number of random samples to create
  #' @param pcent.samp the percentage of samples to sample, per class
  
  # combine x and y and create stratified samples
  
  xy.mat <- as.data.frame(x) %>% tibble::rownames_to_column(var = "analysisID")
  
  xy.mat$class <- y
  
  # sample the dataset the required number of times
  
  xy.samples <- replicate(n.samp, sampleFunction(xy.mat, pcent.samp), simplify=FALSE)
  
  return(xy.samples)
}

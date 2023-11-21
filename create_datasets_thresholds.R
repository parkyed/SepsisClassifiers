# script with functions that create alternative filtered datasets, based on a range of filter criteria
# dataset filters are
# noise: filtering out low expression genes. A threshold of zero, means that all genes with all zeros are filtered out, since filter is > threshold.
# outlier features: filtering out genes that have influential outlier values in at least one sample, based on Cook's distance
# low effect size features: filtering out genes where the log2 fold change between the classes is below a given threshold
# scripts rely on normalisation functions


createNoiseThresholdDatasets <- function(raw.counts, norm.counts, sample.info, t.range){
  
  #' Generate a list of datasets, each with a different low read count (noise) filter, vst transformed and transposed for ml input
  #'
  #' @param raw.counts the matrix of raw read counts
  #' @param norm.counts the matrix of median-ration (DESeq2) normalised read counts
  #' @param sample.info matrix of phenotype / clinical data for each sample
  #' @param t.range vector of noise thresholds to use
  
  # filter normalised counts on alternative thresholds, and apply index to filter the raw counts
  
  print("creating filtered raw count datasets based on alternative noise thresholds")
  datasets.filtered <- vector("list", length = length(t.range))
  for(i in seq_along(t.range)){
    idx.threshold <- which(apply(norm.counts, 1, function(x){sum(x > t.range[i]) >= 1}))
    raw.counts.filtered <-  dplyr::filter(as.data.frame(raw.counts), row_number() %in% idx.threshold)
    datasets.filtered[[i]] <- raw.counts.filtered
  }
  
  ## vst transform the filtered raw counts
  
  print("vst transforming filtered raw count datasets")
  datasets.noise.vst <- map(datasets.filtered, vst_transform, targets=sample.info)
  
  # transpose all the datasets ready for the model input
  
  print("transposing vst transformed datasets for input to ml model")
  datasets.noise.vst.t <- map(datasets.noise.vst, function(x){t(x)})
  
  return(datasets.noise.vst.t)
}

createOutlierThresholdDatasets <- function(raw.counts,norm.counts,sample.info,noise.fil,q.range){

  #' Generate a list of datasets, each with a different outlier read count (outlier) filter, vst transformed and transposed for ml input
  #'
  #' @param raw.counts the matrix of raw read counts
  #' @param norm.counts the matrix of median-ration (DESeq2) normalised read counts
  #' @param sample.info matrix of phenotype / clinical data for each sample
  #' @param noise.fil the noise filter threshold to apply to all datasets
  #' @param t.range vector of noise thresholds to use
  
  ## apply a reasonable baseline counts filter given as noise.fil parameter
  
  idx.threshold <- which(apply(norm.counts, 1, function(x){sum(x > noise.fil) >= 1}))
  counts.data.filtered <-  dplyr::filter(as.data.frame(raw.counts), row_number() %in% idx.threshold) # line could be removed, slice above

  ## run full DESeq function, does the differential expression analysis and calculates the Cooks distance
  
  print("DESeq2 differential expression analysis to calculate Cook's distances")
  dds_cooks <- dds_object(counts.data.filtered, sample.info)
  DEanalysisObject <-  DESeq(dds_cooks)
  
  # extract the calculated cooks distances for each gene in each sample
  
  print("Extract per gene cooks distances and create alternative datasets")
  
  cooks.mat <- assays(DEanalysisObject)[["cooks"]]

  # calculate the cooks outlier threshold based on the F-distribution
  
  m <- ncol(DEanalysisObject)     # number of samples
  p <- 2                          # number of model parameters (in the two condition case)
  
  # set a range of cooks distribution thresholds
  
  threshold.range.cooks <- qf(q.range, p, m - p)
  
  # filter on alternative thresholds
  
  datasets.cooks <- vector("list", length = length(threshold.range.cooks))
  for(i in seq_along(threshold.range.cooks)){
    idx.threshold.cooks <- which(apply(cooks.mat, 1, function(x){(max(x) < threshold.range.cooks[i])}))
    counts.filtered <-  dplyr::filter(as.data.frame(counts.data.filtered), row_number() %in% idx.threshold.cooks) # line could be removed, slice above
    datasets.cooks[[i]] <- counts.filtered
  }
  
  ## vst transform the filtered data set
  
  print("vst transforming outlier filtered raw count datasets")
  datasets.cooks.vst <- map(datasets.cooks, vst_transform, targets=sample.info)
  
  # transpose the datasets
  
  print("transposing vst transformed datasets for input to ml model")
  datasets.cooks.vst.t <- map(datasets.cooks.vst, function(x){t(x)})
  
  return(datasets.cooks.vst.t)
}


createLCOLDatasets <- function(raw.counts, norm.counts, sample.info, noise.fil, quantile){
  
  #' Filter an RNASeq raw count dataset based on a low counts and a cooks quantile filter, vst tranform and transpose
  #'
  #' @param raw.counts the matrix of raw read counts
  #' @param norm.counts the matrix of median-ration (DESeq2) normalised read counts
  #' @param sample.info matrix of phenotype / clinical data for each sample
  #' @param noise.fil the noise filter threshold to apply to all datasets
  #' @param t.range vector of noise thresholds to use
  #' @param quantile the quantile of the cooks distribution to use as a cutoff
  
  ## filter out zero counts
  raw.counts.nozeros <-  raw.counts[which(apply(norm.counts, 1, function(x){sum(x > 0) >= 1})),]
  norm.counts.nozeros <-  norm.counts[which(apply(norm.counts, 1, function(x){sum(x > 0) >= 1})),]
  
  ## get index for genes that pass the low counts filter
  idx.threshold <- apply(norm.counts.nozeros, 1, function(x){sum(x > noise.fil) >= 1})
  
  ## get cooks distances using deseq2
  dds_cooks <- dds_object(raw.counts.nozeros, sample.info)
  DEanalysisObject <-  DESeq(dds_cooks)
  cooks.mat <- assays(DEanalysisObject)[["cooks"]]
  
  # calculate the cooks outlier threshold based on the F-distribution
  m <- ncol(DEanalysisObject)     # number of samples
  p <- 2                          # number of model parameters (in the two condition case)
  threshold.cooks <- qf(quantile, p, m - p)  
  
  # index based on cooks distance filter
  idx.threshold.cooks <- apply(cooks.mat, 1, function(x){(max(x) < threshold.cooks)})
  
  # combine filters
  idx.lcol <- which(as.logical(idx.threshold * idx.threshold.cooks))
  counts.filtered <-  dplyr::filter(as.data.frame(raw.counts.nozeros), row_number() %in% idx.lcol) # line could be removed, slice above
  
  ## vst transform the filtered data set
  counts.filtered.vst <- vst_transform(counts.filtered, targets=sample.info)
  
  # transpose the dataset
  counts.filtered.vst.t <- t(counts.filtered.vst)
  
  return(counts.filtered.vst.t)
}


# sort by log2 fold change and add column to identify 'low effect' size genes below absolute value of 1.5 log2 fold change

createFoldChangeThresholdDatasets <- function(raw.counts, norm.counts, sample.info, noise.fil, fc.range){
  
  #' Generate a list of datasets, each with a different log2 fold change (effect size) filter, vst transformed
  #'
  #' @param raw.counts the matrix of raw read counts
  #' @param norm.counts the matrix of median-ration (DESeq2) normalised read counts
  #' @param sample.info matrix of phenotype / clinical data for each sample
  #' @param noise.fil the noise filter threshold to apply to all datasets
  #' @param fc.range vector of fold change thresholds to use
  
  ## apply a reasonable baseline counts filter given as noise.fil parameter
  idx.threshold <- which(apply(norm.counts, 1, function(x){sum(x > noise.fil) >= 1}))
  counts.data.filtered <-  dplyr::filter(as.data.frame(raw.counts), row_number() %in% idx.threshold)
  
  ## run full DESeq function, does the differential expression analysis and calculates the Log2FoldChanges
  print("DESeq2 differential expression analysis to calculate Log2 Fold Changes")
  dds_cooks <- dds_object(counts.data.filtered, sample.info)
  DEanalysisObject <-  DESeq(dds_cooks)
  
  # extract the results of the deseq2 analysis and save as a dataframe
  (deseq.results.df <- as.data.frame(results(DEanalysisObject)) %>% 
      dplyr::select(log2FoldChange))
  
  # filter on alternative thresholds
  
  print("Filtering dataset on alternative Log2FoldChange Thresholds")
  datasets.fc <- vector("list", length = length(fc.range))
  for(i in seq_along(fc.range)){
    idx.threshold.fc <- which(apply(deseq.results.df, 1, function(x){abs(x) >= fc.range[i]}))
    counts.filtered <-  dplyr::filter(as.data.frame(counts.data.filtered), row_number() %in% idx.threshold.fc)
    datasets.fc[[i]] <- counts.filtered
  }
  
  ## vst transform the filtered data set
  
  print("vst transforming Fold Change filtered raw count datasets")
  datasets.fc.vst <- map(datasets.fc, vst_transform, targets=sample.info)
  
  # transpose the datasets
  
  print("transposing vst transformed datasets for input to ml model")
  datasets.fc.vst.t <- map(datasets.fc.vst, function(x){t(x)})
  
  return(datasets.fc.vst.t)
}

dds_object <- function(rc, targets) {
  
  #' Generate a deseq data set with a given experimental design
  #'
  #' Creates analysis object from table of raw counts and table of gene lengths
  #' Analysis object contains experimental design object which maps sample Ids to sample groups
  #' @param rc The table of raw counts with ensembl_gene_ID as row index and patient sampleID as column name
  #' @param targets the meta data. Rows correspond to the columns of the raw counts
  #' @param design_col the column in the data that defined the experimental design, e.g. control, infected
  #' dds_object()
  
  exptDesign <- data.frame(targets,
                           row.names = 'analysisID')
  
  exptObject <- DESeqDataSetFromMatrix(
    countData = rc,
    colData = exptDesign,
    design = ~ condition)   # note: the design factor is hard coded here, ideally should be passed to the function
  
  return(exptObject)    
}


vst_transform <- function(rc, targets) {
  
  #' Perform variance stabilising transformation on input counts data
  #'
  #' Creates dds object, then performs transformation
  #' @param rc The table of raw counts with ensembl_gene_ID as row index and patient sample ID as column name
  #' @param targets the meta data. Rows correspond to the columns of the raw counts
  #' vst_transform()

  exptObject <- dds_object(rc, targets)
  
  ## perform the transformations, set blind=false as transformation should be aware of the design for downstream analysis
  vsd.analysisObject <- varianceStabilizingTransformation(exptObject, blind=FALSE)
  
  ## pull out the counts data only, use the following:
  vsd <- assay(vsd.analysisObject)
  
  return(vsd)
}


rlog_transform <- function(rc, targets) {
  
  #' Perform rlog transformation on input counts data
  #'
  #' Creates dds object, then performs transformation
  #' @param rc The table of raw counts with ensembl_gene_ID as row index and patient sample ID as column name
  #' @param targets the meta data. Rows correspond to the columns of the raw counts

  exptObject <- dds_object(rc, targets)
  
  ## perform the transformations. set blind=false as transformation should be aware of the design for downstream analysis
  rlog.analysisObject <- rlog(exptObject, blind=FALSE)
  
  rlogd <- assay(rlog.analysisObject)
  
  return(rlogd)
}


dds_add_gl <- function(dds, gl) {
  
  #' Add a column for gene length to a deseq dataset
  #'
  #' Adds an addtional column to the mcols data in a deseq dataset for the gene length
  #' @param dds the deseq data set that we want to add the gene length to
  #' @param gl The table of gene lengths that are added to the analysis object to enable the fpkm normalisation
  #' dds_add_gl()
  
  ## ensure that the ensembl gene ids match first
  geneLengthsTemp <- gl[gl$ensemblGeneID %in% rownames(dds),]
  
  ## add in a $basepairs column to the ads - required for FPKM normalisation based on gene length
  mcols(dds)$basepairs <- geneLengthsTemp[match(rownames(dds), geneLengthsTemp$ensemblGeneID),]$geneLength
  
  return(dds)    
}

## add summary statistics to the raw counts matrix
variance_metrics <- function(x){
  as.data.frame(x) %>%
    mutate( total =       rowSums(x),
            mean  =       rowMeans(x),
            median =      rowMedians(x),
            num_zeros =   rowSums(x == 0),
            var    =      rowVars(x),
            sd    =       rowSds(x),
            cov   =       apply(x, 1, function(x){(sd(x) / mean(x) * 100)}),
    )
}

# plot a simple histogram of counts
histogram_counts <- function(counts_matrix){
  # melt dataframe or normalised counts, add the annotation column
  counts.melt <- as.data.frame(counts_matrix) %>%               
    tibble::rownames_to_column("ensemblID")  %>%      
    melt(id.vars = 'ensemblID', value.name='count', variable.name=c('analysisID'))
  
  # density plot of the distribution of each sample, split out by annotation
  histogram_plot <- ggplot(counts.melt, aes(x = count)) + 
    geom_histogram(bins = 50) + 
    ggtitle("Histogram of counts") + 
    xlab("counts") + 
    ylab("frequency") +
    theme(legend.position="none")
  
  return(histogram_plot)
}


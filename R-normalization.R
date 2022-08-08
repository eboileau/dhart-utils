if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("DESeq2")

library(edgeR)
library(DESeq2)


dge <- calcNormFactors(dge, method = "TMM")
tmm <- cpm(dge)





normalized_tmm_data <- function(input_data) {
  #'
  #' This function normalizes data using edgeR's TMM,
  #' then passes the normalized data back.
  #'
  #' @param input_data The input data passed into the function
  #

  dge <- calcNormFactors(input_data, method = "TMM")
  edger_tmm_normalized_data <- cpm(dge)

  return(edger_tmm_normalized_data)
}


normalized_deseq2_data <- function(input_data) {
  #'
  #' This function normalizes data using DESeq2,
  #' then passes the normalized data back.
  #'
  #' @param input_data The input data passed into the function
  #

  ## Create DESeq2Dataset object
  dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
  


  deseq2_normalized_data <- calcNormFactors(input_data)
  return(deseq2_normalized_data)
}
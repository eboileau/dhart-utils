if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.15")

BiocManager::install("edgeR")
BiocManager::install("DESeq2")

library(edgeR)
library(DESeq2)

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
  deseq2_normalized_data <- DESeqDataSetFromMatrix(countData = input_data, colData = meta, design = ~ sampletype)

  return(deseq2_normalized_data)
}
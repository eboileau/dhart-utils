# TODO: Next two lines can be removed if individual functions are called from libraries
library(edgeR)
library(DESeq2)

normalized_tmm_data <- function(input_data) {
  #'
  #' This function normalizes data using edgeR's TMM,
  #' then passes the normalized data back.
  #'
  #' @param input_data The input data passed into the function
  #
  dgList <- DGEList(counts=input_data)
  dgList <- calcNormFactors(dgList, method="TMM")
  edger_tmm_normalized_data <- cpm(dgList, normalized.lib.sizes=TRUE) # Obtain the absolute TMM values of the data

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
  dds <- DESeqDataSetFromMatrix(countData = input_data)
  dds <- estimateSizeFactors(dds)
  deseq2_normalized_data <- counts(dds, normalized=TRUE)

  return(deseq2_normalized_data)
}
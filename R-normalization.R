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

  meta_data <- data.frame(matrix(NA, nrow = ncol(input_data), ncol = 1))
  rownames(meta_data) <- colnames(input_data)

  dds <- DESeqDataSetFromMatrix(countData = round(input_data), colData = meta_data, design=~1)
  dds <- estimateSizeFactors(dds)
  deseq2_normalized_data <- counts(dds, normalized=TRUE)

  return(deseq2_normalized_data)
}

normalized_tpm_data <- function(input_data, gene_length) {
  #'
  #' This function normalizes data using TPM,
  #' then passes the normalized data back.
  #'
  #' @param input_data The input data passed into the function
  #

  x <- input_data / gene_length
  tpm <- t( t(x) * 1e6 / colSums(x) )

  return(tpm)
}
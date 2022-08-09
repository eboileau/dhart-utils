if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.15")

BiocManager::install("edgeR")
BiocManager::install("DESeq2")
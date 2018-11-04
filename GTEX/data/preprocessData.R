library(data.table)
library('recount')
library(DESeq2)
library(ggplot2)
library(dplyr)
library(magrittr)
library(plyr)
options(stringsAsFactors = FALSE)

# adapt to the data directory used in getData.Rmd
datadir <- "/g/huber/users/bvelten/Documents/GTEx/data/"
load(file.path(datadir, "rse_gene_ext.RData"))

# scale counts from recount
rse_gene_ext <- scale_counts(rse_gene_ext)

# select mRNA
unique(rse_gene_ext$AssemblyName_s)
# Ensembl data from release 75 for GRCh37.
# annotation files retrieved from BioMart
mRNA_file <- "/g/huber/users/bvelten/Documents/BioMart/Hsapiens_genes_BioMart.75.txt"
mRNA = read.csv(file=mRNA_file, header=T, sep="\t", stringsAsFactors=F)
rse_mRNA <- rse_gene_ext[sapply(strsplit(rownames(rse_gene_ext), "\\."), '[[',1) %in% mRNA$ens_id,]
save(rse_mRNA, file=file.path(datadir,"rse_mRNA.RData"))

# normalize using DESeq2
  #DESeq object
  dds <- DESeqDataSet(rse_mRNA, design =~ 1)

  #Filter out genes with low counts
  rs    <- rowSums(counts(dds))
  minrs <- 100
  hist(rs[ rs < 500 ], breaks=seq(0, 500, by = 10), col = "skyblue")
  abline(v = minrs, col = "red")
  dds<-dds[ rs >= minrs, ]

  # Normalize
  dds <- estimateSizeFactors(dds)
  dds_vst <- varianceStabilizingTransformation(dds)

# save output
save(dds, file=file.path(datadir,"dds.RData"))
save(dds_vst, file=file.path(datadir,"dds_vst.RData"))


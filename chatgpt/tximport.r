#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input_files <- strsplit(args[1], ",")[[1]]
output_file <- args[2]

library(tximport)
library(readr)
library(DESeq2)
library(optparse)

# Read Salmon output files
samples <- data.frame(sample = gsub("_quant/quant.sf", "", input_files))
files <- file.path(samples$sample, "quant.sf")
names(files) <- samples$sample

# Import quantification data using tximport
txi <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "lengthScaledTPM")
dds <- DESeqDataSetFromTximport(txi, samples, ~1)

# Generate gene expression matrix
counts <- counts(dds)
write_tsv(counts, output_file)

library(data.table)
library(rtracklayer)
library(parallel)
library(magrittr)
library(plyr)
library(tidyr)
library(dplyr)
library(GenomicFeatures)
library(GenomicRanges)
library(stringr)

args <- commandArgs(T)

gtf <- args[1]
output_gtf <- args[2]

layer_list <- rtracklayer::import(gtf, "gtf")

coding_gene <-
  unique(layer_list[which(layer_list$type == "gene" &
                            layer_list$gene_biotype == "protein_coding"),]$gene_id)


coding_tx_of_codinggene <-
  unique(layer_list[which((layer_list$gene_id %in% coding_gene) &
                            (layer_list$transcript_biotype == "protein_coding")
  ),]$transcript_id)

codinggene_have_codingTx <-
  unique(layer_list[which((layer_list$gene_id %in% coding_gene) &
                            (layer_list$transcript_biotype == "protein_coding")
  ),]$gene_id)
coding_gtf <-
  layer_list[which(((layer_list$type == "gene") &
                      (layer_list$gene_id %in% codinggene_have_codingTx)
  ) |
    (layer_list$transcript_id %in% coding_tx_of_codinggene)), ]

rtracklayer::export(coding_gtf, output_gtf, "gtf")

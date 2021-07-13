#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Date    : 2019-10-11 16:24:10
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$

if (!require("optparse", quietly = TRUE)) {
    install.packages("optparse")

    if (!require("optparse", quietly = TRUE)) {
        stop("optparse not be installed, please install it first.")
    }
}

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("data.table"))

opt_list <- list(
    make_option(
        c("-i", "--indir"),
        type = "character",
        help = "ExonicPart PSI output directory, in which only contains the results form the same species.",
        metavar = "character"
        ),
    make_option(
        c("-o", "--out_prefix"),
        type = "character",
        help = "Output prefix of merged file.",
        metavar = "character"
        ),
    make_option(
        c("-t", "--type"),
        type = "character",
        default = "PSI",
        help = "Which column to be merged, 'inclusion', 'exclusion', or 'PSI'.",
        metavar = "character"
        )
)

opts <- parse_args(OptionParser(
    option_list = opt_list, 
    usage = "usage: %prog [options]", 
    add_help_option = TRUE, 
    prog = "", 
    description = "Merge Exonic parts PSI output files from input directory."))

if (is.null(opts$indir)) {
    stop("-i --indir not set")
} else if (is.null(opts$out_prefix)) {
    stop("-o --out_prefix not set")
} else if (is.null(opts$type)) {
    stop("-o --type not set")
}

if (!(opts$type %in% c("exclusion", "inclusion", "PSI"))) {
    stop("'-t' should be one of 'exclusion', 'inclusion', 'PSI'")
}

MergePSI <- function(path, 
                      type = "PSI"){
  filenames <- sort(list.files(path=path, pattern="exonic_parts.psi", full.names=TRUE))
  datalist <- lapply(filenames, function(x){
    tmp <- fread(x, header=TRUE, sep = "\t")
    tmp <- tmp[, c("exon_ID", type), with = F]
    samplename <- sub("[\\._]*exonic_parts.psi$", "", basename(x))
    setnames(tmp, type, samplename)
    setkey(tmp, exon_ID)
    return(tmp)})
  all <- Reduce(function(x, y) {merge(x, y, all = T, by = "exon_ID")}, datalist)
  outname <- paste0(opts$out_prefix, "_ExonicPart_", type, ".tsv")
  fwrite(all, file = outname, sep = "\t", na = "NA", quote = FALSE)
}


MergePSI(opts$indir, opts$type)

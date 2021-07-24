library(tximport)
args <- commandArgs(T)
parser_loc <- args[1]
files <-
  list.files(parser_loc, pattern = 'genes.results$', full.names = T)
txi.rsem <-
  tximport(files,
           type = "rsem",
           txIn = FALSE,
           txOut = FALSE)
sample_label <- gsub('.genes.results', '', basename(files))
colnames(txi.rsem$counts) <- sample_label
colnames(txi.rsem$abundance) <- sample_label

saveRDS(txi.rsem, file = glue::glue('{parser_loc}/RSEM.Rds'))

library("GenomicFeatures")
library("GenomicAlignments")
options(mc.cores = 4)

args = commandArgs(T)
config = args[1]
gtf = args[2]
output = args[3]

sampleList <- read.table(config,
                         stringsAsFactors = F,
                         sep = '\t',
                         header = T)

bamFile <- paste('./data/STAR/',
                 sampleList[['names']],
                 '.Aligned.sortedByCoord.out.bam',
                 sep = '')
bamFile <- './data/STAR/B7M1.Aligned.sortedByCoord.out.bam'
bamFileL <- BamFileList(bamFile,
                        yieldSize = 200000)
bamName <- do.call(rbind, strsplit(names(bamFileL), split = "[.]"))[, 1]
names(bamFile) <- bamName

txdb <-
  makeTxDbFromGFF(gtf,
                  format = 'gtf')

ebg <- exonsBy(txdb, by = "gene")

quantFunc <- function(bamFile, ebg = ebg) {
  return(tryCatch(
    summarizeOverlaps(
      features = ebg,
      reads = bamFile,
      mode = "Union",
      singleEnd = FALSE,
      ignore.strand = TRUE,
      fragments = TRUE
    ),
    error = function(e)
      NULL
  ))
}

res <- parallel::mclapply(names(bamFile),
                          function(bamName) {
                            bamList <- BamFileList(bamFile[[bamName]], yieldSize = 200000)
                            res <- quantFunc(bamList, ebg = ebg)
                            saveRDS(res, file = glue::glue("data/quant/{bamName}.Rds"))
                            return(res)
                          }, mc.cores = 4)


# colnames(geneReadsCount.pe.nonsp) <-
#   sampleList[['names']]


# saveRDS(geneReadsCount.pe.nonsp,
#         file = output)


B7M1 <- summarizeOverlaps(
      features = ebg,
      reads = bamFileL,
      mode = "Union",
      singleEnd = FALSE,
      ignore.strand = TRUE,
      fragments = TRUE
    )
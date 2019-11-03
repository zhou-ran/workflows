library("GenomicFeatures")
library("GenomicAlignments")
options(mc.cores = 10)

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

txdb <-
  makeTxDbFromGFF(gtf,
                  format = 'gtf')

ebg <- exonsBy(txdb, by = "gene")



bamFile <- BamFileList(bamFile,
                       yieldSize = 2000000)


geneReadsCount.pe.nonsp <-
  summarizeOverlaps(
    features = ebg,
    reads = bamFile,
    mode = "Union",
    singleEnd = FALSE,
    ignore.strand = TRUE,
    fragments = TRUE
  )


colnames(geneReadsCount.pe.nonsp) <-
  sampleList[['names']]


saveRDS(geneReadsCount.pe.nonsp,
        file = output)

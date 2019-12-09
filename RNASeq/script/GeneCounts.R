library("GenomicFeatures")
library("GenomicAlignments")
options(mc.cores = 10)

args = commandArgs(T)
config = args[1]
gtf = args[2]
output = args[3]

txdb <-
  makeTxDbFromGFF(gtf,
                  format = 'gtf')

ebg <- exonsBy(txdb, by = "gene")

strandMode_ <- c(
  "none" = 0,
  "forward" = 1,
  "reverse" = 2
)

## bam parser information

flag <- scanBamFlag(isSecondaryAlignment = FALSE,
                    isNotPassingQualityControls = FALSE,
                    isUnmappedQuery = FALSE,
                    isDuplicate = FALSE)
sbp <- ScanBamParam(flag=flag, mapqFilter = 255)


## read sampleinformation
sampleList <- read.table(config,
                         stringsAsFactors = F,
                         sep = '\t',
                         header = T)

## here to split into different strand specific mode
fileNames <- split(sampleList$names, sampleList$Strand)
se <- lapply(names(fileNames),function(x){
  
  name <- fileNames[x]
  bamFile <- paste('./data/STAR/',
                 name,
                 '.Aligned.sortedByCoord.out.bam',
                 sep = '')

  bamFile <- BamFileList(bamFile,
                        yieldSize = 2000000)


  se<- summarizeOverlaps(features = ebg, 
                                reads = bamFile, 
                                mode = "IntersectionStrict",
                                ignore.strand = FALSE, 
                                inter.feature = FALSE, 
                                singleEnd = FALSE,
                                fragments = FALSE, 
                                strandMode = unname(strandMode_[x]), 
                                param = sbp, 
                                preprocess.reads = NULL)
  se
})


saveRDS(se,
        file = output)

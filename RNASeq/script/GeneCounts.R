library("GenomicFeatures")
library("GenomicAlignments")
options(mc.cores = 4)

args = commandArgs(T)
config = args[1]
gtf = args[2]
output = args[3]


txdb <-
  makeTxDbFromGFF(gtf,
                  format = 'gtf')

ebg <- exonsBy(txdb, by = "gene")

strandMode_ <- c("none" = 0,
                 "forward" = 1,
                 "reverse" = 2)

## bam parser information

flag <- scanBamFlag(
  isSecondaryAlignment = FALSE,
  isNotPassingQualityControls = FALSE,
  isUnmappedQuery = FALSE,
  isDuplicate = FALSE
)
sbp <- ScanBamParam(flag = flag, mapqFilter = 255)


## read sampleinformation
sampleList <- read.table(config,
                         stringsAsFactors = F,
                         sep = '\t',
                         header = T)

## here to split into different strand specific mode
fileNames <-
  split(sampleList[, c('type', 'Strand')], sampleList$names)
se <- lapply(names(fileNames), function(x) {
  name <- x
  strandInfo <- unname(unlist(fileNames[[x]])['Strand'])
  fragType <- unname(unlist(fileNames[[x]])['type'])
  # Must change the search dir or pattern based your output
  bamFile <- paste('./data/STAR/',
                   name,
                   '.Aligned.sortedByCoord.out.bam',
                   sep = '')
  
  bamFile <- BamFileList(bamFile,
                         yieldSize = 200000)
  
  if (fragType == 'SE') {
    se <- summarizeOverlaps(
      features = ebg,
      reads = bamFile,
      mode = "Union",
      ignore.strand = ifelse(strandInfo == 'none', TRUE, FALSE),
      inter.feature = FALSE,
      singleEnd = TRUE,
      fragments = FALSE,
      # strandMode = unname(strandMode_[strandInfo]),
      param = sbp,
      preprocess.reads = NULL
    )
  } else {
    se <- summarizeOverlaps(
      features = ebg,
      reads = bamFile,
      mode = "Union",
      ignore.strand = FALSE,
      inter.feature = FALSE,
      singleEnd = FALSE,
      fragments = FALSE,
      strandMode = unname(strandMode_[strandInfo]),
      param = sbp,
      preprocess.reads = NULL
    )
    
  }
  se
})

names(se) <- names(fileNames)

saveRDS(se,
        file = output)
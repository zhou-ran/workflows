args = commandArgs(T)
config = args[1]
output = args[2]
output_dir = args[3]


## read sampleinformation
sampleList <- read.table(config,
                         stringsAsFactors = F,
                         sep = "\t",
                         header = T)



processReadsFile <- function(fullFileName) {
  dat <- read.table(fullFileName, stringsAsFactors = T)
  dat <-
    data.frame(dat[, -1], row.names = dat[, 1])
  return(dat[-c(1:4),])
}

getGivenCol <- function(listOfDf, colId) {
  do.call(cbind, lapply(listOfDf, "[", colId)) -> res
  colnames(res) <- names(listOfDf)
  res
}


fullFileName <- file.path(output_dir,
                          paste0(sampleList$names, ".ReadsPerGene.out.tab"))
names(fullFileName) <- sampleList$names

lapply(fullFileName, processReadsFile) -> res

dat <- list(
  "none" = getGivenCol(res, "V2"),
  "forward" = getGivenCol(res, "V3"),
  "reverse" = getGivenCol(res, "V4")
)

saveRDS(dat,
        file = output)

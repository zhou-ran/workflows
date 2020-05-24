library(DNAcopy)
args <- commandArgs(T)
files <- unlist(strsplit(args[1], split=','))
labels <- unlist(strsplit(args[2], split=','))
outdir <- args[3]
sd <- 2.5

forMarkers <- lapply(1:length(files), function(inds) {
  file <- files[inds]
  dir_name <- dirname(file)
  x <- labels[inds]
  cn <-
    read.table(file, header = TRUE)
  cn$chromosome <- as.character(cn$chromosome)
  cn$end <- as.numeric(as.character(cn$end))
  return(cn)
})

forMarkers <- do.call(rbind, forMarkers)

CNA.object <-
  CNA(
    genomdat = forMarkers$log2,
    chrom = forMarkers$chromosome,
    maploc = forMarkers$end,
    data.type = 'logratio',
    sampleid = 'forMarkers'
  )

smoothed.CNA.object <- smooth.CNA(CNA.object)
  
segment.smoothed.CNA.object <-
  segment(
    smoothed.CNA.object,
    undo.splits = "sdundo",
    undo.SD = sd,
    verbose = 1
  )

markers <-
  data.frame(
    paste(
      segment.smoothed.CNA.object$data$chrom,
      segment.smoothed.CNA.object$data$maploc,
      sep = ":"
    ),
    segment.smoothed.CNA.object$data$chrom,
    segment.smoothed.CNA.object$data$maploc
  )
colnames(markers) <-
  c("Marker Name", "Chromosome", "Marker Position")

write.table(
  unique(markers),
  file = glue::glue("{outdir}/markers.tsv"),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# sd <- 2.5
res <- parallel::mclapply(1:length(files), function(inds) {
  file <- files[inds]
  dir_name <- dirname(file)
  x <- labels[inds]
  cn <-
  read.table(file, header = TRUE)
  cn$chromosome <- as.character(cn$chromosome)
  cn$end <- as.numeric(as.character(cn$end))
  
  CNA.object <-
    CNA(
      genomdat = cn$log2,
      chrom = cn$chromosome,
      maploc = cn$end,
      data.type = 'logratio',
      sampleid = x
    )
  
  smoothed.CNA.object <- smooth.CNA(CNA.object)
  
  segment.smoothed.CNA.object <-
    segment(
      smoothed.CNA.object,
      undo.splits = "sdundo",
      undo.SD = sd,
      verbose = 1
    )
  segmentation.values <- segment.smoothed.CNA.object$output
  colnames(segmentation.values) <-
    c("Sample",
      "Chromosome",
      "Start Position",
      "End Position",
      "Num markers",
      "Seg.CN")
  write.table(
    segmentation.values,
    file = glue::glue("{dir_name}/{x}.seg.tsv"),
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
  
  
  markers <-
    data.frame(
      paste(
        segment.smoothed.CNA.object$data$chrom,
        segment.smoothed.CNA.object$data$maploc,
        sep = ":"
      ),
      segment.smoothed.CNA.object$data$chrom,
      segment.smoothed.CNA.object$data$maploc
    )
  colnames(markers) <-
    c("Marker Name", "Chromosome", "Marker Position")
  write.table(
    markers,
    file = glue::glue("{dir_name}/{x}.markers.tsv"),
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
  return(segmentation.values)
}, mc.cores = 4)



res <- do.call(rbind,res)

write.table(
  res,
  file = glue::glue("{outdir}/allsample.seg.tsv"),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

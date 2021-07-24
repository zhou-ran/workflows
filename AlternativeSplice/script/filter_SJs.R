library(data.table)
library(parallel)
library(stringr)
args <- commandArgs(T)

multi_merge <-
  function(path,
           pattern = 'SJ.out.tab',
           minSJ = 10,
           minSJs = 100,
           minSamps = 2,
           uniqueMapOnly = TRUE,
           SJtype = "allAnnotatedAndCanonicallyNovel",
           cores = 16) {
    file_list <- list.files(path, pattern = pattern, full.names = F)
    # print('processing file list')
    tab_list <-
      mclapply(file_list,
               function(f) {
                 tab <- fread(file.path(path, f), sep = '\t', header = F)
                 tab[, junctions := paste(V1, V2, V3, V4, V5, V6, sep = ":")]
                 
                 if (SJtype == "all") {
                   tab <- tab
                 } else if (SJtype == "allAnnotated") {
                   tab <- tab[V6 == 1, ]
                 } else if (SJtype == "allAnnotatedAndCanonicallyNovel") {
                   tab <- tab[(V6 == 1) | (V5 != 0), ]
                 } else if (SJtype == "allCanonical") {
                   tab <- tab[V5 != 0, ]
                 } else {
                   stop("'SJtype' error")
                 }
                 
                 if (uniqueMapOnly) {
                   tab <- tab[, .(junctions, V7)]
                 } else {
                   tab <- tab[, .(junctions, V7 + V8)]
                 }
                 
                 ### assign the unique colnames for each cell
                 colnames(tab)[2] <-
                   sub("[^0-9a-zA-Z]{0,}SJ.out.tab$", "", f)
                 
                 setkey(tab, junctions)
                 return(tab)
               },
               mc.cores = cores)
    
    # print('start merging')
    merged_tab <-
      Reduce(function(x, y) {
        merge(x, y, by = 'junctions', all = T)
      },
      tab_list)
    
    info <-
      data.frame(
        stringr::str_match(merged_tab$junctions, "(.*):(.*):(.*):(.*):(.*):(.*)"),
        stringsAsFactors = F
      )
    info <- info[, 2:7]
    colnames(info) <-
      c("seqname", "start", "end", "strand", "motif", "annotation")
    info <- data.table(info)
    
    merged_tab <- merged_tab[, -1]
    
    SJ_2_keep <-
      merged_tab[, .(
        Total_cnt = apply(.SD, 1, sum, na.rm = TRUE),
        Total_sample = apply(.SD, 1, function(x) {
          sum(sum(x >= minSJ, na.rm = TRUE))
        })
      )]
    
    info <-
      info[SJ_2_keep[, Total_cnt >= minSJs &
                       Total_sample >= minSamps], ]
    
    return(info)
  }

raw_SJs_path <-
  args[1] ## STAR output directory for a species, STAROUT
filtered_SJs <- args[2] ## SJs for calculate PSI

res <-
  multi_merge(
    path = raw_SJs_path,
    pattern = ".SJ.out.tab",
    minSJ = 10,
    minSJs = 100,
    minSamps = 2,
    uniqueMapOnly = TRUE,
    SJtype = "allAnnotatedAndCanonicallyNovel",
    cores = 20
  )

fwrite(
  res,
  filtered_SJs,
  sep = "\t",
  row.names = F,
  col.names = F
)

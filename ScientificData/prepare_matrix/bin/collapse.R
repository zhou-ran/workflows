arg <- commandArgs(T)
all_files <- arg[-length(arg)]
out_files <- arg[length(arg)]

# all_files <- c("data/fastqc_raw/A1C1_R1.txt","data/fastqc_raw/B7M1_1_R1.txt","data/fastqc_raw/D6M1_1_R1.txt")
label <- gsub('.txt', '', basename(all_files))

col_name <- c("Base","Mean","Median","Lower_Quartile","Upper_Quartile","10th_Percentile","90th_Percentile")

dfs <- lapply(all_files, function(file){
	df <- read.table(file, sep = '\t', stringsAsFactor=F)
	colnames(df) <- col_name
	df
	})

names(dfs) <- label
saveRDS(dfs, file = out_files)
#! /usr/bin/env r

options(warn=-1)
suppressMessages(source("/home/brad/src/songanalysis/song_util.R"))

suppressMessages(library(parallel))
suppressMessages(library(docopt)) 
suppressMessages(library(caret))
suppressMessages(library(gbm))
suppressMessages(library(randomForest))

doc <- "Usage: findSubSongs.r [-d directory] [-f no_filter_by_size] [-t date] [-c cores] [-n nsongs]

-d --dir DIRECTORY  directory to be processed, default = '.'
-f --no_filter_by_size skip filtering files by size (0.5 to 5 MB)
-t --date DATE only process files written by this date
-c --cores CORES number of cores to use, default = 4
-n --nsongs NSONGS copy not_songs to new directory 
-h --help           show this help text"

opt <- docopt(doc)
#print(opt)
#print(length(opt$date))
if(opt$dir == "NA")
  opt$dir = "."
if(length(opt$cores) == 0)
  opt$cores = 4
if(length(opt$no_filter_by_size) > 0) {
  opt$no_filter_by_size = TRUE
} else {
  opt$no_filter_by_size = FALSE
}
  
doMC::registerDoMC(cores=opt$cores)

print("Getting file list...")
files = list.files(opt$dir, recursive = T, pattern=".wav", full.names = T)
ind = grep("songs|sized", files)
if (length(ind) > 0) files = files[-ind]
print(paste("Number of files to process: ", length(files), sep=""))

#### Filter by size ####  
if (!opt$no_filter_by_size) {
  print("Filtering by size...")
  info = file.info(files)
  if (length(opt$date) > 0) {
    info = info[as.Date(info$mtime) >= as.Date(opt$date),]
  }
  info = info[order(-info$size),]
  files = rownames(info)
  sizes = info$size
  
  sized_range = c(3E5, 2E6)
  sized_files = files[sizes>sized_range[1] & sizes<sized_range[2]]
  print(paste("Number of files with sizes between ", sized_range[1] / 1E6, " and ", sized_range[2] / 1E6, " Mbytes: ", length(sized_files), sep=""))
} else {
  sized_files = files
}
#### Find songs ####
print("Finding songs...")
chunkSize = 500
nchunks = ceil(length(sized_files) / chunkSize)
i = 1
#mod = readRDS("~/src/songanalysis/models/gbm_model.rds")
mod = readRDS("~/src/songanalysis/models/rf_model.rds")
pre = readRDS("~/src/songanalysis/models/preprocess.rds")
res = foreach(chunk=isplitVector(sized_files, chunkSize=chunkSize), .combine="rbind") %do% {
  rm(wavs)
  rm(data)
  gc()
  print(paste("Processing chunk ", i, " of ", nchunks, " chunks", sep=""))
  i = i + 1
  
  #chunk = sized_files[1:100]
  #chunk = sized_files[1000:length(sized_files)]
  #chunk = sized_files
  wavs = foreach(x=chunk, .inorder = T) %dopar% readWave(x)
  data = foreach(wav=wavs, .inorder = T) %dopar% {
    songfinder_features_clean(wav, 
                              rmnoise=F, 
                              min_duration = 30, 
                              max_gap = 10, 
                              max_duration=200,
                              min_num_peaks=4, 
                              max_num_peaks=50)
    
  }
  data = do.call(rbind, data)
  #data = as.data.frame(data)
  data = predict(pre, data)
  rownames(data) = chunk
  data = na.omit(data)
  if (nrow(data) == 0) 
    return(NULL)
  data_mc = predict(mod, data, verbose=F)
  #data_mc = Mclust(data, G=2)
  #data$predict = data_mc
  #data_summary = data.frame(data) %>% group_by(predict) %>% summarize(mean_durations = mean(mean_durations))
  #subsong = data_summary$predict[which.max(data_summary$mean_durations)]
  
  #return(data.frame(name=rownames(data), song=data$predict==subsong))
  return(data.frame(name=rownames(data), song=as.logical(data_mc)))
}
if (is.null(res)) {
  print("No songs found.")
  stop()
}
print(summary(res[,2]))
print(paste("Number of identified songs: ", nrow(res[res[,2],]), sep=""))
songs_dir = paste(opt$dir, "subsongs", sep="/")
dir.create(songs_dir, showWarnings = F)
file.copy(as.character(res[res[,2],1]), songs_dir, copy.date=T, overwrite=F)

if (length(opt$nsongs) > 0) {
  nsongs_dir = paste(opt$dir, "not_subsongs", sep="/")
  dir.create(nsongs_dir, showWarnings = F)
  file.copy(as.character(res[!res[,2],1]), nsongs_dir, copy.date=T, overwrite=F)
}
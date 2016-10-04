#! /usr/bin/env r

suppressMessages(library(parallel))
suppressMessages(library(docopt))      

doc <- "Usage: findSongs.r [-d directory] [-f no_filter_by_size] [-t date] [-m time] [-c cores] [-n nsongs] [-p param] [-l high_noise]

-d --dir DIRECTORY  directory to be processed, default = '.'
-f --no_filter_by_size skip filtering files by size (0.5 to 5 MB)
-t --date DATE only process files written by this date
-m --time TIME only process files written by this time
-c --cores CORES number of cores to use, default = 4
-n --nsongs NSONGS copy not_songs to new directory 
-r --remove delete non song files
-l --high_noise BOOLEAN wav recorded using relatively high noise setup (e.g. not powered amp, Presonus). Default is false.
-p --param PARAM
-h --help           show this help text"

opt <- docopt(doc)

options(warn=-1)
suppressMessages(source("~/src/songanalysis/song_util.R"))
suppressMessages(source("~/src/songanalysis/threshold.r"))
library(stringr)
#opt = list(dir="/mnt/bengal_home/song/bk89", date="2016-08-14", cores=7, param="/mnt/bengal_home/song/bk89/findSongs_param.config")
#opt = list(dir="/mnt/bengal_home/song/wh96pk55_2", no_filter_by_size=1, cores=10, )

if (length(opt$date) == 0 & length(opt$time) > 0)
  stop("Must specify date with time")


if(opt$dir == "NA")
  opt$dir = "."
if(length(opt$thresh) == 0)
  opt$thresh = .4
if(length(opt$cores) == 0)
  opt$cores = 4
if(length(opt$time) == 0)
  opt$time = "0:00:00"
if(length(opt$high_noise)==0) {
  opt$high_noise = FALSE
} else {
  opt$high_noise = TRUE
}
  

if(length(opt$no_filter_by_size) > 0) {
  opt$no_filter_by_size = TRUE
} else {
  opt$no_filter_by_size = FALSE
}


files = list.files(opt$dir, recursive = T, pattern=".wav", full.names = T)
ind = grep("songs", files)
if (length(ind) > 0) files = files[-ind]
print(paste("Number of files to process: ", length(files), sep=""))

## Filter by date -------------------------------------------------
file_times = parse_timestamp(files)
date = NULL
if (length(opt$date) > 0) {
  date = paste(opt$date, opt$time)
  date = as.POSIXct(date)
  print(paste("Filtering for recordings after:", date, sep=" "))
  files = files[file_times > date]
} 

#### Filter by size ####
sized_files = NULL
if (!opt$no_filter_by_size) {
  info = file.info(files)
  info = info[order(-info$size),]
  files = rownames(info)
  sizes = info$size
  
  sized_range = c(6E5, 5E6)
  sized_files = files[sizes>sized_range[1] & sizes<sized_range[2]]
  print(paste("Number of files with sizes between ", sized_range[1] / 1E6, " and ", sized_range[2] / 1E6, " Mbytes: ", length(sized_files), sep=""))
} else {
  sized_files = files
}
#### Find songs ####
print("Finding songs...")
#sized_files = sized_files[1:10]
#print(opt$low_noise)
if (!opt$high_noise) {
param = list( min_duration=15,
              max_gap=10,
              max_duration=400,
              min_num_peaks=10, 
              max_num_peaks=NULL,
              amp_ratio_max_gap=120,
              amp_ratio_min=0,
              amp_ratio_max=.2,
              rate_min=3,
              rate_max=12)
} else {
  param = list( min_duration=15,
                max_gap=10,
                max_duration=400,
                min_num_peaks=10, 
                max_num_peaks=NULL,
                amp_ratio_max_gap=120,
                amp_ratio_min=0,
                amp_ratio_max=.4,
                rate_min=3,
                rate_max=12)
}

if (length(opt$param) > 0 ) {
  param_tmp = read.delim(opt$param, sep="=", header=F, as.is = T)
  param_tmp[,1] = str_trim(param_tmp[,1])
  param_tmp[,2] = str_trim(param_tmp[,2])
  print(param_tmp)
  for (i in 1:nrow(param_tmp)) {
    param[[param_tmp[i,1]]] = param_tmp[i,2]
  }
} 

print("---- Params ----")
ind = which(unlist(lapply(param, function(x) !is.null(x))))
for (i in ind) {
  param[[i]] = as.numeric(param[[i]])
}
for (i in 1:length(param)) {
  print(paste(names(param)[i], param[[i]],  sep = " = "))
}

#res = unlist(lapply(sized_files, function(file) {
res = unlist(mclapply(sized_files, function(file) {
  w = readWave(file)
  songfinder2(w, 
              min_duration = param[["min_duration"]], 
              max_gap = param[["max_gap"]], 
              max_duration= param[["max_duration"]],
              min_num_peaks = param[["min_num_peaks"]], 
              max_num_peaks = param[["max_num_peaks"]], 
              amp_ratio_max_gap = param[["amp_ratio_max_gap"]], 
              amp_ratio_min = param[["amp_ratio_min"]],
              amp_ratio_max = param[["amp_ratio_max"]],
              rate_min = param[["rate_min"]],
              rate_max = param[["rate_max"]],
              low_noise=!opt$high_noise)
}
#))
, mc.cores=opt$cores))

res = data.frame(name=sized_files, song=res) 
res1 = res
res1$song = as.character(res$song)
song_ind = res1$song=="TRUE"

#print(summary(res))
print(paste("Number of identified songs: ", nrow(res[song_ind,]), sep=""))
songs_dir = paste(opt$dir, "songs", sep="/")
dir.create(songs_dir, showWarnings = F)
file.copy(as.character(res[song_ind,1]), songs_dir, copy.date=T, overwrite=F)

if (length(opt$nsongs) > 0) {
  nsongs_dir = paste(opt$dir, "not_songs", sep="/")
  dir.create(nsongs_dir, showWarnings = F)
  file.copy(as.character(res[!song_ind,1]), nsongs_dir, copy.date=T, overwrite=F)
}

if (length(opt$remove) > 0) {
  file.remove(as.character(res[!song_ind,1]))
}

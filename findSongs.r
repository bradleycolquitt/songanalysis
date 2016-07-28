#! /usr/bin/env r

suppressMessages(library(parallel))
suppressMessages(library(docopt))      

doc <- "Usage: findSongs.r [-d directory] [-f no_filter_by_size] [-t date] [-m time] [-c cores] [-n nsongs] [-l low_noise]

-d --dir DIRECTORY  directory to be processed, default = '.'
-f --no_filter_by_size skip filtering files by size (0.5 to 5 MB)
-t --date DATE only process files written by this date
-m --time TIME only process files written by this time
-c --cores CORES number of cores to use, default = 4
-n --nsongs NSONGS copy not_songs to new directory 
-r --remove delete non song files
-l --low_noise BOOLEAN wav recorded using low noise setup (e.g. Presonus). Default is true.
-h --help           show this help text"

opt <- docopt(doc)

options(warn=-1)
suppressMessages(source("~/src/songanalysis/song_util.R"))
suppressMessages(source("~/src/songanalysis/threshold.r"))

#opt = list(dir="/mnt/bengal_home/song/bk47/2016-07-25", cores=10)
#opt = list(dir="/mnt/bengal_home/song/gr44gr48", no_filter_by_size=1, cores=10)

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
if(length(opt$low_noise)==0)
  opt$low_noise = TRUE

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
res = unlist(mclapply(sized_files, function(file) {
  w = readWave(file)
  songfinder2(w, min_duration = 15, max_gap = 10, max_duration=200,
                 min_num_peaks=10, max_num_peaks=NULL, amp_ratio_max_gap=120, low_noise=opt$low_noise)
}
#))
, mc.cores=opt$cores))

res = data.frame(name=sized_files, song=res) 
#print(summary(res))
print(paste("Number of identified songs: ", nrow(res[res[,2],]), sep=""))
songs_dir = paste(opt$dir, "songs", sep="/")
dir.create(songs_dir, showWarnings = F)
file.copy(as.character(res[res[,2],1]), songs_dir, copy.date=T, overwrite=F)

if (length(opt$nsongs) > 0) {
  nsongs_dir = paste(opt$dir, "not_songs", sep="/")
  dir.create(nsongs_dir, showWarnings = F)
  file.copy(as.character(res[!res[,2],1]), nsongs_dir, copy.date=T, overwrite=F)
}

if (length(opt$remove) > 0) {
  file.remove(as.character(res[!res[,2],1]))
}

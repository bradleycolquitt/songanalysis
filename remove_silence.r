#! /usr/bin/env r

suppressMessages(library(parallel))
suppressMessages(library(docopt))      

doc <- "Usage: remove_silence.r [-d directory] [-t date] [-c cores]

-d --dir DIRECTORY  directory to be processed, default = '.'
-t --date DATE only process files written by this date
-c --cores CORES number of cores to use, default = 4

-h --help           show this help text"

opt <- docopt(doc)

options(warn=-1)
suppressMessages(source("~/src/songanalysis/song_util.R"))
suppressMessages(source("~/src/songanalysis/threshold.r"))

opt = list(dir="/mnt/tutor_home/data/or46pu77/recordings/2016-03-02", cores=6)
#opt = list(dir="/mnt/bengal_home/song/gr44gr48", no_filter_by_size=1, cores=10)
if(opt$dir == "NA")
  opt$dir = "."
if(length(opt$thresh) == 0)
  opt$thresh = .4
if(length(opt$cores) == 0)
  opt$cores = 4
if(length(opt$no_filter_by_size) > 0) {
  opt$no_filter_by_size = TRUE
} else {
  opt$no_filter_by_size = FALSE
}

files = list.files(opt$dir, recursive = T, pattern=".wav", full.names = T)
ind = grep("songs", files)
if (length(ind) > 0) files = files[-ind]
print(paste("Number of files to process: ", length(files), sep=""))

#### Calc amp params ####
print("Finding songs...")
#files1 = files[201:400]
#l = profvis({
res = unlist(mclapply(files, function(file) {
  #w = filtersong(readWave(file))
  w = highpass_filter(readWave(file), from=500, wl=1024, ovlp=0)
  rat = calc_max_mean_ratio(w, subsamp=10)
  rat <= 100
  #songfinder2(w, min_duration = 15, max_gap = 10, max_duration=200,
  #               min_num_peaks=5, max_num_peaks=NULL, amp_ratio_max_gap=120)
}
#))
, mc.cores=opt$cores))
#})
res = data.frame(name=files, song=res) 
#print(summary(res))
print(paste("Number of silence files: ", nrow(res[res[,2],]), sep=""))
silence_dir = paste(opt$dir, "silence", sep="/")
dir.create(silence_dir, showWarnings = F)
to_move = as.character(res[res[,2],1])
file.rename(to_move, paste(silence_dir, basename(to_move), sep="/"))

# if (length(opt$nsongs) > 0) {
#   nsongs_dir = paste(opt$dir, "not_songs", sep="/")
#   dir.create(nsongs_dir, showWarnings = F)
#   file.copy(as.character(res[!res[,2],1]), nsongs_dir, copy.date=T, overwrite=F)
# }

# if (length(opt$remove) > 0) {
#   file.remove(as.character(res[!res[,2],1]))
# }

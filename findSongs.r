#! /usr/bin/env r

options(warn=-1)
suppressMessages(source("/home/brad/src/songanalysis/song_util.R"))
suppressMessages(library(parallel))
suppressMessages(library(docopt))      

doc <- "Usage: cleanSongDir.r [-d directory] [-t date] [-c cores] [-n nsongs]

-d --dir DIRECTORY  directory to be processed, default = '.'
-t --date DATE only process files written by this date
-c --cores CORES number of cores to use, default = 4
-n --nsongs NSONGS copy not_songs to new directory 
-h --help           show this help text"

opt <- docopt(doc)
#print(opt)
#print(length(opt$date))
if(opt$dir == "NA")
  opt$dir = "."
if(length(opt$thresh) == 0)
  opt$thresh = .4
if(length(opt$cores) == 0)
  opt$cores = 4

files = list.files(opt$dir, recursive = T, pattern=".wav", full.names = T)
ind = grep("songs", files)
if (length(ind) > 0) files = files[-ind]
print(paste("Number of files to process: ", length(files), sep=""))

#### Filter by size ####

info = file.info(files)
if (length(opt$date) > 0) {
  info = info[as.Date(info$mtime) >= as.Date(opt$date),]
}
info = info[order(-info$size),]
files = rownames(info)
sizes = info$size

sized_range = c(5E5, 5E6)
sized_files = files[sizes>sized_range[1] & sizes<sized_range[2]]
print(paste("Number of files with sizes between ", sized_range[1] / 1E6, " and ", sized_range[2] / 1E6, " Mbytes: ", length(sized_files), sep=""))

#### Define amp thresh ####
print("Defining amplitude threshold...")
num_files_for_amp = 100
max_ind = ifelse(length(sized_files)>num_files_for_amp,  num_files_for_amp, length(sized_files))
amps = unlist(mclapply(sized_files[1:max_ind], function(x) {
  wav = filtersong(readWave(x))
  env(wav, envt="abs", plot=F)
}, mc.cores=opt$cores))
amp_thresh = median(amps) * 10
print(paste("Amplitude threshold: ", amp_thresh, sep=""))
rm(amps)
gc()

#### Find songs ####
print("Finding songs...")
res = unlist(mclapply(sized_files, function(file) {
#res = unlist(lapply(sized_files, function(file) {
  #print(file)
  w = readWave(file)
  #songfinder(w, band = c(5000,7000), wein_thresh = opt$thresh)
  songfinder2(w, min_duration = 15, max_gap = 5, min_num_peaks=5, thresh=amp_thresh)
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
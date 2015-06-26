#! /usr/bin/env r

suppressMessages(source("/media/data2/rstudio/birds/song/song_util.R"))

suppressMessages(library(docopt))      

doc <- "Usage: cleanSongDir.r [-d directory]

-d --dir DIRECTORY  directory to be processed
-t --thresh THRESH Weiner entropy threshold, 0 < t < 1
-h --help           show this help text"

opt <- docopt(doc)
if(opt$dir == "NA")
  opt$dir = "."
if(opt$thresh == "NA")
  opt$thresh = .4

files = list.files(opt$dir, pattern=".wav", full.names = T)

#### Filter by size ####
sized_thresh = 400000
sized_dir = paste(opt$dir, "sized", sep="/")
print(paste("Moving files greater than ", sized_thresh / 1000, " Mbytes to sized", sep=""))
sizes = sapply(files, function(x) file.info(x)$size)
dir.create(sized_dir, showWarnings = F)
sapply(files[sizes>sized_thresh], function(x) file.copy(x, paste(sized_dir, x, sep="/")))

#### Filter by spectral info ####
print("Loading WAVs...")
sized_files = list.files(sized_dir, pattern=".wav", full.names = T)
wav = lapply(sized_files, readWave)
print("Calcuating PSDs...")
psd = lapply(wav, function(x) spec(filtersong(x), wl=256, PSD=T, plot=F))
print("Finding songs...")
res = sapply(psd, function(x) {
  songfinder_psd(x, band=c(5,7), sum_thresh=0, var_thresh=0, wein_thresh=.4)
})
res = data.frame(name=sized_files, song=res) 
songs_dir = paste(opt$dir, "songs", sep="/")
dir.create(songs_dir, showWarnings = F)

#file.copy(as.character(res$name[res[,2]]), songs_dir)
for (ind in which(res[,2])) {
  print(as.character(res[ind,1]))
  file.rename(as.character(res[ind,1]), paste(songs_dir, basename(as.character(res[ind,1])), sep="/"))
  #file.rename(as.character(res[ind,1]), songs_dir)
}

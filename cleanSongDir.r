#! /usr/bin/env r

suppressMessages(source("/media/data2/rstudio/birds/song/song_util.R"))
suppressMessages(library(parallel))
suppressMessages(library(docopt))      

doc <- "Usage: cleanSongDir.r [-d directory] [-t thresh] [-c cores]

-d --dir DIRECTORY  directory to be processed, default = '.'
-t --thresh THRESH Weiner entropy threshold, 0 < t < 1, default = 0.4
-c --cores CORES number of cores to use, default = 4
-h --help           show this help text"

opt <- docopt(doc)
if(opt$dir == "NA")
  opt$dir = "."
if(length(opt$thresh) == 0)
  opt$thresh = .4
if(length(opt$cores) == 0)
  opt$cores = 4

files = list.files(opt$dir, pattern=".wav", full.names = T, )

#### Filter by size ####
sized_thresh = 400000
sized_dir = paste(opt$dir, "sized", sep="/")
print(paste("Moving files greater than ", sized_thresh / 1000, " Mbytes to sized", sep=""))
sizes = sapply(files, function(x) file.info(x)$size)
dir.create(sized_dir, showWarnings = F)
file.copy(files[sizes>sized_thresh], sized_dir, copy.date=T)

#### Filter by spectral info ####
print("Loading WAVs...")
sized_files = list.files(sized_dir, pattern=".wav", full.names = T)
wav = mclapply(sized_files, readWave, mc.cores=opt$cores)

print("Finding songs...")
res = unlist(mclapply(wav, function(w) {
  songfinder(w, band = c(5,7), wein_thresh = opt$thresh)
}, mc.cores=opt$cores))
res = data.frame(name=sized_files, song=res) 

songs_dir = paste(opt$dir, "songs", sep="/")
dir.create(songs_dir, showWarnings = F)
file.copy(as.character(res[,1]), songs_dir, copy.date=T)

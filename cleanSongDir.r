#! /usr/bin/env r

options(warn=-1)
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
print(paste("Number of files to process: ", length(files), sep=""))

#### Filter by size ####
sized_range = c(2E6, 10E6)
sized_dir = paste(opt$dir, "sized", sep="/")
#print(paste("Copying files between ", sized_range[1] / 1000, " and ", sized_range[2] / 1000, " Mbytes to sized", sep=""))
sizes = sapply(files, function(x) file.info(x)$size)
#dir.create(sized_dir, showWarnings = F)
#file.copy(files[sizes>sized_thresh], sized_dir, copy.date=T)

#### Filter by spectral info ####
sized_files = files[sizes>sized_range[1] & sizes<sized_range[2]]
print(paste("Number of files with sizes between ", sized_range[1] / 1000, " and ", sized_range[2] / 1000, " Mbytes: ", length(sized_files), sep=""))
print("Loading WAVs...")
#sized_files = list.files(sized_dir, pattern=".wav", full.names = T)


#wav = mclapply(sized_files, readWave, mc.cores=opt$cores)

print("Finding songs...")
res = unlist(mclapply(sized_files, function(file) {
  w = readWave(file)
  #songfinder(w, band = c(5,7), wein_thresh = opt$thresh)
  songfinder2(w)
}, mc.cores=opt$cores))
res = data.frame(name=sized_files, song=res) 

songs_dir = paste(opt$dir, "songs2", sep="/")
dir.create(songs_dir, showWarnings = F)
res = file.copy(as.character(res[res[,2],1]), songs_dir, copy.date=T)

if (prod(res)) unlink(sized_dir, recursive=TRUE)

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

#opt = list(dir="/mnt/tutor_home/data/bk69wh82/recordings/2016-04-23", cores=7)
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

files = list.files(opt$dir, recursive = F, pattern=".wav", full.names = T)
#ind = grep("songs", files)
#if (length(ind) > 0) files = files[-ind]
print(paste("Number of files to process: ", length(files), sep=""))

#### Get size info ####
print("Getting sizes...")
info = file.info(files)
info = info[order(-info$size),]
info$fname = rownames(info)
#files = rownames(info)
#sizes = info$size

min_size = 1000
sized_files = files[sizes>sized_range[1] & sizes<sized_range[2]]
#print(paste("Number of files with sizes between ", sized_range[1] / 1E6, " and ", sized_range[2] / 1E6, " Mbytes: ", length(sized_files), sep=""))
#### Calc amp params ####
print("Finding the silence...")
#files1 = files[10610:10614]
#l = profvis({
info1 =info[1:10,]
res = unlist(foreach(row=isplitRows(info, chunkSize=1)) %dopar% {
#res = unlist(apply(info[,c("info", "size")], 1, function(row) {
  #w = filtersong(readWave(file))
  file = as.character(row["fname"])
  size = as.numeric(row["size"])
  print(file)
  if (size < min_size)
    return(TRUE)
  
  w = highpass_filter(readWave(file), from=500, wl=1024, ovlp=0)
  if (is.null(w))
    return(TRUE)
  rat = calc_max_mean_ratio(w, subsamp=10)
  rat <= 100
})
#))
#, mc.cores=opt$cores))
#})
res = data.frame(name=files, song=res) 
#print(summary(res))
print(paste("Number of silence files: ", nrow(res[res[,2],]), sep=""))
silence_dir = paste(opt$dir, "silence", sep="/")
print(silence_dir)
dir.create(silence_dir, showWarnings = T)
to_move = as.character(res[res[,2],1])
if (length(to_move) == nrow(info)) {
  stop("Error in process.")
}

copy_res = file.copy(to_move, silence_dir, copy.date=T, overwrite = T)

if (all(copy_res)) {
  file.remove(as.character(to_move))
  print("Files moved successfully.", quote = F)
} else {
  print("Error moving files", quote=F)
}
# if (length(opt$nsongs) > 0) {
#   nsongs_dir = paste(opt$dir, "not_songs", sep="/")
#   dir.create(nsongs_dir, showWarnings = F)
#   file.copy(as.character(res[!res[,2],1]), nsongs_dir, copy.date=T, overwrite=F)
# }

# if (length(opt$remove) > 0) {
#   file.remove(as.character(res[!res[,2],1]))
# }

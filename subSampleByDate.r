#! /usr/bin/env r

options(warn=-1)
suppressMessages(source("/media/data2/rstudio/birds/song/song_util.R"))
suppressMessages(library(parallel))
suppressMessages(library(docopt))      
suppressMessages(library(dplyr))

doc <- "Usage: subSampleByDate.r [-d directory] [-n nsongs] [-t date]

Subsample wav files in specified directory. Selects 'nsongs' per date. 

-d --dir DIRECTORY  directory to be processed, default = '.'
-n --nsongs NSONGS number of songs to subsample per day
-t --date DATE only 
-h --help           show this help text"

opt <- docopt(doc, strict=F)

suppressMessages(source("~/src/songanalysis/song_util.R"))
#print(opt)
if(opt$dir == "NA")
  opt$dir = "."
if(length(opt$nsongs) == 0 )
  opt$nsongs = 20

opt$nsongs = as.numeric(opt$nsongs)

files = list.files(opt$dir, recursive = T, pattern=".wav$", full.names = T)
ind = grep("song", files)
if (length(ind) > 0) files = files[-ind]

info = file.info(files)
# if (length(opt$date) > 0) {
#   info = info[as.Date(info$mtime) >= as.Date(opt$date),]
# }

file_times = parse_timestamp(files)
date = NULL
print(dim(info))
if (length(opt$date) > 0) {
  date = opt$date
  date = as.POSIXct(date)
  print(paste("Filtering for recordings after:", date, sep=" "))

  info = info[file_times > date,]
  print(dim(info))
}

info$fname = rownames(info)
info$date = as.character(strftime(info$mtime, '%Y-%m-%d'))

song_by_date = table(info$date)

set.seed(1)

info1 = info %>% dplyr::filter(date %in% names(song_by_date)) %>% dplyr::group_by(date) %>% do({
   if (nrow(.) < opt$nsongs) {
     return(.)
   } else {
     sample_n(. , opt$nsongs)
   }
})

select_dir = paste(opt$dir, "select", sep="/")
dir.create(select_dir, showWarnings = T)
file.copy(info1$fname, select_dir, copy.date=T)

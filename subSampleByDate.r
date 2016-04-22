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
#print(opt)
if(opt$dir == "NA")
  opt$dir = "."
if(length(opt$nsongs) == 0 )
  opt$nsongs = 20

opt$nsongs = as.numeric(opt$nsongs)

files = list.files(opt$dir, recursive = T, pattern=".wav$", full.names = T)
ind = grep("songs", files)
if (length(ind) > 0) files = files[-ind]

info = file.info(files)
if (length(opt$date) > 0) {
  info = info[as.Date(info$mtime) >= as.Date(opt$date),]
}

info$fname = rownames(info)
info$date = strftime(info$mtime, '%Y-%m-%d')

song_by_date = table(info$date)
#if (song_by_date > opt$nsongs)
#song_by_date = song_by_date[song_by_date>opt$nsongs]

set.seed(1)

info1 = info %>% filter(info$date %in% names(song_by_date)) %>% group_by(date) %>% do({
  if (nrow(.) < opt$nsongs) {
    return(.)
  } else {
    sample_n(. , opt$nsongs)
  }
})

select_dir = paste(opt$dir, "select", sep="/")
dir.create(select_dir, showWarnings = F)
file.copy(info1$fname, select_dir, copy.date=T)

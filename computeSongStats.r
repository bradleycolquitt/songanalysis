#! /usr/bin/env r

###################################################################
# Compute basics song statistics from directory of wav files      #   
# Stats will be inserted into specified MySQL database            #
###################################################################


suppressMessages(source("~/src/songanalysis/song_util.R"))
suppressMessages(library(parallel))
suppressMessages(library(docopt))  
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(DBI))

doc <- "Usage: computeSongStats.r [-d directory] [-o db] [-i bird_info] [-c cores]

Compute basics song statistics from directory of wav files
Stats will be inserted into specified MySQL database

-d --dir DIRECTORY  directory to be processed, default = '.'
-i --bird_info INFO csv file containing information about given bird
-o --database DATABASE database to write to
-c --cores CORES number of cores to use, default = 4
-h --help           show this help text"

opt <- docopt(doc)
if(opt$dir == "NA")
  opt$dir = "."
if(length(opt$cores) == 0)
  opt$cores = 4
if(length(opt$bird_info) == 0) {
  opt$bird_info = paste(opt$dir, "bird_info.csv", sep="/")
}
opt$dir = normalizePath(opt$dir)

#### Read in data, establish DB connection ####
conn = dbConnect(RMySQL::MySQL(), group="songdb")
bird_df = read.csv(opt$bird_info, header=T)

#### Check presence of data ####
dirs = str_split(opt$dir, "/")[[1]]
tags = dirs[length(dirs)-1]
bird = bird_df[bird_df$tags==tags,]
res = dbSendQuery(conn, paste("SELECT count(*) FROM birds WHERE tags = '", tags, "'", sep=""))
res1 = dbFetch(res)
dbClearResult(res)

if (res1 > 0) {   
  print(paste("Song data is present for ", bird[2], ". Overwrite[y/n]?", sep=""))
  x = readLines(n=1)
  if (x == "n") {
      stop("Exiting")
  } else if (x == "y") {
      print("Overwriting")
      dbSendQuery(conn, paste("DELETE FROM birds WHERE birds.tags = '", tags, "'", sep=""))
  }
}
dbWriteTable(conn, "birds", bird, append=TRUE, row.names=FALSE)

#setwd("/media/data2/song/bl63bl44/songs")
#opt = list(dir=".", thresh=.4, cores=10, bird_info="/media/data2/song/bl63bl44/bird_info.csv")

#### Process wav files ####
files = list.files(opt$dir, pattern=".wav", full.names = T)
file_info = file.info(files)
file_info$date = strftime(file_info$mtime, '%Y-%m-%d')
file_info$time = strftime(file_info$mtime, '%H:%M:%S')
results = mclapply(1:length(files), function(ind) {
  result = list() 
  wav = readWave(files[[ind]])
  
  ### Song-level analysis ###
  print(files[[ind]])
  result$songs = findpeaks(wav, min_duration=1000, max_gap=300, thresh=0.05)
  result$songs$song_number = ind
  
  ### Motif level analysis ###
  result$motifs = findpeaks(wav, min_duration=200, max_gap=50, thresh=0.05)
  if (is.null(result$motifs)) return(result)
  result$motifs$song_number = ind
  result$motifs$motif_number = 1:nrow(result$motifs)
  colnames(result$motifs)[1:2] = c("motif_start", "motif_end")
  
  #### Syllable level analysis####
  result$syllables = findpeaks(wav, min_duration=20, max_gap=10, thresh=0.05)
  if (is.null(result$syllables)) return(result)
  
  result$syllables$midpoint = rowSums(result$syllables[,1:2]) / 2
  motif_breaks = c(0, interpeak_midpoint(result$motifs), length(wav@left) / wav@samp.rate)
  result$syllables$motif_number = cut(result$syllables$midpoint, 
                        breaks=motif_breaks, 
                        include.lowest=F, 
                        labels=F)
  result$syllables$song_number = ind
  result$syllables = result$syllables %>% group_by(motif_number) %>% mutate(syllable_number=1:n())
  
  #### Weiner entropies ####
  result$syllables$ent = apply(result$syllables, 1, function(syl) {
    weiner_entropy(wav, band=c(5,10), subregion=syl)
  })
  return(result)
}
#)
,mc.cores=opt$cores)
names(results) = files

#### Songs ####
print("Inserting songs...")
songs_timings = lapply(results, function(result) {
                      result[['songs']]
                    })
songs_df = do.call("rbind", songs_timings)
file_info = cbind(file_info, songs_df)
file_info$bird_id = bird$bird_id

song_formatted = file_info[,c("bird_id", "mtime", "ins", "outs", "song_number")]
#song_formatted = file_info[,c("bird_id", "date", "time", "ins", "outs", "song_number")]
colnames(song_formatted) = c("bird_id", "song_datetime", "song_start", "song_end", "song_number")
#colnames(song_formatted) = c("bird_id", "song_date", "song_time", "song_start", "song_end", "song_number")
dbWriteTable(conn, "songs", song_formatted, append=TRUE, row.names=F)

#### Motifs ####
print("Inserting motifs...")
motifs = lapply(results, function(result) {
  song_ids = apply(result[['motifs']], 1, function(motif) {
    res = dbSendQuery(conn, paste("SELECT song_id FROM songs WHERE song_number = ", motif['song_number'], 
                            " AND bird_id = ", bird$bird_id, sep=""))
    df = dbFetch(res)
    dbClearResult(res)
    return(df)
  })
  result[['motifs']] = data.frame(result[['motifs']], song_id=unlist(song_ids))
  result[['motifs']]
})
motifs_df = do.call("rbind", motifs)
dbWriteTable(conn, 'motifs', motifs_df[,c("song_id", "motif_number", "motif_start", "motif_end")], append=T, row.names=F)

#### Syllables ####
print("Inserting syllables...")
syllables = lapply(results, function(result) {
   motif_ids = apply(result[['syllables']], 1, function(syllables) {
    res = dbSendQuery(conn, paste("SELECT motif_id FROM songs JOIN motifs ON motifs.song_id=songs.song_id WHERE",
                                  " bird_id = ", bird$bird_id,
                                  " AND songs.song_number = ", syllables['song_number'],
                                  " AND motifs.motif_number = ", syllables['motif_number'], 
                                  sep=""))
    df = dbFetch(res)
    dbClearResult(res)
    return(df)
  })
  result[['syllables']] = data.frame(result[['syllables']], motif_id=unlist(motif_ids))
  result[['syllables']]
})
syllables_df = do.call("rbind", syllables)
syllables_formatted = syllables_df[,c("motif_id", "syllable_number", "ins", "outs", "ent")]
to_select = colnames(syllables_formatted) = c("motif_id", "syllable_number", "syllable_start", "syllable_end", "syllable_weinent")
dbWriteTable(conn, 'syllables', syllables_formatted[,to_select], append=T, row.names=F)


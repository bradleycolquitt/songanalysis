#! /usr/bin/env r

options(warn=-1)
suppressMessages(library(docopt))
doc <- "Usage: findSubSongs.r [-d directory] [-f no_filter_by_size] [-t date] [-c cores] [-n nsongs]

-d --dir DIRECTORY  directory to be processed, default = '.'
-f --no_filter_by_size skip filtering files by size (0.5 to 5 MB)
-t --date DATE only process files written by this date
-c --cores CORES number of cores to use, default = 4
-n --nsongs NSONGS copy not_songs to new directory 
-h --help           show this help text"

opt <- docopt(doc)

suppressMessages(source("~/src/songanalysis/song_util.R"))
suppressMessages(source("~/src/songanalysis/song_features.R"))
suppressMessages(library(parallel))
 
#suppressMessages(library(caret))
#suppressMessages(library(gbm))
#suppressMessages(library(randomForest))


#print(opt)
#print(length(opt$date))
#opt = list(dir="/mnt/tutor_home/data/bk69wh82/recordings/2016-04-21", cores=4)
if(opt$dir == "NA")
  opt$dir = "."
if(length(opt$cores) == 0)
  opt$cores = 4
if(length(opt$no_filter_by_size) > 0) {
  opt$no_filter_by_size = TRUE
} else {
  opt$no_filter_by_size = FALSE
}
  
doMC::registerDoMC(cores=opt$cores)

print("Getting file list...")
files = list.files(opt$dir, recursive = T, pattern=".wav", full.names = T)
ind = grep("songs|sized", files)
if (length(ind) > 0) files = files[-ind]
dirs = unique(unlist(lapply(str_split(files, "/"), function(x) x[1])))
print(dirs)
print(paste("Number of files to process: ", length(files), sep=""))

#### Filter by size ####  
if (!opt$no_filter_by_size) {
  print("Filtering by size...")
  info = file.info(files)
  if (length(opt$date) > 0) {
    info = info[as.Date(info$mtime) >= as.Date(opt$date),]
  }
  info = info[order(-info$size),]
  files = rownames(info)
  sizes = info$size
  
  sized_range = c(3E5, 2E6)
  sized_files = files[sizes>sized_range[1] & sizes<sized_range[2]]
  print(paste("Number of files with sizes between ", sized_range[1] / 1E6, " and ", sized_range[2] / 1E6, " Mbytes: ", length(sized_files), sep=""))
} else {
  sized_files = files
}
#### Find songs ####
print("Finding songs...")
chunkSize = 1000
nchunks = ceil(length(sized_files) / chunkSize)
i = 1

## Params
min_duration = 30
max_duration = 500
max_gap = 10
min_num_peaks = 5
min_mean_freq = 1.5 # kHz
max_q95_freq = 8 # kHz
min_total_duration = 0.5
max_amp_ratio = .5
subsamp = 10
max_mean_gap_durations = 0.6
max_cv_gap_durations = 1

songs_dir = paste(opt$dir, "subsongs", sep="/")
dir.create(songs_dir, showWarnings = T)

total_songs = 0

#sized_files = sized_files[1:10]
foreach(chunk=isplitVector(sized_files, chunkSize=chunkSize)) %do% {
  cat(paste("chunk ", i , " of ", nchunks, ": ", sep=""))
  i = i +1
  res = foreach(fname=chunk) %dopar% {
    wav = readWave(fname)
    wavf = highpass_filter(wav, from = 500, wl = 1024, ovlp = 25)
    peak_info = findpeaks_range(wavf, 
                                min_duration = min_duration, 
                                max_gap = max_gap, 
                                max_duration = max_duration, 
                                absolute = T,
                                thresh_range=seq(1, 10, .5))
    peaks = peak_info$peaks
    #    print(nrow(peaks))
    if (is.null(peaks) || nrow(peaks)<=2)
      return(FALSE)
    wl = 1024
    frate = wav@samp.rate
    window = wl / frate
    peaks = peaks %>% dplyr::filter(ins>window, outs<(seewave::duration(wav)-window))
    #    print(nrow(peaks))
    if (nrow(peaks)<=2)
      return(FALSE)
    pdists = peak_dist(peaks)
    if(mean(pdists) > max_mean_gap_durations)
      return(FALSE)
    if (sd(pdists)/mean(pdists) > max_cv_gap_durations)
      return(FALSE)
    
    ### Filter large gaps
    peaks = filter_by_gaps(peaks, 0.5)
    if (is.null(peaks))
      return(FALSE)
    num_peaks = nrow(peaks)
    #    print(num_peaks)
    if (num_peaks < min_num_peaks) 
      return(FALSE)
    
    total_duration = peaks[num_peaks,2] - peaks[1,1]
    #  print(total_duration)
    if (total_duration < min_total_duration) 
      return(FALSE)
    
    #pdists = peak_dist(peaks)
    #if (sd(pdists) / )
    amp = abs(wavf@left[seq(1,length(wavf@left), subsamp)])
    amp_ratio = calc_amp_ratio(amp, frate, peaks, max_gap=500, subsamp=subsamp)
    if (amp_ratio>max_amp_ratio) 
      return(FALSE)
    
    syl_ents = apply(peaks[,1:2], 1, function(peak) {
      #calc_freq_stats(wavf, q =c(.5,.95), subregion = c(peak[1], peak[2]), wl=256, overlap=0 )
     
     calc_mean_freq(wavf, subregion = c(peak[1], peak[2]), wl=256)
     #calc_spectral_features(wavf, subregion = c(peak[1], peak[2]), overlap=0, wl=wl, freq)
    })
   # syl_ents = do.call("rbind", syl_ents)
    mean_freq = mean(syl_ents, na.rm=T)
    #freq_stats = apply(syl_ents, 1, mean, na.rm=T)
    #if ((freq_stats[1]<min_mean_freq) || (freq_stats[2]>max_q95_freq))
    if (mean_freq < min_mean_freq)
      return(FALSE)
    
    return(TRUE)
  }
    
  res = data.frame(fname=chunk, isSong=unlist(res))
  if (is.null(res)) {
    print("No songs found.")
    #stop()
  }
  #print(summary(res[,2]))
  print(paste("Number of identified songs: ", nrow(res[res[,2],]), sep=""))
  total_songs = total_songs + sum(res[,2])
  
  file.copy(as.character(res[res[,2],1]), songs_dir, copy.date=T, overwrite=F)
  
  if (length(opt$nsongs) > 0) {
    nsongs_dir = paste(opt$dir, "not_subsongs", sep="/")
    dir.create(nsongs_dir, showWarnings = F)
    file.copy(as.character(res[!res[,2],1]), nsongs_dir, copy.date=T, overwrite=F)}
  
}
print(paste("Total songs: ", total_songs, " of ", length(sized_files), sep=""))


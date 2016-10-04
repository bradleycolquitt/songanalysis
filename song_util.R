library(stats)
library(seewave)
library(tuneR)
library(signal)
library(reshape2)
library(ggplot2)
library(doMC) 
library(foreach)
library(itertools)
library(R.matlab)
library(stringr)
library(plyr)
library(dplyr)
library(gridExtra)
library(parallel)
library(matlab)
library(lubridate)
library(accelerometry)
#library(smoother)
#registerDoMC(cores=8)

source("~/src/songanalysis/threshold.r")
#source("~/src/songanalysis/clustering.R")
source("~/src/songanalysis/song_features.R")

plot_2dspec = function(wav, peaks=NULL) {
  theme_set(theme_classic())
  manipulate({
    gg = ggspectro(wav, wl=256, ovlp=75) + geom_tile(aes(fill=amplitude)) 
    gg = gg + scale_fill_gradientn(colours=spectro.colors(30), limits=c(-90, 0), na.value="transparent")
    gg = gg + xlim(xmin, xmax) + ylim(2, 15)
    gg
  }, 
  xmax=slider(0, length(wav) / wav@samp.rate, initial=length(wav)/wav@samp.rate),
  xmin=slider(0, length(wav) / wav@samp.rate, initial=0))
}

##### FILE UTIL #####
parse_timestamp = function(fname, prefix="output_") {
  times = str_replace(basename(fname), prefix, "")
  times = str_replace(times, ".wav", "")
  times = as.numeric(times)
  times = as.POSIXct(times,  origin="1970-01-01")
  return(times)
}

######## FILTERING/EDITING ########

#' taken from seewave::ffilter 
highpass_filter = function(wav, from=1800, wl=1024, ovlp=75, wn="hanning", fftw=T) {
  input = inputw(wave = wav, f = f)
  wave = input$w
  f = input$f
  bit = input$bit
  #rm(input)
  n = nrow(wave)
  if (n < wl) 
    return(NULL)
  step = seq(1, n - wl, wl - (ovlp * wl/100))
  z = stft(wave = wave, f = f, wl = wl, zp = 0, step = step, 
           wn = wn, fftw = fftw, complex = TRUE, scale = FALSE)
  
  to = f/2
  F = ceiling(wl * (from/f))
  T = floor(wl * (to/f))
  z[-c(F:T), ] = 0
  
  res = istft(z, wl = wl, ovlp = ovlp, wn = wn, output = "Wave", f = f)
  return(res)
}
#' High pass filter Wave object
#' @param wav, the input Wave object
filtersong = function(wav, band=500, span=50) {
  if (length(band) > 1) {
    for (i in band) {
      ffilter(wav, from=(i-span), to=(i+span), fftw=T, output="Wave")
    }
  } else {
    wav = ffilter(wav, from=band, output="Wave", fftw=T)
  }
  return(wav)
}

filtersong_signal = function(wav, band=500) {
  require(signal)
  w = band / (wav@samp.rate / 2)
  f = fir1(n=256, type = "high", w = w)
  out = filtfilt(f, wav@left)
  wav2 = wav
  wav2@left = out
}

region_on_edge = function(wav_duration, fs, wl, start, stop) {
  buffer = wl / (2*fs)
  
  min_val = buffer
  max_val = wav_duration - buffer
  return(start<=min_val || stop>=max_val)
  
}
#' Moving average smooth Wave object
#' @param wav, the input Wave object
#' @param window, the smoothing window in milliseconds
#' @return Wave object containing smoothed data in left slot 
smoothsong = function(wav, window=0.2) {
  window_length = round(wav@samp.rate * window / 1000)
  data = wav@left
  filt = rep(1, times=window_length) / window_length
  data1 = convolve(data, filt, type="filter")
  return(Wave(data1, samp.rate=wav@samp.rate, bit=wav@bit))
}

#' Split wav into chunks
#' @param wav, the input Wave object
#' @param chunk_size, the size of each chunk in seconds
#' @return list of Wave chunks
chunk_wav = function(wav, chunk_size=2) {
  dur = duration(wav)
  starts = seq(0, dur, chunk_size)
  ends = starts + chunk_size
  ends[length(ends)] = dur
  
  wavs = lapply(1:length(starts), function(i) {
    seewave::cutw(wav, from = starts[i], to = ends[i], plot = F, output = "Wave")
  })
  return(wavs)
}

remove_peaks_at_edge = function(peaks, fs, wl, wav_duration) {
  # Filter out peaks on edge of recording
  buffer = wl / fs
  
  ind = peaks[,1]< buffer | peaks[,2] > (wav_duration - buffer)
  peaks = peaks[!ind,]
  return(peaks)
}

######## PROCESSING ########

parse_song = function(file, 
                      fast_filter=F,
                      cluster=F, 
                      plot=F, 
                      log_transform=F, 
                      absolute=F,
                      thresh_method=mean_sd, 
                      thresh_factor=0.5, 
                      thresh_range=seq(0,10,.25),
                      min_dur=15, 
                      max_gap=10, 
                      max_duration=300,
                      subsamp=1,
                      include_wavs = TRUE,
                      peak_source = "new") {
  wav = NULL
  fname = NULL
  fname_base = NULL
  
  if (class(file) == "character") {
    wav = readWave(file)
    #fname = basename(str_replace(file, ".wav", ""))
    fname = file
    fname_base = basename(str_replace(file, ".wav", ""))
  } else if (class(file) == "Wave") {
    wav = file
    fname = "dummy"
    fname_base = "dummy"
  }
  if (fast_filter) {
    wavf = highpass_filter(wav, from=600, ovlp=50, wl=2048)
  } else {
    wavf = filtersong(wav, band=600)
  }
  
  
  #thresh = threshold_auto(wavf, otsu)
  syls = NULL
  thresh = NULL
  if (peak_source %in% c("mat", "automat") && file.exists(paste(file, ".not.mat", sep=""))) {
    syls = load_mat(paste(file, ".not.mat", sep=""))
    syls = as.data.frame(syls %>% rename(ins=onsets, outs=offsets))
  } else if (peak_source == "new"){
    if (thresh_method=="range") {
      peak_info = findpeaks_range(wavf, 
                                  min_duration = min_dur, 
                                  max_duration = max_duration, 
                                  max_gap = max_gap, 
                                  absolute=absolute,
                                  thresh_range=thresh_range,
                                  subsamp=subsamp,
                                  floor=F)
      syls = peak_info$peaks
      thresh = peak_info$thresh
    } else {
      thresh = threshold_auto(wavf, get(thresh_method), log=log_transform, factor=thresh_factor)
      syls = findpeaks_abs(wavf, min_duration = min_dur, max_duration = max_duration, max_gap = max_gap, thresh=thresh, subsamp=subsamp)
    }
  }

  if (is.null(syls))  return(list(fname = fname, 
                                  wav=wavf, 
                                  syllable=NULL, 
                                  motifs=NULL,
                                  thresh=thresh, 
                                  min_dur = min_dur, 
                                  max_gap = max_gap))
  
  if (nrow(syls) < 2) return(list(fname = fname, 
                                  wav=wavf, 
                                  syllable=NULL, 
                                  motifs=NULL,
                                  thresh=thresh, 
                                  min_dur = min_dur, 
                                  max_gap = max_gap))
  
  ### Motif level analysis ###
  motifs = find_motifs(syls, max_gap=1000)
  
  if (plot) {
    par(mfrow=c(1,1))
    dur = seewave::duration(wavf)
    if (absolute) {
      env_data = env(wavf, envt="abs", from=0, to=ifelse(dur<4, dur, 4), plot=T)
      
    } else {
      max_plot = 8
      env_data = wavf@left^2
      xmin = 1
      xmax = ifelse(dur<max_plot, dur * wavf@samp.rate, max_plot*wavf@samp.rate) 
      xind = seq(0, xmax / wavf@samp.rate, 1 / wavf@samp.rate)
      plot(xind, env_data[xmin:(xmax+1)], type="l", 
           ylab="Amplitude power", xlab="seconds",
           ylim=c(-1 * max(env_data)/10, max(env_data)/2))
      plot_peaks(syls, col=4, yval = -1 * max(env_data) / 20)
    }
    abline(h=thresh, col=2)
  }
  
  syls$id = paste(fname_base, 1:nrow(syls), sep="-")
  if (include_wavs) {
    song = list(fname = fname, wav=wavf, syllable=syls, motifs=motifs, thresh=thresh, min_dur = min_dur, max_gap = max_gap)
  } else {
    song = list(fname = fname, syllable=syls, motifs=motifs, thresh=thresh, min_dur = min_dur, max_gap = max_gap)
  }
  if (cluster) {
    syl_data = extract_features2(song)
    mod = mclust_syllable_data(syl_data, plot=T)
    song$syllable$called = mod$classification
    return(list(song=song, model=mod))
  }
  return(song)
}

parse_song_batch2 = function(input, num_of_songs=NULL, peak_source="new", thresh_range = seq(0,10,.25)) {
  #wdir = input
  cat("   Loading data\n")
  if (length(input) == 1) {
    ## Input is character vector of working dir
    wdir = input
    files = list.files(wdir, pattern=".wav$", full.names = T)
  } else {
    ## Input is vector of file names
    files = input
  }
  if (!is.null(num_of_songs) && length(files)>num_of_songs) {
    files = files[sample(1:length(files), num_of_songs, replace = F)]
  }
  
  #songd = foreach(i=1:length(files)) %dopar% {
  songd = mclapply(1:length(files), function(i) {
  #songd = lapply(1:length(files), function(i) {
    file = files[i]
    #print(i)
    song = parse_song(file, 
                      cluster = F, 
                      plot = F, 
                      log_transform = F, 
                      absolute=F, 
                      thresh_method="range", 
                      thresh_range=thresh_range,
                      peak_source=peak_source)
    return(song)
  #}
  }, mc.cores=10)
  #})
  return(list(songs=songd))
}

process_song_batch = function(songs,
                              models = NULL,
                              feature_set="mid3_fixed",
                              wl=512,
                              smoothed=TRUE,
                              threshed=FALSE,
                              weighted=FALSE,
                              freq_limits=c(0,12),
                              range_to_test=c(4:9), 
                              distmethod="cor", 
                              power=2, 
                              sample_size=120,
                              cluster=TRUE, 
                              return_data=FALSE) 
{
  songd = NULL
  if (length(songs) == 1) {
    ### Expects input to be list of length 1 containing list of songs
    songs1 = songs[[1]]
    songd = foreach(i=1:length(songs1)) %dopar% {
      song = songs1[[i]]
      #print(i)
      if (is.null(song$syllable)) return(list(song=song, syl_data=NULL))
      syl_data = extract_features2(song, wl=wl, feature_set=feature_set, smoothed=smoothed, threshed=threshed, weighted=weighted, freq_limits=freq_limits)
      return(list(song=song, syl_data=syl_data))
    }
    fnames = unlist(lapply(songs1, function(x) x$fname))
    names(songd) = fnames
  } else if (length(songs) == 2) {
    ### Expects list returned by parse_song_batch:cluster=TRUE
    songd = songs$songs
  } 
  cat("   Rearranging data...\n")
  # remove non-song/no peak files
  ind = unlist(lapply(songd, function(x) is.null(x$syl_data)))
  if (!is.null(ind))
    songd = songd[!ind]
  wavs = lapply(songd, function(x) x$song$wav)
  
  syl_data = lapply(songd, function(x) x$syl_data)
  syl_data = do.call("rbind", syl_data)
  
  syls = lapply(songd, function(song) song$song$syllable)
  syls = do.call("rbind", syls)
  
  
  if (return_data) {
    songd2 = list(wavs=wavs, syl_data=syl_data, syls=syls)
    return(songd2)
  }
  
  
  if (cluster) {
    print("Clustering syllables...")
    syl.mc = mclust_syllable_data(syls, syl_data, models=models, sample_size=sample_size, distmethod=distmethod, power=power, range_to_test = range_to_test, plot=T)
    if (is.null(syl.mc[[1]])) 
      return(NULL)
    
    to_out = list(models=syl.mc$mc, ref_mat=syl.mc$ref_mat, inds=syl.mc$inds, train_data=syl.mc$train_data)
    syls = syl.mc$syls
    syls$wav = paste(wdir, str_replace(syls$id, "-[0-9]+$", ".wav"), sep="/")
    syls.list = split(syls, syls$wav)
    
    for(i in names(songd)) {
      songd[[i]]$song$syllable = syls.list[[i]]
    }
    return(list(songs=songd, models=to_out))
  } else {
    return(syl_data)
  }
}  
process_song_batch_chunk = function(songs,
                                    model = NULL,
                                    feature_set="mid3_fixed",
                                    wl=512,
                                    smoothed=TRUE,
                                    threshed=FALSE,
                                    weighted=FALSE,
                                    range_to_test=c(4:9), 
                                    distmethod="cor", 
                                    power=2, 
                                    sample_size=120,
                                    cluster=TRUE, 
                                    return_data=FALSE) 
{
  songd = NULL
  batch_size = 500
  if (length(songs) == 1) {
    ### Expects input to be list of length 1 containing list of songs
    songs1 = songs[[1]]
    songd = foreach(i=1:length(songs1)) %dopar% {
      # songd = mclapply(1:length(files), function(i) {
      song = songs1[[i]]
      #print(i)
      #song = parse_song(file, cluster = F, plot = F, log_transform = T, thresh_method="range")
      if (is.null(song$syllable)) return(list(song=song, syl_data=NULL))
      syl_data = extract_features2(song, wl=wl, feature_set=feature_set, smoothed=smoothed, threshed=threshed, weighted=weighted)
      return(list(song=song, syl_data=syl_data))
    }
    #}, mc.cores=6)  
    
    #return(songd)
    fnames = unlist(lapply(songs1, function(x) x$fname))
    names(songd) = fnames
  } else if (length(songs) == 2) {
    ### Expects list returned by parse_song_batch:cluster=TRUE
    songd = songs$songs
  } 
  print("Rearranging data...")
  
  # remove non-song/no peak files
  ind = unlist(lapply(songd, function(x) is.null(x$syl_data)))
  songd = songd[!ind]
  wavs = lapply(songd, function(x) x$song$wav)
  
  syl_data = lapply(songd, function(x) x$syl_data)
  syl_data = do.call("rbind", syl_data)
  
  syls = lapply(songd, function(song) song$song$syllable)
  syls = do.call("rbind", syls)
  
  songd2 = list(wavs=wavs, syl_data=syl_data, syls=syls)
  if (return_data) 
    return(songd2)
  
  if (cluster) {
    print("Clustering syllables...")
    syl.mc = mclust_syllable_data(syls, syl_data, model=model, sample_size=sample_size, distmethod=distmethod, power=power, range_to_test = range_to_test, plot=T)
    
    if (is.null(syl.mc[[1]])) 
      return(NULL)
    syls = syl.mc$syls
    syls$wav = paste(wdir, str_replace(syls$id, "-[0-9]+$", ".wav"), sep="/")
    syls.list = split(syls, syls$wav)
    
    for(i in names(songd)) {
      songd[[i]]$song$syllable = syls.list[[i]]
    }
    return(list(songs=songd, model=syl.mc$mc))
    #return(list(wav=wavs, syls=syls.list, syl_data=syl_data, model=syl.mc$mc))
    #ssyls$called = mod$classification
    #return(list(syls=syl.mc$syls, model=syl.mc$mc))
  } else {
    return(syl_data)
  }
  
  
}
parse_song_batch = function(input, 
                            feature_set="mid3_fixed",
                            wl=512,
                            smoothed=TRUE,
                            threshed=FALSE,
                            weighted=FALSE,
                            log_transform=FALSE,
                            num_of_songs=NULL, 
                            range_to_test=c(4:9), 
                            distmethod="cor", 
                            power=2, 
                            cluster=TRUE, 
                            return_data=FALSE,
                            model_outdir = NULL) {
  if (class(input) == "character") {
    ### Expects input to be working directory of song files
    wdir = input
    print("Loading data")
    files = list.files(wdir, pattern=".wav$", full.names = T)
    
    if (!is.null(num_of_songs)) {
      files = files[1:num_of_songs]
    }
    #lapply(1:length(files), function(i) {
    songd = foreach(i=1:length(files)) %dopar% {
      #songd = mclapply(1:length(files), function(i) {
      file = files[i]
      print(i)
      song = parse_song(file, cluster = F, plot = F, thresh_method="range", )
      if (is.null(song$syllable)) return(list(song=song, syl_data=NULL))
      syl_data = extract_features2(song, 
                                   wl=wl, 
                                   feature_set=feature_set, 
                                   smoothed=smoothed, 
                                   threshed=threshed, 
                                   weighted=weighted)  
      return(list(song=song, syl_data=syl_data))
    }
    #)
    #}, mc.cores=6)  
    
    #return(songd)
    names(songd) = files
  } else if (class(input) == "list") {
    ### Expects list returned by parse_song_batch:cluster=TRUE
    songd = input
  } 
  print("Rearranging data...")
  
  # remove non-song/no peak files
  ind = unlist(lapply(songd, function(x) is.null(x$syl_data)))
  songd = songd[!ind]
  
  songs = lapply(songd, function(x) x$song)
  #wavs = lapply(songd, function(x) x$song$wav)
  
  syl_data = lapply(songd, function(x) x$syl_data)
  syl_data = do.call("rbind", syl_data)
  
  syls = lapply(songd, function(song) song$song$syllable)
  syls = do.call("rbind", syls)
  
  songd2 = list(songs=songs, syl_data=syl_data, syls=syls)
  if (return_data) 
    return(songd2)
  
  if (cluster) {
    print("Clustering syllables...")
    syl.mc = mclust_syllable_data(syls, syl_data, distmethod=distmethod, power=power, range_to_test = range_to_test, plot=T)
    #return(syl.mc)
    if (is.null(syl.mc[[1]])) 
      return(NULL)
    syls = syl.mc$syls
    syls$wav = paste(wdir, str_replace(syls$id, "-[0-9]+$", ".wav"), sep="/")
    syls.list = split(syls, syls$wav)
    
    for(i in names(songd)) {
      songd[[i]]$song$syllable = syls.list[[i]]
    }
    
    #if ()
    return(list(songs=songd, model=syl.mc$mc))
    #return(list(wav=wavs, syls=syls.list, syl_data=syl_data, model=syl.mc$mc))
    #ssyls$called = mod$classification
    #return(list(syls=syl.mc$syls, model=syl.mc$mc))
  } else {
    return(syl_data)
  }
  
}


# parse_song_batch2 = function(info, cluster=TRUE) {
#   print("Loading data")
#   #songd = mclapply(1:nrow(info), function(i) {
#   songd = foreach(i=1:nrow(info)) %dopar% {
#    # file = row[1,"wav"]
#     file = info[i,"wav"]
#     print(i)
#     wav = readWave(file)
#     fname = basename(str_replace(file, ".wav", ""))
#     wavf = filtersong(wav)
#     
#     thresh = threshold_auto(wavf, otsu)
#     syls = findpeaks_abs(wavf, min_duration = 15, max_gap = 5, thresh=thresh)
#     
#     syls$id = paste(fname, 1:nrow(syls), sep="-")
#     song = list(fname = fname, wav=wavf, syllable=syls)
#     syl_data = extract_features2(song)  
#     return(list(wav = wavf, syls = syls, syl_data=syl_data))
#   #})
#   #}, mc.cores=10)
#   }
#   return(songd)
#   print("Rearranging data...")
#   syl_data = lapply(songd, function(x) x$syl_data)
#   syl_data = do.call("rbind", syl_data)
#   
#   syls = lapply(songd, function(x) x$song.syllable)
#   syls = do.call("rbind", syls)
#   par(mfrow=c(1,1))
# 
#   if (cluster) {
#     print("Clustering syllables...")
#     mod = mclust_syllable_data(syl_data, range_to_test = c(4:12), plot=T)
#     songd$syls$called = mod$classification
#     return(list(song=songd, model=mod))
#   } else {
#     return(syl_data)
#   }
#    
# }

#' Read in info for wav and processed mat files
#' @param wdir, list or character of working directory containing .wav and .not.mat
#' @param file_ex, pattern searched for by list.files
#' @return data.frame containing file info 
load_mat_info = function(wdir, file_ex = ".not.mat$") {
  files = NULL
  if (length(wdir) > 1) {
    files = unlist(lapply(wdir, list.files, full.names = TRUE, pattern = file_ex))
  } else {
    files = list.files(wdir, full.names = TRUE, pattern = file_ex)
  }  
  
  files.wav = str_replace(files, ".not.mat", "")
  info = file.info(files.wav)
  info$date = strftime(info$mtime, '%Y-%m-%d')
  info$time = strftime(info$mtime, '%H:%M:%S')
  info$time_h = strftime(info$mtime, '%H')
  info$mat = files
  info$wav = files.wav
  return(info)
}

#' Read in info for wav from an sqlite database, 
#'   generated by insert_song_info_to_db
#' @param db_name, full name of database
#' @param birds, birds to select for, default NULL will return all
#' @param local, load data into local data.frame
#' @return if local, data.frame containing file info. 
#'   if not local, sql_tbl
load_song_info_from_db = function(db_name, birds=NULL, local=TRUE) {
  db = src_sqlite(db_name)
  song_tbl = tbl(db, "songs")
  
  if (!is.null(birds)) {
    if (length(birds) > 1) {
      filter_q = translate_sql(bird %in% birds, vars=c("bird"))
    } else {
      filter_q = translate_sql(bird == birds, vars=c("bird"))
    }
    query = sql(paste("SELECT * FROM songs WHERE ", filter_q, sep=""))
    song_tbl = tbl(db, query)
  } 
  if (local) {
    info = collect(song_tbl, n=Inf)
    info$mtime = as.POSIXct(info$mtime, tz = "PST", origin="1970-01-01")
  } else {
    info = song_tbl
  }
  return(info)
}

update_song_info_db = function(db_name, wdir="/mnt/bengal_home/song", subdir="songs", file_ex="wav$") {
  files = list.files(wdir)
  files = files[grep("[a-z]{2,2}[1-9]{2,3}", files)]
  res = lapply(files, function(f) {
    cat(sprintf("%s\n", f))
    insert_song_info_to_db(f, wdir=wdir, subdir=subdir, file_ex=file_ex, db_name=db_name)
  })
  unlist(res)
}
insert_song_info_to_db = function(bird, wdir="/mnt/bengal_home/song", subdir="songs", file_ex = ".wav$", db_name) {
  dir = paste(wdir, bird, subdir, sep="/")
  if (!dir.exists(dir))
    return("No song directory.")
  info = load_mat_info(dir, file_ex)
  info$bird = bird
  
  db = NULL
  if (!file.exists(db_name)) {
    db = src_sqlite(db_name, create=T)    
  } else {
    db = src_sqlite(db_name)
  }
  
  if (db_has_table(db$con, "songs")) {
    sql = sprintf("DELETE FROM songs WHERE bird='%s'", bird)
    dbSendQuery(db$con, sql)
    db_insert_into(db$con, "songs", info, temporary=F)
  } else {
    copy_to(db, info, "songs", temporary=F)    
  }

}


load_mat = function(file, parse_labels=TRUE) {
  mat = readMat(file)
  if (parse_labels) {
    if (is.null(mat$labels)) {
      labels = "n"
    } else {
      labels = unlist(str_split(mat$labels, ""))
    }
  } else {
    labels = "n"
  }
  out = data.frame(onsets=mat$onsets, offsets=mat$offsets, labels=labels)
  out$id = paste(basename(str_replace(file, ".wav.not.mat", "")), 1:nrow(out), sep="-")
  out[,1:2] = out[,1:2] / 1000
  out
}

load_mat_batch = function(info, parse_labels=TRUE) {
  mat = foreach(row=isplitRows(info, chunkSize=1), .combine="rbind") %do% {
    a = load_mat(row["mat"], parse_labels=parse_labels)
    a
  }
  #mat$id = 1:nrow(mat)
  return(mat)
}

process_mat = function(mat) {
  mat = mat[(mat[,2]>mat[,1]),]
  mat
}

song_to_mat = function(song, subdir=NULL) {
  fname = song$fname
  wav = song$wav
  # syls = song$syllable
  ### TODO change NA syllableto '-'
  print(fname)
  song$syllable$called = as.character(song$syllable$called)
  song$syllable$called[is.na(song$syllable$called)] = "n"
  mat = list(Fs=wav@samp.rate, 
             fname=song$fname,
             labels=paste(as.character(song$syllable$called), collapse=""),
             onsets=song$syllable$ins * 1E3,
             offsets=song$syllable$outs * 1E3,
             min.int=song$max_gap,
             min.dur=song$min_dur,
             threshold=song$thresh,
             sm.win=2)
  if (!is.null(subdir)) {
    file = basename(fname)
    prefix = str_replace(fname, file, "")
    prefix = paste(prefix, subdir, sep="/")
    if (!dir.exists(prefix)) dir.create(prefix)
    fname = paste(prefix, file, sep="/")
  }
  out = paste(fname, ".not.mat", sep="")
  write_mat(mat, out)
}
write_mat = function(mat, outfile) {
  writeMat(Fs=mat$Fs, fname=mat$fname, labels=mat$labels, onsets=mat$onsets, offsets=mat$offsets,
           min_int=mat$min.int, min_dur=mat$min.dur, threshold=mat$threshold, sm_win=mat$sm.win, con=outfile)
  
}

######## SYLLABLE CLASSIFICATION ########

#' Extract song features from wav file and peaks file generated by evsonganaly
#' @param wav, unfiltered Wave file
#' @param peaks, as generated by load_mat function
#' @return data.frame with song features that classify syllables well 
#'         (durations, wiener entropy, time-to-max amplitude, etc.)
extract_features = function(wav, peaks) {
  wavf = filtersong(wav)
  peaks = peaks[as.numeric(peaks$onsets)>.005,]
  peaks = peaks %>% filter(!labels %in% c(" ", "n"))
  if (nrow(peaks) == 0) return(NULL)
  #peaks$ind = 1:nrow(peaks)
  peaks$onsets.frames = peaks[,"onsets"] * wav@samp.rate
  peaks$offsets.frames = peaks[,"offsets"] * wav@samp.rate
  
  durations = peaks[,"offsets"] - peaks[,"onsets"]
  
  amp = env(wavf, envt="abs", plot=F)
  data = apply(peaks, 1, function(row) {
    amp_window = amp[row["onsets.frames"]:row["offsets.frames"]]
    amplitude = mean(amp_window)
    max_amp = max(amp_window)
    time_to_max = which.max(amp_window)
    time_from_max = length(amp_window) - which.max(amp_window)
    went = wiener_entropy(wavf, band=c(2000, 10000), subregion=c(as.numeric(row["onsets"]), as.numeric(row["offsets"])))
    psd1 = foreach(step=c(.005, .007, .009, .011), .combine="cbind") %do% {
      psd = spec(wavf,  fftw = T, PSD = T, plot=F, at=as.numeric(row["onsets"]) + step)
      psd
    }
    psd1 = psd1[,c(1, seq(2, ncol(psd1), 2))]
    psd = cbind(psd1[,1], apply(psd1[,2:ncol(psd1)], 1, mean))
    mid = round(nrow(psd)/2)
    psd_ratio = sum(psd[mid:nrow(psd),2]) / sum(psd[1:mid,2])
    max_freq = psd[which.max(psd[,2]),1]
    return(c(amplitude, time_to_max, time_from_max, went, max_freq, psd_ratio, psd[1:128,2]))
  })
  data = t(data)
  colnames(data) = c("amp", "time_to_max", "time_from_max", "went", "max_freq", "psd_ratio", paste("freq", 1:(ncol(data)-6), sep=""))
  data = data.frame(id=peaks$id, labels=peaks$labels, durations, data)
  data$labels = factor(data$labels)
  return(data)
}

extract_features2 = function(song, wl=512, feature_set="mid5_fixed", weighted=FALSE,
                             threshed=FALSE, threshed_val=10,
                             smoothed=TRUE, smoothed_window=4,
                             freq_limits=c(0,12)) {
  suppressMessages(require(smoother))
  feature_sets = c("mid3_fixed", "mid5_fixed", "mid3_rel", "mid5_fixed_5ms", "mid_mean", "mid2_mean")
  stopifnot(feature_set %in% feature_sets)
  
  wav = song$wav
  peaks = song$syllable
  if (is.null(peaks)) return(NULL)
  #peaks = peaks[as.numeric(peaks$ins)>.005,]
  
  ## Define frequency range (in kHz)
  res = wl / wav@samp.rate
  freq_range = seq(0, wav@samp.rate / 2000, length.out=wl/2)
  ind_range = 1:length(freq_range)
  if(!is.null(freq_limits))
    ind_range = which(freq_range > freq_limits[1] & freq_range < freq_limits[2])
  
  ## Define frequency weights
  weights_init_x = c(1, length(ind_range))
  weights_init_y = c(1, 10)
  weights_target = 1:length(ind_range)
  weights = approx(weights_init_x, weights_init_y, weights_target)$y
  
  
  offset = .010
  peaks$mids = (peaks[,"outs"] + peaks[,"ins"]) / 2
  
  points = NULL
  buffer = wl / (2*wav@samp.rate)
  if (feature_set == "mid3_fixed") {
    points = matrix(0, nrow=nrow(peaks), ncol=3)   
    points[,1] = peaks[,"mids"] - offset
    points[,2] = peaks[,"mids"]
    points[,3] = peaks[,"mids"] + offset
  } else if (feature_set == "mid5_fixed") {
    points = matrix(0, nrow=nrow(peaks), ncol=5)  
    points[,1] = peaks[,"mids"] - offset * 2
    points[,2] = peaks[,"mids"] - offset
    points[,3] = peaks[,"mids"]
    points[,4] = peaks[,"mids"] + offset
    points[,5] = peaks[,"mids"] + offset * 2 
  } else if (feature_set == "mid5_fixed_5ms") {
    offset= .005
    points = matrix(0, nrow=nrow(peaks), ncol=5)  
    points[,1] = peaks[,"mids"] - offset * 2
    points[,2] = peaks[,"mids"] - offset
    points[,3] = peaks[,"mids"]
    points[,4] = peaks[,"mids"] + offset
    points[,5] = peaks[,"mids"] + offset * 2 
  } else if (feature_set %in% c("mid3_rel")) {
    points = matrix(0, nrow=nrow(peaks), ncol=3)  
    points[,1] = peaks[,"ins"] + buffer
    points[,2] = peaks[,"mids"]
    points[,3] = peaks[,"outs"] - buffer
  } else if (feature_set == "mid_mean") {
    points = matrix(0, nrow=nrow(peaks), ncol=2)
    points[,1] = peaks[,"ins"] - buffer
    points[,2] = peaks[,"outs"] + buffer
  } else if (feature_set == "mid2_mean") {
    points = matrix(0, nrow=nrow(peaks), ncol=4)  
    points[,1] = peaks[,"ins"] + buffer
    points[,2] = peaks[,"mids"]
    points[,3] = peaks[,"mids"]
    points[,4] = peaks[,"outs"] - buffer
    ind = which((points[,2] - points[,1]) < (res))
    points[ind,2] = points[ind,1] + (res)
    ind = which((points[,4] - points[,3]) < (res))
    points[ind,3] = points[ind,4] - (res)
  }
  dur = seewave::duration(wav)
  
  # Filter out peaks on edge of recording
  if (feature_set == "mid2_mean") {
    ind = points[,2] > (dur - (2*buffer))  | points[,3] < (2*buffer)
    points = points[!ind,]
    peaks = peaks[!ind,]
  }
  ind = points[,1]<(2*buffer) | points[,ncol(points)] > (dur- (2*buffer))
  points = points[!ind,]
  peaks = peaks[!ind,]
  
  if(nrow(peaks) == 0) 
    return(NULL)
  # points[points>(dur-buffer)] = (dur-buffer)
  #points = points[(points[,2] - points[,1]) > 2*buffer,]
  peaks$ins.frames = peaks[,"ins"] * wav@samp.rate
  peaks$outs.frames = peaks[,"outs"] * wav@samp.rate
  # peaks[,1:2] = peaks[,1:2] + offset
  
  psd1 = NULL
  data = NULL
  data = foreach(row=isplitRows(points, chunkSize=1), .inorder = T, .combine="rbind") %do% {
    
    d = NULL
    psd1 = NULL
    if (feature_set %in% c("mid_mean")) {
      psd1 = seewave::meanspec(wav, wl=wl, fftw=T, PSD=T, norm=T, plot=F, from=row[1], to=row[2])[ind_range,"y"]
    } else if (feature_set == "mid2_mean") {
      p1 = seewave::meanspec(wav, wl=wl, fftw=T, PSD=T, norm=T, plot=F, from=row[1], to=row[2])[ind_range,"y"]
      p2 = seewave::meanspec(wav, wl=wl, fftw=T, PSD=T, norm=T, plot=F, from=row[3], to=row[4])[ind_range,"y"]
      #p1 = 0
      psd1 = c(p1, p2)
    } else {
      psd1 = foreach(x= isplitCols(row, chunkSize=1), .combine="c") %do% {
        d = seewave::spec(wav, wl=wl, fftw=T, PSD=T, norm=T,  plot=F, at=x)[ind_range,"y"]
        if (weighted) d = d * weights
        d
      }
    } 
    
    if (threshed) {
      psd1 = psd1 - song$thresh * threshed_val
      psd1[psd1 < 1] = 1
    }
    if (smoothed) {
      psd1 = na.omit(smth.gaussian(psd1, window = smoothed_window ))
    }
    return(psd1)
  }
  rownames(data) = peaks$id
  return(data)
}

extract_features_batch = function(info) {
  data = foreach(row=isplitRows(info, chunkSize=1), .combine="rbind") %dopar% {
    wav = readWave(row[1, "wav"])
    mat = load_mat(row[1, "mat"])
    ex = extract_features(wav, mat)
    ex
  }
  syl_table = table(data$labels)
  data = droplevels(data[!data$labels %in% names(syl_table)[syl_table<(nrow(data) * .05)],])
  data
}

extract_features_syl = function(wav, mat, subsamp=1, wl=1024) {
  fs = wav@samp.rate
  window = wl / fs
  if (is.null(mat[[3]]))
    return(list(syls=NULL,
                gaps=NULL))
  
  peaks = data.frame(ins=mat$onsets, outs=mat$offsets, labels=parse_labels(mat))
  peaks[,1:2] = peaks[,1:2] / 1000
  peaks = peaks %>% filter(!(labels %in% c(" ", "-")))
  peaks = peaks %>% filter(ins>window, outs<(seewave::duration(wav)-window))
  if (nrow(peaks) < 2)
    return(list(syls=NULL, 
                gaps=NULL))
  
  durations = peaks[,2] - peaks[,1]
  # Gap durations
  gap_durations = vector("numeric", length = nrow(peaks) - 1)
  for (i in 1:(length(gap_durations)-1)) {
    gap_durations[i] = peaks[i+1,1] - peaks[i,2]
  }
  
  # amplitude properties
  amp = abs(wav@left[seq(1, length(wav@left), subsamp)])
  syl_amps = peaks %>% rowwise() %>% do({
    amps = amp[(round(as.numeric(.[1])*fs / subsamp)):(round(as.numeric(.[2])*fs / subsamp))]
    time_to_max = which.max(amps) * subsamp / fs
    time_from_max = (length(amps) - which.max(amps)) * subsamp / fs
    half_max_1 = time_to_max / 2
    half_max_2 = time_from_max + (time_from_max / 2)
    half_max_duration = half_max_2 - half_max_1
    
    spec_feat = calc_spectral_features(wav, wl=wl, subregion=c(as.numeric(.[1]), as.numeric(.[2])), overlap=0)
    data.frame(mean_syl_amp=mean(amps),
               sd_syl_amp=sd(amps),
               time_to_max = time_to_max,
               spec_feat
    )
  })
  
  out = cbind(peaks,
              durations,
              syl_amps)
  return(list(syls=out, gaps=gap_durations))
}

#' Train random forest classifier on syllable labels and song features
#' @param data, song feature file as generated by extract_features
#' @return randomForest object
rf_syllables = function(data) {
  require(randomForest)
  #ind = sample(1:nrow(data), replace = F, size = round(nrow(data)*.80))
  #train = data[ind,]
  #test = data[-ind,]
  #fit.rf = randomForest(labels~., train, xtest = test[,2:ncol(test)], ytest = test[,1], ntree=5000)
  fit.rf = randomForest(labels~., data, ntree=5000)
  return(fit.rf)
}

#' Train random forest classifier on syllable labels and song features
#' @param data, song feature file as generated by extract_features
#' @return randomForest object
rf_syllables2 = function(data) {
  #ind = sample(1:nrow(data), replace = F, size = round(nrow(data)*.80))
  #train = data[ind,]
  #test = data[-ind,]
  #fit.rf = randomForest(labels~., train, xtest = test[,2:ncol(test)], ytest = test[,1], ntree=5000)
  fit.rf = randomForest(labels~., data, ntree=5000)
  return(fit.rf)
}


predict_syllables = function(rf, data) {
  pred = predict(rf, data)  
}

test_rf_sample_size = function(data) {
  ind = sample(1:nrow(data), round(nrow(data)*.8), replace=F)
  train = data[ind,]
  test = data[-ind,]
  test_ind = seq(100, nrow(train), 100)
  rs = mclapply(test_ind, function(x) {
    rf = lapply(1:10, function(x) {
      
      rf = randomForest(labels~., train[sample(1:nrow(train), test_ind, replace=F),], xtest = test[,2:ncol(test)], ytest = test[,1], ntree=5000)
      confusion = rf$confusion
      n = confusion
      off_ind = lower.tri(confusion) | upper.tri(confusion)
      off = sum(confusion[off_ind])
      on = sum(confusion[!off_ind])
      return(off / on)
    })
    mean(unlist(rf))
  }, mc.cores=10)
  rs = unlist(rs)
  return(data.frame(ind = test_ind, error= rs))
}

plot_clustering_3d = function(data) {
  require(plotly)
  plot_ly(data, x = PC1, 
          y = PC2, 
          z = PC3, 
          color = sc,
          #size = PC4,
          type = "scatter3d", 
          mode = "markers", 
          filename="r-docs/3d-scatter")
}

# plot_classified_syllable_multi = function(songs, syls, num_syls=5) {
#   
#   syls_sub = syls %>% group_by(called) %>%$ do({sample_n(., num_syls)})
#   song_names = str_replace(syls_sub$id, "-[0-9]", ".wav")
#   
#   
# }

plot_syllable_range = function(song, duration=2) {
  wav = song$wav
  peaks = song$syllable
  dots <- lapply(label_name, as.symbol)
  
  peaks_s = peaks %>% group_by_(.dots=dots) %>% do({
    if (nrow(.)<=num_syls) {
      return(.)
    } else {
      return(sample_n(., num_syls))  
    }
  })
  peaks_s1 = split(peaks_s, peaks_s[,label_name])
  #par(mfrow)
  
  
  
}
plot_classified_syllable = function(song, label_name="called", num_syls=5) {
  
  wav = song$wav
  peaks = song$syllable
  dots <- lapply(label_name, as.symbol)
  
  peaks_s = peaks %>% group_by_(.dots=dots) %>% do({
    if (nrow(.)<=num_syls) {
      return(.)
    } else if (nrow(.) == 0) {
      return(NA)
    } else {
      return(sample_n(., num_syls))  
    }
  })
  peaks_s = droplevels(peaks_s)
  peaks_s1 = split(peaks_s, peaks_s[,label_name])
  #par(mfrow)
  plots = lapply(1:length(peaks_s1), function(d)  {
    suppressWarnings(plot_multiple_syl(wav, peaks_s1[[d]][,1:2], title=names(peaks_s1)[d] ))
  })
  
  ncol = 3
  nrow = ceiling(length(plots) / ncol)
  plot.list = c(plots, list(ncol=ncol, nrow=nrow))
  do.call(grid.arrange, plot.list)
}

plot_spectro = function(wav, subregion=NULL, wl=512, overlap=50, labels=NULL, label_name = "called") {
  require(matlab)
  options(stringsAsFactors=T) # kludge to get stat_contour to work
  spectrogram = spectro(wav, plot = FALSE, tlim=subregion, fftw=T)
  frequency = rep(spectrogram$freq, times = ncol(spectrogram$amp))
  time = rep(spectrogram$time, each = nrow(spectrogram$amp))
  amplitude = as.vector(spectrogram$amp)
  df = data.frame(time, frequency, amplitude)
  #df = df %>% mutate(density = MASS::kde2d(time, frequency))
  gg = ggplot(df, aes_string(x = "time", 
                             y = "frequency",
                             z = "amplitude"))
  gg = gg + stat_contour(aes(fill=..level..), geom="polygon", binwidth=1)
  # gg = gg + geom_raster( interpolate=T)
  #gg = gg + stat_density(aes(fill = ..density..), geom = "raster", position = "identity")
  gg = gg + scale_fill_gradientn("Amplitude, (db)",
                                 colours=jet.colors(7), 
                                 breaks=c(-40, -15, 0),
                                 limits=c(-40, -15),
                                 na.value="transparent")
  gg = gg + labs(x="Time", y="Frequency") + ylim(0,10) + xlim(subregion)
  
  if (!is.null(labels)) {
    labels$time = (labels[,2] + labels[,1]) / 2
    labels = cbind(labels, data.frame(frequency=0, amplitude=0))
    gg = gg + geom_text(data=labels, aes_string(label=label_name))
  }
  #print(gg)
  return(gg)
}

plot_single_syl = function(wav, subregion, wl=512, overlap=50) {
  input = inputw(wav)
  wave = input$w
  f = input$f
  
  buffer = .010
  subregion1 = subregion
  subregion1[,1] = subregion1[,1] - buffer
  subregion1[,2] = subregion1[,2] + buffer
  sf = subregion1 * f
  steps = seq(sf[1,1], sf[1,2], wl - (overlap * wl/100))
    
   # step = steps[[row]]
    data = stft(wave = wave, f = f, wl = wl, zp = 0, step = steps, 
             wn = "hanning", fftw = TRUE, scale = T, complex = FALSE)
    freqs = data.frame(freq_ind=1:nrow(data), freq=seq(0, (wav@samp.rate-1) / 2, wav@samp.rate / wl)/1000)
    times = data.frame(time_ind=1:ncol(data), time=steps / f)
    #z = 20 * log10(z)
  
  datam = melt(data)
  colnames(datam) = c("freq_ind", "time_ind", "value")
  datam = datam %>% left_join(freqs) %>% left_join(times)
  maxz = max(datam$value)
  theme_set(theme_classic())
  gg = ggplot(datam, aes(time, freq, z=value)) + stat_contour(aes(fill=..level..), geom="polygon", binwidth=.1)
  #gg = gg + scale_fill_gradientn(colours=rainbow(100), limits=c(50, maxz), breaks = seq(maxz-10, maxz), na.value="white" )
  #gg = gg + scale_fill_gradientn(colours=rainbow(7), guide=T, values = seq(maxz-100, maxz), rescaler = function(x, ...) x, oob = identity, na.value="white" )
  gg = gg + scale_fill_gradientn(colours=jet.colors(7), breaks=c(.1, 1),
                                 limits=c(0, 1), na.value="transparent")
  gg = gg + geom_vline(xintercept=as.numeric(subregion[,1]), linetype=3)
  gg = gg + geom_vline(xintercept=as.numeric(subregion[,2]), linetype=3)
  gg = gg + scale_x_continuous(breaks=c(round(as.numeric(subregion[,1]), 3), round(as.numeric(subregion[,2]), 3)))
  gg = gg + labs(x="seconds", y="kHz", title="") + theme(axis.line=element_blank(),
                                               #axis.text.x=element_blank(),
                                                #  axis.text.y=element_blank(),
                                               #axis.ticks=element_blank(),
                                                  #axis.title.x=element_blank(),
                                                  #axis.title.y=element_blank(),
                                               plot.margin=unit(c(-1,0,-1,0), "lines"),
                                                  legend.position="none",
                                                  panel.background=element_blank(),
                                               panel.border=element_blank(),
                                               panel.grid.major=element_blank(),
                                                  panel.grid.minor=element_blank(),
                                               plot.background=element_blank())
  gg = gg + ylim(0,10)
  return(gg)
}

plot_multiple_syl = function(wav, subregions, wl=512, overlap=50, title="") {
  input = inputw(wav)
  wave = input$w
  f = input$f
  subregions$ind = 1:nrow(subregions)
  
  steps =  lapply(1:nrow(subregions), function(i) {
    sf = subregions[i,] * f
    seq(sf[1,1], sf[1,2], wl - (overlap * wl/100))
  })
  steps_num = sum(unlist(lapply(steps, length)))
  steps_max = max(unlist(lapply(steps, length)))
  data = matrix(0, nrow=(wl/2), ncol=(steps_max + 5) * nrow(subregions))
  for(row in 1:nrow(subregions)) {
    
    step = steps[[row]]
    z = stft(wave = wave, f = f, wl = wl, zp = 0, step = step, 
             wn = "hanning", fftw = TRUE, scale = T, complex = FALSE)
    #z = 20 * log10(z)
    data[,(row -1) * (steps_max + 5) + (1 : ncol(z))] = z
  }  
  
  datam = melt(data)
  maxz = max(datam$value)
  theme_set(theme_classic())
  gg = ggplot(datam, aes(Var2, Var1, z=value)) + stat_contour(aes(fill=..level..), geom="polygon", binwidth=.1)
  #gg = gg + scale_fill_gradientn(colours=rainbow(100), limits=c(50, maxz), breaks = seq(maxz-10, maxz), na.value="white" )
  #gg = gg + scale_fill_gradientn(colours=rainbow(7), guide=T, values = seq(maxz-100, maxz), rescaler = function(x, ...) x, oob = identity, na.value="white" )
  gg = gg + scale_fill_gradientn(colours=jet.colors(7), breaks=c(.1, 1),
                                 limits=c(0, 1), na.value="transparent")
  gg = gg + labs(x="", y="", title=title) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                  axis.title.x=element_blank(),
                                                  axis.title.y=element_blank(),legend.position="none",
                                                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                                  panel.grid.minor=element_blank(),plot.background=element_blank())
  gg = gg + ylim(0,128)
  return(gg)
}

calc_syllable_number_from_mat = function(info, thresh=.01) {
  if(!all(c("bird", "labels") %in% colnames(info)))
    stop("'bird' and 'labels' must be columns")
  info %>% 
    dplyr::filter(!(labels %in% c("-", " "))) %>% 
    group_by(bird, labels) %>% 
    summarize(nlabel=n()) %>% 
    ungroup(.) %>% group_by(bird) %>%
    mutate(nlabel_norm = nlabel / sum(nlabel)) %>% 
    filter(nlabel_norm > thresh) %>% 
    summarize(length(unique(labels)))
}

######## PEAK ANALYSIS ########

#' Find amplitude peaks in wav file
#' @param wav, the input Wave object
#' @param min_duration, minimum duration required for peak calling
#' @param max_gap, maximum distance (in ms) between peaks before merging
#' @param thresh, threshold for peak definition, defined as fraction of max amplitude
#' @return data.frame containing peak starts (col1) and stops (col2) defined in seconds
findpeaks = function(wav, min_duration=50, max_gap=75, thresh=0.1) {
  samp_rate = wav@samp.rate
  min_duration = min_duration * samp_rate / 1000 
  max_gap = max_gap * samp_rate / 1000
  wav_env = env(wav, plot=F)
  max_val = max(wav_env)
  wav_env = wav_env / max_val
  above_thresh_ind = which(wav_env[,1]>thresh) 
  thresh_diff = which(diff(above_thresh_ind)>max_gap)
  df = data.frame(ins=c(above_thresh_ind[1], 
                        above_thresh_ind[thresh_diff+1]),
                  outs=c(above_thresh_ind[thresh_diff], 
                         above_thresh_ind[length(above_thresh_ind)]))
  return(df[(df[,2] - df[,1])>min_duration,] / samp_rate)
  if ((df[,2] - df[,1]) == 0) return(NULL)
  return (df / samp_rate)
}

#' Find amplitude peaks in wav file, using absolute amplitude threshold
#' @param wav, the input Wave object
#' @param min_duration, minimum duration required for peak calling
#' @param max_gap, maximum distance (in ms) between peaks before merging
#' @param thresh, threshold for peak definition, defined as fraction of max amplitude
#' @return data.frame containing peak starts (col1) and stops (col2) defined in seconds
findpeaks_abs = function(wav, min_duration=50, max_gap=75, max_duration=300, thresh=1E7) {
  samp_rate = wav@samp.rate
  min_duration = min_duration * samp_rate / 1000 
  max_duration = max_duration * samp_rate / 1000 
  max_gap = max_gap * samp_rate / 1000
  #wav_env = seewave::env(wav, envt = "abs", plot=F)
  wav_env = as.matrix(wav@left^2)
  #max_val = max(wav_env)
  #wav_env = wav_env / max_val
  above_thresh_ind = which(wav_env[,1]>thresh) 
  if (length(above_thresh_ind) == 0) return(data.frame())
  thresh_diff = which(diff(above_thresh_ind)>max_gap)
  df = data.frame(ins=c(above_thresh_ind[1], 
                        above_thresh_ind[thresh_diff+1]),
                  outs=c(above_thresh_ind[thresh_diff], 
                         above_thresh_ind[length(above_thresh_ind)]))
  durations = df[,2] - df[,1]
  ind = durations>min_duration & durations<max_duration
  return(df[ind,] / samp_rate)
  #if ((df[,2] - df[,1]) == 0) return(data.frame())
  #return (df / samp_rate)
}

#' Find amplitude peaks in wav file, using absolute amplitude threshold
#' @param wav, the input Wave object
#' @param min_duration, minimum duration required for peak calling
#' @param max_gap, maximum distance (in ms) between peaks before merging
#' @param thresh, threshold for peak definition, defined as fraction of max amplitude
#' @return data.frame containing peak starts (col1) and stops (col2) defined in seconds
findpeaks_abs_env = function(wav_env, samp_rate, min_duration=50, max_gap=75, max_duration=300, thresh=1E7, subsamp=1) {
  wav_env = as.matrix(wav_env[seq(1, nrow(wav_env), subsamp),])
  
  min_duration = min_duration * samp_rate / (1000 * subsamp) 
  max_duration = max_duration * samp_rate / (1000 * subsamp) 
  max_gap = max_gap * samp_rate / (1000 * subsamp)
  above_thresh_ind = which(wav_env[,1]>thresh) 
  if (length(above_thresh_ind) == 0) return(data.frame())
  thresh_diff = which(diff(above_thresh_ind)>max_gap)
  df = data.frame(ins=c(above_thresh_ind[1], 
                        above_thresh_ind[thresh_diff+1]),
                  outs=c(above_thresh_ind[thresh_diff], 
                         above_thresh_ind[length(above_thresh_ind)]))
  durations = df[,2] - df[,1]
  ind = durations>min_duration & durations<max_duration
  return(subsamp * df[ind,] / samp_rate)
  #if ((df[,2] - df[,1]) == 0) return(data.frame())
  #return (df / samp_rate)
}

findpeaks_range = function(wav, min_duration=50, max_gap=75, max_duration=300, thresh_range=seq(-1,1,.1), absolute=F, subsamp=1, floor=F) {
  if (absolute) {
    thresh = threshold_auto(wav, mean_sd2, log=F, abs=T, factor=thresh_range, floor=floor)
    wav_env = seewave::env(wav, envt = "abs", plot=F)
  } else {
    thresh = threshold_auto(wav, mean_sd2, log=F, abs=F, factor=thresh_range, floor=floor)
    wav_env = as.matrix(wav@left^2)
  }
  #peaks = foreach(th=1:length(thresh_range)) %do% {
  
  
  samp_rate = wav@samp.rate
  peaks = lapply(1:length(thresh_range), function(th){
    peak = findpeaks_abs_env(wav_env, samp_rate, min_duration=min_duration, max_gap=max_gap, max_duration=max_duration, thresh=thresh[th], subsamp=subsamp)
    peak
  })
  peak_info = data.frame(num_peak = unlist(lapply(peaks, nrow)), thresh_factor=thresh_range)
  peak_info$num_peak_diff = c(0, diff(peak_info$num_peak))
  
  ind1 = which.max(peak_info$num_peak_diff) 
  #flat_value = median(peak_info$num_peak) / 2
  flat_value = 2
  peak_info$num_peak_diff = ifelse(abs(peak_info$num_peak_diff) <= flat_value, 0, peak_info$num_peak_diff )
  diff_rle = rle2(peak_info[ind1:nrow(peak_info),"num_peak_diff"], indices = T)
  
  plateau = round((1 / diff(thresh_range)[1]) / 2)
  ind = diff_rle[which(diff_rle[,1]==0 & diff_rle[,4] > plateau)[1],2]
  ind = ind + plateau + ind1 - 1
  #ind = which.max(peak_info$num_peak_diff) + 2
  ind = ifelse(ind > length(peaks), length(peaks), ind)
  #print(peak_info$thresh_factor[ind])
  selected = peaks[[ind]]
  return(list(peaks=peaks[[ind]], thresh_factor=thresh_range[ind], thresh=thresh[ind]))
}


#' Find frequency peaks in PSD
#' @param mat, the input PSD matrix, col1 is frequncy, col2 is power
#' @param min_value, minimum frequency value, in kHz
#' @param min_size, minimum band size of peak
#' @param max_gap, maximum distance (in kHz) between peaks before merging
#' @param thresh, threshold for peak definition, defined as fraction of max amplitude
#' @return data.frame containing peak starts (col1) and stops (col2) defined in seconds
findpeaks_freq = function(mat, min_value=2, min_size=50, max_gap=100, thresh=0.1, absolute=FALSE) {
  
  mat = mat[mat[,1]>min_value,]
  if (!absolute) {
    max_val = max(mat[,2])
    mat[,2] = mat[,2] / max_val
  }
  scale_val = (mat[2,1] - mat[1,1])
  above_thresh_ind = which(mat[,2]>thresh) 
  if (length(above_thresh_ind) == 0) return(data.frame())
  thresh_diff = which(diff(above_thresh_ind)>(max_gap/scale_val))
  df = data.frame(ins=c(mat[above_thresh_ind[1],1], 
                        mat[above_thresh_ind[thresh_diff+1],1]),
                  outs=c(mat[above_thresh_ind[thresh_diff],1], 
                         mat[above_thresh_ind[length(above_thresh_ind)],1]))
  return(df[(df[,2] - df[,1])>min_size,])
}

find_motifs = function(peaks, max_gap=500, min_duration=500) {
  max_gap = max_gap / 1000
  min_duration = min_duration / 1000
  peak_mids = apply(peaks[,1:2], 1, mean)
  thresh_diff = which(diff(peak_mids) > max_gap)
  df = data.frame(ins=c(peaks[1,1], 
                        peaks[thresh_diff+1,1]),
                  outs=c(peaks[thresh_diff,2], 
                         peaks[nrow(peaks),2]))
  durations = df[,2] - df[,1]
  ind = durations>min_duration 
  return(df[ind,])
  
}
#' Find temporal distance between peaks
#' @param peaks, data frame as outputted by findpeaks, col1 is peak starts, col2 is peak stops
#' @return vector of interpeak distances
peak_dist = function(peaks) {
  if (nrow(peaks)<2) return(NA)
  dists = vector("numeric", nrow(peaks)-1)
  for (i in 1:(nrow(peaks)-1)) {
    dists[i] = peaks[i+1,1] - peaks[i,2]
  }
  return(dists)
}

filter_by_gaps = function(peaks, max_gap_size) {
  pdists = peak_dist(peaks)
  ind = which(pdists<max_gap_size)
  if (length(ind)>0) {
    ind1 = unique(unlist(lapply(ind, function(x) c(x, x+1))))
    peaks[ind1,]
  } else {
    return(NULL)
  }
  
}


interpeak_midpoint = function(peaks) {
  if (nrow(peaks)<2) return(NA)
  dists = vector("numeric", nrow(peaks)-1)
  for (i in 1:(nrow(peaks)-1)) {
    dists[i] = sum(peaks[i+1,1] + peaks[i,2]) / 2
  }
  return(dists)
}

match_peaks = function(d1, d2, multisong=NULL, tol=.01) {
  d1_song = unlist(lapply(str_split(d1$id, "-"), function(x) x[1]))
  d2_song = unlist(lapply(str_split(d2$id, "-"), function(x) x[1]))
  d1_s = split(d1, d1_song)
  d2_s = split(d2, d2_song)
  
  if (length(d1_s) < length(d2_s)) {
    d2_s = d2_s[names(d2_s) %in% names(d1_s)]
  } else {
    d1_s = d1_s[names(d1_s) %in% names(d2_s)]
  }
  
  for(i in 1:length(d2_s)) {
    d1_curr = d1_s[[i]]
    d2_curr = d2_s[[i]]
    d1_mids = (d1_curr[,2] + d1_curr[,1]) / 2
    d2_mids = (d2_curr[,2] + d2_curr[,1]) / 2  
    mid_diffs = dist(d1_mids, d2_mids)
    
    if (length(d1_mids) < length(d2_mids)) {
      margin = 1
    } else {
      margin = 2
    }
    mins = apply(mid_diffs, margin, which.min)
    
    if (margin == 2) {
      d2_curr$id = paste(names(d2_s)[i], mins, sep="-")
      d2_s[[i]] = d2_curr
    } else {
      d1_curr$id = paste(names(d1_s)[i], mins, sep="-")
      d1_s[[i]] = d1_curr
    }
  }
  d1 = do.call(rbind, d1_s)
  d2 = do.call(rbind, d2_s)
  
  return(list(d1=d1, d2=d2))
}


plot_peaks = function(peaks, yval=-1E7, ...) {
  apply(peaks, 1, function(peak) {
    id = NULL
    
    if (length(grep('-', peak[3]))>0) {
      id = str_split(peak[3], "-")[[1]][2]
    } else {
      id = peak[3]
    }
    segments(as.numeric(peak[1]), yval, as.numeric(peak[2]), yval, ...)
    text(as.numeric((as.numeric(peak[2]) + as.numeric(peak[1]))/2), y=yval + yval*.5, labels = id)
  })
}



######### SONG FINDERS ########
#' Classify Wave object as song or not song
#'     First filters for small Weiner entropy values within given band
#'     Next filters for small interpeak distances
#' @param wav, Wave object
#' @param band, selected frequency band
#' @param wein_thresh, maximum allowed Wiener entropy in frequency band
#' @param dist_thresh, maximum interpeak distance allowed
#' @return logical, song or not song
songfinder = function(wav, band=c(5000,7000), 
                      wein_thresh=0.4, 
                      dist_thresh = 0.06,
                      min_duration = 1000) {
  wav = filtersong(wav)
  wein = wiener_entropy(wav, band=band)
  if (wein > wein_thresh) {
    return(FALSE)
  } 
  peaks = findpeaks(wav, min_duration=min_duration, max_gap=200, thresh=0.05) 
  return(nrow(peaks) > 0)
}


songfinder2 = function(wav, 
                       max_gap = 10, 
                       min_duration = 15,
                       max_duration = 300 , 
                       thresh=0.2, 
                       min_num_peaks = 5, 
                       max_num_peaks=10,
                       amp_ratio_max_gap=100,
                       amp_ratio_min=0,
                       amp_ratio_max=.4,
                       rate_min=5,
                       rate_max=15,
                       low_noise=T) {
  wavf = highpass_filter(wav, from = 500, wl = 1024, ovlp = 25)
  subsamp = 10
  peak_info = findpeaks_range(wavf, min_duration = min_duration, 
                              max_duration = max_duration, 
                              max_gap = max_gap, 
                              absolute = F, 
                              thresh_range=seq(0, 5, .2),
                              subsamp=subsamp)
  
  # peak_info = findpeaks_range(wavf, 
  #                             min_duration = min_dur, 
  #                             max_duration = max_duration, 
  #                             max_gap = max_gap, 
  #                             absolute=absolute,
  #                             thresh_range=thresh_range,
  #                             subsamp=subsamp,
  #                             floor=F)
  # syls = peak_info$peaks
  peaks = peak_info$peaks
  thresh = peak_info$thresh
  
  if (is.null(peaks))
    return(FALSE)
  
  num_peaks = nrow(peaks)
  thresh = peak_info$thresh
  if (num_peaks < min_num_peaks) return(FALSE)
  
  # Filter for poor amplitude definition of peaks
  amp = abs(wavf@left[seq(1,length(wavf@left), subsamp)])
  amp_ratio = calc_amp_ratio(amp, wavf@samp.rate, peaks, max_gap=amp_ratio_max_gap, subsamp=subsamp)
  
  rate = num_peaks / (peaks[num_peaks,2] - peaks[1,1])
  total_duration = peaks[num_peaks,2] - peaks[1,1]
  peak_dists = peak_dist(peaks)
  peak_dists = peak_dists[peak_dists<0.5]
  cv_peak_dist = sd(peak_dists) / mean(peak_dists)
  
  if (low_noise) {
    res = (total_duration > 0.5) & 
      (amp_ratio < amp_ratio_max) & (amp_ratio > amp_ratio_min) &
      (rate > rate_min) & (rate < rate_max)
  } else {
    res = (total_duration > 0.5) & 
      (amp_ratio < amp_ratio_max) & (amp_ratio > amp_ratio_min) &
      (rate > rate_min) & (rate < rate_max)
  }
  #(cv_peak_dist < 1.1) & (cv_peak_dist > 0)
  
  #  (num_peaks >= min_num_peaks)
  # (rate > min_num_peaks ) & 
  #  (rate < max_num_peaks ) & 
  #      nrow(peaks) > 12 &
  
  return(res)
}

songfinder3 = function(wav, 
                       to_return=NULL, 
                       max_gap = 10, 
                       min_duration = 15, 
                       max_duration = 300, 
                       min_num_peaks = 5, 
                       max_num_peaks=NULL,
                       thresh_method=mean_sd,
                       thresh_factor=1.5) {
  wavf = filtersong(wav)
  wavf = rmnoise(wavf, output="Wave")
  fs = wav@samp.rate
  
  # Identify peak regions
  thresh = threshold_auto(wavf, thresh_method, log=T, factor=thresh_factor)
  peaks = findpeaks_abs(wavf, min_duration=min_duration, max_gap=max_gap, max_duration=max_duration, thresh=thresh)
  
  #return(peaks)
  if (nrow(peaks) < 8) return(FALSE)
  
  # Durations
  durations = peaks[,2] - peaks[,1]
  
  # Filter by gap durations. Subsong has regular gaps
  gap_durations = vector("numeric", length = nrow(peaks) - 1)
  for (i in 1:(length(gap_durations)-1)) {
    gap_durations[i] = peaks[i+1,1] - peaks[i,2]
  }
  if (to_return=="gap_durations")
    return(gap_durations)
  
  return((total_duration > 1) & (rate > min_num_peaks ) & (rate < max_num_peaks ))
}

songfinder_features = function(wav, 
                               filtersong=TRUE, 
                               rmnoise=FALSE, 
                               max_gap = 10, 
                               min_duration = 15, 
                               max_duration = 300 , 
                               thresh_method=mean_sd, 
                               thresh_factor=.25, 
                               min_num_peaks = 5, 
                               max_num_peaks=10,
                               subsamp=1) {
  wavf = NULL
  wl = 1024
  if (filtersong) { 
    #wavf = filtersong(wav)
    wavf = highpass_filter(wav, from = 500, wl = 1024, ovlp = 25)
  } else {
    wavf = wav
  }
  if (rmnoise)
    wavf = rmnoise(wavf, output="Wave")
  frate = wav@samp.rate
  #thresh = threshold_auto(wavf, mean_sd, sd_factor=sd_factor)
  
  # Identify peak regions
  #thresh = threshold_auto(wavf, thresh_method, log=T, factor=thresh_factor)
  #peaks = findpeaks_abs(wavf, min_duration=min_duration, max_gap=max_gap, max_duration=max_duration, thresh=thresh)
  peak_info = findpeaks_range(wavf, 
                              min_duration = min_duration, 
                              max_gap = max_gap, 
                              max_duration = max_duration, 
                              absolute = F,
                              thresh_range=seq(0, 2, .2))
  peaks = peak_info$peaks
  # fs = wav@samp.rate
  window = wl / frate
  peaks = peaks %>% filter(ins>window, outs<(seewave::duration(wav)-window))
  
  peak_dists = peak_dist(peaks)
  peak_dists = peak_dists[peak_dists<0.5]
  cv_peak_dist = sd(peak_dists) / mean(peak_dists)
  num_peaks = nrow(peaks)
  thresh_factor = peak_info$thresh_factor
  thresh = peak_info$thresh
  # syllable rate, or number of syllables in peak region
  if (nrow(peaks) > 0) {
    rate = nrow(peaks) / (peaks[nrow(peaks),2] - peaks[1,1])
  } else {
    rate = 0
  }
  
  # total duration of peak region.
  if (nrow(peaks) > 0) {
    total_duration = peaks[nrow(peaks),2] - peaks[1,1]
  } else {
    total_duration = 0
  }
  
  # Durations
  if (nrow(peaks) > 1) {
    durations = peaks[,2] - peaks[,1]
    mean_durations = mean(durations)
    med_durations = median(durations)
    q75_durations = quantile(durations, .75)
    sd_durations = sd(durations)
  } else if (nrow(peaks) == 1) {
    durations = unlist(peaks[2] - peaks[1])
    mean_durations = durations
    med_durations = durations
    q75_durations = durations
    sd_durations = 0
  } else {
    mean_durations = 0
    med_durations = 0
    q75_durations = 0
    sd_durations = 0
  }
  
  # Gap durations
  if (nrow(peaks) > 2) {
    gap_durations = vector("numeric", length = nrow(peaks) - 1)
    for (i in 1:(length(gap_durations)-1)) {
      gap_durations[i] = peaks[i+1,1] - peaks[i,2]
    }
    mean_gap_durations = mean(gap_durations)
    med_gap_durations = median(gap_durations)
    sd_gap_durations = sd(gap_durations)
  } else if (nrow(peaks) == 2) {
    mean_gap_durations = peaks[2,1] - peaks[1,2]
    med_gap_durations = peaks[2,1] - peaks[1,2]
    sd_gap_durations = 0
  } else {
    mean_gap_durations = 0
    med_gap_durations = 0
    sd_gap_durations = 0
  }
  
  # amplitude properties
  amp = abs(wavf@left[seq(1,length(wavf@left), subsamp)])
  
  #amp = seewave::env(wavf, envt="abs", plot=F)
  mean_amp = mean(amp)
  sd_amp = sd(amp)
  
  if (nrow(peaks) > 1) { 
    syl_amps = apply(peaks[,1:2], 1, function(peak) {
      amps = amp[round(peak[1]*frate/subsamp):round(peak[2]*frate/subsamp)]
      time_to_max = which.max(amps)
      time_from_max = length(amps) - which.max(amps)
      half_max_1 = time_to_max / 2
      half_max_2 = time_from_max + (time_from_max / 2)
      half_max_duration = half_max_2 - half_max_1
      c(mean_syl_amp=mean(amps), 
        sd_syl_amp=sd(amps),
        time_to_max = time_to_max,
        time_from_max = time_from_max,
        half_max_duration = half_max_duration)
    })
    mean_syl_amps = mean(syl_amps[1,])
    sd_syl_amps = mean(syl_amps[2,])
    mean_time_to_max = mean(syl_amps[3,])
    med_time_to_max = median(syl_amps[3,])
    min_time_to_max = min(syl_amps[3,])
    max_time_to_max =max(syl_amps[3,])
    q25_time_to_max = quantile(syl_amps[3,], .25)
    q75_time_to_max = quantile(syl_amps[3,], .75)
    sd_time_to_max = sd(syl_amps[3,])
    mean_time_from_max = mean(syl_amps[4,])
    med_time_from_max = median(syl_amps[4,])
    sd_time_from_max = sd(syl_amps[4,])
    #mean_syl_amps = mean(syl_amps[1,])
    #sd_syl_amps = mean(syl_amps[2,])
    mean_half_max_duration = mean(syl_amps[5,])
    med_half_max_duration = median(syl_amps[5,])
    min_half_max_duration = min(syl_amps[5,])
    max_half_max_duration =max(syl_amps[5,])
    q25_half_max_duration = quantile(syl_amps[5,], .25)
    q75_half_max_duration = quantile(syl_amps[5,], .75)
  } else {
    mean_syl_amps = 0
    sd_syl_amps = 0
    mean_time_to_max = 0
    med_time_to_max = 0
    sd_time_to_max = 0
    mean_time_from_max = 0
    med_time_from_max = 0
    sd_time_from_max = 0
    min_time_to_max =0
    max_time_to_max =0
    q25_time_to_max = 0
    q75_time_to_max =0
    mean_half_max_duration = 0
    med_half_max_duration = 0
    min_half_max_duration = 0
    max_half_max_duration = 0
    q25_half_max_duration = 0
    q75_half_max_duration = 0
  }
  
  # spectral
  if (nrow(peaks) > 1) { 
    syl_ents = apply(peaks[,1:2], 1, function(peak) {
      calc_spectral_features(wavf, subregion = c(peak[1], peak[2]), overlap=0, wl=1024 )
      # wiener_entropy_var2(wavf, subregion = c(peak[1], peak[2]), overlap = 75)
    })
    freq_stats = apply(peaks[,1:2], 1, function(peak) {
      calc_freq_stats(wavf, q=c(.1,.5,.9, .95), subregion = c(peak[1], peak[2]), overlap=0, wl=1024 )
      # wiener_entropy_var2(wavf, subregion = c(peak[1], peak[2]), overlap = 75)
    })
    syl_ents = do.call(rbind, syl_ents)
    mean_freq = mean(syl_ents[,1])
    mean_freq_cv = mean(syl_ents[,2], na.rm=T)
    mean_syl_went = mean(syl_ents[,3])
    sd_syl_went = sd(syl_ents[,3], na.rm=T)
    mean_syl_wev = mean(syl_ents[,4])
    sd_syl_wev = sd(syl_ents[,4])
    
    freq_stats = apply(freq_stats, 1, mean, na.rm=T)
    freq_1 = freq_stats[1]
    freq_5 = freq_stats[2]
    freq_9 = freq_stats[3]
    freq_95 = freq_stats[4]
    
  } else {
    mean_freq = 0
    mean_freq_cv = 0
    mean_syl_went = 0
    sd_syl_went = 0
    mean_syl_wev = 0
    sd_syl_wev = 0
  }
  
  amp_ratio = calc_amp_ratio(amp, frate, peaks, max_gap=500, subsamp=subsamp)
  
  thresh_log=log(thresh+1)
  mean_amp_log=log(mean_amp+1)
  sd_amp_log=log(sd_amp+1)
  mean_syl_amps_log=log(mean_syl_amps+1)
  sd_syl_amps_log=log(sd_syl_amps+1)
  cv_gap_durations = unlist(sd_gap_durations / mean_gap_durations)
  cv_durations = unlist(sd_durations / mean_durations)
  cv_syl_amps = unlist(sd_syl_amps_log / mean_syl_amps_log)
  cv_gap_durations[is.nan(cv_gap_durations)] = 0
  cv_durations[is.nan(cv_durations)] = 0
  cv_syl_amps[is.nan(cv_syl_amps)] = 0
  
  cnames = c("thresh",
             "thresh_factor",
             "num_peaks",
             "rate",
             "total_duration",
             "mean_durations",
             "med_durations",
             "q75_durations",
             "sd_durations",
             "mean_gap_durations",
             "med_gap_durations",
             "sd_gap_durations",
             "mean_amp",
             "sd_amp",
             "mean_syl_amps",
             "sd_syl_amps",
             "mean_time_to_max",
             "med_time_to_max",
             "sd_time_to_max",
             "min_time_to_max",
             "max_time_to_max",
             "q25_time_to_max",
             "q75_time_to_max",
             "mean_time_from_max",
             "med_time_from_max",
             "sd_time_from_max",
             "mean_half_max_duration",
             "med_half_max_duration",
             "min_half_max_duration",
             "max_half_max_duration",
             "q25_half_max_duration",
             "q75_half_max_duration",
             "mean_syl_went",
             "sd_syl_went",
             "mean_syl_wev",
             "sd_syl_wev", 
             "thresh_log",
             "mean_amp_log",
             "sd_amp_log",
             "mean_syl_amps_log",
             "sd_syl_amps_log",
             "cv_gap_durations",
             "cv_durations",
             "cv_syl_amps",
             "amp_ratio",
             "mean_freq", 
             "mean_freq_cv",
             "freq_1",
             "freq_5",
             "freq_9",
             "freq_95")
  
  data = data.frame(mget(cnames), check.names = F)
  #   data = c(thresh=thresh, 
  #            thresh_factor=thresh_factor,
  #            num_peaks=nrow(peaks),
  #            rate=rate,
  #            total_duration=total_duration,
  #            mean_durations=mean_durations,
  #            med_durations=med_durations,
  #            q75_durations=q75_durations,
  #            sd_durations=sd_durations,
  #            mean_gap_durations=mean_gap_durations,
  #            med_gap_durations=med_gap_durations,
  #            sd_gap_durations=sd_gap_durations,
  #            mean_amp=mean_amp,
  #            sd_amp=sd_amp,
  #            mean_syl_amps=mean_syl_amps,
  #            sd_syl_amps=sd_syl_amps,
  #            mean_time_to_max=mean_time_to_max,
  #            med_time_to_max=med_time_to_max,
  #            sd_time_to_max=sd_time_to_max,
  #            min_time_to_max=min_time_to_max,
  #            max_time_to_max=max_time_to_max,
  #            q25_time_to_max=q25_time_to_max,
  #            q75_time_to_max=q75_time_to_max,
  #            mean_time_from_max=mean_time_from_max,
  #            med_time_from_max=med_time_from_max,
  #            sd_time_from_max=sd_time_from_max,
  #            mean_half_max_duration =  mean_half_max_duration,
  #            med_half_max_duration = med_half_max_duration,
  #            min_half_max_duration =min_half_max_duration,
  #            max_half_max_duration =max_half_max_duration,
  #            q25_half_max_duration = q25_half_max_duration,
  #            q75_half_max_duration =  q75_half_max_duration,
  #            mean_syl_went=mean_syl_went,
  #            sd_syl_went=sd_syl_went,
  #            mean_syl_wev=mean_syl_wev,
  #            sd_syl_wev=sd_syl_wev)
  return(data)
}

songfinder_features_clean = function(wav, rmnoise=TRUE, max_gap = 10, min_duration = 15, max_duration = 300 , thresh_method=mean_sd, thresh_factor=.25, min_num_peaks = 5, max_num_peaks=10) {
  #wavf = filtersong(wav)
  wavf = highpass_filter(wav, from = 500, wl = 1024, ovlp = 25)
  if (rmnoise)
    wavf = rmnoise(wavf, output="Wave")
  fs = wav@samp.rate
  
  # Identify peak regions
  peak_info = findpeaks_range(wavf, 
                              min_duration = min_duration, 
                              max_gap = max_gap, 
                              max_duration = max_duration, 
                              thresh_range = seq(0,2,.2), 
                              absolute = T)
  peaks = peak_info$peaks
  thresh_factor = peak_info$thresh_factor
  thresh = peak_info$thresh
  thresh_log = log(thresh + 1)
  num_peaks = nrow(peaks)
  
  #   output_cnames = c("mean_time_to_max",    
  #                      "q75_time_to_max.75.", 
  #                      "med_time_to_max",     
  #                      "med_durations",       
  #                      "cv_syl_amps",         
  #                      "mean_durations",      
  #                      "cv_gap_durations",   
  #                      "mean_syl_went",       
  #                      "thresh",              
  #                      "q25_time_to_max.25.", 
  #                      "sd_durations",        
  #                      "rate",                
  #                      "mean_gap_durations",  
  #                      "total_duration",     
  #                      "mean_syl_wev",        
  #                      "sd_time_to_max",      
  #                      "cv_durations",        
  #                      "sd_gap_durations")
  
  #   output_cnames = c("q75_time_to_max.75.",
  #                     "med_durations",
  #                     "mean_time_to_max",
  #                     "med_time_to_max",
  #                     "sd_gap_durations",
  #                     "mean_durations",     
  #                     "cv_syl_amps",
  #                     "cv_gap_durations",
  #                     "med_gap_durations",
  #                     "q25_time_to_max.25.",
  #                     "mean_syl_wev",
  #                     "mean_gap_durations", 
  #                     "min_time_to_max",
  #                     "med_time_from_max",
  #                     "sd_time_to_max",
  #                     "sd_amp",
  #                     "thresh",
  #                     "mean_syl_went")
  
  output_cnames = c("med_time_to_max", 
                    "q75_time_to_max", 
                    "q25_time_to_max",
                    "cv_gap_durations", 
                    "cv_durations",
                    "med_durations",
                    "med_gap_durations",
                    "mean_syl_wev", 
                    "mean_syl_went",
                    "num_peaks",
                    "cv_syl_amps",
                    "rate",
                    "thresh_log",
                    "mean_amp_log")  
  
  if (is.null(min_num_peaks)) min_num_peaks = 0
  if (is.null(max_num_peaks)) max_num_peaks = 1E5
  if ((nrow(peaks) < min_num_peaks) || (nrow(peaks) > max_num_peaks)) 
    return(data.frame(matrix(NA, ncol=length(output_cnames), dimnames=list("", output_cnames))))
  
  # number of peaks / second
  num_peaks = nrow(peaks)
  rate = num_peaks / (peaks[num_peaks,2] - peaks[1,1])
  
  # total duration of peak region.
  total_duration = peaks[nrow(peaks),2] - peaks[1,1]
  
  # Durations
  durations = peaks[,2] - peaks[,1]
  med_durations = median(durations)
  mean_durations = mean(durations)
  sd_durations = sd(durations)
  cv_durations = sd_durations / mean_durations
  
  # Gap durations
  gap_durations = vector("numeric", length = nrow(peaks) - 1)
  for (i in 1:(length(gap_durations)-1)) {
    gap_durations[i] = peaks[i+1,1] - peaks[i,2]
  }
  mean_gap_durations = mean(gap_durations)
  med_gap_durations = median(gap_durations)
  sd_gap_durations = sd(gap_durations)
  cv_gap_durations = sd_gap_durations / mean_gap_durations
  
  # amplitude properties
  amp = seewave::env(wavf, envt="abs", plot=F)
  mean_amp = mean(amp)
  mean_amp_log = log(mean_amp + 1)
  sd_amp = sd(amp)
  #max_amp = max(amp)
  #mean_max_amp = mean(amp) / max_amp
  syl_amps = apply(peaks[,1:2], 1, function(peak) {
    amps = amp[(peak[1]*fs):(peak[2]*fs)]
    time_to_max = which.max(amps)
    time_from_max = length(amps) - which.max(amps)
    c(mean_syl_amp=mean(amps), 
      sd_syl_amp=sd(amps),
      time_to_max = time_to_max,
      time_from_max = time_from_max)
  })
  mean_syl_amps = mean(syl_amps[1,])
  sd_syl_amps = mean(syl_amps[2,])
  cv_syl_amps = log(sd_syl_amps+1) / log(mean_syl_amps+1)
  mean_time_to_max = mean(syl_amps[3,])
  med_time_to_max = median(syl_amps[3,])
  q75_time_to_max = quantile(syl_amps[3,], .75)
  q25_time_to_max = quantile(syl_amps[3,], .25)
  min_time_to_max = min(syl_amps[3,])
  sd_time_to_max = sd(syl_amps[3,])
  #mean_time_from_max = mean(syl_amps[4,])
  #med_time_from_max = median(syl_amps[4,])
  #sd_time_from_max = sd(syl_amps[4,])
  
  # Spectral
  syl_ents = apply(peaks[,1:2], 1, function(peak) {
    wiener_entropy_var2(wavf, subregion = c(peak[1], peak[2]), overlap = 75)
  })
  syl_ents = do.call(rbind, syl_ents)
  mean_syl_went = mean(syl_ents[,1])
  #sd_syl_went = sd(syl_ents[,1])
  mean_syl_wev = mean(syl_ents[,2])
  #sd_syl_wev = sd(syl_ents[,2])
  
  data = as.data.frame(mget(output_cnames))
  colnames(data) = output_cnames
  return(data)
}


songfinder_psd = function(psd, 
                          band=c(5,8), 
                          wein_thresh=0.45) {
  wein = weiner_entropy_psd(psd, region=band)
  return(wein<wein_thresh)
}

compute_class = function(res, trues) {
  tp = vector("numeric", 4)
  for(i in 1:nrow(trues)) {
    m = match(trues[i,1], res[,1])
    if (trues[i,2] & res[m,2]) { #TP
      tp[1] = tp[1] + 1
    } else if (trues[i,2] & !res[m,2]) { #FN
      tp[2] = tp[2] + 1
    } else if (!trues[i,2] & res[m,2]) { #FP
      tp[3] = tp[3] + 1
    } else if (!trues[i,2] & !res[m,2]) { #TN
      tp[4] = tp[4] +1
    }
  }
  names(tp) = c("TP", "FN", "FP", "TN")
  return(tp)
}

######## SEQUENCING ########
parse_labels = function(mat, select=NULL) {
  d.labels = as.character(mat$labels)
  label_list = unlist(str_split(d.labels, ""))
  return(label_list)
  #d.labels = unlist(lapply(mat, function(x) as.character(x$labels)))
  #label_list = lapply(d.labels, function(x) unlist(str_split(x, "")))
}

syllable_transition = function(mat, select=NULL, deselect=NULL) {
  label_list = lapply(mat, parse_labels)
  
  return(syl_trans(label_list, select=select, deselect=deselect))
}
syllable_transition_df = function(df_list, select=NULL, deselect=NULL, max_gap_size) {
  
  df_list = lapply(df_list, function(x) {
    x$labels = as.character(x$labels)
    x
  })
  return(syl_trans_df(df_list, select=select, deselect=deselect, max_gap_size))
}

syl_trans = function(label_list, select=NULL, deselect=NULL) {
  label_list_cat = unlist(label_list)
  
  totals = table(label_list_cat)
  labels = NULL
  if (is.null(select)) {
    labels = names(totals)
  } else {
    labels = select
  }
  
  labels2 = c("start", labels, "end")
  counts = matrix(0, nrow=length(labels2), ncol=length(labels2), dimnames=list(labels2, labels2))
  
  for (i in 1:length(label_list)) {
    #print(paste("i", i, sep=""))
    curr = label_list[[i]]
    curr = curr[curr %in% labels]
    if (length(curr) <= 1) next
    counts["start", curr[1]] = counts["start", curr[1]] + 1
    for (j in 1:(length(curr)-1)) {
      counts[curr[j],curr[j+1]] = counts[curr[j],curr[j+1]] + 1
    }
    counts[curr[length(curr)], "end"] = counts[curr[length(curr)], "end"] + 1
  }
  return(counts)  
}

syl_trans_df = function(label_list, select=NULL, deselect=NULL, max_gap_size=NULL) {
  label_list_cat = bind_rows(label_list)
  
  totals = table(label_list_cat$labels)
  labels = NULL
  #if (is.null(select) | is.null(deselect)) {
  labels = names(totals)
  #} else {
  if (!is.null(select)) {
    labels = labels[labels %in% select]
  }
  
  if (!is.null(deselect)) {
    labels = labels[!(labels %in% deselect)]
  }
  
  
  
  labels2 = c("start", labels, "end")
  counts = matrix(0, nrow=length(labels2), ncol=length(labels2), dimnames=list(labels2, labels2))
  
  for (i in 1:length(label_list)) {
    #print(paste("i", i, sep=""))
    curr = label_list[[i]]
    curr = curr[curr$labels %in% labels,]
    if (nrow(curr) <= 1) next
    counts["start", curr$labels[1]] = counts["start", curr$labels[1]] + 1
    for (j in 1:(nrow(curr)-1)) {
      if (!is.null(max_gap_size)) 
        if ((curr$onsets[j+1] - curr$offsets[j]) > max_gap_size)
          next
      counts[curr$labels[j],curr$labels[j+1]] = counts[curr$labels[j],curr$labels[j+1]] + 1
    }
    counts[curr$labels[nrow(curr)], "end"] = counts[curr$labels[nrow(curr)], "end"] + 1
  }
  return(counts)  
}

process_syllable_matrix_mat = function(mats, select=NULL, norm="total") {
  syl = syllable_transition(mats, select)
  if (!is.null(select)) {
    syl = syl[rownames(syl) %in% select, colnames(syl) %in% select]
  }
  
  if (norm=="total") {
    syl = syl / sum(c(syl))
  } else if (norm=="row") {
    syl = sweep(syl, MARGIN = 1, STATS = rowSums(syl), "/")
  } else if (norm=="none") {
    syl = syl
  }
  m = melt(syl)
  colnames(m)[1:2] = c("From", "To")
  m = m %>% filter(From!="end") %>% filter(To!="start") 
  #print(m)
  return(m)
}

process_syllable_matrix = function(info, select=NULL, deselect=NULL, norm="total", max_gap_size=NULL, randomize=F) {
  # print(info$mat)
  d = lapply(info$mat, function(x) {
    m = readMat(x)
    if (is.null(m[[3]])) 
      return(NULL)
    return(m)})
  if (is.null(d[[1]]))
    return(data.frame(From="-", To="-", value=0))
  
  
  d = lapply(d, function(m) {
    d1 = data.frame(onsets=m$onsets, offsets=m$offsets, labels=parse_labels(m))
    d1 = d1 %>% filter(!(labels %in% deselect))
    if (randomize)
      d1$labels = sample(d1$labels, length(d1$labels), replace=F)
    d1
  })
  
  
  
  syl = syllable_transition_df(d, select, deselect, max_gap_size=max_gap_size)
  
  #if (!is.null(max_gap_size)) {
  
  #} else {
  #  syl = syllable_transition(d, select, deselect)
  #}  
  
  if (!is.null(select)) {
    syl = syl[rownames(syl) %in% select, colnames(syl) %in% select]
  }
  
  if (norm=="total") {
    syl = syl / sum(c(syl))
  } else if (norm=="row") {
    syl = sweep(syl, MARGIN = 1, STATS = rowSums(syl), "/")
  } else if (norm=="none") {
    syl = syl
  }
  m = melt(syl)
  colnames(m)[1:2] = c("From", "To")
  m = m %>% filter(From!="end") %>% filter(To!="start") 
  #print(m)
  return(m)
}



#' Calculate repeat lengths by bout
#' @param data, info file as returned from load_mat_info
#' @return data.frame
#'   \item{lengths}{length as from rle}
#'   \item{values}{repeated syllable}
#'   \item{mat}{original .not.mat file}
repeat_lengths = function(data, select=NULL) {
  d = foreach(row=isplitRows(data, chunkSize=1)) %dopar% {
    mat = readMat(row["mat"])
    labels = as.character(mat$labels)
    labels1 = unlist(str_split(labels, ""))
    rles = rle(labels1)
    rd = data.frame(lengths=rles$lengths, values=rles$values, mat=row["mat"])
    if (is.null(select))
      select = unique(labels1)
    rd = rd[rd$values %in% select,]
    rd
  }
  do.call(rbind, d)
}

mean_transition_entropy_long = function(long_data) {
  mean_transition_entropy(acast(long_data, From~To, fun.aggregate = sum))
}

mean_transition_entropy = function(mat) {
  #mat = mat / sum(c(mat))
  #res = apply(mat, 1, calc_entropy)
  #data.frame(transition_entropy = sum(unlist(res)))
  data.frame(transition_entropy = calc_entropy(c(mat)))
}




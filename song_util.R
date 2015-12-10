library(seewave)
library(tuneR)
library(signal)
library(manipulate)
library(reshape2)
library(ggplot2)
library(doMC)
library(foreach)
library(itertools)
library(R.matlab)
library(stringr)
library(mclust)
library(dplyr)
library(gridExtra)
library(parallel)
library(matlab)
registerDoMC(cores=2)

source("/home/brad/src/songanalysis/threshold.r")
source("/home/brad/src/songanalysis/clustering.R")

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

######## FILTERING/EDITING ########

#' taken from seewave::ffilter 
highpass_filter = function(wav, from=1800, wl=1024, ovlp=75, wn="hanning", fftw=T) {
  input = inputw(wave = wav, f = f)
  wave = input$w
  f = input$f
  bit = input$bit
  #rm(input)
  n = nrow(wave)
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
#' High pass filter Wave object at 1800 Hz
#' @param wav, the input Wave object
filtersong = function(wav, band=500) {
  return(ffilter(wav, from=band, output="Wave", fftw=T, ))
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

######## PROCESSING ########

parse_song = function(file, cluster=F, plot=F) {
  min_dur = 15
  max_gap = 10
  wav = readWave(file)
  fname = basename(str_replace(file, ".wav", ""))
  wavf = filtersong(wav)
  
  #thresh = threshold_auto(wavf, otsu)
  thresh = threshold_auto(wavf, mean_sd)
 
  syls = findpeaks_abs(wavf, min_duration = min_dur, max_gap = max_gap, thresh=thresh)
  if (plot) {
  par(mfrow=c(1,2))
  env(wavf, envt="abs", from=1, to=2)
  abline(h=thresh, col=2)
  plot_peaks(syls, col=4)
  }
  syls$id = paste(fname, 1:nrow(syls), sep="-")
  song = list(fname = fname, wav=wavf, syllable=syls, thresh=thresh, min_dur = min_dur, max_gap = max_gap)
  
  if (cluster) {
    syl_data = extract_features2(song)
    mod = mclust_syllable_data(syl_data, plot=T)
    song$syllable$called = mod$classification
    return(list(song=song, model=mod))
  }
  
  return(song)
}

parse_song_batch = function(wdir, range_to_test=c(4:9), cluster=TRUE) {
  print("Loading data")
  files = list.files(wdir, pattern=".wav$", full.names = T)

  songd = foreach(i=1:length(files)) %dopar% {
    # file = row[1,"wav"]
    file = files[i]
    print(i)
    song = parse_song(file)
    syl_data = extract_features2(song)  
    return(list(song=song, syl_data=syl_data))
  }
  
  #return(songd)
  names(songd) = files
  print("Rearranging data...")
  
  wavs = lapply(songd, function(x) x$song$wav)
  
  syl_data = lapply(songd, function(x) x$syl_data)
  syl_data = do.call("rbind", syl_data)
 
  syls = lapply(songd, function(x) x$song$syllable)
  syls = do.call("rbind", syls)
  #return(list(syl_data=syl_data, syls=syls))
  if (cluster) {
    print("Clustering syllables...")
    syl.mc = mclust_syllable_data(syls, syl_data, range_to_test = range_to_test, plot=T)
    #return(syl.mc)
    
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


parse_song_batch2 = function(info, cluster=TRUE) {
  print("Loading data")
  #songd = mclapply(1:nrow(info), function(i) {
  songd = foreach(i=1:nrow(info)) %dopar% {
   # file = row[1,"wav"]
    file = info[i,"wav"]
    print(i)
    wav = readWave(file)
    fname = basename(str_replace(file, ".wav", ""))
    wavf = filtersong(wav)
    
    thresh = threshold_auto(wavf, otsu)
    syls = findpeaks_abs(wavf, min_duration = 15, max_gap = 5, thresh=thresh)
    
    syls$id = paste(fname, 1:nrow(syls), sep="-")
    song = list(fname = fname, wav=wavf, syllable=syls)
    syl_data = extract_features2(song)  
    return(list(wav = wavf, syls = syls, syl_data=syl_data))
  #})
  #}, mc.cores=10)
  }
  return(songd)
  print("Rearranging data...")
  syl_data = lapply(songd, function(x) x$syl_data)
  syl_data = do.call("rbind", syl_data)
  
  syls = lapply(songd, function(x) x$song.syllable)
  syls = do.call("rbind", syls)
  par(mfrow=c(1,1))

  if (cluster) {
    print("Clustering syllables...")
    mod = mclust_syllable_data(syl_data, range_to_test = c(4:12), plot=T)
    songd$syls$called = mod$classification
    return(list(song=songd, model=mod))
  } else {
    return(syl_data)
  }
   
}

#' Read in info for wav and processed mat files
#' @param wdir, list or character of working directory containing .wav and .not.mat
#' @param file_ex, pattern searched for by list.files
#' @return data.frame containing file info 
load_mat_info = function(wdir, file_ex = ".mat") {
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

load_mat = function(file) {
  mat = readMat(file)
  labels = unlist(str_split(mat$labels, ""))
  out = data.frame(onsets=mat$onsets, offsets=mat$offsets, labels=labels)
  out$id = paste(basename(str_replace(file, ".wav.not.mat", "")), 1:nrow(out), sep="-")
  out[,1:2] = out[,1:2] / 1000
  out
}

load_mat_batch = function(info) {
  mat = foreach(row=isplitRows(info, chunkSize=1), .combine="rbind") %do% {
    a = load_mat(row["mat"])
    #a$id = paste(basename(str_replace(row["wav"], ".wav", "")), 1:nrow(a), sep="-")
    a
  }
  #mat$id = 1:nrow(mat)
  return(mat)
}

song_to_mat = function(song) {
  #fname = song$fname
  wav = song$wav
 # syls = song$syllable
  ### TODO change NA syllableto '-'
  mat = list(Fs=wav@samp.rate, 
             fname=song$fname,
             labels=song$syllable$called,
             onsets=song$syllable$ins,
             offsets=song$syllable$outs,
             min.int=song$max_gap,
             min.dur=song$min_dur,
             threshold=song$thresh,
             sm.win=2)
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
      psd = spec(wavf,  fftw = T, PSD = T, plot=F, at=as.numeric(row["onsets"]) + step,)
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

extract_features2 = function(song) {
  wav = song$wav
  peaks = song$syllable
  peaks = peaks[as.numeric(peaks$ins)>.005,]
  
  #offset = .010
  peaks$mids = (peaks[,"outs"] + peaks[,"ins"]) / 2
  peaks$ins = (peaks[,"ins"] + peaks[,"mids"]) / 2
  peaks$outs = (peaks[,"outs"] + peaks[,"mids"]) / 2
  
  points = matrix(0, nrow=nrow(peaks), ncol=3)
#   points[,1] = peaks[,"mids"] - 0.02
#   points[,2] = peaks[,"mids"] - 0.01
#   points[,3] = peaks[,"mids"]
#   points[,4] = peaks[,"mids"] + 0.01
#   points[,5] = peaks[,"mids"] + 0.02

   points[,1] = peaks[,"mids"] - 0.01
   points[,2] = peaks[,"mids"]
   points[,3] = peaks[,"mids"] + 0.01
  
  peaks$ins.frames = peaks[,"ins"] * wav@samp.rate
  peaks$outs.frames = peaks[,"outs"] * wav@samp.rate
 # peaks[,1:2] = peaks[,1:2] + offset

#   amp = env(wavf, envt="abs", plot=F, norm = T)
#   amp = amp / (mean(amp) + sd(amp))
#   amp[amp>1] = 1

  
  res = 512 / wav@samp.rate
  freq_limits = c(0, 12)
  freq_range = seq(0, wav@samp.rate / 2000, length.out=256)
  ind_range = which(freq_range > freq_limits[1] & freq_range < freq_limits[2])
  psd1 = NULL
  
  weights_init_x = c(1, length(ind_range))
  weights_init_y = c(1, 10)
  weights_target = 1:length(ind_range)
  weights = approx(weights_init_x, weights_init_y, weights_target)$y
  
  
  data = foreach(row=isplitRows(points, chunkSize=1), .inorder = T, .combine="rbind") %do% {
    #tomid = row[1,"mids"] - row[1,"ins"]
    
    #if (tomid < res) {
    #  psd1 = c(meanspec(wav, wl=512, fftw=T, PSD=T, norm = F, plot=F,  at=as.numeric(row["in"]))[ind_range,"y"],
    #           meanspec(wav, wl=512, fftw=T, PSD=T, norm = F, plot=F,  at=as.numeric(row["mids"]))[ind_range,"y"])
               #meanspec(wav, wl=512, fftw=T, PSD=T, norm = F, plot=F,  at=as.numeric(row["outs"]))[ind_range,"y"])
    #} else {
     # psd1 = c(meanspec(wav, wl=512, fftw=T, PSD=T, norm=F,  plot=F, from=as.numeric(row["ins"]), to=as.numeric(row["mids"]))[ind_range,"y"],
    #         meanspec(wav, wl=512, fftw=T, PSD=T, norm=F, plot=F, from=as.numeric(row["mids"]), to=as.numeric(row["outs"]))[ind_range,"y"])
    psd1 = foreach(x= isplitCols(row, chunkSize=1), .combine="c") %do% {
      d = seewave::spec(wav, wl=512, fftw=T, PSD=F, norm=F,  plot=F, at=x)[ind_range,"y"]
      d * weights
      #d
    } 
      #c(spec(wav, wl=512, fftw=T, PSD=T, norm=F,  plot=F, at=as.numeric(row["ins"]), to=as.numeric(row["mids"]))[ind_range,"y"],
      #               spec(wav, wl=512, fftw=T, PSD=T, norm=F, plot=F, at=as.numeric(row["mids"]), to=as.numeric(row["outs"]))[ind_range,"y"])
      #psd1 = c(spec(wav, wl=512, fftw=T, PSD=T, norm = F, plot=F, at=as.numeric(row["ins"]))[ind_range,"y"],
      #         spec(wav, wl=512, fftw=T, PSD=T, norm = F, plot=F, at=as.numeric(row["mids"]))[ind_range,"y"])
               #spec(wav, wl=512, fftw=T, PSD=T, norm = F, plot=F, at=as.numeric(row["outs"]))[ind_range,"y"])
      #psd1 = spectro(wav, wl=512, fftw=T, PSD=T, norm = F, plot=F, from=as.numeric(row["ins"]), to=as.numeric(row("outs")))
    #}
#     amp1 = amp[as.numeric(row["ins.frames"]):as.numeric(row["outs.frames"])]
#     a = round(seq(1, length(amp1), length.out=200))
#     cuts = matrix(0, nrow=length(a), ncol=2)
#     cuts[,1] = a
#     for (i in 1:(nrow(cuts)-1)) {
#       cuts[i,2] = cuts[i+1, 1] - 1
#     }
#   
#     amp1 = apply(cuts, 1, function(x) {
#       mean(amp1[x[1]:x[2]])
#     }) 
    #amp1 = amp1 / max(amp1)
    
    #psd1 = c(psd1, amp1)
    #psd1 = predict(loess(psd1~c(1:length(psd1)), span=.05))
    #psd1 = na.omit(stats::filter(psd1, rep(1, times=20) / 20))
    #max_psd = max(psd1)
    #print(max_psd)
    #perc = .1
    #psd1 = psd1 - (max_psd * perc)
    #psd1[psd1<1] = 1
    #psd1[psd1<1] = 1
    #return()
    return(psd1)
    #psd1 = meanspec(wavf, fftw = T, from = as.numeric(row["ins"]), to = as.numeric(row["outs"]))
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
    plot_multiple_syl(wav, peaks_s1[[d]][,1:2], title=names(peaks_s1)[d] )
  })
  
  ncol = 3
  nrow = ceiling(length(plots) / ncol)
  plot.list = c(plots, list(ncol=ncol, nrow=nrow))
  do.call(grid.arrange, plot.list)
}

plot_spectro = function(wav, subregion=NULL, wl=512, overlap=50, labels=NULL, label_name = "called") {
  require(matlab)
  spectrogram = spectro(wav, plot = FALSE, tlim=subregion, fftw=T)
  frequency = rep(spectrogram$freq, times = ncol(spectrogram$amp))
  time = rep(spectrogram$time, each = nrow(spectrogram$amp))
  amplitude = as.vector(spectrogram$amp)
  df = data.frame(time, frequency, amplitude)
  gg = ggplot(df, aes_string(x = "time", 
                            y = "frequency",
                            z = "amplitude"))
  gg = gg + stat_contour(aes(fill=..level..), geom="polygon", binwidth=1)
  gg = gg + scale_fill_gradientn("Amplitude, (db)",
                                 colours=jet.colors(7), 
                                 breaks=c(-40, -15, 0),
                                 limits=c(-40, -15),
                                 na.value="transparent")
  gg = gg + labs(x="Time", y="Frequency") + ylim(0,10) + xlim(subregion)
  
  if (!is.null(labels)) {
    labels$time = (labels$outs + labels$ins) / 2
    labels = cbind(labels, data.frame(frequency=0, amplitude=0))
    gg = gg + geom_text(data=labels, aes_string(label=label_name))
  }
  #print(gg)
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
  return(gg)
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
  wav_env = env(wav, envt = "abs", plot=F)
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


#' Find frequency peaks in PSD
#' @param wav, the input Wave object
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

#' Find temporal distance between peaks
#' @param peaks, data frame as outputted by findeaks, col1 is peak starts, col2 is peak stops
#' @return vector of interpeak distances
peak_dist = function(peaks) {
  if (nrow(peaks)<2) return(NA)
  dists = vector("numeric", nrow(peaks)-1)
  for (i in 1:(nrow(peaks)-1)) {
    dists[i] = peaks[i+1,1] - peaks[i,2]
  }
  return(dists)
}

interpeak_midpoint = function(peaks) {
  if (nrow(peaks)<2) return(NA)
  dists = vector("numeric", nrow(peaks)-1)
  for (i in 1:(nrow(peaks)-1)) {
    dists[i] = sum(peaks[i+1,1] + peaks[i,2]) / 2
  }
  return(dists)
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

######## WIENER ENTROPY ########
#' Calculate Wiener entropy of given frequency band in Wave object
#' @param wav, Wave object
#' @param band, selected frequency band
#' @param subregion, vector giving start and end of wav subregion to analyze, in seconds
#' @references https://en.wikipedia.org/wiki/Spectral_flatness
#' @return scalar, calculated Weiner entropy
wiener_entropy = function(wav, band=NULL, subregion=NULL) {
  psd = NULL
  if (!is.null(subregion)) {
    psd = seewave::spec(wav, wl=512, plot=F, PSD=T, from=subregion[1], to=subregion[2], fftw=T)
  } else {
    psd = spec(wav, wl=512, plot=F, PSD=T)
  }
  
  psd[,1] = 1000 * psd[,1] # convert to Hertz from kHertz
  if (is.null(band)) band = c(1, nrow(psd))
  psd1 = psd[psd[,1]>band[1] & psd[,1]<band[2],]
  res = exp(mean(log(psd1[,2]))) / mean(psd1[,2])
  return(res)
}

wiener_entropy_var = function(wav, window=16, overlap=50, band=NULL, subregion=NULL) {
  subregion = subregion * 1000
  steps = seq(0, subregion[2]-subregion[1], window * overlap / 100)
  wav_duration = duration(wav)
  res = sapply(steps, function(step) {
    subregion2 = subregion 
    subregion2[1] = subregion[1] + step
    subregion2[2] = subregion2[1] + window
    subregion2 = subregion2 / 1000
    if (subregion2[2] > wav_duration) return(NA)
    wiener_entropy(wav, band=band, subregion=subregion2)
  })
  res = na.omit(res)
  var(res)
}

wiener_entropy_var2 = function(wav,  window=.016, overlap=50, band=NULL, subregion=NULL) {
  
  steps = seq(0, subregion[2]-subregion[1], window * overlap / 100)
  #print(steps)
  wav_duration = seewave::duration(wav)
  
  
  #subregion2 = matrix(0, nrow=length(steps), ncol=2)
  #subregion2[,1] = subregion[1] + step
  #subregion2[,2] = subregion2[,1] + window
  subregion2 = subregion
  #subregion2 = subregion2[subregion2[,2]<wav_duration,]
  res = sapply(steps, function(step) {
    subregion2[1] = subregion[1] + step
    subregion2[2] = subregion2[1] + window
    if (subregion2[2] > wav_duration) return(NA)
    wiener_entropy(wav, band=band, subregion=subregion2)
  })
  res = na.omit(res)
  data.frame(went_mean = log2(mean(res)), went_var = var(res))
}

wiener_entropy_var2.1 = function(wav,  window=.016, overlap=50, band=NULL, subregion=NULL) {
  
  steps = seq(0, subregion[2]-subregion[1], window * overlap / 100)
  #print(steps)
  wav_duration = duration(wav)
  #subregion2 = subregion
  subregion2 = matrix(0, nrow=length(steps), ncol=2)
  subregion2[,1] = subregion[1] + steps
  subregion2[,2] = subregion2[,1] + window
  #print(subregion2)
  #subregion2 = subregion2 / 1000
  subregion2 = subregion2[subregion2[,2]<wav_duration,]
  res = apply(subregion2, 1, function(row) {
    
    wiener_entropy(wav, band=band, subregion=row)
    #subregion2[1] = subregion[1] + step
    #subregion2[2] = subregion2[1] + window
    #if (subregion2[2] > wav_duration) return(NA)
    
  })
  res = na.omit(res)
  data.frame(went_mean = log2(mean(res)), went_var = var(res))
}

#' Calculate Weiner entropy of given frequency band in PSD
#' @param psd, PSD as output from spec(PSD=T)
#' @param band, selected frequency band
#' @references https://en.wikipedia.org/wiki/Spectral_flatness
#' @return scalar, calculated Weiner entropy
weiner_entropy_psd = function(psd, region=c(6,8)) {
  psd1 = psd[psd[,1]>region[1] & psd[,1]<region[2],]
  res = exp(mean(log(psd1[,2]))) / mean(psd1[,2])
  return(res)
}

#' Calculate Weiner entropy and Wiener entropy variance given info file
#' as returned by load_mat_file
#' @param info, info file as returned by load_mat_file
#' @param band, selected frequency band
#' @return data.frame, given info file supplmented with calculated values
#'   \item{went}{wiener entropy}
#'   \item{wev}{variance of wiener entropy across syllable}
wiener_stats = function(info, band=c(2000, 10000)) {
  d = foreach(row=isplitRows(info, chunkSize=1)) %dopar% {
   # print(row["mat"])
    out = load_mat(row["mat"])
    wav = readWave(row[1, "wav"])
    #out = transform(out, onsets=onsets * 1000, offsets = offsets * 1000)
    out = out %>% rowwise() %>% do({
      cbind(., wiener_entropy_var2(wav, band=band, subregion=c(.$onsets, .$offsets)))
    })
    out$mat = row[1, "mat"]
    out
  }
  do.call(rbind, d)
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

songfinder2 = function(wav, max_gap = 10, min_duration = 15, max_duration = 300 , thresh=0.2, min_num_peaks = 5, max_num_peaks=10) {
  wavf = filtersong(wav)
  #wavf = highpass_filter(wav)
  peaks = findpeaks_abs(wavf, min_duration=min_duration, max_gap=max_gap, max_duration=max_duration, thresh=thresh)
  #return(peaks)
  if (nrow(peaks) == 0) return(FALSE)
  rate = nrow(peaks) / (peaks[nrow(peaks),2] - peaks[1,1])
  total_duration = peaks[nrow(peaks),2] - peaks[1,1]
  #print(rate)
  return((total_duration > 1) & (rate > min_num_peaks ) & (rate < max_num_peaks ))
}

songfinder_psd = function(psd, 
                          band=c(5,8), 
                          wein_thresh=0.45) {
  wein = weiner_entropy_psd(psd, region=band)
  return(wein<wein_thresh)
}

######## Cepstrum ########
format_cepstrum = function(data) {
  wav_cep_form = melt(data$amp)
  colnames(wav_cep_form) = c("time", "quef", "amp")
  wav_cep_form$time = data$time
  wav_cep_form$quef = rep(data$quef, each=length(data$time))
  return(wav_cep_form)
}

#' Calcualte dynamic cepstrum of wav file 
#' @param Wave object
#' @return cepstrum as dataframe: 1:time, 2:quef, 3:amp
compute_cepstrum_time = function(wav) {
  wav_cep = cepstro(wav, wl=256, ovlp=50)
  return(format_cepstrum(wav_cep))
}

#' Calculate derivative of dynamic cepstrum using least-square approximation
#' @param Wave object
#' @param w, window size
#' @return delta cepstrum as data.frame: 1:time, 2:quef, 3: local derivative
delta_cepstrum = function(wav, w=10, region=NULL) {
  wav_cep = NULL
  if (is.null(region)) {
    wav_cep = cepstro(wav, wl=512, ovlp=50, collevels = seq(0, .01, .001), ylim=c(0, 1) )
  } else {
    wav_cep = cepstro(wav, wl=512, ovlp=50, from=region[1], to=region[2])
  }
  wav_cepf = format_cepstrum(wav_cep)
  wav_cep_mat = wav_cep$amp
  delta_cep = matrix(nrow=(nrow(wav_cep_mat) - (2*w)), ncol=ncol(wav_cep_mat))
  inds = seq(-1*w, w)
  
  #### loop through rows of wav_cep
  norm_factor = inds %*% inds
  for (i in 1:nrow(delta_cep)) {
      delta_cep[i,] = t(inds) %*% wav_cep_mat[(i+w + inds),]
  }
  delta_cep = delta_cep / norm_factor[1,1]
  
  delta_cep_form = melt(delta_cep)
  colnames(delta_cep_form) = c("time", "quef", "amp")
  delta_cep_form$time = wav_cep$time[(w+1):(nrow(wav_cep_mat) - w)]
  delta_cep_form$quef = rep(wav_cep$quef, each=nrow(delta_cep))
  return(delta_cep_form)
}

ggplot_cepstrum = function(data) {
  limits = c(min(data$amp), quantile(data$amp, probs=.99))
  gg = ggplot(data, aes(time, quef, fill=amp)) +  geom_tile() + scale_fill_gradient(limits=limits)
  gg = gg + ylim(0, .7)
  return(gg) 
}

######## Fundamental frequency ########
calc_ff = function(wav, peaks) {
  
  psd = apply(peaks, 1, function(x) spec(wav, wl = 512, PSD=T, from=x[1], to=x[2], plot=F))
  psd_peaks = lapply(psd, function(data) {
    total_power = sum(data[,2])
    d = findpeaks_freq(data, min_size=0)
    d %>% rowwise() %>% summarize(ins=ins, 
                                  outs=outs, 
                                  prop=sum(data[data[,1]>ins & data[,1]<outs,2])/total_power)
  })
  return(psd_peaks)
}

calc_ff = function(wav, mat) {
  
  #psd = apply(mat, 1, function(x) spec(wav, PSD=T, from=as.numeric(x["onsets"]), to=as.numeric(x["offsets"]), plot=F))
  #if (class(psd) == "matrix") psd = list(cbind(1:nrow(psd), psd))
  freqs = unlist(lapply(psd, function(data) {
    d = findpeaks_freq(data, min_value=2, max_gap=1000, min_size=0, thresh=.1)
    maxes = apply(d, 1, function(x) max(psd[x[1]:x[2]]))
    d = d[which.max(maxes),]
    #d_mean = apply(d, 1, mean)
    #d_mean = data.frame(ind=1:length(d_mean), d_mean=d_mean)
    #d_freq = apply(d_mean, 1, function(x) x[2]/x[1])
    #mean(d_freq)
  }))
  
  #if (nrow(mat) == 1) {
  #  return(c(mat, freqs))
  #} else {
  return(cbind(mat, freqs))
  #}
}

#' Calculate fundamental frequency by averaging harmonics
#' @param wav, input Wave file
#' @param mat, data.frame containing onsets, offsets, and label as generated by load_mat
#' @return mat file supplemented with calculated fundamental frequency
calc_ff2 = function(wav, mat) {
  
  psd = apply(mat, 1, function(x) spec(wav, PSD=T, from=as.numeric(x["onsets"]), to=as.numeric(x["offsets"]), plot=F))
  if (class(psd) == "matrix") psd = list(cbind(1:nrow(psd), psd))
  freqs = unlist(lapply(psd, function(data) {
    
    d = findpeaks_freq(data, min_value=2, max_gap=1000, min_size=0, thresh=.1)
    d_mean = apply(d, 1, mean)
    d_mean = data.frame(ind=1:length(d_mean), d_mean=d_mean)
    d_freq = apply(d_mean, 1, function(x) x[2]/x[1])
    mean(d_freq)
  }))
  
  return(cbind(mat, freqs))

}

calc_ff2_batch = function(info, label="a") {
  d = foreach(row=isplitRows(info, chunkSize=1), .combine="rbind" ) %dopar% {
    wav = readWave(row[1, "wav"])
    wavf = filtersong(wav)
    mat = load_mat(row[1, "mat"])
    
    mat = mat %>% filter(labels==label)
    if (nrow(mat) > 0) {
      return(cbind(calc_ff2(wavf, mat), mat = row[1,"mat"]))
    } else {
      return(NULL)
    }
  }
  d
}
calc_freq_rolling_batch = function(info, label="a", offset=8) {
  d = foreach(row=isplitRows(info, chunkSize=1), .combine="rbind" ) %do% {
    print(row[1,"wav"])
    wav = readWave(row[1, "wav"])
    wavf = filtersong(wav)
    mat = load_mat(row[1, "mat"])
    
    
    mat = mat %>% filter(labels==label)
    #mat = mat[mat$offsets - mat$onsets]
    if (nrow(mat) > 0) {
      res = mat %>% rowwise() %>% do({calc_freq_rolling(wavf, offset=offset, subregion=unlist(.[1:2]))})
      return(cbind(mat, res, mat = row[1,"mat"]))
    } else {
      return(NULL)
    }
  }
  d
}
calc_freq_rolling = function(wav, offset=4, step = 4, duration = 8, subregion=NULL, fraction_max=.05) {
  wl = 512
  res = 512 / wav@samp.rate
  offset = offset / 1000
  null_df = data.frame(freq = NA)
  if (!is.null(subregion)) {
    from = subregion[1] + offset
    to = subregion[2] - offset
    if ((to - from) < res ) return(null_df)
    psd = spectro(wav, wl=512, norm=F, fftw = T, dB = NULL, plot=F, PSD=T, tlim=c(from, to), ovlp=50)
  } else {
    stop("Please specify subregion!")
  }
  
  max_overall = max(psd[[3]])
  
  peaks = apply(psd[[3]], 2, function(x) {
    f = findpeaks_freq(cbind(psd[[2]], x), 
                       min_value=2, 
                       max_gap=1, 
                       min_size=0, 
                       thresh=max_overall*fraction_max, 
                       absolute=T)
    if (nrow(f)==0) return(NA)
  
    power = apply(f, 1, function(peak) {
      inds = c(which(psd[[2]]==peak[1]), which(psd[[2]]==peak[2])) # select frequency corresponding to power peak
      mean(x[inds])
    })
    
    fmean = cbind(freq = rowMeans(f), power) 
    return(fmean[which.min(fmean[,1]),])
    })
  
  if (length(na.omit(peaks)) == 0) return(null_df)
  if (class(peaks) == "list") {
    peaks = na.omit(do.call(rbind, peaks))
  } else {
    peaks = t(peaks)
    #print("here")
  }
  total_power = sum(peaks[,"power"])
  return( data.frame(freq = as.numeric(t(peaks[,"freq"]) %*% peaks[,"power"] / total_power )))
}

calc_freq = function(wav, offset=10, duration=8, band=NULL, subregion=NULL) {
  psd = NULL
  if (!is.null(subregion)) {
    from = subregion[1] + offset
    to = from + duration
    psd = spec(wav, wl=512, norm=F, plot=F, PSD=T, from=(from / 1000), to=(to / 1000))
  } else {
    stop("Please specify subregion!")
  }
 
  psd[,1] = 1000 * psd[,1] # convert to Hertz from kHertz
  if (is.null(band)) band = c(1, nrow(psd))
  psd1 = psd[psd[,1]>band[1] & psd[,1]<band[2],]
  colnames(psd1) = c("freq", "power")
  psd1 = as.data.frame(psd1)
  psd1$onsets = subregion[1]
  psd1$offsets = subregion[2]
  #return(as.data.frame(psd1))
  return(psd1[which.max(psd1[,2]),1])
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
syllable_transition = function(mat, select=NULL) {
  d.labels = unlist(lapply(mat, function(x) as.character(x$labels)))
  label_list = lapply(d.labels, function(x) unlist(str_split(x, "")))
  return(syl_trans(label_list, select=select))
}
syllable_transition_df = function(df_list) {
    label_list = lapply(df_list, function(x) as.character(x$labels))
  return(syl_trans(label_list))
#   label_data = df$labels
#   totals = table(label_data)
#   labels = names(totals)
#   
#   labels2 = c("start", labels, "end")
#   counts = matrix(0, nrow=length(labels2), ncol=length(labels2), dimnames=list(labels2, labels2))
#   
#   for (i in 1:length(label_list)) {
#     curr = label_list[[i]]
#     counts["start", curr[1]] = counts["start", curr[1]] + 1
#     for (j in 1:(length(curr)-1)) {
#       counts[curr[j],curr[j+1]] = counts[curr[j],curr[j+1]] + 1
#     }
#     counts[curr[length(curr)], "end"] = counts[curr[length(curr)], "end"] + 1
#   }
#   return(counts)
  
  
}

syl_trans = function(label_list, select) {
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


process_syllable_matrix = function(info, select=NULL, norm="total") {
    d = lapply(info$mat, readMat)
    syl = syllable_transition(d, select)
    if (!is.null(select)) {
      syl = syl[rownames(syl) %in% select, colnames(syl) %in% select]
    }
    if (norm=="total") {
      syl = syl / sum(c(syl))
    } else if (norm=="row") {
      syl = sweep(syl, MARGIN = 1, STATS = rowSums(syl), "/")
    }
    m = melt(syl)
    colnames(m)[1:2] = c("From", "To")
    m = m %>% filter(From!="end") %>% filter(To!="start") 
    return(m)
}



#' Calculate repeat lengths by bout
#' @param data, info file as returned from load_mat_info
#' @return data.frame
#'   \item{lengths}{length as from rle}
#'   \item{values}{repeated syllable}
#'   \item{mat}{original .not.mat file}
repeat_lengths = function(data) {
  require(foreach)
  require(itertools)
  d = foreach(row=isplitRows(data, chunkSize=1)) %dopar% {
    mat = readMat(row["mat"])
    labels = as.character(mat$labels)
    labels1 = unlist(str_split(labels, ""))
    rles = rle(labels1)
    rd = data.frame(lengths=rles$lengths, values=rles$values, mat=row["mat"])
    rd = rd[rd$values %in% select,]
    rd
  }
  do.call(rbind, d)
}

mean_transition_entropy_long = function(long_data) {
  mean_transition_entropy(acast(long_data, From~To, fun.aggregate = sum))
}

mean_transition_entropy = function(mat) {
  mat = mat / sum(c(mat))
  res = apply(mat, 1, function(x) {
    y = x[x>0]
    sum(y * log2(y))
  })
  data.frame(transition_entropy = -1 * sum(res))
}



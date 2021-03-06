
find_closest = function(target, range) {
  
  if (is.unsorted(range))
    stop("range must be sorted")
  
  if (length(range)==1) {
    return(range)
  }
  val = median(range)
  mid = round(length(range) / 2)
  if (target <= val) {
    find_closest(target, range[1:mid])
  } else {
    find_closest(target, range[(mid+1):length(range)])
  }
}
# 
# find_freqband = function(psd, band) {
#   ind1 = which.min(dist(band[1], psd[,1]))
#   ind2 = which.min(dist(band[2], psd[,1]))
#   return(ind1:ind2)
#   #c(ind1, ind2)
# }

find_closest_inds = function(targets, range) {
  ind1 = which(range==find_closest(targets[1], range))
  ind2 = which(range==find_closest(targets[2], range))
  return(ind1:ind2)
  #c(ind1, ind2)
}


######### SPECTRAL ##########
calc_psd = function(wav, wl=512, subregion=NULL) {
  if (is.null(subregion)) {
    spec(wav, wl=wl, fftw=T, PSD=T, plot=F)
  } else {
    spec(wav, wl=wl, fftw=T, PSD=T, from=subregion[1], to=subregion[2], plot=F)
  }
}


calc_fm = function(wav, wl=512, subregion=NULL) {
  
}

#' stripped down version of seewave::spec
calc_psd_region = function(wav_matrix, f, wl = 512, wn = "hanning", fftw = TRUE, norm = TRUE, 
                           PSD = TRUE, from = NULL, to = NULL, at = NULL, ovlp=0) {
  
  # input <- inputw(wave = wave, f = f)
  #   wave <- input$w
  #   f <- input$f
  #wave = as.matrix(wav@left)
  wave = wav_matrix
  #f = wav@samp.rate
  #rm(input)
  
  if (!is.null(from) && !is.null(to)) {
    
    a <- round(from * f)
    b <- round(to * f)
    #    }
    wl <- (b - a) + 1
    # wave = wave[a:b,]
    wave <- as.matrix(wave[a:b, ])
  }
  if (!is.null(at)) {
    c <- round(at * f)
    wl2 <- wl%/%2
    wave <- as.matrix(wave[(c - wl2):(c + wl2), ])
  }
  n <- nrow(wave)
  W <- ftwindow(n, wn = wn)
  wave <- wave * W
  p <- fftw::planFFT(n)
  y <- Mod(fftw::FFT(wave[, 1], plan = p))
  
  y <- y[1:(n%/%2)]
  if (norm) {
    y <- y/max(y)
  }
  y <- ifelse(y == 0, yes = 1e-06, no = y)
  x <- seq(0, (f/2) - (f/wl), length.out = n%/%2)/1000
  if (PSD) 
    y <- y^2
  spec <- cbind(x, y)
  return(spec)
  
}


######### AMPLITUDE #########
calc_amp = function(wav, subsamp=1, band=NULL) {
  if (is.null(band)) {
    (wav@left[seq(1,length(wav@left), subsamp)])^2
  } else {
    psd = spec(wav, wl=512, fftw=T, PSD=T)
    calc_power(psd, band)
  }
}

calc_amp_subregion = function(wav, subsamp=1, band=NULL, subregion=NULL) {
  if (is.null(band)) {
    if (is.null(subregion)) {
      (wav@left[seq(1,length(wav@left), subsamp)])^2
    } else {
      (wav@left[seq(subregion[1] * wav@samp.rate, subregion[2] * wav@samp.rate, subsamp)])^2
    }
  } else {
    if (is.null(subregion)) {
      psd = spec(wav, wl=512, fftw=T, PSD=T)
    calc_power(psd, band)
    } else {
      psd = spec(wav, wl=512, fftw=T, PSD=T, from=subregion[1], to=subregion[2])
      calc_power(psd, band)
    }
  }
}

calc_mean_amp = function(wav, subsamp=1, band=NULL, subregion=NULL) {
  mean(calc_amp_subregion(wav, subsamp, band, subregion))
}

calc_amp_stats = function(wav, subsamp=1, band=NULL, subregion=NULL) {
  d = calc_amp_subregion(wav, subsamp, band, subregion)
  fs = wav@samp.rate
  
  time_to_max = which.max(d) / fs
  return(c(amp_mean = mean(d), 
           amp_sd = sd(d),
           amp_cv = sd(d)/mean(d),
           time_to_max = time_to_max))
}

calc_max_mean_ratio = function(wav, subsamp=10, band=NULL) {
  amp_env = calc_amp(wav, subsamp=subsamp)
  return(max(amp_env) / mean(amp_env))
}

calc_max_median_ratio = function(wav, subsamp=10, band=NULL) {
  amp_env = calc_amp(wav, subsamp=subsamp)
  return(max(amp_env) / median(amp_env))
}

calc_power = function(psd, band=NULL) {
  if (is.null(band)) {
    sqrt(sum(psd[,2]))
  } else {
    sqrt(sum(psd[find_freqband(psd, band),2]))    
  }
}

calc_amp_ratio = function(amp, samp.rate, peaks, max_gap=50, subsamp=1, min_num_gaps=5) {
  
  if (nrow(peaks) <= 2) 
    return(1)
  max_gap = max_gap / 1000
  
  gaps = matrix(0, nrow=nrow(peaks) - 1, ncol=2)
  gaps[,1] = peaks[1:(nrow(peaks) - 1), 2]
  gaps[,2] = peaks[2:(nrow(peaks)), 1]
  
  ind = (gaps[,2] - gaps[,1])<max_gap
  gaps = matrix(gaps[ind,], ncol=2)
  
  if (nrow(gaps)<min_num_gaps) {
    return(2)
  } else {
    peak_ind = sapply(which(ind), function(x) c(x, x+1))
    peaks = peaks[unique(c(peak_ind)),]
    peak_amp = mean(apply(peaks, 1, function(x) {
      d = amp[round((x[1]*samp.rate)/subsamp):(round(x[2]*samp.rate)/subsamp)]
      mean(d)
    }))
    gap_amp = mean(apply(gaps, 1, function(x) {
      d = amp[(round(x[1]*samp.rate)/subsamp):(round(x[2]*samp.rate)/subsamp)]
      mean(d)
    }))
  }
  amp_ratio = gap_amp / peak_amp
}

######### ENTROPY #########
calc_entropy = function(x) {
  y = x[x>0]
  -1 * sum(y * log2(y))
}

calc_entropy_psd = function(psd, band=NULL) {
  if (is.null(band)) {
    calc_entropy(psd[,2])
  } else {
    calc_entropy(psd[find_freqband(psd, band),2])    
  }
  
}
######### WIENER ENTROPY ########


#' Calculate Wiener entropy of given frequency band in Wave object
#' @param wav, Wave object
#' @param band, selected frequency band
#' @param subregion, vector giving start and end of wav subregion to analyze, in seconds
#' @references https://en.wikipedia.org/wiki/Spectral_flatness
#' @return scalar, calculated Weiner entropy
wiener_entropy = function(wav, band=NULL, subregion=NULL) {
  psd = NULL
  wl = 512
  if (!is.null(subregion)) {
    # psd = seewave::spec(wav, wl=wl, plot=F, PSD=T, from=subregion[1], to=subregion[2], fftw=T)
    psd = calc_psd_region(as.matrix(wav@left), f=wav@samp.rate, wl=wl, from=subregion[1], to=subregion[2])
  } else {
    psd = spec(wav, wl=wl, plot=F, PSD=T)
  }
  
  psd[,1] = 1000 * psd[,1] # convert to Hertz from kHertz
  if (is.null(band)) band = c(1, nrow(psd))
  psd1 = psd[psd[,1]>band[1] & psd[,1]<band[2],]
  res = calc_wiener_entropy(psd1[,2])
  return(res)
}

#' Calculate Weiner entropy of given frequency band in PSD
#' @param psd, PSD as output from spec(PSD=T)
#' @param band, selected frequency band
#' @references https://en.wikipedia.org/wiki/Spectral_flatness
#' @return scalar, calculated Weiner entropy
calc_wiener_entropy_psd = function(psd, band=NULL, freq_ind=NULL) {
  if (is.null(band)) {
    calc_wiener_entropy(psd[,2])
  } else {
    if (!is.null(freq_ind)) {
      calc_wiener_entropy(psd[freq_ind,2])       
    } else {
      calc_wiener_entropy(psd[find_freqband(psd, band),2]) 
    }

  }
}

calc_wiener_entropy = function(psd) {
  exp(mean(log(psd))) / mean(psd)
  #prod(psd)^(1/length(psd)) / mean(psd)
}

#' Calculate Wiener entropy of given frequency band in Wave object
#' @param wav_matrix, wav signal as matrix
#' @param band, selected frequency band
#' @param subregion, vector giving start and end of wav subregion to analyze, in seconds
#' @references https://en.wikipedia.org/wiki/Spectral_flatness
#' @return scalar, calculated Weiner entropy
wiener_entropy_mat = function(wav_matrix, f, wl = 512, band=NULL, subregion=NULL) {
  psd = NULL
  if (!is.null(subregion)) {
    # psd = seewave::spec(wav, wl=wl, plot=F, PSD=T, from=subregion[1], to=subregion[2], fftw=T)
    psd = calc_psd_region(wav_matrix, f=f, wl=wl, from=subregion[1], to=subregion[2])
  } else {
    psd = spec(wav, wl=wl, plot=F, PSD=T)
  }
  
  psd[,1] = 1000 * psd[,1] # convert to Hertz from kHertz
  if (is.null(band)) band = c(1, nrow(psd))
  psd1 = psd[psd[,1]>band[1] & psd[,1]<band[2],]
  res = exp(mean(log(psd1[,2]))) / mean(psd1[,2])
  return(res)
}

#' Calculate spectral features iniven frequency band in Wave object
#' @param wav_matrix, wav signal as matrix
#' @param band, selected frequency band
#' @param subregion, vector giving start and end of wav subregion to analyze, in seconds
#' @references https://en.wikipedia.org/wiki/Spectral_flatness
#' @return scalar, calculated Weiner entropy
spectral_features_mat = function(wav_matrix, f, wl = 512, band=NULL, subregion=NULL) {
  psd = NULL
  if (!is.null(subregion)) {
    # psd = seewave::spec(wav, wl=wl, plot=F, PSD=T, from=subregion[1], to=subregion[2], fftw=T)
    psd = calc_psd_region(wav_matrix, f=f, wl=wl, from=subregion[1], to=subregion[2])
  } else {
    psd = spec(wav, wl=wl, plot=F, PSD=T)
  }
  
  psd[,1] = 1000 * psd[,1] # convert to Hertz from kHertz
  if (is.null(band)) band = c(1, nrow(psd))
  psd1 = psd[psd[,1]>band[1] & psd[,1]<band[2],]
  res = exp(mean(log(psd1[,2]))) / mean(psd1[,2])
  return(res)
}

######### AGGREGATE ########
calc_spectral_features_batch = function(info, wl=512, band=NULL, subregion=NULL, overlap=50, thresh=0) {
  #data = apply(info, 1, function(row) {
 # info = info[9:10,]
  data = foreach(row=isplitRows(info, chunkSize=1)) %dopar% {
    #print(as.data.frame(row))
    mat = load_mat(row$mat)
   # mat = mat[!(mat$labels %in% c("-", "")),]
    mat$wav =row$wav
    if(nrow(mat) < 10) {(NULL)
      return(NULL)
    }
    wav = filtersong(readWave(row$wav))
    d = foreach(syl=isplitRows(mat, chunkSize=1), .combine="rbind") %do% {
      calc_spectral_features(wav, wl=wl, band=band, subregion=c(as.numeric(syl$onsets), as.numeric(syl$offsets)), overlap=overlap, thresh=thresh)
    }
    data.frame(mat, d)
  }
  data = data[unlist(lapply(data, function(x) !is.null(x)))]
  #})
  #data1 = bind_rows(data)
  return(bind_rows(data))
}
calc_spectral_features = function(wav, wl=512, band=NULL, subregion=NULL, overlap=50, thresh=0) {
  window = wl / wav@samp.rate
  fs = wav@samp.rate
  
  steps = NULL
  if ((subregion[2] - subregion[1]) < window) {
    steps = c(0)
  } else {
    steps = seq(0, (subregion[2]-subregion[1] - window), window * (100-overlap) / 100)
  }
  
  wav_duration = length(wav@left) / fs
  subregion2 = subregion
  wav_mat = as.matrix(wav@left)
  freqs = seq(0, (fs-1) / 2, fs / (wl)) / 1000
  res = lapply(steps, function(step) {
    subregion2[1] = subregion[1] + step
    subregion2[2] = subregion2[1] + window
    #print(subregion2)
    if (subregion2[1] <= window |  subregion2[2] >= (wav_duration-window)) return(NA)
    a = calc_psd_region(wav_mat, 
                        f = fs, 
                        wl = wl, 
                        fftw = T, 
                        PSD=T, 
                        at = subregion2[1],
                        norm=F)[,2]
    ##anorm = a / max(a)
    #a[anorm<thresh] = min(a)
    a
  })
  res = do.call(cbind, res)
  res = t(na.omit(t(res)))
  if(nrow(res)< (wl / 2)) {
    freq_mean = NA
    freq_mean_cv = NA
    freq_median = NA
    freq_median_cv = NA
    freq_max = NA
    freq_max_cv = NA
    
    went_mean = NA
    went_var = NA
    good_ff = NA
    goodness = NA
  } else {
    freq_ind = find_closest_inds(band, freqs)
    spec = apply(res, 2, function(x) {
      psd = cbind(freqs, x)
      c(calc_wiener_entropy_psd(psd, freq_ind=freq_ind),
        calc_mean_freq_psd(psd),
        calc_median_freq_psd(psd),
        calc_max_freq_psd(psd))
      }) 
    #went_mean = mean(spec[1,])
    
    
    #res1 = data.frame(ind=1:nrow(res), val=rowSums(res))
    #res1 = res1[order(-res1$val),]
    #ress = cumsum(res1$val)
    #ress = ress / max(ress)
    #inds = res1$ind[ress>.9]
    #res1$val[res1$ind %in% inds] = min(res1$val)
    #res1 = res1$val[order(res1$ind)]
    res1 = rowSums(res)
    went_mean = calc_wiener_entropy_psd(cbind(freqs, res), freq_ind=freq_ind)
    went_var = var(spec[1,])
    
    freq_mean = mean(spec[2,])
    freq_mean_cv = sd(spec[2,]) / freq_mean
    
    freq_median = mean(spec[3,])
    freq_median_cv = sd(spec[3,], freq_median)
    
    freq_max = mean(spec[4,])
    freq_max_cv = sd(spec[4,], freq_median)
    
    good = calc_goodness(wav, wl=wl, wn="hanning", subregion=subregion, freq_range=band)
    good_ff = good$ff
    goodness = good$goodness
  }
  cnames = c("freq_mean", "freq_mean_cv", 
             "freq_median", "freq_median_cv",
             "freq_max", "freq_max_cv",
             "went_mean", "went_var", "good_ff", "goodness")
  data = data.frame(mget(cnames), check.names = F)
  return(data)
}

calc_amplitude_features_batch = function(info, band=NULL,  subsamp=1) {
  data = foreach(row=isplitRows(info, chunkSize=1)) %dopar% {
    mat = load_mat(row$mat)
    mat$wav =row$wav
    if(nrow(mat) < 10) {
      return(NULL)
    }
    
    mat = process_mat(mat)
    wav = filtersong(readWave(row$wav))
    d = foreach(syl=isplitRows(mat, chunkSize=1), .combine="rbind") %do% {
      calc_amp_stats(wav,  band=band, subregion=c(as.numeric(syl$onsets), as.numeric(syl$offsets)), subsamp=subsamp)
    }
    data.frame(mat, d)
  }
  #data1 = bind_rows(data)
  return(bind_rows(data))
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

######### FREQUENCY #########
calc_mean_freq = function(wav, wl=512, subregion=NULL, band=NULL) {
  psd = calc_psd(wav, wl=wl, subregion=subregion)
  calc_mean_freq_psd(psd, band=band)
}

calc_mean_freq_psd = function(psd, band=NULL) {
  psd[,2] = psd[,2] / sum(psd[,2])
  if (is.null(band)) {
    return(t(psd[,1]) %*% psd[,2])  
  } else {
    inds = find_freqband(psd, band)
    return(t(psd[inds,1]) %*% psd[inds,2])
  }
}

calc_sd_freq = function(psd, band=NULL, mean_val=NULL) {
  if(is.null(mean_val)) 
    mean_val = calc_mean_freq_psd(psd, band=band)
  if (is.null(band)) {
    sqrt(sum((t(psd[,1] - mean_val) %*% psd[,2])^2))
  } else {
    inds = find_freqband(psd, band)
    sqrt(sum((t(psd[inds,1] - mean_val) %*% psd[inds,2])^2)) 
  }
}

calc_max_freq_psd = function(psd) {
  psd[which.max(psd[,2]),1]
}

calc_median_freq_psd = function(psd) {
  cumamp = cumsum(psd[,2])
  psd[length(cumamp[cumamp <= 0.5]) + 1]
}

calc_freq_stats = function(wav, q=c(0.5), wl=512, band=NULL, subregion=NULL, overlap=50) {
  ms =  seewave::meanspec(wav, wl = wl, ovlp = overlap, fftw = T, PSD = T, from = subregion[1], to = subregion[2], plot = F)
  c1 = cumsum(ms[,2])
  c1 = c1/max(c1)
  #r = seq(min(ms[,1]), max(ms[,1]), .1)
  #c2 = interp1(ms[,1], c1, r)
  res = sapply(q, function(x) {
    ms[which((c1>(x-.05)) & (c1 < (x+.05)))[1],1]
  })
  names(res) = q
  return(res)
}

######## CEPSTRUM ########
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
calc_cepstrum = function(wav, wl=512, overlap=0, subregion=NULL) {
  cp = NULL
  if (is.null(subregion)) {
    cp = ceps(wav, wl=wl, plot = F)
  } else {
    cp = ceps(wav, wl=wl, at=mean(subregion), plot=F)
  }
  #return(format_cepstrum(wav_cep))
}

calc_goodness = function(wav, wl=512, wn="hanning", subregion=NULL, freq_range=c(1, 4)) {
  freq_range = freq_range * 1000
  cp = NULL
  if (!is.null(subregion)) {
    cp = ceps(wav, wl=wl, at=mean(subregion), plot=F)
  } else { 
    cp = ceps(wav, wl=wl, plot=F)
  }
  #cp_ind1 = find_closest_inds(1/freq_range[2], cp[,1])
  #cp_ind2 = find_closest(1/freq_range[1], cp[,1])
  #cp_ind1 = which.min(dist(cp[,1], 1/freq_range[2]))
  #cp_ind2 = which.min(dist(cp[,1], 1/freq_range[1]))
  inds = find_closest_inds(c(1/freq_range[2], 1/freq_range[1]), cp[,1])
  #cp = cp[cp_ind1:cp_ind2,]
  cp = cp[inds,]
  cp_max_ind = which.max(cp[,2])
  cp_max_power = cp[cp_max_ind,2]
  cp_max_freq = 1/cp[cp_max_ind,1]
  return(list(ff=cp_max_freq / 1000, goodness=cp_max_power))
}

calc_cepstrum = function(wav, wl=512, wn="hanning", subregion=NULL) {
  spec_floor = 0.05
  input <- inputw(wave = wav, f = f)
  wave <- input$w
  f <- input$f
  if (!is.null(subregion)) {
    a <- round(subregion[1] * f)
    b <- round(subregion[2] * f)
    wave = as.matrix(wave[a:b,])
    #ind = round(mean(subregion) * f)
    #wl2 = wl%/%2
    #wave = as.matrix(wave[(ind - wl2):(ind + wl2), ])
  }
  n <- nrow(wave)
  N = round(n/2)
  #W <- ftwindow(n, wn = wn)
  #wave <- wave * W
  p = fftw::planFFT(nrow(wave))
  #y <- abs(fftw::FFT(wave[,1], plan = p)^2)
  y = abs(fft(wave[,1]))
  y = pmax(y, spec_floor)
  yl = log(y)
  #yl = c(0, diff(yl))
  #p2 = fftw::planFFT()
  #cp = Re(fftw::IFFT(yl, plan=p))
  cp = Re(fft(yl, inverse=T))
  cp = cp[1:N]
  x = seq(0, N/f, length.out=N)
  return(cbind(x, cp))
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

######## TEMPO ########
calc_tempo = function(peaks, gap_thresh=.5) {
  return(nrow(peaks) / (peaks[nrow(peaks),2]-peaks[1,1]))
}

calc_rolling_tempo = function(peaks, chunk_size=5) {
  npeaks = nrow(peaks)
  if (npeaks<5) 
    return(list(tempo_mean=calc_tempo(peaks), tempo_sd = NA))
  
  start = 1
  stop = chunk_size
  samples = npeaks - chunk_size + 1
  tempos = vector(mode="numeric", length = samples )
  for (i in 1:samples) {
    tempos[i] = calc_tempo(peaks[start:stop,])
    start = start + 1
    stop = stop + 1
  }
  
  return(list(tempo_mean = mean(tempos), tempo_sd = sd(tempos)))
  
}


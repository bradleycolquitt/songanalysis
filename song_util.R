library(seewave)
library(tuneR)
library(signal)
library(manipulate)

plot_2dspec = function(wav) {
  theme_set(theme_classic())
  manipulate({
    gg = ggspectro(wav, wl=256) + geom_tile(aes(fill=amplitude)) 
    gg = gg + scale_fill_gradientn(colours=spectro.colors(30), limits=c(-75, 0), na.value="transparent")
    gg = gg + xlim(xmin, xmax)
    gg
  }, 
  xmax=slider(0, length(wav) / wav@samp.rate, initial=length(wav)/wav@samp.rate),
  xmin=slider(0, length(wav) / wav@samp.rate, initial=0))
}

#' High pass filter Wave object at 1800 Hz
#' @param wav, the input Wave object
filtersong = function(wav) {
  return(ffilter(wav, from=1800, output="Wave"))
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

#' Find amplitude peaks in wav file
#' @param wav, the input Wave object
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

plot_peaks = function(peaks) {
  apply(peaks, 1, function(peak) {
    segments(peak[1], -1E7, peak[2], -1E7, col=2)
  })
}
#' Calculate Weiner entropy of given frequency band in Wave object
#' @param wav, Wave object
#' @param band, selected frequency band
#' @references https://en.wikipedia.org/wiki/Spectral_flatness
#' @return scalar, calculated Weiner entropy
weiner_entropy = function(wav, band=NULL, subregion=NULL) {
  psd = NULL
  if (!is.null(subregion)) {
    psd = spec(wav, wl=256, plot=F, PSD=T, from=subregion[1], to=subregion[2])
  } else {
    psd = spec(wav, wl=256, plot=F, PSD=T)
  }
  
  if (is.null(band)) band = c(1, nrow(psd))
  psd1 = psd[psd[,1]>band[1] & psd[,1]<band[2],]
  res = exp(mean(log(psd1[,2]))) / mean(psd1[,2])
  return(res)
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

#' Classify Wave object as song or not song
#'     First filters for small Weiner entropy values within given band
#'     Next filters for small interpeak distances
#' @param wav, Wave object
#' @param band, selected frequency band
#' @param wein_thresh, maximum allowed Weiner entropy in frequency band
#' @param dist_thresh, maximum interpeak distance allowed
#' @return logical, song or not song
songfinder = function(wav, band=c(5,7), 
                      wein_thresh=0.4, 
                      dist_thresh = 0.06,
                      min_duration = 1000) {
  wav = filtersong(wav)
  wein = weiner_entropy(wav, band=band)
  if (wein > wein_thresh) {
    return(FALSE)
  } 
  peaks = findpeaks(wav, min_duration=min_duration, max_gap=200, thresh=0.05) 
  #median_peaks_dist = median(peak_dist(peaks))
  #return(median_peaks_dist < dist_thresh)
  return(nrow(peaks) > 0)
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

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

filtersong = function(wav) {
  return(ffilter(wav, from=1800, output="Wave"))
}

#' Moving average smooth Wave object
#' @param wav the input Wave object
#' @param window the smoothing window in milliseconds
#' @return Wave object containing smoothed data in left slot 
# FIX
smoothsong = function(wav, window=2) {
  # window is in units of ms
  window_length = round(wav@samp.rate * window / 1000)
  data = wav@left
  filt = rep(1, times=window_length) / window_length
  data1 = conv(data, filt)
  offset = round((length(data1) - length(data))/2)
  return(Wave(data1[(1+offset):(length(data)+offset)], samp.rate=wav@samp.rate, bit=wav@bit))
}

findpeaks = function(wav, max_gap=75, thresh=0.5) {
  samp_rate = wav@samp.rate
  #min_width = min_width * samp_rate / 1000 
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
  return (df / samp_rate)
}

peak_dist = function(peaks) {
  dists = vector("numeric", nrow(peaks)-1)
  for (i in 1:(nrow(peaks)-1)) {
    dists[i] = peaks[i+1,1] - peaks[i,2]
  }
  return(dists)
}


weiner_entropy = function(wav, band=c(6,8)) {
  psd = spec(wav, wl=256, plot=F, PSD=T)
  psd1 = psd[psd[,1]>band[1] & psd[,1]<band[2],]
  res = exp(mean(log(psd1[,2]))) / mean(psd1[,2])
  return(res)
}


weiner_entropy_psd = function(psd, region=c(6,8)) {
  psd1 = psd[psd[,1]>region[1] & psd[,1]<region[2],]
  res = exp(mean(log(psd1[,2]))) / mean(psd1[,2])
  return(res)
}


songfinder = function(wav, band=c(5,7), 
                      wein_thresh=0.4, 
                      dist_thresh = 0.06) {
  wav = filtersong(wav)
  wein = weiner_entropy(wav, band=band)
  if (wein > wein_thresh) {
    return(FALSE)
  } 
  peaks = findpeaks(wav, max_gap=25, thresh=0.1) 
  median_peaks_dist = median(peak_dist(peaks))
  return(median_peaks_dist < dist_thresh)
}


songfinder_psd = function(psd, band=c(5,8), 
                      sum_thresh=0.4, var_thresh=.002,
                      wein_thresh=0.45) {
  #psd_sum = sum(psd[psd[,1]>band[1] & psd[,1]<band[2],2]) / sum(psd[,2])
  #psd_var = var(psd[psd[,1]>band[1] & psd[,1]<band[2],2]) 
  wein = weiner_entropy_psd(psd, region=band)
  return(wein<wein_thresh)
  #return(psd_sum>sum_thresh & psd_var>var_thresh & wein<wein_thresh)
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

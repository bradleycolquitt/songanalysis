suppressWarnings(library(tuneR))
suppressWarnings(library(graphics))
suppressWarnings(library(seewave))
INT_MIN = 1

threshold_auto = function(wav, method, log=F, absolute=T, floor=F, ...) {
  amp = NULL
  if (absolute) {
    amp = seewave::env(wav, envt="abs", plot=F)
  } else {
    amp = wav@left^2
  }
  
  median_amp = median(amp)
  if (floor) 
    amp =  amp[amp>median_amp]
  
  if (log) {
    amp = log(amp+1)
  }
  #amp = amp[amp[,1]<=5E6,]
  val = method(amp, ...)
  if (log) {
    val = exp(val)
  }
  #if(floor)
  #  val = val + median_amp
  val
}

# currently takes too long, will need dynamic programming approach
renyi = function(data) {
  ih = 0 
  it = 0
  threshold = 0
  opt_threshold = 0
  first_bin = 0
  last_bin = 0
  tmp_var = 0
  t_star1 = 0
  t_star2 = 0
  t_star3 = 0
  beta1 = 0 
  beta2 = 0 
  beta3 = 0
  alpha = 0			# alpha parameter of the method 
  term= 0
  tot_ent = 0		# total entropy 
  max_ent = 0		# max entropy 
  ent_back = 0		# entropy of the background pixels at a given threshold 
  ent_obj = 0		# entropy of the object pixels at a given threshold 
  omega = 0
 
  #data			 normalized histogram data 
  #*P1			 cumulative normalized histogram 
  #*P2			 see below 
  #Histo *norm_histo		 normalized histogram 
  
  #  Calculate the normalized histogram 
  breaks =seq(min(data), max(data)+100, 1000)
  data = data[data<max(breaks)]
  histogram = hist(data, breaks=breaks, plot=F)
  #return(histogram)
  data = data.frame(breaks=histogram$breaks[1:(length(histogram$breaks)-1)], counts=histogram$counts / sum(histogram$counts))
 
  # Calculate the cumulative normalized histogram 
  data$counts.cumsum = cumsum(data$counts)
  data$counts.cumsumRev = 1 - data$counts.cumsum
  #return(data)
  # Determine the first non-zero bin starting from the first bin  
  first_bin = min(which(data$counts > 0))
  # Determine the first non-one bin starting from the last bin
  last_bin = max(which(data$counts > 0))
  print(last_bin)
  ##  Maximum Entropy Thresholding - BEGIN 
  #   ALPHA = 1.0 
  #   
  #  Calculate the total entropy at each gray-level and find the threshold that maximizes it 
  print("alpha 1")
  for ( it in first_bin:last_bin) {
    print(it)
    #if (it %% 100 == 0) print(it)
    #   Entropy of the background pixels 
    ent_back = 0
    for (ih in 1:it) {
      if ( data$counts[ih] > 0 ) {
        ent_back = ent_back - ( data$counts[ih] / data$counts.cumsum[it] ) * log ( data$counts[ih] / data$counts.cumsum[it] )
      }
    }

    # Entropy of the object pixels 
    ent_obj = 0
    for (ih in c((it):nrow(data))) {
      if ( data$counts[ih] > 0 ) {
        ent_obj = ent_obj - ( data$counts[ih] / data$counts.cumsumRev[it] ) * log ( data$counts[ih] / data$counts.cumsumRev[it] )
      }
    }
    #  Total entropy 
    tot_ent = ent_back + ent_obj
  
  if ( tot_ent > max_ent )
  {
    max_ent = tot_ent
    threshold = data$breaks[it]
  }
  }
  t_star2 = threshold
  
  ### alpha 0.5
  threshold = INT_MIN
  max_ent = 0
  alpha = 0.5
  term = 1.0 / ( 1.0 - alpha )
  print("alpha 0.5")
  for ( it in first_bin:last_bin) {
    # Entropy of the background pixels 
    ent_back = 0
    for ( ih in 1:it) {
      ent_back = ent_back +  sqrt ( data$counts[ih] / data$counts.cumsum[it] )
    }
    
    # Entropy of the object pixels 
    ent_obj = 0
    for ( ih  in (it+1):nrow(data)) {
      ent_obj = ent_obj + sqrt( data$counts[ih] / data$counts.cumsumRev[it] )
    }
    
    # Total entropy 
    if (ent_back * ent_obj > 0) {
      tot_ent = term * log( ent_back * ent_obj )
    } else {
      tot_ent = 0
    }
    
    if ( tot_ent > max_ent ) {
      max_ent = tot_ent
      threshold = data$breaks[it]
    }
  }
  t_star1 = threshold
  
  print("alpha 2")
  threshold = INT_MIN
  max_ent = 0
  alpha = 2
  term = 1.0 / ( 1.0 - alpha )
  for ( it in  first_bin:last_bin) {
    # Entropy of the background pixels 
    ent_back = 0
    for (ih in 1:it) {
      ent_back = ent_back + ( data$counts[ih]^2 ) / ( data$counts.cumsum[it]^2)
    }
    
    # Entropy of the object pixels 
    ent_obj = 0
    for ( ih in (it+1):nrow(data)) {
      ent_obj = ent_obj + ( data$counts[ih]^2 ) / ( data$counts.cumsumRev[it]^2)
    }
    
    # Total entropy 
    if (ent_back * ent_obj > 0) {
      tot_ent = term * log( ent_back * ent_obj )
    } else {
      tot_ent = 0
    }
    
    if ( tot_ent > max_ent )
    {
      max_ent = tot_ent
      threshold = data$breaks[it]
    }
  }
  t_star3 = threshold
  
  #* Sort t_star values 
  if ( t_star2 < t_star1 ) {
    tmp_var = t_star2
    t_star2 = t_star1
    t_star1 = tmp_var
    }
  if ( t_star3 < t_star2 )
  {
    tmp_var = t_star3
    t_star3 = t_star2
    t_star2 = tmp_var
  }
  if ( t_star2 < t_star1 )
  {
    tmp_var = t_star2
    t_star2 = t_star1
    t_star1 = tmp_var
  }
  
  # Adjust beta values 
  if ( abs ( t_star1 - t_star2 ) <= 5 )
    {
      if ( abs ( t_star2 - t_star3 ) <= 5 )
      {
        beta1 = 1
        beta2 = 2
        beta3 = 1
      }
      else
      {
        beta1 = 0
        beta2 = 1
        beta3 = 3
      }
    }
  else
  {
    if ( abs ( t_star2 - t_star3 ) <= 5 )
    {
      beta1 = 3
      beta2 = 1
      beta3 = 0
    }
    else
    {
      beta1 = 1
      beta2 = 2
      beta3 = 1
    }
  }
  
  # Determine the optimal threshold value 
  omega = data$counts.cumsum[t_star3] - data$counts.cumsum[t_star1]
  opt_threshold = t_star1 * ( data$counts.cumsum[t_star1] + 0.25 * omega * beta1 )
  + 0.25 * t_star2 * omega * beta2
  + t_star3 * ( data$counts.cumsumRev[t_star3] + 0.25 * omega * beta3 )
  
  return (opt_threshold)
}

huang = function(data, ...) {
  threshold = 0
  
  breaks =seq(min(data), max(data), length.out = 1000)
  data = data[data<max(breaks)]
  histogram = hist(data, breaks=breaks, plot=F)
  data = data.frame(breaks=histogram$breaks[1:(length(histogram$breaks)-1)], 
                    counts=histogram$counts / sum(histogram$counts))
  
  # Determine the first non-zero bin starting from the first bin  
  first_bin = min(which(data$counts > 0))
  # Determine the first non-one bin starting from the last bin
  last_bin = max(which(data$counts > 0))
  
  # 
  #  Equation (4) in Ref. 1
  #  C = g_max - g_min
  #  This is the term ( 1 / C ) 
  #
  term = 1 / ( last_bin - first_bin )
  
  # Equation (2) in Ref. 1 
  mu_0 = vector("numeric", nrow(data))
  sum_pix = num_pix = 0
  for (ih in first_bin:nrow(data)) {
    sum_pix = sum_pix + data$breaks[ih] * data$counts[ih]
    num_pix = num_pix + data$counts[ih]
    
    mu_0[ih] = sum_pix / num_pix
  }
  
  # Equation (3) in Ref. 1 
  mu_1 = vector("numeric", nrow(data))
  sum_pix = num_pix = 0
  for (ih in last_bin:1) {
    sum_pix = sum_pix + data$breaks[ih] * data$counts[ih]
    num_pix = num_pix + data$counts[ih]
    
    # NUM_PIX cannot be zero ! 
    mu_1[ih] = sum_pix / num_pix
  }
  
  # Determine the threshold that minimizes the fuzzy entropy 
  threshold = INT_MIN
  min_ent = .Machine$double.xmax

  ent = 0
  
  #mu_x = matrix(nrow=nrow(data), ncol=nrow(data))
  print("Prefilling entropy matrices")
  mu_x = 0
  ent0 = matrix(0, nrow=nrow(data), ncol=nrow(data))
  ent1 = matrix(0, nrow=nrow(data), ncol=nrow(data))
  
  for (it in 1:nrow(ent0)) {
      if (it == round(nrow(ent0)/2)) {
        print("here")
      }
        
      mu_x0 = 1.0 / ( 1.0 + term * abs ( data$breaks - mu_0[it] ) )
      

      ind = mu_x0 > 0 & mu_x0 < 1
      mu_x0_ind = mu_x0[ind]
    
      ent0[it,ind] = data$counts[ind] * ( -mu_x0_ind * log ( mu_x0_ind ) -
                                    ( 1.0 - mu_x0_ind ) * log ( 1.0 - mu_x0_ind ) )
     
      mu_x1 = 1.0 / ( 1.0 + term * abs ( data$breaks - mu_1[it] ) )
      ind = mu_x1 > 0 & mu_x1 < 1
      mu_x1_ind = mu_x1[ind]
 
      ent1[it,ind] = data$counts[ind] * ( -mu_x1_ind * log ( mu_x1_ind ) -
                                            ( 1.0 - mu_x1_ind ) * log ( 1.0 - mu_x1_ind ) )
  
  }
  

  print("Calculating cumulative entropies")
  for (it in 1:nrow(data)) {
    #ent = 0
    ent = sum(ent0[it,1:(it-1)], ent1[it,(it):ncol(ent1)])
    #for (ih in 1:it) {
      
      # Equation (4) in Ref. 1 
      #mu_x = 1.0 / ( 1.0 + term * abs ( data$breaks[ih] - mu_0[it] ) )
      #if ( mu_x > 0 & mu_x < 1 ) {
        # Equation (6) & (8) in Ref. 1 
        #ent = ent + data$counts[ih] * ( -mu_x * log ( mu_x ) -
        #                 ( 1.0 - mu_x ) * log ( 1.0 - mu_x ) )
      #  ent = ent + ent0[it, ih]
      #}

    
    #for ( ih in it:nrow(data) ) {
      # Equation (4) in Ref. 1 
#       mu_x = 1.0 / ( 1.0 + term * abs ( data$breaks[ih] - mu_1[it] ) )
#       if ( mu_x > 0 & mu_x < 1 ) {
#         # Equation (6) & (8) in Ref. 1 
#         ent = ent + 
#         data$counts[ih] * ( -mu_x * log ( mu_x ) -
#                          ( 1.0 - mu_x ) * log ( 1.0 - mu_x ) )
     #   ent = ent + ent1[it, ih]
    #}

    if ( ent < min_ent ) {
      
      min_ent = ent
      threshold = data$breaks[it]
    }
  }
  
  return (threshold)
}

#' Implementation of Otsu automatic thresholding algo
otsu = function(data, factor=1, ...) {
  breaks = seq(min(data), max(data), length.out = 1000)
  data = data[data<max(breaks)]
  total = length(data)
  
  histogram = hist(data, breaks=breaks, plot=F)
  histogram = data.frame(breaks=histogram$breaks[1:(length(histogram$breaks)-1)], counts=histogram$counts)
  sumB = 0
  wB = 0
  maximum = 0
  threshold1 = 0
  threshold2 = 0
  sum1 = sum(as.numeric(apply(histogram, 1, prod)))
  for (i in 1:nrow(histogram)) {
    wB = wB + histogram[i,2]
    if (wB == 0) next
    
    wF = total - wB
    if (wF == 0) break
    
    sumB = sumB +  histogram[i,1] * histogram[i,2]
    mB = sumB / wB
    mF = (sum1 - sumB) / wF
    
    between = wB * wF * (mB - mF) * (mB - mF)
    if ( between >= maximum ) {
      threshold1 = histogram[i,1]
      if ( between > maximum ) {
        threshold2 = histogram[i,1]
      }
      maximum = between
    }
  }
  (threshold1 + threshold2 )/2 * factor
}

mean_sd = function(data, factor = 0) {
  return(mean(data) + sd(data) * factor)
}

median_sd = function(data, factor = 0) {
  return(median(data) + sd(data) * factor)
}

max_amp = function(data, factor) {
  return(max(data) * factor)
}
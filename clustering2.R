library(mclust)
library(foreach)
library(doMC)
library(proxy)
library(cluster)
library(caret)
library(dplyr)
library(gridExtra)
library(mvtnorm)
library(stringr)
#suppressMessages(library(WGCNA))
#registerDoMC(cores=10) 

pam_syllable = function(syls, data, range_to_test=c(4:15)) {
  print("Computing distance...")
  data_d = as.matrix(dist(data, method="kullback")) 
  print("Running PAM...")
  d = lapply(range_to_test, function(x) {
    ppam(data_d, k=x, is_dist=T, cluster.only = )
  })
  d_class = lapply(d, function(x) x$cluster)
  d_class = t(do.call(rbind, d_class))
  
  colnames(d_class) = paste("pam", range_to_test, sep="-")
  syls = cbind(syls, d_class)
  return(list(syls=syls, mc=d))
}

mclust_par = function(data, G, modelNames="VVI", plot=T, parallel=F) {
  gc() 
  mcs = NULL
  if (parallel) {
    mcs = foreach(g=G) %dopar% {
      Mclust(data, G=g, modelNames=modelNames)
    }
  } else {
    mcs = foreach(g=G) %do% {
      Mclust(data, G=g, modelNames=modelNames)
    }
  }

  mcs = mcs[!unlist(lapply(mcs, is.null))]
  
  recovered_G = unlist(lapply(mcs, function(x) x$G))
  bics = unlist(lapply(mcs, function(x) x$BIC))
  if (length(bics) <= 1) 
    return(list(mcs=NULL, stats=NULL))
  ll = unlist(lapply(mcs, function(x) x$loglik))
  lr = vector("numeric", length(ll))
  lr[1] = 1
  for (i in 2:length(lr)) {
    lr[i] = bics[i] / bics[i-1]
  }
  ps = pchisq(lr, df = 1)
  
  sils = rep(0, length(ll))
  d = data.frame(G=recovered_G, stat=rep(c("BIC", "logLik", "Sil"), each=length(recovered_G)), 
                 value=c(bics, ll, sils))
  
  if (plot) {
    g = ggplot(d, aes(G, value, color=stat, group=stat)) + geom_point() + geom_line() + facet_grid(.~stat, scales="free_y")
    print(g)
  }
  bic_max_ind = which.max(bics)
  return(list(mcs=mcs, stats=d))
}

dist_syllable_data = function(data, compare_data=NULL, distmethod="cor", power=2, selectCols=NULL) {
  #data = sweep(data, 1, apply(data, 1, max), FUN = "/")
  
  data_s = NULL
  if (is.null(selectCols)) 
    selectCols = 1:nrow(data)
  if (is.null(compare_data))
    compare_data = data
  
  if (distmethod == "cor") {
    data_s = as.matrix(proxy::dist(data, y=compare_data[selectCols,], method="cor"))
    #data_s = adjacency(t(data), power=1, type="signed", corFnc="cor", selectCols=selectCols)
  } else if (distmethod == "bicor") {
    data_s = adjacency(t(data), power=1, type="signed", corFnc="bicor", selectCols=selectCols)
  }  else if (distmethod == "spearman") {
    #data_s = as.matrix(dist(x=))
    data_s = adjacency(t(data), power=1, type="signed", corOptions="use='p', method='spearman'", selectCols=selectCols)
  } else if (distmethod == "euclidean") {
    data_s = proxy::dist(x=data, y=compare_data[selectCols,], method="euclidean", convert_similarities = F )
  } else if (distmethod == "KL") {
    data_s = as.matrix(proxy::dist(x=data, y=compare_data[selectCols,], method="kullback", convert_similarities = F))
  }
  data_s = data_s ^ power
  return (data_s)
}

plot_cluster_stats = function(data_stats) {
  data_stats_mean = data_stats %>% group_by(G, stat) %>% dplyr::summarize(value=mean(value), rep=1) 
  theme_set(theme_bw())
  g = ggplot(data_stats, aes(G, value, group=factor(rep))) + geom_line(alpha=I(1/2))
  g = g + geom_line(data=data_stats_mean, color="red") + geom_point(data=data_stats_mean, color="red")  + facet_wrap(~stat, scales="free_y")
  g = g + scale_x_discrete(breaks=c(unique(data_stats$G)))
  
  data_stats_mean = ungroup(data_stats_mean) %>% group_by(stat) %>% dplyr::mutate(rel_diff_max = abs(value - max(value)) / max(value),
                                                                           diffs = c(NA, diff(value)) / max(value),
                                                                           diffs2 = c(NA, diff(diffs)))
  
  g1 = ggplot(data_stats_mean, aes(G, diffs)) + geom_point() + facet_wrap(~stat)
 # g1 = g1 + scale_x_discrete(breaks=c(unique(data_stats$G)))
  
  g2 = ggplot(data_stats_mean, aes(G, diffs2)) + geom_point() + facet_wrap(~stat)
  #g2 = g2 + scale_x_discrete(breaks=c(unique(data_stats$G)))
  
  grid.arrange(g, g1, g2, nrow=3)
  return(data_stats_mean)
}

determine_num_clusters = function(data, sample_size, feature_size, reps=20, range_to_test=4:12, distmethod, power) {
  #reps = 20
  inds = lapply(1:reps, function(x) sample(1:nrow(data), sample_size))
  feature_inds = lapply(inds, function(ind) sample(1:length(ind), feature_size))
  #feature_inds = lapply(1:reps, function(x) sample(1:nrow(data), feature_size))
  ## Train several models on subsets of data
  mcs_stats = lapply(1:reps, function(x) {
    # Calculate inter-syllable distances
    print(paste("rep", x, sep=""))
    #data1_s = dist_syllable_data(data[inds[[x]],], distmethod=distmethod, power=power)
    data1_s = dist_syllable_data(data[inds[[x]],], selectCols=feature_inds[[x]], distmethod=distmethod, power=power)
    data1_s = matrix(data1_s, nrow=nrow(data1_s), ncol=ncol(data1_s))
    colnames(data1_s) = rownames(data)[feature_inds[[x]]]
    #data2_s = as.matrix(dist(data1_s))
    mc_stats = mclust_par(data1_s, G=range_to_test, modelNames="VVI", plot=F)
    if (is.null(mc_stats[[1]])) 
      return(list(syls = NULL, mc=NULL))
    mc_stats$train_syls = rownames(data)[inds[[x]]][feature_inds[[x]]]
    mc_stats$stats$rep = x
    
    mc_stats
  })
  
  data_stats = lapply(mcs_stats, function(x) x$stats)
  data_stats = do.call("rbind", data_stats)
  data_stats_mean = plot_cluster_stats(data_stats)
  
  max_bic_ind = data_stats_mean%>% filter(diffs2==max(diffs2)) 
  max_bic_ind = round(mean(max_bic_ind$G))
  models = lapply(mcs_stats, function(x) x$mcs[[which(range_to_test == max_bic_ind)]])
  train_syls = lapply(mcs_stats, function(x) x$train_syls)
  return(list(G=max_bic_ind, mc = models, train_syls = train_syls))
}

predict_syllables = function(data, models, train_syls, G, distmethod, power) {
  preds = lapply(1:length(models), function(i) {
    print(paste("model", i, sep=""))
    # inds = which(rownames(data) %in% train_syls[[i]])
    data1_s = dist_syllable_data(data, selectCols=train_syls[[i]], distmethod=distmethod, power=power)
    data1_s = matrix(data1_s, nrow=nrow(data1_s), ncol=ncol(data1_s))
    rownames(data1_s) = rownames(data)
    #data2_s = as.matrix(dist(data1_s))
    mc.pred = predict(models[[i]], data1_s)
    classification = unlist(mc.pred$classification)
    #names(classification) = rownames(data)
    ## Compute average/median spectral profiles for classes
    data1 = data[match(names(classification), rownames(data)),]
    data_split = split.data.frame(data1, classification)
    data_split1 = do.call("rbind", lapply(data_split, function(d) apply(d, 2, median)))
    data_split1 = sweep(data_split1, 1, apply(data_split1, 1, max), "/" )
    list(classification=classification, avg_syls=data_split1)
    #unlist(mc.pred$classification)
  })
  
  preds1 = lapply(preds, function(x) x$classification)
  preds1 = preds1[unlist(lapply(preds1, function(x) !is.null(x)))]
  avg_syls = lapply(preds, function(x) x$avg_syls)
  avg_syls = avg_syls[unlist(lapply(avg_syls, function(x) !is.null(x)))]
  
  avg_syls_m = melt(t(avg_syls[[1]]))
  print(ggplot(avg_syls_m, aes(Var1, value, group=Var2)) + geom_line() + facet_wrap(~Var2))
  ## Populate matrix with model predictions
  num_models = length(preds1)
  num_data = nrow(data)
  num_syls = G
  pred_mat = matrix(nrow=num_data, ncol=num_models, dimnames=list(rownames(data)))
  for (i in 1:ncol(pred_mat)) {
    pred_mat[,i] = preds1[[i]][match(rownames(pred_mat), names(preds1[[i]]))]
  }
  
  predictions = letters[pred_mat[,1]]
  names(predictions) = rownames(pred_mat)
  return(predictions)
}


update_syllables = function(syls, predicted) {
  syls$called = predicted[match(syls$id, names(predicted))]
  syls$called[is.na(syls$called)] = "-"
  return(syls)
}

write_song_batch = function(songs, subdir="autolabel") {
  songs$syls$wav = unlist(lapply(str_split(songs$syls$id, "-"), function(x) x[1]))
  syls_s = split(songs$syls, songs$syls$wav)
  #  syls_s = syls_s[match(names(songs[[1]]), names(syls_s))]
  
  for (i in 1:length(songs[[1]])) {
    songs$songs[[i]]$syllable = syls_s[[i]]
    song_to_mat(songs$songs[[i]], subdir=subdir)
  }
}

calc_lik_from_ref = function(bird, 
                             reference_date, 
                             peak_source="mat",
                             
                             ## Filtering
                             deselect = c(" ", "-", "i"),
                             
                             ##PSD params
                             wl=512,
                             freq_limits=c(1,10),
                             feature_set = "mid_mean",
                             
                             # Distance calculation
                             distmethod = "euclidean",
                             power = 2,
                             reference_size = 20,
                             
                             # Modeling params
                             test_fraction = .1,
                             nreps = 10,
                             
                             # Plotting params
                             plot_dir = NULL
) {
  require(mvtnorm)
  
  dir = paste("/mnt/bengal_home/song", bird, "songs/select", sep="/")
  info = load_mat_info(dir, file_ex="wav.not.mat")
  songs = parse_song_batch2(info$wav, peak_source=peak_source)
  
  ## Get song features
  song_data = process_song_batch(songs, 
                                 feature_set=feature_set, 
                                 wl=wl, 
                                 smoothed=F, 
                                 cluster=F, 
                                 return_data=T, 
                                 freq_limits=freq_limits)
  
  song_syls = song_data$syls
  song_syls = song_syls %>% 
    mutate(wav_strip =str_replace(id, "-[0-9]+", "") ) %>% 
    filter(!(labels %in% deselect))
  
  song_psds = song_data$syl_data
  song_psds = song_psds[rownames(song_psds) %in% song_syls$id,]
  song_psds_wav = str_replace(rownames(song_psds), "-[0-9]+", "")
  
  info = info %>% mutate(wav_base = basename(wav), 
                         wav_strip =str_replace(wav_base, ".wav", ""))
  base_wavs  = info %>% filter(date<=reference_date)
  
  ## Run replicates
  print("Fitting models, calculating LL...")
  #res2 = foreach(i=1:nreps) %do% {
  res2 = mclapply(1:nreps, function(i) { 
    #res2 = lapply(1:nreps, function(i) {
    print(paste("Rep ", i, sep=""))
    set.seed(i)
    ## Define data sets
    train_wavs = base_wavs %>% sample_frac(size = (1-test_fraction))
    suppressMessages(test_base_wavs = anti_join(base_wavs, train_wavs))
    suppressMessages(test_rest_wavs = info %>% anti_join(base_wavs))
    
    train_syls = song_syls[song_syls$wav_strip %in% train_wavs$wav_strip,]
    train_psds = song_psds[song_psds_wav %in% train_wavs$wav_strip,]
    
    test_base_syls = song_syls[song_syls$wav_strip %in% test_base_wavs$wav_strip,]
    test_base_psds = song_psds[song_psds_wav %in% test_base_wavs$wav_strip,]
    
    test_rest_syls = song_syls[song_syls$wav_strip %in% test_rest_wavs$wav_strip,]
    test_rest_psds = song_psds[song_psds_wav %in% test_rest_wavs$wav_strip,]
    
    psds_list = list(train=train_psds, test_base=test_base_psds, test_rest=test_rest_psds)
    
    ## Define syllable number
    syls_tab = table(train_syls$labels) 
    syls_tab = as.matrix(syls_tab[!(names(syls_tab) %in% deselect)])
    syls_tab = syls_tab[,1] / sum(syls_tab[,1])
    syls_tab = syls_tab[syls_tab>.01]
    num_syls = length(syls_tab)
    
    ## Calculate distances
    print("..Calculating distances..")
    ref_inds = sample(1:nrow(train_psds), reference_size)
    ref_data = train_psds[ref_inds,]
    train_dist = dist_syllable_data(train_psds, compare_data = train_psds, 
                                    distmethod=distmethod, power=power, selectCols=ref_inds)
    
    test_base_dist = dist_syllable_data(test_base_psds, compare_data = train_psds, 
                                        distmethod=distmethod, power=power, selectCols=ref_inds)
    
    test_rest_dist = dist_syllable_data(test_rest_psds, compare_data = train_psds,
                                        distmethod=distmethod, power=power, selectCols=ref_inds)
    
    dist_list = list(train=train_dist, test_base=test_base_dist, test_rest=test_rest_dist)
    
    if (!is.null(plot_dir) & i==1) {
      train_dist1 = adjacency(t(train_dist), type="signed")  
      test_base_dist1 = adjacency(t(test_base_dist), type = "signed")
      test_rest_dist1 = adjacency(t(test_rest_dist), type="signed")
      
      suppressMessages(require(NMF))
      
      ind = 50
      cairo_pdf(paste(plot_dir, "heatmaps.pdf", sep="/"), width=12, height=5)
      par(mfrow=c(1,3), mar=c(0, 1, 1, 1))
      aheatmap(train_dist1[1:ind,1:ind], hclustfun = "ward.D2", color="topo", main="Training")
      aheatmap(test_base_dist1[1:ind,1:ind], hclustfun = "ward.D2", color="topo", main="Baseline - held out")
      aheatmap(test_rest_dist1[1:ind,1:ind], hclustfun = "ward.D2", color="topo", main="Post baseline")
      dev.off()
    }
    
    ## Fit GMM
    print("..Fitting GMM..")
    mod = Mclust(train_dist, G=num_syls, modelNames=c("VII"))
    
    ## Calculate loglikelihoods
    print("..Calculating likelihoods..")
    params = mod$parameters
    res1 = lapply(dist_list[2:3], function(p) {
      calc_likelihood(p, params)
    })
    res_df = bind_rows(res1, .id="set")
    res_df = left_join(res_df, song_syls, by="id")
    res_df = left_join(res_df, info, by="wav_strip")
    res_df$rep = i
    res_df
    #}
  }, mc.cores = 10)
  #})
  res2_df = bind_rows(res2)
  res2_df
}

calc_lik_auto = function(bird, 
                         reference_date, 
                         peak_source = "auto",
                         max_songs_per_day = 50,
                         recalculate_psds = T,
                         recalculate_models = T,
                         
                         # PSD params
                         wl=512,
                         freq_limits=c(1,10),
                         feature_set = "mid_mean",
                         
                         # Distance calculation
                         distmethod = "euclidean",
                         power = 2,
                         reference_size = 20,
                        
                         # Modeling params
                         range_to_test = 5:10,
                         test_fraction = .1,
                         train_number = NULL,
                         nreps = 10,
                         
                         # Plotting params
                         plot_dir = NULL,
                         
                         # Computation params
                         ncores = 6,
                         
                         # Output
                         output_db_name = "data.db"
) {
  require(mvtnorm)
  set.seed(1)
  dir = paste("/mnt/bengal_home/song", bird, "songs", sep="/")
  dir1 = paste(dir, "psd", sep="/")
  #info = load_mat_info_from(dir, file_ex="wav$")
  info = load_song_info_from_db("/mnt/bengal_home/song/song_files.db", bird, local=T)
  info = info %>% mutate(wav_base = basename(wav), 
                         wav_strip =str_replace(wav_base, ".wav", ""))
  
  if (!dir.exists(dir1))
    dir.create(dir1)
  

  fname = paste(paste(feature_set, 
                      distmethod, 
                      sprintf("power%s", power), 
                      sprintf("refsize%s", reference_size), sep="-"), ".rds", sep="")
  
  model_fname = paste(plot_dir, "models.rds", sep="/")
  if (recalculate_models | !file.exists(model_fname)) {
  song_data = NULL
  psd_fname = paste(dir1, fname, sep="/")
  if ((recalculate_psds | !file.exists(psd_fname))) {
    cat("  Calculating PSDs...\n")
    song_data = calc_song_psds(bird,
                               peak_source = peak_source,
                               max_songs_per_day=max_songs_per_day,
                               wl=wl,
                               freq_limits=freq_limits,
                               feature_set=feature_set,
                               distmethod=distmethod,
                               power=power,
                               reference_size=reference_size)
    saveRDS(song_data, psd_fname)
  } else {
    cat("  Loadings PSDs...\n")
    song_data = readRDS(psd_fname)
  }
  
  song_syls = song_data$syls
  song_syls = song_syls %>% 
    mutate(wav_strip =str_replace(id, "-[0-9]+", ""))
  
  song_psds = song_data$syl_data
  song_psds = song_psds[rownames(song_psds) %in% song_syls$id,]
  song_psds_wav = str_replace(rownames(song_psds), "-[0-9]+", "")
  
  base_wavs  = info %>% filter(date<=reference_date)
  base_syls_ids = song_syls %>% filter(wav_strip %in% base_wavs$wav_strip)
  cat(sprintf("Plot directory: %s\n", plot_dir))
  #print(dim(base_wavs))
  cat(sprintf("   Number of baseline syllables: %s\n", nrow(base_syls_ids)))
  ## Find best number of mixtures
  cat("  Calculating distances/Fitting models..\n")
    #mcs_stats = lapply(1:nreps, function(i) {
    mcs_stats = mclapply(1:nreps, function(i) { 
      # Calculate inter-syllable distances
      cat(sprintf("   Rep %s\n", i))
      set.seed(i)
      
      ## Define data sets
      if (is.null(train_number)) {
        train_wavs = base_wavs %>% sample_frac(size = (1-test_fraction))
      
        suppressMessages({test_base_wavs = anti_join(base_wavs, train_wavs)})
        suppressMessages({test_rest_wavs = info %>% anti_join(base_wavs)})
        #print(dim(test_base_wavs))
        #print("post segment")
        train_syls = song_syls[song_syls$wav_strip %in% train_wavs$wav_strip,]
        train_psds = song_psds[song_psds_wav %in% train_wavs$wav_strip,]
        #   print(dim(test_syls))
        #print("post filter")
        test_base_syls = song_syls[song_syls$wav_strip %in% test_base_wavs$wav_strip,]
        test_base_psds = song_psds[song_psds_wav %in% test_base_wavs$wav_strip,]
        #print(dim(test_base_syls))
        #print("post filter1")
        test_rest_syls = song_syls[song_syls$wav_strip %in% test_rest_wavs$wav_strip,]
        test_rest_psds = song_psds[song_psds_wav %in% test_rest_wavs$wav_strip,]
        #print("post fitler2")
      } else {
        if (nrow(base_syls_ids) >= train_number) {
          train_ids = base_syls_ids %>% sample_n(size = train_number)
        } else {
          train_ids = base_syls_ids
        }

        suppressMessages({test_base_ids = anti_join(base_syls_ids, train_ids)})
        suppressMessages({test_rest_ids = song_syls %>% anti_join(base_syls_ids)})
        #print(dim(test_base_wavs))
        #print("post segment")
        train_syls = song_syls[song_syls$id %in% train_ids$id,]
        train_psds = song_psds[rownames(song_psds) %in% train_ids$id,]
        #   print(dim(test_syls))
        #print("post filter")
        test_base_syls = song_syls[song_syls$id %in% test_base_ids$id,]
        test_base_psds = song_psds[rownames(song_psds) %in% test_base_ids$id,]
        #print(dim(test_base_syls))
        #print("post filter1")
        test_rest_syls = song_syls[song_syls$id%in% test_rest_ids$id,]
        test_rest_psds = song_psds[rownames(song_psds) %in% test_rest_ids$id,]
      }
      psds_list = list(train=train_psds, test_base=test_base_psds, test_rest=test_rest_psds)
      
      ## Calculate distances
      #print("distances")
      ref_inds = sample(1:nrow(train_psds), reference_size)
      ref_data = train_psds[ref_inds,]
      train_dist = dist_syllable_data(train_psds, compare_data = train_psds, 
                                      distmethod=distmethod, power=power, selectCols=ref_inds)
      
      test_base_dist = dist_syllable_data(test_base_psds, compare_data = train_psds, 
                                          distmethod=distmethod, power=power, selectCols=ref_inds)
      
      test_rest_dist = dist_syllable_data(test_rest_psds, compare_data = train_psds,
                                          distmethod=distmethod, power=power, selectCols=ref_inds)
      
      dist_list = list(train=train_dist, test_base=test_base_dist, test_rest=test_rest_dist)
      
      if (i == 1) {
        cat(sprintf("    Number training syllables: %s\n", nrow(train_dist)))
        cat(sprintf("    Number of held training syllables: %s\n", nrow(test_base_dist)))
        cat(sprintf("    Number of test syllables: %s\n", nrow(test_rest_dist)))
      }
      
      ## Plot out syllable x syllable distance heatmaps
      if (!is.null(plot_dir) & i==1) {
        plot_sylsyl_heatmaps3(dist_list, 
                              c("Baseline - Training", 
                                "Baseline - Held out",
                                "Post"),
                              plot_dir)
      }

      ## Test range of mixtures
      mc_stats = NULL
      mc_stats = mclust_par(train_dist, G=range_to_test, modelNames="VVI", plot=F, parallel=F)
      
      if (!is.null(plot_dir) & i == 1) {
        lapply(mc_stats$mcs, function(x) {
          plot_fname = paste(plot_dir, sprintf("psds_G%s.pdf", x$G), sep="/")
          plot_example_psds(x, train_psds, plot_fname)
          
          plot_fname3d = paste(plot_dir, sprintf("3dpsds_G%s.pdf", x$G), sep="/")
          plot_example_3d_psds(x, song_data, plot_fname3d)
        })
      }
      
      if (is.null(mc_stats[[1]])) 
        return(list(mc_stats = NULL, dist_list=NULL))
      
      mc_stats$stats$rep = i
      return(list(mc_stats=mc_stats, dist_list=dist_list))
    }, mc.cores=ncores)
    #})
    saveRDS(mcs_stats, model_fname)
  } else {
    mcs_stats = readRDS(model_fname)
  }
  
  mcs_stats = mcs_stats[unlist(lapply(mcs_stats, function(x) !is.null(x$mc_stats)))]
  data_stats = lapply(mcs_stats, function(x) x$mc_stats$stats)
  data_stats = do.call("rbind", data_stats)
  
  if (!is.null(plot_dir)) 
    cairo_pdf(paste(plot_dir, "component_estimation.pdf", sep="/"), width=9, height=9)
  
  data_stats_mean = plot_cluster_stats(data_stats)
  
  if (!is.null(plot_dir)) 
    dev.off()
  
  max_bic_ind = na.omit(data_stats_mean) %>% filter(stat=="BIC") %>% filter(diffs2==min(diffs2)) 
  max_bic_ind = round(mean(max_bic_ind$G))
  selected_model = which(range_to_test == max_bic_ind)
  cat(sprintf("  Selected number of mixtures: %s\n", range_to_test[selected_model]))
  
  ## Run replicates
  cat("  Calculating LL...\n")
  #res2 = foreach(i=1:nreps) %do% {
  res2 = mclapply(1:nreps, function(i) { 
  #res2 = lapply(1:nreps, function(i) {
    
    ## Calculate loglikelihoods
    dist_list = mcs_stats[[i]]$dist_list
    mod = mcs_stats[[i]]$mc_stats$mcs[[selected_model]]
    params = mod$parameters
    
    res1 = lapply(dist_list[2:3], function(p) {
      d = calc_likelihood(p, params)
    })
    res_df = bind_rows(res1, .id="set")
    res_df$rep = i
    res_df
    #}
  }, mc.cores = ncores)
  #})
  res2_df = bind_rows(res2)
  res2_df = res2_df %>% mutate(wav_strip = str_replace(id, "-[0-9]+", ""))
  #res2_df = left_join(res2_df, song_syls, by="id")
  res2_df = left_join(res2_df, info %>% select(wav_strip, bird, wav))
  output_db = src_sqlite(output_db_name, create=T)
  replace_by_key(output_db, table_name = "like", res2_df, "bird")
  res2_df
}

calc_lik_auto2 = function(bird, 
                         reference_dates1,
                         reference_dates2,
                         peak_source = "auto",
                         max_songs_per_day = 50,
                         recalculate_psds = T,
                         recalculate_models = T,
                         
                         # PSD params
                         wl=512,
                         freq_limits=c(1,10),
                         feature_set = "mid_mean",
                         
                         # Distance calculation
                         distmethod = "euclidean",
                         power = 2,
                         reference_size = 20,
                         
                         # Modeling params
                         range_to_test = 5:10,
                         test_fraction = .1,
                         train_number = NULL,
                         nreps = 10,
                         
                         # Plotting params
                         plot_dir = NULL,
                         
                         # Computation params
                         ncores = 6,
                         
                         # Output
                         output_db_name = "data.db"
) {
  require(mvtnorm)
  set.seed(1)
  dir = paste("/mnt/bengal_home/song", bird, "songs", sep="/")
  dir1 = paste(dir, "psd", sep="/")
  #info = load_mat_info_from(dir, file_ex="wav$")
  info = load_song_info_from_db("/mnt/bengal_home/song/song_files.db", bird, local=T)
  info = info %>% mutate(wav_base = basename(wav), 
                         wav_strip =str_replace(wav_base, ".wav", ""))
  
  if (!dir.exists(dir1))
    dir.create(dir1)
  
  
  fname = paste(paste(feature_set, 
                      distmethod, 
                      sprintf("power%s", power), 
                      sprintf("refsize%s", reference_size), sep="-"), ".rds", sep="")
  
  model_fname = paste(plot_dir, "base_models.rds", sep="/")
  
  
  if (recalculate_models | !file.exists(model_fname)) {
    song_data = NULL
    
    # Calculate PSDS ----------------------------------------------------------------------
    psd_fname = paste(dir1, fname, sep="/")
    if ((recalculate_psds | !file.exists(psd_fname))) {
      cat("  Calculating PSDs...\n")
      song_data = calc_song_psds(bird,
                                 peak_source = peak_source,
                                 max_songs_per_day=max_songs_per_day,
                                 wl=wl,
                                 freq_limits=freq_limits,
                                 feature_set=feature_set,
                                 distmethod=distmethod,
                                 power=power,
                                 reference_size=reference_size)
      saveRDS(song_data, psd_fname)
    } else {
      cat("  Loadings PSDs...\n")
      song_data = readRDS(psd_fname)
    }
    
    song_syls = song_data$syls
    song_syls = song_syls %>% 
      mutate(wav_strip =str_replace(id, "-[0-9]+", ""))
    
    song_psds = song_data$syl_data
    song_psds = song_psds[rownames(song_psds) %in% song_syls$id,]
    song_psds_wav = str_replace(rownames(song_psds), "-[0-9]+", "")
    
    base_wavs  = info %>% filter(date>=reference_dates1[1], date<=reference_dates1[2])
    base_syls_ids = song_syls %>% filter(wav_strip %in% base_wavs$wav_strip)
    
    compare_wavs = info %>% filter(date>=reference_dates2[1], date<=reference_dates2[2])
    compare_syls_ids = song_syls %>% filter(wav_strip %in% compare_wavs$wav_strip)
    
    cat(sprintf("  Plot directory: %s\n", plot_dir))

    cat(sprintf("   Number of baseline syllables: %s\n", nrow(base_syls_ids)))
    cat(sprintf("   Number of compare syllables: %s\n", nrow(compare_syls_ids)))
    

    # Fit baseline models ------------------------------------------------------------
    cat("  Calculating distances. Fitting baseline models..\n")
    mcs_stats = mclapply(1:nreps, function(i) { 
    #mcs_stats = lapply(1:nreps, function(i) {
      cat(sprintf("   Rep %s\n", i))
      set.seed(i)
      
      ## Define data sets -------------------------------------------------------------
      if (is.null(train_number)) {
        base_train_wavs = base_wavs %>% sample_frac(size = (1-test_fraction))
        compare_train_wavs = compare_wavs %>% sample_frac(size = (1-test_fraction))
        
        suppressMessages({base_test_wavs = anti_join(base_wavs, base_train_wavs)})
        suppressMessages({compare_test_wavs = anti_join(compare_wavs, compare_train_wavs)})
      
        base_train_syls = song_syls[song_syls$wav_strip %in% train_wavs$wav_strip,]
        base_train_psds = song_psds[song_psds_wav %in% train_wavs$wav_strip,]
        base_test_syls = song_syls[song_syls$wav_strip %in% base_test_wavs$wav_strip,]
        base_test_psds = song_psds[song_psds_wav %in% base_test_wavs$wav_strip,]
      
        compare_train_syls = song_syls[song_syls$wav_strip %in% compare_train_wavs$wav_strip,]
        compare_train_psds = song_psds[song_psds_wav %in% compare_train_wavs$wav_strip,]
        compare_test_syls = song_syls[song_syls$wav_strip %in% compare_test_wavs$wav_strip,]
        compare_test_psds = song_psds[song_psds_wav %in% compare_test_wavs$wav_strip,]
        #print("post fitler2")
      } else {
        if (nrow(base_syls_ids) >= train_number) {
          base_train_ids = base_syls_ids %>% sample_n(size = train_number)
        } else {
          base_train_ids = base_syls_ids
        }
        
        if (nrow(compare_syls_ids) >= train_number) {
          compare_train_ids = compare_syls_ids %>% sample_n(size = train_number)
        } else {
          compare_train_ids = compare_syls_ids
        }
        
        suppressMessages({base_test_ids = anti_join(base_syls_ids, base_train_ids)})
        suppressMessages({compare_test_ids = anti_join(compare_syls_ids, compare_train_ids)})
     
        base_train_syls = song_syls[song_syls$id %in% base_train_ids$id,]
        base_train_psds = song_psds[rownames(song_psds) %in% base_train_ids$id,]
        base_test_syls = song_syls[song_syls$id %in% base_test_ids$id,]
        base_test_psds = song_psds[rownames(song_psds) %in% base_test_ids$id,]
      
        compare_train_syls = song_syls[song_syls$id %in% compare_train_ids$id,]
        compare_train_psds = song_psds[rownames(song_psds) %in% compare_train_ids$id,]
        compare_test_syls = song_syls[song_syls$id%in% compare_test_ids$id,]
        compare_test_psds = song_psds[rownames(song_psds) %in% compare_test_ids$id,]
      }
      psds_list = list(base_train=base_train_psds, 
                       base_test=base_test_psds, 
                       compare_train=compare_train_psds,
                       compare_test=compare_test_psds)
      
      ## Calculate distances -------------------------------------------------------------------
      ref_inds = sample(1:nrow(base_train_psds), reference_size)
      ref_data = base_train_psds[ref_inds,]
      
      dist_list = lapply(psds_list, function(p) {
        dist_syllable_data(p, 
                           compare_data = base_train_psds, 
                           distmethod=distmethod, 
                           power=power, 
                           selectCols=ref_inds)
      })
      
      if (i == 1) {
        cat(sprintf("    Number baseline training syllables: %s\n", nrow(dist_list$base_train)))
        cat(sprintf("    Number baseline held syllables: %s\n", nrow(dist_list$base_test)))
        cat(sprintf("    Number compare training syllables: %s\n", nrow(dist_list$compare_train)))
        cat(sprintf("    Number compare held syllables: %s\n", nrow(dist_list$compare_test)))
      }
      
      ## Plot out syllable X syllable distance heatmaps ---------------------------------------------------
      if (!is.null(plot_dir) & i==1) {
        plot_sylsyl_heatmaps4(dist_list, 
                              c("Baseline - Training", 
                                "Baseline - Held out",
                                "Compare - Training",
                                "Compare - Held out"),
                              plot_dir)
      }
      
      ## Test range of model components ------------------------------------------------------------------
      mc_stats = NULL
      mc_stats = mclust_par(dist_list$base_train, G=range_to_test, modelNames="VVI", plot=F, parallel=F)
      
      ## Plot PSDs ----------------------------------------------------------------------------------------
      if (!is.null(plot_dir) & i == 1) {
        lapply(mc_stats$mcs, function(x) {
          plot_fname = paste(plot_dir, sprintf("psds_G%s.pdf", x$G), sep="/")
          plot_example_psds(x, base_train_psds, plot_fname)
          
          plot_fname3d = paste(plot_dir, sprintf("3dpsds_G%s.pdf", x$G), sep="/")
          plot_example_3d_psds(x, song_data, plot_fname3d)
        })
      }
      
      if (is.null(mc_stats[[1]])) 
        return(list(mc_stats = NULL, dist_list=NULL))
      
      mc_stats$stats$rep = i
      return(list(mc_stats=mc_stats, dist_list=dist_list))
    }, mc.cores=ncores)
    #})
    saveRDS(mcs_stats, model_fname)
  } else {
    mcs_stats = readRDS(model_fname)
  }
  
  # Pick number of components ---------------------------------------------------------------------
  mcs_stats = mcs_stats[unlist(lapply(mcs_stats, function(x) !is.null(x$mc_stats)))]
  data_stats = lapply(mcs_stats, function(x) x$mc_stats$stats)
  data_stats = do.call("rbind", data_stats)
  
  if (!is.null(plot_dir)) 
    cairo_pdf(paste(plot_dir, "component_estimation.pdf", sep="/"), width=9, height=9)
  
  data_stats_mean = plot_cluster_stats(data_stats)
  
  if (!is.null(plot_dir)) 
    dev.off()
  
  max_bic_ind = na.omit(data_stats_mean) %>% filter(stat=="BIC") %>% filter(diffs2==min(diffs2)) 
  max_bic_ind = round(mean(max_bic_ind$G))
  selected_model = which(range_to_test == max_bic_ind)
  cat(sprintf("  Selected number of mixtures: %s\n", range_to_test[selected_model]))
  
  # Train compare models --------------------------------------------------------------------------
  model2_fname = paste(plot_dir, "compare_models.rds", sep="/")
  if (recalculate_models | !file.exists(model2_fname)) {
    cat("  Training compare model...\n")
    plot_compare_dir = paste(plot_dir, "compare", sep="/")
    if (!dir.exists(plot_compare_dir))
      dir.create(plot_compare_dir)
    
    #compare_mcs = lapply(1:nreps, function(i) {
    compare_mcs = mclapply(1:nreps, function(i) { 
      cat(sprintf("   Rep %s\n", i))
      set.seed(i)
      
      cur_dist_list = mcs_stats[[i]]$dist_list
      
      mc_stats = NULL
      mc = Mclust(cur_dist_list$compare_train, G=range_to_test[selected_model], modelNames="VVI")
      
      if (!is.null(plot_dir) & i == 1) {
        plot_fname = paste(plot_compare_dir, sprintf("psds_G%s.pdf", mc$G), sep="/")
        plot_example_psds(mc, song_psds, plot_fname)
          
        plot_fname3d = paste(plot_compare_dir, sprintf("3dpsds_G%s.pdf", mc$G), sep="/")
        plot_example_3d_psds(mc, song_data, plot_fname3d)
      }
      
      #if (is.null(mc_stats[[1]])) 
      #  return(list(mc_stats = NULL))

      return(list(mc_stats=mc))
    }, mc.cores=ncores)
    #})
    saveRDS(compare_mcs, model2_fname)
  } else {
    compare_mcs = readRDS(model2_fname)
  }
  
  # Calc log-likelihoods -----------------------------------------------------------
  cat("  Calculating LL...\n")
  #res2 = mclapply(1:nreps, function(i) { 
  res2 = lapply(1:nreps, function(i) {
    dist_list = mcs_stats[[i]]$dist_list
    mod1 = mcs_stats[[i]]$mc_stats$mcs[[selected_model]]
    mod2 = compare_mcs[[i]]$mc_stats
    
    params1 = mod1$parameters
    params2 = mod2$parameters
    
    res1 = lapply(dist_list[grep("test", names(dist_list))], function(p) {
      if (nrow(p) > 1) {
        d1 = calc_likelihood(p, params1) %>% rename(lik1 = lik)
        d2 = calc_likelihood(p, params2) %>% rename(lik2 = lik)
        left_join(d1, d2, by="id")
      } else{
        data.frame(id="dummy", lik1=NA, lik2=NA)
      }
      
    })
    res_df = bind_rows(res1, .id="set")
    res_df$rep = i
    res_df
    #}
  #}, mc.cores = ncores)
  })
  res2_df = bind_rows(res2)
  
  scores = res2_df %>% group_by(set, id) %>% summarize(kl12=calc_kl_divergence(lik1, lik2))
  scores = scores %>% mutate(wav_strip = str_replace(id, "-[0-9]+", ""))
  scores = left_join(scores, info %>% select(wav_strip, bird, wav), by="wav_strip")
  scores$bird[is.na(scores$bird)] = bird
  
  # Write out data -----------------------------------------------------------------------------
  output_db = src_sqlite(output_db_name, create=T)
  replace_by_key(output_db, table_name = "like", scores, "bird")
  
  return(scores)
}

calc_song_psds = function(bird,
                               max_songs_per_day,
                               peak_source = "mat",
                               wl=512,
                               freq_limits=c(1,10),
                               feature_set="mid_mean",
                               distmethod="euclidean",
                               power=2,
                               reference_size=20) {
  dir = paste("/mnt/bengal_home/song", bird, "songs", sep="/")
  info = load_mat_info(dir, file_ex="wav$")
  info = info %>% group_by(date) %>% do ({
    d = . 
    nsize =  ifelse(length(d$date)<=max_songs_per_day, length(d$date), max_songs_per_day)
    ind = sample(1:nrow(d), nsize, replace=F)
    d[ind,]
  })
  
  cat(sprintf("   Total number of songs: %s\n", nrow(info)))
  songs = parse_song_batch2(info$wav, peak_source=peak_source, thresh_range = seq(-5, 5, .5))
  
  l = unlist(sapply(songs[[1]], function(x) nrow(x$syllable)))
  cat(sprintf("   Average number of syllables / song: %s\n", round(mean(l))))
  ## Get song features
  process_song_batch(songs, 
                     feature_set=feature_set, 
                     wl=wl, 
                     smoothed=F, 
                     cluster=F, 
                     return_data=T, 
                     freq_limits=freq_limits)
}

plot_sylsyl_heatmaps3 = function(dist_list, dist_names, plot_dir=NULL) {
  dists = lapply(1:length(dist_list), function(i) {
    list(data=adjacency(t(dist_list[[i]]), type="signed"), 
         title=dist_names[i])
    })
  
  suppressMessages(require(NMF))
  ind = 40
  
  if (!is.null(plot_dir))
  cairo_pdf(paste(plot_dir, "auto_heatmaps.pdf", sep="/"), width=12, height=5)
  
  par(mfrow=c(1,3), mar=c(0, 1, 1, 1))
  for (i in 1:length(dists)) {
    suppressMessages(aheatmap(dists[[i]]$data[1:ind,1:ind], 
                              hclustfun = "ward", 
                              color="topo", 
                              main=dists[[i]]$title))
  }

  if (!is.null(plot_dir))
    dev.off()
}

plot_sylsyl_heatmaps4 = function(dist_list, dist_names, plot_dir=NULL) {
  dists = lapply(1:length(dist_list), function(i) {
  
    d = NULL
    if (nrow(dist_list[[i]]) > 1) {
    d = adjacency(t(dist_list[[i]]), type="signed")
    } 
    list(data=d, 
         title=dist_names[i])
  })
  
  suppressMessages(require(NMF))
  ind = 40
  
  if (!is.null(plot_dir))
    cairo_pdf(paste(plot_dir, "auto_heatmaps.pdf", sep="/"), width=12, height=12)
  
  par(mfrow=c(2,2), mar=c(0, 1, 1, 1))
  for (i in 1:length(dists)) {
    if (!is.null(dists[[i]]$data)) {
      suppressMessages(aheatmap(dists[[i]]$data[1:ind,1:ind], 
                                hclustfun = "ward", 
                                color="topo", 
                                main=dists[[i]]$title))
    }
  }
  
  if (!is.null(plot_dir))
    dev.off()
}


plot_example_psds = function(mc, psds, plot_fname) {
  cairo_pdf(plot_fname, width=12, height=5)
 
  classification = mc$classification
  classification = na.omit(classification[match(rownames(psds), names(classification))])
  psds1 = psds[match(names(classification), rownames(psds)),]
  psds_m = melt(psds1)
  psds_m$class = classification[match(psds_m$Var1, names(classification))]
  gg = ggplot(psds_m, aes(Var2, value))
  gg = gg + geom_line(aes(group=Var1), alpha=I(1/10))
  gg = gg + stat_summary(geom="line", fun.y="mean")
  gg = gg + facet_wrap(~class)
  gg = gg + labs(x="", y ="")
  print(gg)
  dev.off()
}

plot_example_3d_psds = function(mc, song_data, plot_fname, num_syls=5) {
  wav = song_data$wav
  wav_names = str_replace(basename(names(wav)), ".wav", "")
  #selected_wav = sample(1:length(wav_names), 1)
  syls = song_data$syls
  syls = syls %>% mutate(wav_strip = str_replace(id, "-[0-9]+", ""))
  classification = mc$classification
  syls$called = classification[match(syls$id, names(classification))]
  syls = na.omit(syls)

  nrow= num_syls
  ncol = length(unique(classification))
  selected_syls = na.omit(syls %>% group_by(called) %>% sample_n(num_syls)) %>% mutate(ind=1:num_syls) %>% arrange(ind)
  plots = lapply(1:nrow(selected_syls), function(i) {
    subregion = selected_syls[i,1:2]
    wav = song_data$wav[[which(wav_names==selected_syls$wav_strip[i])]]
    try_default({
    suppressMessages(plot_single_syl(wav, subregion))
    }, ggplot())
  })
# 
#   x_title = plots[[1]]$labels$x
#   y_title = plots[[1]]$labels$y
#   
  plots = lapply(plots, function(p) {
    p + labs(x="", y="") + theme(axis.text.y=element_blank())
  })
  
  # a = do.call(grid.arrange, c(c(plots, left=tex, ncol=ncol, nrow=nrow))#, left=textGrob(y_title, rot=90), bottom=textGrob(x_title)))
  # grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), 
  #                          p2 + theme(legend.position="none"),
  #                          p3 + theme(legend.position="none"),
  #                          p4 + theme(legend.position="none"), 
  #                          nrow = 2,
  #                          top = textGrob("Main Title", vjust = 1, gp = gpar(fontface = "bold", cex = 1.5)),
  #                          left = textGrob("Global Y-axis Label", rot = 90, vjust = 1)), 
  #              legend, 
  cairo_pdf(plot_fname, width=12, height=8)
  do.call(grid.arrange, c(plots, ncol=ncol, nrow=nrow))
  dev.off()
}

find_clusters_by_hclust = function(dist_mat) {
  require(WGCNA)
  hc = hclust(as.dist(dist_mat), method = "ward.D2")
  dynamicMods = cutreeDynamic(dendro = hc, 
                              distM = dist_mat, 
                              method = "hybrid",
                              deepSplit = 2, 
                              pamRespectsDendro = TRUE,
                              minClusterSize = 5)
  return(length(unique(dynamicMods)))
}

calc_likelihood = function(data, params) {
  norms = sapply(1:ncol(params$mean), function(i) {
    dmvnorm(x=data, mean=params$mean[,i], sigma=params$variance$sigma[,,i])
  })
  norms1 = rowSums(sweep(norms, 1, params$pro, "*"))
  res = data.frame(id=rownames(norms), lik=log(norms1))
  res$lik[is.infinite(res$lik)] = min(res$lik[res$lik>0])
  res
}

calc_kl_divergence = function(lik1, lik2) {
  score = log2(exp(1))*(mean(lik1)-mean(lik2))                                                                                                                                                     
  #score = score/len1                                                                                                                                                                                                  
  score                                                                                                                                                                                              
}

mclust_syllable_data = function(syls, data, models=NULL, range_to_test=c(4:15),sample_size=120, feature_size=20,
                                pca=FALSE, distmethod="cor", power=2, plot=FALSE) {
  
  mc = NULL
  ref_mat = NULL
  inds = NULL
  mcs_stats = NULL
  
  ## Train models
  if (nrow(data) > sample_size) {
    if (is.null(models)) 
      ## Train new models
    {
      ## Subsample data
      reps = 20
      inds = lapply(1:reps, function(x) sample(1:nrow(data), sample_size))
      feature_inds = lapply(inds, function(ind) sample(1:length(ind), feature_size))
      #feature_inds = lapply(1:reps, function(x) sample(1:nrow(data), feature_size))
      ## Train several models on subsets of data
      mcs_stats = lapply(1:reps, function(x) {
        # Calculate inter-syllable distances
        print(paste("rep", x, sep=""))
        #data1_s = dist_syllable_data(data[inds[[x]],], distmethod=distmethod, power=power)
        data1_s = dist_syllable_data(data[inds[[x]],], selectCols=feature_inds[[x]], distmethod=distmethod, power=power)
        data1_s = matrix(data1_s, nrow=nrow(data1_s), ncol=ncol(data1_s))
        data2_s = as.matrix(dist(data1_s))
        mc_stats = NULL
        if (length(range_to_test) > 1) {
          mc_stats = mclust_par(data2_s, G=range_to_test, modelNames="VVI", plot=F)
          if (is.null(mc_stats[[1]])) 
            return(list(syls = NULL, mc=NULL))
          mc_stats$stats$rep = x
        } else {
          #mc_stats = Mclust(data1_s, G=range_to_test, modelNames="VVI")
          mc_stats = Mclust(data2_s, G=range_to_test)
          mc_stats$stats$rep = x
        }
        mc_stats
      })
      
      ### FIXME: don't return, continue into next
      if (length(range_to_test) > 1) {
        data_stats = lapply(mcs_stats, function(x) x$stats)
        data_stats = do.call("rbind", data_stats)
        data_stats_mean = plot_cluster_stats(data_stats)
        
        max_bic_ind = data_stats_mean%>% filter(diffs2==min(diffs2)) 
        max_bic_ind = round(mean(max_bic_ind$G))
        models = lapply(mcs_stats, function(x) x$mcs[[which(range_to_test == max_bic_ind)]])
        return(list(syls = syls, mc = models))
      } 
      names(mcs_stats) = paste("rep", 1:reps, sep="-")
      #ind1 = 1:reps
      #ind = NULL
      #bics = unlist(lapply(mcs_stats, function(mcs_stat) mcs_stat$bic))
      #mc = mcs_stats[[which.max(bics)]]
      #ind = inds[[which.max(bics)]]
      train_data = data
    } 
    else 
    {
      mcs_stats = models$models
      ref_mat = models$ref_mat
      inds = models$inds
      train_data = models$train_data
    }
    ## Predict classification from each model
    preds = foreach (i=1:length(mcs_stats)) %dopar% {
      #  for (i in 1:length(mcs_stats)) {
      if(is.null(mcs_stats[[i]]))
        return(NULL)
      data_s = dist_syllable_data(data, compare_data=train_data, distmethod=distmethod, power=power, selectCols=inds[[i]])
      mc.pred = predict(mcs_stats[[i]], data_s)
      
      classification = unlist(mc.pred$classification)
      
      ## Compute average/median spectral profiles for classes
      data1 = data[match(names(classification), rownames(data)),]
      data_split = split.data.frame(data1, classification)
      data_split1 = do.call("rbind", lapply(data_split, function(d) apply(d, 2, median)))
      data_split1 = sweep(data_split1, 1, apply(data_split1, 1, max), "/" )
      list(classification=classification, avg_syls=data_split1)
      #unlist(mc.pred$classification)
    }  
    preds = lapply(preds, function(x) x$classification)
    preds = preds[unlist(lapply(preds, function(x) !is.null(x)))]
    avg_syls = lapply(preds, function(x) x$avg_syls)
    avg_syls = avg_syls[unlist(lapply(avg_syls, function(x) !is.null(x)))]
    
    num_models = length(preds)
    
    ## Populate matrix with model predictions
    pred_mat = matrix(nrow=nrow(data), ncol=num_models, dimnames=list(rownames(data)))
    for (i in 1:ncol(pred_mat)) {
      pred_mat[,i] = preds[[i]][match(rownames(pred_mat), names(preds[[i]]))]
    }
    
    ## Generate matrix aligning cluster labels
    if (is.null(ref_mat)) {
      ref_mats = list()
      #ref_mats = foreach (i=1:(num_models-1), .combine="rbind") %do% {
      #ref_mats = foreach (i=1:num_models, .combine="rbind") %do% {
      ref_mats = foreach(i=1:num_models) %do% {
        # Loop through models
        ref_mat = matrix(nrow=range_to_test, ncol=num_models)
        ref_mat[,i] = 1:range_to_test # define selected model as reference
        for (j in setdiff(1:ncol(ref_mat), i)) {
          tab = table(pred_mat[,i], pred_mat[,j])
          ref_mat[,j] = apply(tab, 1, which.max)
        }
        ref_mat
      }
      uniques = lapply(ref_mats, function(d) apply(d, 2, function(e) length(unique(e))))
      uniques_sum = unlist(lapply(uniques, sum))
      ref_mat = ref_mats[[which.max(uniques_sum)]]
      
      #        foreach (j=(i+1):(num_models), .combine="rbind") %do% {
      # #         foreach (j=1:num_models, .combine="rbind") %do% {
      # #           dists = dist(avg_syls[[i]], avg_syls[[j]], method="kullback")
      # #           dists = matrix(dists, ncol=range_to_test, nrow=range_to_test)
      # #           dists_m = melt(as.data.frame(dists))
      # #           colnames(dists_m) = c("set1", "value")
      # #           dists_m$set1 = str_replace(dists_m$set1, "V", "")
      # #           dists_m$set2 = 1:range_to_test
      # #           dists_m$model_num1 = i
      # #           dists_m$model_num2 = j
      # #           dists_m$syl1 = paste("model", i, "-syl", dists_m$set1, sep="")
      # #           dists_m$syl2 = paste("model", j, "-syl", dists_m$set2, sep="")
      # #           dists_m
      #           #conf_mat = table(pred_mat[,1], pred_mat[,i])
      #         #  ref_mat[,j] = apply(sims, 1, which.min) # align cluster assignments by maximum tabulation
      #         }
      #        # ref_mats[[i]] = ref_mat
      #       }
      #       ref_mats1 = ref_mats %>% filter(value<.2, value>0)
      #       ref_mats2 = ref_mats1 %>% group_by(syl1) %>% summarize()
      #  ref_mats2 = acast(ref_mats, syl1~syl2, value.var="value")
      #       mc = mclust(ref_mats2, G=(range_to_test-2):(range_to_test+2))
      
      #       avg_syls_comb = do.call("rbind", avg_syls)
      #        avg_syls_comb = data.frame(avg_syls_comb)
      #       info = expand.grid( syl=paste("syl", 1:range_to_test, sep=""), model=paste("model", 1:num_models, sep=""))
      #       info$modelsyl = paste(info$model, info$syl, sep="-")
      
      #       avg_syls_comb$cl = mc$classification[match(avg_syls_comb$modelsyl, names(mc$classification))]
      #       avg_syls_comb_m = melt(avg_syls_comb, id.vars=c("cl", "modelsyl"))
      #       avg_syls_comb_m$variable = as.numeric(str_replace(avg_syls_comb_m$variable, "X", ""))
      #       ggplot(avg_syls_comb_m, aes(variable, value, group=cl)) + geom_line(alpha=I(1/10)) + facet_wrap(~cl)
      
      #       rownames(avg_syls_comb) = info$modelsyl
      #       mc1 = Mclust(avg_syls_comb, G=4:10)
      
      #       avg_syls_comb$modelsyl = info$modelsyl
      #       avg_syls_comb_m = melt(avg_syls_comb, id.vars="modelsyl" )
      #       avg_syls_comb_m$cl = mc1$classification[match(avg_syls_comb_m$modelsyl, names(mc1$classification))]
      #      # avg_syls_comb_m$variable = as.numeric(str_replace(avg_syls_comb_m$variable, "X", ""))
      #       ggplot(avg_syls_comb_m, aes(variable, value, group=cl)) + geom_line(alpha=I(1/10)) + facet_wrap(~cl)
      #       ggplot(avg_syls_comb_m %>% filter(cl==6), aes(variable, value)) + geom_line() + facet_wrap(~modelsyl)
      #       
      #       km1 = kmeans(avg_syls_comb, centers=8)
      #       avg_syls_comb$modelsyl = info$modelsyl
      #       avg_syls_comb_m = melt(avg_syls_comb, id.vars=c("modelsyl"))
      #       avg_syls_comb_m$cl = km1$cluster
      #       ggplot(avg_syls_comb_m, aes(variable, value, group=cl)) + geom_line(alpha=I(1/10)) + facet_wrap(~cl)
      #       ggplot(avg_syls_comb_m %>% filter(cl==6), aes(variable, value)) + geom_line() + facet_wrap(~modelsyl)
      #       ref_mat = matrix(km1$cluster, nrow=range_to_test, ncol=num_models)
      #       ## Doesn't work well
      #       km1 = kmeans(ref_mats2, centers=8)
      #       avg_syls_comb$modelsyl = info$modelsyl
      #       avg_syls_comb_m = melt(avg_syls_comb, id.vars=c("modelsyl"))
      #       avg_syls_comb_m$cl = km1$cluster
      #       ggplot(avg_syls_comb_m, aes(variable, value, group=cl)) + geom_line(alpha=I(1/10)) + facet_wrap(~cl)
      #       ggplot(avg_syls_comb_m %>% filter(cl==6), aes(variable, value)) + geom_line() + facet_wrap(~modelsyl)
      
      #       avg_syls_comb = do.call("rbind", avg_syls)
      #       avg_syls_comb = data.frame(avg_syls_comb)
      #       avg_syls_comb_dist = dist_syllable_data(avg_syls_comb, distmethod="cor", power=2)
      #       hc1 = hclust(dist(avg_syls_comb_dist), method="ward")
      #        plot(hc1)
      #        hc1_cut = cutree(hc1, k = range_to_test)
      #       avg_syls_comb$modelsyl = info$modelsyl
      #       avg_syls_comb_m = melt(avg_syls_comb, id.vars=c("modelsyl"))
      #       avg_syls_comb_m$cl = hc1_cut
      #       ggplot(avg_syls_comb_m, aes(variable, value, group=cl)) + geom_line(alpha=I(1/10)) + facet_wrap(~cl)
      #       ggplot(avg_syls_comb_m %>% filter(cl==6), aes(variable, value)) + geom_line() + facet_wrap(~modelsyl)
      #       
      #       ref_mat = matrix(hc1_cut, nrow=range_to_test, ncol=num_models)
      
    }
    
    ## Standardize labelling
    pred_mat1 = pred_mat
    for (i in 1:ncol(pred_mat1)) {
      pred_mat1[,i] = ref_mat[match(pred_mat1[,i], ref_mat[,i]),1]
    }
    #     for (i in 1:ncol(pred_mat1)) {
    #       pred_mat1[,i] = ref_mat[pred_mat[,i],i]
    #     }
    pred_mat1[is.na(pred_mat1)] = 0
    consen_mat = matrix(0, nrow=nrow(pred_mat), ncol=range_to_test)
    for (i in 1:nrow(consen_mat)) {
      tab = table(pred_mat1[i,])
      consen_mat[i,] = tab[match(1:range_to_test, names(tab))]
    }
    consen_mat[is.na(consen_mat)] = 0
    consen_norm = sweep(consen_mat, 1, rowSums(consen_mat), "/")
    consen_norm[is.na(consen_norm)] = 0
    plot(density(apply(consen_norm, 1, max), adjust=1/2))
    cl = apply(consen_norm, 1, which.max)
    uncertain = apply(consen_norm, 1, function(x) x[which.max(x)]) < 0.6
    
    names(cl) = rownames(pred_mat1)
    cl[uncertain] = NA
    cl = unlist(cl)
    syls$called = factor(cl[match(syls$id, names(cl))], labels=(letters[1:length(table(cl))]))
    return(list(syls = syls, mc=mcs_stats, ref_mat=ref_mat, inds=inds, train_data=train_data))
  } else {
    data= sweep(data, 1, apply(data, 1, max), FUN = "/")
    data_s = adjacency(t(data), power = 2, type =  "signed", corFnc="cor")
    mc = Mclust(data_s, G=range_to_test, modelNames="VVI") 
    if (plot) {
      par(mfrow=c(1,1))
      #print(logLik(mc, data = data_s))
      plot(mc, what="BIC")
    }
    unique_class = unique(mc$classification)
    #mc$classification[abs(mc$uncertainty) > 0 ] = max(mc$classification) + 1
    
    mc$classification = factor(mc$classification, levels=sort(unique_class), labels=c(letters[1:(length(unique_class))]))
    syls$called = mc$classification
    return(list(syls = syls, mc=mc))
  }
}

det_slope = function(data, val, window=1) {
  res = vector("numeric", length=nrow(data))
  res[1:window] = 1
  res[length(res) - (window - 1)] = 1
  for (i in (window+1):(nrow(data) - window)) {
    res[i] = (data[i+window,val] - data[i-window,val]) / (window * 2 + 1)
  }
  data[,val] = unlist(res)
  return(data)
  #return(unlist(res))
}
spectral = function(data, centers=c(4:10)) {
  scs = lapply(centers, function(center) {
    sc = specc(data, centers=center, iterations=5000)
    data.k = kernelMatrix(sc@kernelf, data)
    sil = silhouette(sc, dmatrix=(1-data.k))
    list(sc = sc, sil = sil)
  })
  
  avg.widths = unlist(lapply(scs, function(x) {
    sil = x$sil
    sil$
      s = summary(x$sil)
    s$avg.width
  }))
  plot(centers, avg.widths)
  return(scs)
}

hclust_range = function(data, group=c(4:12)) {
  data.d = dist()
}

kmeans_range = function(data, range=c(2:10)) {
  d = lapply(range, function(x) {
    kmeans(data, centers = x)$cluster
  })
  d = do.call(cbind, d)
  colnames(d) = paste("k", range, sep="-")
  d
}
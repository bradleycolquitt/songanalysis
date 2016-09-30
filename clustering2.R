library(mclust)
library(foreach)
library(doMC)
library(proxy)
library(cluster)
library(caret)

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

mclust_par = function(data, G, modelNames="VVI", plot=T) {
  gc() 
  mcs = foreach(g=G) %do% {
    Mclust(data, G=g, modelNames=modelNames)
  }
  mcs = mcs[!unlist(lapply(mcs, is.null))]
  
  recovered_G = unlist(lapply(mcs, function(x) x$G))
  bics = unlist(lapply(mcs, function(x) x$BIC))
  #print(bics)
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
  suppressMessages(require(WGCNA))
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
  data_stats_mean = data_stats %>% group_by(G, stat) %>% summarize(value=mean(value), rep=1) 
  theme_set(theme_bw())
  g = ggplot(data_stats, aes(G, value, group=factor(rep))) + geom_line(alpha=I(1/2))
  g = g + geom_line(data=data_stats_mean, color="red") + geom_point(data=data_stats_mean, color="red")  + facet_wrap(~stat, scales="free_y")
  g = g + scale_x_discrete(breaks=c(unique(data_stats$G)))
  
  data_stats_mean = ungroup(data_stats_mean) %>% group_by(stat) %>% mutate(rel_diff_max = abs(value - max(value)) / max(value),
                                                                           diffs = c(0, diff(value)) / max(value),
                                                                           diffs2 = c(0, diff(diffs)))
  
  g1 = ggplot(data_stats_mean, aes(G, diffs)) + geom_point() + facet_wrap(~stat)
  g1 = g1 + scale_x_discrete(breaks=c(unique(data_stats$G)))
  
  g2 = ggplot(data_stats_mean, aes(G, diffs2)) + geom_point() + facet_wrap(~stat)
  g2 = g2 + scale_x_discrete(breaks=c(unique(data_stats$G)))
  
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
                             freq_limits=c(1,8),
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
#  source("/home/brad/src/songanalysis/song_util.R")
#  source("/media/data2/rstudio/birds/utils/stats.R")
#  source("~/data2/rstudio/birds/utils/db.R")
#  source("~/src/songanalysis/clustering2.R")
#  source("/home/brad/src/songanalysis/deafen_plot.R")
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
    test_base_wavs = anti_join(base_wavs, train_wavs)
    test_rest_wavs = info %>% anti_join(base_wavs)
    
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
      
      require(NMF)
      par(mfrow=c(1,3))
      ind = 50
      cairo_pdf(paste(plot_dir, "heatmaps.pdf"), width=12, height=8)
      aheatmap(train_dist1[1:ind,1:ind], hclustfun = "ward", color="topo", main="Training")
      aheatmap(test_base_dist1[1:ind,1:ind], hclustfun = "ward", color="topo", main="Baseline - held out")
      aheatmap(test_rest_dist1[1:ind,1:ind], hclustfun = "ward", color="topo", main="Post baseline")
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

calc_likelihood = function(data, params) {
  norms = sapply(1:ncol(params$mean), function(i) {
    dmvnorm(x=data, mean=params$mean[,i], sigma=params$variance$sigma[,,i])
  })
  norms1 = rowSums(sweep(norms, 1, params$pro, "*"))
  res = data.frame(id=rownames(norms), lik=log(norms1))
  res$lik[is.infinite(res$lik)] = min(res$lik[res$lik>0])
  res
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
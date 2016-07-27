source("~/src/songanalysis/song_util.R")
source("~/src/songanalysis/song_features.R")
source("~/src/songanalysis/song_ff.R")
source("/media/data2/rstudio/birds/tutor_monitor/tutor_util.R")

library(googlesheets)
library(tidyr)
options(stringsAsFactors=F)
servers = data.frame(server=c("tutor", "egret", "osprey"),
                     dir=c("/mnt/tutor_home/data", "/mnt/egret/brad/tutor", "/mnt/osprey/tutor"))

# MAIN -----------------------------------------------------------------------------------------------------------
load_tutoring_info = function() {
  sheet = gs_title("Foster Cage Status")  
  info = sheet %>% gs_read("tutoring")
  colnames(info) = tolower(colnames(info))
  colnames(info) = make.names(colnames(info))
  info = info %>% mutate(date.tutoring.begin = as.Date(date.tutoring.begin, format="%Y-%m-%d"))
  return(info)
}

run_subsong_analysis_batch = function(num_per_day=NULL, take_db_action=FALSE, type="deep") {
  info = load_tutoring_info()
  info = info %>% dplyr::filter(db_action!="skip")
  
  data = info %>% group_by(bird) %>% do({
    print(.$bird)
    if (type=="deep") {
      run_subsong_analysis(., .$bird, num_per_day=num_per_day, take_db_action=take_db_action)
    } else if (type=="shallow") {
      run_subsong_shallow_analysis(. , .$bird)
    }
    
  })
}

run_subsong_analysis = function(info, curr_bird, num_per_day = NULL, take_db_action=FALSE){
  info1 = info %>% dplyr::filter(bird==curr_bird)
  dir = paste(servers[servers[,1]==info1$server,"dir"], curr_bird, "subsongs", sep="/")
  files = list.files(dir, pattern = ".wav", full.names = T)
  
  ### Parse fname
  times = parse_timestamp(files)
  times = data.frame(fname=files, datetime=times)
  times = times %>% mutate(date=floor_date(datetime, unit = "day"),
                           rel_mtime = as.numeric(difftime(datetime, as.POSIXct(info1$date.tutoring.begin), units="days")))
  times$tutored = times$date >= info1$date.tutoring.begin
  
  if (!is.null(num_per_day)) {
    times = times %>% group_by(date) %>% sample_n(num_per_day, replace=T) %>% distinct(fname)
  }
  ### Parse songs
  #files1 = files[4838:4839]
  files = times$fname
  #songs = profvis({lapply(files1, function(file) {
  #songs = lapply(files, function(file) {
  songs = foreach(file=files) %dopar% {
    print(file)
    song = parse_song(file, log_transform=TRUE,  thresh_method="range", include_wavs=T, thresh_range=seq(0,10, .5), subsamp=10)
    if (is.null(song$syllable)) {
      return(NULL)
    } else {
      return(calc_subsong_features(song))
    }
  }
  #})
  #})})
  npeaks = unlist(lapply(songs, function(x) is.null(x)))
  if (length(npeaks)>0)
    songs = songs[!npeaks]
  songs = bind_rows(songs)
  songs = inner_join(times, songs, by="fname")
  songs = ungroup(songs)
  if (take_db_action) {
    insert_subsong_data(songs, curr_bird, info1$db_action)
  }
  return(songs)
}

run_subsong_shallow_analysis = function(info, curr_bird) {
  info1 = info %>% dplyr::filter(bird==curr_bird)
  dir = paste(servers[servers[,1]==info1$server,"dir"], curr_bird, "subsongs", sep="/")
  files = list.files(dir, pattern = ".wav", full.names = T)
  
  ### Parse fname
  times = parse_timestamp(files)
  times = data.frame(fname=files, datetime=times)
  times = times %>% mutate(date=floor_date(datetime, unit = "day"),
                           rel_mtime = as.numeric(difftime(datetime, as.POSIXct(info1$date.tutoring.begin), units="days")))
  times$tutored = times$date >= info1$date.tutoring.begin
  
  ### Parse songs
  #files1 = files[1:10]
  files = times$fname
  #songs = profvis({lapply(files1, function(file) {
  #songs = lapply(files1, function(file) {
  songs = foreach(file=files) %dopar% {
    print(file)
    song = parse_song(file, fast_filter=T, log_transform=TRUE,  thresh_method="range", include_wavs=F, thresh_range=seq(0,10, .5), subsamp=10)
    #return(song)
    if (is.null(song$syllable)) {
      return(NULL)
    } else {
      return(calc_duration(song))
    }
  }
  #})
  #})})
  npeaks = unlist(lapply(songs, function(x) is.null(x)))
  #songs1 = NULL
  #files1 = NULL
  
  songs1 = songs[!npeaks]
  files1 = files[!npeaks]
  
  names(songs1) = files1
  songs1 = bind_rows(songs1, .id="fname")
  out = data.frame(fname=files1, duration = unlist(songs))
  songs1 = inner_join(times, songs1, by="fname")
  
  insert_subsong_shallow_data(songs1, curr_bird)
  
  return(songs1)
}

# DB ----------------------------------------------------------------------------------------------------------
insert_subsong_data = function(data, bird, action, db_name="rstudio/birds/tutoring/song_analysis/aggregate_data2.db") {
  db_prefix = NULL
  nn = Sys.info()["nodename"]
  if (nn == "lyre") {
    db_prefix = "/media/data2"
  } else if (nn == "osprey") {
    db_prefix = "/mnt/lyre_data2"
  }
  
  data$bird = bird
  data$date_analyzed = Sys.Date()
  table_name = "syllable_data"
  
  db_name = paste(db_prefix, db_name, sep="/")
  db = NULL
  if (!file.exists(db_name)) {
    db = src_sqlite(db_name, create=T)
  } else {
    db = src_sqlite(db_name, create=F)
  }
  
  #db_create_table(db$con, "syllable_data", types=unlist(sapply(data1, class)))
  #db_insert_into(db$con, "syllable_data", data1)
  #replace=F
  data = data %>% mutate(datetime = as.character(datetime),
                           date = as.character(date))
  data = data.frame(data)
  if (db_has_table(db$con, table_name)) {
    if (action=="replace") {
      sql = sprintf("DELETE FROM %s WHERE bird='%s'", table_name, bird)
      dbSendQuery(db$con, sql)
      db_insert_into(db$con, table_name, data, overwrite=F)
    } else if (action=="append") {
      old_data = collect(tbl(db, table_name))
      new_data =  data[!(data$id %in% old_data$id),]
      db_insert_into(db$con, table_name, data[!(data$id %in% old_data$id),], overwrite=F)
    } else if (action=="add") {
      db_insert_into(db$con, table_name, data, overwrite=F)
    }
  } else {
    #if (action=="insert") 
      copy_to(db, data, table_name, temporary=F)
  }
}

insert_subsong_shallow_data = function(data, bird,  db_name="rstudio/birds/tutoring/song_analysis/aggregate_data2.db") {
  data$bird = bird
  data$date_analyzed = Sys.Date()
  table_name = "duration_data"

  db = load_syllable_db(db_name=db_name)
  #db_create_table(db$con, "syllable_data", types=unlist(sapply(data1, class)))
  #db_insert_into(db$con, "syllable_data", data1)
  #replace=F
  #data = data %>% mutate(datetime = as.character(datetime),
  #                       date = as.character(date))
  data = data.frame(data)
  if (db_has_table(db$con, table_name)) {
      sql = sprintf("DELETE FROM %s WHERE bird='%s'", table_name, bird)
      dbSendQuery(db$con, sql)
      db_insert_into(db$con, table_name, data, overwrite=F)
  } else {
    copy_to(db, data, table_name, temporary=F)
  }
}


insert_log_data = function(info) {
  info = info %>% dplyr::filter(!is.na(server))
  dirs = paste(servers[match(info$server, servers[,1]),"dir"], sep="/")
  logs = lapply(1:nrow(info), function(i) {
    load_data(servers[match(info[i,"server"], servers[,1]), "dir"], as.character(info[i,"bird"]))
  })
  logs = bind_rows(logs)
  
  log_table = "logs"
  db = load_syllable_db()
  if (db_has_table(db$con, log_table))
    db_drop_table(db$con, log_table)
  
  copy_to(db, logs, log_table, temporary=F)
  
 # logs = lapply(dirs, parse_log)
}

load_syllable_db =  function(db_name="rstudio/birds/tutoring/song_analysis/aggregate_data2.db") {
  nn = Sys.info()["nodename"]
  if (nn == "lyre") {
    db_prefix = "/media/data2"
  } else if (nn == "osprey") {
    db_prefix = "/mnt/lyre_data2"
  }
  db_name = paste(db_prefix, db_name, sep="/")
  db = src_sqlite(db_name, create=F)
  return(db)
}
load_syllable_data = function(db_name="rstudio/birds/tutoring/song_analysis/aggregate_data2.db") {
  db = load_syllable_db(db_name=db_name)
  table_name = "syllable_data"
  data = collect(tbl(db, table_name))
  data = data %>% mutate(datetime = as.POSIXct(datetime, tz="PST"),
                         date = as.Date(date))
  return(data) 
}

load_log_data = function(db_name="rstudio/birds/tutoring/song_analysis/aggregate_data2.db") {
  db = load_syllable_db(db_name=db_name)
  table_name = "logs"
  data = collect(tbl(db, table_name))
  data = data %>% mutate(time = as.POSIXct(time,  origin="1970-01-01"),
                         date = as.Date(floor_date(time, unit="day")))
  return(data) 
}

load_duration_data = function(db_name="rstudio/birds/tutoring/song_analysis/aggregate_data2.db") {
  db = load_syllable_db(db_name=db_name)
  table_name = "duration_data"
  data = collect(tbl(db, table_name))
  data = data %>% mutate(datetime = as.POSIXct(datetime, origin="1970-01-01", tz="PST"),
                         date = as.Date(floor_date(datetime)),
                         date_analyzed = as.Date(date_analyzed, origin="1970-01-01", tz="PST"))
  return(data) 
}


# SUBSONG ANALYSIS ------------------------------------------------------------------------------------------
calc_subsong_features = function(song) {
  
  ## Parameters
  wl = 1024
  
  ## INITIAL PROCESSING 
  peaks = song$syllable
  peaks = remove_peaks_at_edge(peaks, song$wav@samp.rate, wl=wl, length(song$wav@left) / song$wav@samp.rate)
  peaks$fname = song$fname
  
  
  ## SONG LEVEL 
  song_level_data = data.frame(fname = song$fname)
  
  ### Motifs
  #motifs = find_motifs(peaks, max_gap=1000)
  motifs = song$motifs
  
  if (nrow(motifs)==0) {
    motifs = data.frame(ins=peaks[1,1], outs=peaks[nrow(peaks),2])
  }
  motifs = motifs %>% mutate(duration=outs-ins)
  
  song_level_data$motif_duration_mean = mean(motifs$duration)
  
  motif_peaks = lapply(1:nrow(motifs), function(i) {
    peaks[peaks[,1]>motifs[i,1] & peaks[,2] < motifs[i,2],]
  })
  
  ### Song duration
  song_level_data$song_duration = motifs[nrow(motifs),2] - motifs[1,1]
  
  ### Tempo
  tempo = calc_rolling_tempo(peaks)
  song_level_data$tempo_mean = tempo$tempo_mean
  song_level_data$tempo_sd = tempo$tempo_sd

  ## SYLLABLE LEVEL ------------------------------------------------------------------------------------------
  ### Syllable durations
  peaks = peaks %>% mutate(duration = outs - ins)
  
  ### Amplitude
  amps = apply(peaks, 1, function(peak) {
    calc_amp_stats(song$wav, subsamp=10, subregion=c(as.numeric(peak[1]), as.numeric(peak[2])))
  })
  peaks = cbind(peaks, t(amps))
  
  ### Spectral 
  spec_stats = apply(peaks,1 , function(syl) {
  #spec_stats = foreach(syl=isplitRows(peaks, chunkSize=1)) %do% {
    tmp = calc_spectral_features(song$wav, subregion=c(as.numeric(syl[1]), as.numeric(syl[2])), wl=wl, band=c(1.5, 10))
    tmp["id"] = syl[3]
    # rownames(tmp) = tmp$id[,1]
    tmp
  #}
  })
  spec_stats = bind_rows(spec_stats)
  spec_stats = na.omit(spec_stats)
  peaks = inner_join(peaks, spec_stats)
  
  
  #peaks2 = peaks[,1:2]
  #colnames(peaks2)[1:2] = c("onsets", "offsets")
  if (nrow(peaks) > 1) {
    peaks = calc_ff2(song$wav, peaks, min_freq=1.5, error = .05)
  }
  ## GAP LEVEL --------------------------------------------------------------------------------------------------
  ### Gap durations
  gap_level_data = data.frame(fname=song$fname)
  
  gaps = lapply(motif_peaks, function(p) {
    gaps = peak_dist(p)
    gaps_mean = mean(gaps)
    gaps_sd = sd(gaps)
    gaps_cv = gaps_sd / gaps_mean
    return(c(gaps_mean=gaps_mean,
                gaps_sd=gaps_sd,
                gaps_cv=gaps_cv
                ))
  })
  
  gaps1 = do.call("rbind", gaps)
  gaps1 = matrix(apply(gaps1, 2, mean, na.rm=T), nrow=1)
  colnames(gaps1) = c("gaps_mean", "gaps_sd", "gaps_cv")
  gap_level_data = cbind(gap_level_data, gaps1)
  
  data = inner_join(song_level_data, peaks, by="fname")
  data = inner_join(data, gap_level_data, by="fname")
  return(data)
}

calc_duration_batch = function(files) {
  songs = foreach(file=files) %dopar% {
    print(file)
    song = parse_song(file, fast_filter=T, log_transform=TRUE,  thresh_method="range", include_wavs=F, thresh_range=seq(0,10, .5), subsamp=10)
    #return(song)
    if (is.null(song$syllable)) {
      return(NULL)
    } else {
      return(calc_duration(song))
    }
  }
  names(songs) = files
  return(bind_rows(songs, .id="wav"))
}

calc_duration = function(song) {
  motifs = song$motifs

  if (nrow(motifs)==0) {
    peaks = song$syllable
    motifs = data.frame(ins=peaks[1,1], outs=peaks[nrow(peaks),2])
  }
  song_duration = motifs %>% mutate(duration=outs-ins) %>% summarize(duration=sum(duration))
  return(song_duration)
}

summarize_data = function(d, logs) {
  dm = d %>% gather(measure, value, -fname,-datetime,-date,-rel_mtime,-tutored,-bird,-date_analyzed,-id)
  
  dm1 = dm %>% group_by(date, tutored, bird, measure) %>%
    summarize(measure_mean = mean(value),
              measure_sd = sd(value),
              measure_cv = measure_sd/measure_mean,
              measure_upper = measure_mean + measure_sd,
              measure_lower = measure_mean - measure_sd)
  
  l1 = logs %>% group_by(bird, date, event) %>% 
    summarize(event_num = n())
  
  dm2 = left_join(dm1, l1, by=c("bird", "date"))
  return(dm2)
}
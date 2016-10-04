theme_set(theme_classic())
colors = scale_color_manual("", values=c("black", "red"))
fills = scale_fill_manual("", values=c("black", "red"))


prep_wav_info = function(bird, deafened_date, manipulation, subdir="songs/select", file_ex="wav.not.mat") {
  wdir = paste("/mnt/bengal_home/song", bird, subdir, sep="/")
  info = load_mat_info(wdir, file_ex = file_ex)
  if (nrow(info) == 0) 
    stop(paste("Info not found for: ", bird, sep=""))
  info$bird = bird
  info$manipulation = manipulation
  
  info$date_noon = paste(info$date, "12:00:00", sep=" ")
  info$date_noon = parse_date_time(info$date_noon, orders="ymd hms", tz = "PST")
  if (manipulation=="deaf") {
    info$deafened = info$date>deafened_date
  } else {
    info$deafened = FALSE
  }
  
  info = info %>% mutate(rel_mtime = as.numeric(difftime(mtime, as.POSIXct(deafened_date), units="days")),
                         rel_date = as.numeric(difftime(as.Date(date), as.Date(deafened_date)), units="days"),
                         rel_date_noon = as.numeric(difftime(as.Date(date_noon), as.Date(paste(deafened_date, "12:00:00 PM", sep=" ")), units="days")))
  info$deafened_pretty = factor(ifelse(info$deafened, "deaf", "intact"), levels=c("intact", "deaf"))
  info
}

plot_syl_transition_matrices_by_date = function(info, labels=NULL, rel_date=T) {
 
  absent = NULL
  if (rel_date) {
    d = info %>% group_by(deafened, rel_date) %>% do(process_syllable_matrix(., labels)) %>% ungroup(.)
    date_seq = seq(min(d$rel_date), max(d$rel_date))
    absent = setdiff(date_seq, unique((d$rel_date)))
    colnames(d)[grep("rel_date", colnames(d))] = "date"
  } else {
    d = info %>% group_by(deafened, date) %>% do(process_syllable_matrix(., labels))
    date_seq = seq.Date(min(as.Date(d$date)), max(as.Date(d$date)), by = 1)
    absent = as.Date(setdiff(date_seq, unique(as.Date(d$date))), origin="1970-01-01")
  }
  if (length(absent) > 0 ) {
    mats = expand.grid(date=as.character(absent), From=labels, To=labels, value=0)
    mats = mats %>% filter(From!="end") %>% filter(To!="start")
    mats$deafened = TRUE
    d = rbind(d, mats)
    if (rel_date)
      d$date = as.numeric(d$date)
  }
  
  d$value[is.nan(d$value)] = 0
  gg = ggplot(d) + geom_tile(aes(To, From, fill=value)) + facet_wrap(~date)
  gg = gg + scale_fill_continuous("Total transition probability")
  print(gg)
  return(d)
}
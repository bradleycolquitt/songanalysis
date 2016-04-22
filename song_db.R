library(DBI)
library(zoo)

connectToSongDB = function() {
  return(dbConnect(RMySQL::MySQL(), group="songdb"))
}

readSongData = function(conn, bird_id=NULL) {
  res = NULL
  to_select = c("birds.bird_id", "tags", "manipulation", "time_of_death", 
                "songs.song_id", "song_datetime","song_start", "song_end", "song_number",
                "motifs.motif_id", "motif_number", "motif_start", "motif_end", 
                "syllable_id", "syllable_number", "syllable_start", "syllable_end", "syllable_weinent")
  to_select = paste(to_select, collapse=",")
  if (!is.null(bird_id)) {
    bird_id_form = paste(bird_id, collapse=",")
    sql = paste("SELECT ", to_select, " FROM birds, songs, motifs, syllables WHERE 
          birds.bird_id=songs.bird_id AND
          songs.song_id=motifs.song_id AND
          motifs.motif_id=syllables.motif_id AND
          birds.bird_id in (", bird_id_form, ")", sep="")
  } else {
    sql = paste("SELECT ", to_select, " FROM birds, songs, motifs, syllables WHERE 
          birds.bird_id=songs.bird_id AND
          songs.song_id=motifs.song_id AND
          motifs.motif_id=syllables.motif_id", sep="")
  }
  res = dbSendQuery(conn, sql)
  df = dbFetch(res)
  dbClearResult(res)
  return(df)
}

plot_song_per_hour = function(df) {
  #df = readSongData(conn, bird_id)
  df$song_time1 = strptime(df$song_time, "%H:%M:%S")
  df$song_time_h = df$song_time1$h
  ggplot(df[,-grep("song_time1", colnames(df))] 
         %>% group_by(song_date, song_time_h) 
         %>% summarize(nsongs = n_distinct(song_id)), aes(song_time_h, nsongs)) + geom_point() + facet_wrap(~song_date)
}

bin_song_time = function(df, min_time="09:00:00", max_time="23:59:59", unit=1800) {
  #x = as.xts(as.POSIXct(df$song_datetime))
  df$song_datetime1 = as.numeric(as.POSIXct(df$song_datetime))
  ct = as.POSIXct(df$song_datetime)
  df$song_date = strftime(ct, "%Y-%m-%d")
  calc_time_bounds = df %>% group_by(song_date) %>% do({
    time_breaks = seq(find_time_within_date(.$song_date[1], min_time), 
                      time_max=find_time_within_date(.$song_date[1], max_time),
                      by = unit)
    return(data.frame(., song_time_cut = cut(song_datetime1, breaks=time_breaks)))
    return(data.frame(time_breaks))
  #min_time_ct = as.numeric(as.POSIXct(min_time))
  #ep = endpoints(x, 'minutes', 30)
  #x1 = period.apply(x, ep, function(x) nrow(x))
  #z = aggregate(z, time(z) - as.numeric(time(x)) %% 600, mean)
}

find_time_within_date = function(date, time) {
  as.numeric(as.POSIXct(paste(date, time, sep=" "), tz="PDT"))
  #, format="%m/%d/%Y %H:%M:%S", ))
}
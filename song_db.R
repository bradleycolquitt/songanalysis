library(DBI)

createSongDB(fname) {
  con = dbConnect(RSQLite::SQLite(), fname)
  
  createBirdTable(con)
  createSongTable(con)
  createMotifTable(con)
  createSyllableTable(con)
  createGapsTable(con)
}

#' Columns
#'   bird_id
#'   age
#'   manipulation
#'   time_of_death
writeBirdTable = function(con, df) {
  dbWriteTable(con, "birds", df, append=TRUE)
}

#' Columns
#'   bird_id
#'   song_id
#'   date
#'   time
#'   song_number
#'   song_duration
writeSongTable = function(con, df) {
  dbWriteTable(con, "songs", df, append=TRUE)
}

#' Columns
#'   song_id
#'   motif_number
#'   motif_duration
#'   number_of_syllables
writeMotifTable = function(con, df) {
  dbWriteTable(con, "motifs", df, append=TRUE)
}

#' Columns
#'   motif_id
#'   syllable_id
#'   syllable_number
#'   syllable_duration
#'   syllable_weinent
writeSyllableTable = function(con, df) {
  dbWriteTable(con, "syllables", df, append=TRUE)
}

#' Columns
#'   motif_id
#'   gap_id
#'   gap_number
#'   gap_duration

writeGapTable = function(con, df) {
  dbWriteTable(con, "gaps", df, append=TRUE)
}






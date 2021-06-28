##' Data.
##'
##' @param dat Data
##'
##' @return Data frame.
##'
##' @import raster
##' @export
coarsen <- function(dat, fact = 4){

  ## Okay, now coarsen the data.
  library(raster)
  r <- dat
  coordinates(r) <- ~ lon + lat
  gridded(r) <- TRUE

  # Stack each variable as a band in the raster:
  r <- raster(r, layer = 1)

  # Aggregate to 1 degree resolution (4 times as coarse):
  r <- aggregate(r, fact = fact, fun = mean)

  # Get it back to a dataframe with (lat, lon, val) columns
  coords <- as.matrix(coordinates(r))
  df <- tibble(lon = coords[, 1],
               lat = coords[, 2],
               val = getValues(r))
                   ## as.data.frame(r)) ##%>% rename(val = chl)
 return(df)
}

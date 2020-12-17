##' Gets the data (in \url{https://github.com/brorfred/ocean_clustering}).
##' Resolution of (1/2, 1, 2, 4) degrees.
##'
##' @param resolution
##' @param datadir Directory for data.
##'
##' @return Data that contains month, lat, lon, and some other data columns.
##'
##' @export
read_data <- function(resolution = c("4", "2", "1", "0.5"),
                      type = c("Darwin", "Real"),
                      datadir = "/home/shyun/Dropbox/research/usc/ocean-provinces/omd/data"){

  ## Match resolution and type
  resolution = match.arg(resolution) %>% as.numeric()
  ## resolution = which(resolution == c(0.5, 1, 2, 4))
  resolution = which(resolution == c(4, 2, 1, 0.5))
  type = match.arg(type)
  type = which(type == c("Darwin", "Real"))
  ind = resolution + (type - 1) * 4
  ## print(resolution)
  ## if(resolution == 3) browser()

  ## File names
  filenames = c(
                ## Darwin data
                "tabulated_darwin_montly_clim_045_090_ver_0_2_6.csv",
                "tabulated_darwin_montly_clim_090_180_ver_0_2_6.csv",
                "tabulated_darwin_montly_clim_180_360_ver_0_2_6.csv",
                "tabulated_darwin_montly_clim_360_720_ver_0_2_6.csv",
                ## Real data, new version
                "tabulated_geospatial_montly_clim_045_090_ver_0_2_5.csv",
                "tabulated_geospatial_montly_clim_090_180_ver_0_2_5.csv",
                "tabulated_geospatial_montly_clim_180_360_ver_0_2_5.csv",
                "tabulated_geospatial_montly_clim_360_720_ver_0_2_5.csv",
                ## Real data, old version
                "tabulated_geospatial_montly_clim_045_090_ver_0_2.csv",
                "tabulated_geospatial_montly_clim_090_180_ver_0_2.csv",
                "tabulated_geospatial_montly_clim_180_360_ver_0_2.csv",
                "tabulated_geospatial_montly_clim_360_720_ver_0_2.csv")
  filename = filenames[ind]
  ## print(filename)
  dat = read.csv(file.path(datadir, filename))
  return(dat)
}


##' Month and lon/lat range; also standardizes by dividing by the entire data
##' sum (so that all data sums to one).
##'
##' @param dat Data, obtained from \code{read_data()}.
##' @param month Month.
##' @param lonrange Longitude range.
##' @param latrange Latitude range.
##'
##' @return Processed data
##'
##' @export
process_data <- function(dat, month, lonrange, latrange){

  ## Basic checks.
  stopifnot(month %in% c(1:12))

  ## Month
  mydat = dat[which(dat[,"month"] == month),]
  one.dat = mydat[which(lonrange[1] < mydat[,"lon"] & mydat[, "lon"] < lonrange[2] &
                     latrange[1] < mydat[,"lat"] & mydat[, "lat"] < latrange[2]),]

  ## Formatting the data
  one.mat = make_mat(one.dat)
  if(any(is.na(one.mat))) one.mat = one.mat %>% fill_na()
  one.mat = one.mat / sum(one.mat, na.rm = TRUE)

  ## Return
  stopifnot(abs(sum(one.mat, na.rm = TRUE) - 1) < 1E-8)
  return(one.mat)
}


##' Helper function to take dataset and make matrix that contains Chl, lat, and
##' lon. DONT USE IF data is not on a grid.
##'
##' @param dat A matrix that contains Chl, lat and lon; places them on a grid of
##'   values based on Chl, lat and lon.
##'
##' @return Matrix
##'
##' @export
make_mat <- function(dat){
  ## dat = mydat[which(lonrange[1] < mydat[,"lon"] & mydat[, "lon"] < lonrange[2] &
  ##                    latrange[1] < mydat[,"lat"] & mydat[, "lat"] < latrange[2]),]
  lat.ordered = sort(unique(dat[,"lat"]), decreasing=TRUE)
  lon.ordered = sort(unique(dat[,"lon"]))

  x = sapply(dat[,"lon"], function(lon) which(lon == lon.ordered))
  y = sapply(dat[,"lat"], function(lat) which(lat == lat.ordered))

  xy = cbind(x,y)

  nlon = length(unique(dat[,"lon"]))
  nlat = length(unique(dat[,"lat"]))

  mat = matrix(NA, ncol=nlon, nrow=nlat)
  for(ii in 1:nrow(xy)){
    mat[xy[ii,2], xy[ii,1]] = dat[ii, "Chl"]
  }
  rownames(mat) = lat.ordered
  colnames(mat) = lon.ordered
  return(mat)
}

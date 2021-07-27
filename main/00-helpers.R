##' Given some indices (say, \code{time_range = 1:250}) which means the 1st
##' through the 250th day out of the dates in a separate file e.g. \code{
##' datfile = file.path(datadir, "Justin_pacific_box_halfdegree_Darwin.nc")} },
##' produce a daily dataset processed from the raw daily NC data in
##' e.g. \code{datfile = file.path(datadir,
##' "Justin_pacific_box_halfdegree_Darwin.nc")}.
##'
##' Also removes the data near the coastline, using the \code{overreact=TRUE}
##' option in \code{mark_land()} function.
##'
##' @param type Data type
##' @param datadir Directory of data
##' @param time_range Range of time (vector of two dates). If NULL, attempts to
##'   process /all/ months.
##' @param latrange Latitude range (vector of two numbers)
##' @param lonrange Longitude range (vector of two numbers)
##' @param verbose If TRUE, becomes loud but informative.
##'
##' @return A data frame containing columns: lat, lon, val, time, day, mo, year.
get_nonclim_dat <- function(type = c("real", "darwin"),
                            datadir = "/home/sangwonh/Dropbox/research/usc/ocean-provinces/data",
                            time_range = NULL,
                            latrange,
                            lonrange,
                            verbose = TRUE){


  ## Restrict data to a box
  restrictbox <- . %>% filter(lat >= latrange[1],
                              lat <= latrange[2],
                              lon >= lonrange[1],
                              lon <= lonrange[2])

  ## Setup
  library(ncdf4)
  type = match.arg(type)

  if(type=="real"){
    datfile = file.path(datadir, "Justin_pacific_box_halfdegree.nc")
    timefile = file.path(datadir, "Justin_pacific_box_halfdegree_time.csv")
  }
  if (type == "darwin"){
    datfile = file.path(datadir, "Justin_pacific_box_halfdegree_Darwin.nc")
    timefile = file.path(datadir, "Justin_pacific_box_halfdegree_Darwin_time.csv")
  }


  ## Read the data
  library(tidync)
  dat = tidync(datfile)
  nc_data <- nc_open(datfile)
  all_times = ncvar_get(nc_data, "time") %>% unique() %>% sort()
  if(is.null(time_range)){ times = all_times }
  if(!is.null(time_range)){ times = all_times[time_range] }

  ## Subset data into time frame of interest
  if(verbose){
    print("handling between")
    print(all_times[min(time_range)])
    print("and")
    print(all_times[max(time_range)])
  }
  dat = dat %>% hyper_filter(time = time %in% all_times[time_range])
  dat = dat %>% hyper_tibble()
  if(verbose){
    print("Size of data is")
    object.size(dat) %>% format("Mb") %>% print()
  }
  browser()
  ## stopifnot((times %>% length()) ==
  ##           (read.csv(timefile) %>% .[,2] %>% length()))

  ## Get time strings from external file
  realtimes = read.csv(timefile) %>% .[,2] %>% .[time_range]
  conversion = tibble::tibble(time = all_times[time_range], realtime = realtimes)
  dat = dat %>% dplyr::filter(time %in% (conversion$time))
  ## object.size(dat) %>% format("Mb") %>% print()

  ## Add proper times
  dat = dat %>%
    dplyr::full_join(conversion, "time") %>%
    dplyr::filter(!is.na(realtime)) %>%
    dplyr::filter(!is.na(lat)) %>%
    dplyr::filter(!is.na(lon)) %>%
    dplyr::select(lat, lon, val = chl, time = realtime) %>%
    dplyr::mutate(time = lubridate::as_date(time)) %>%
    dplyr::mutate(day = lubridate::day(time),
                  mo = lubridate::month(time),
                  year = lubridate::year(time)) %>%
    ## remove_landlock() %>%
    dplyr::mutate(land = mark_land(lon, lat, overreact = TRUE)) %>%
    dplyr::filter(land == FALSE)  %>%
    dplyr::select(-land) %>%
    restrictbox()

  ## yy = 1998
  ## mm = 9
  ## this_month_dat = dat %>%
  ##   dplyr::filter(mo == mm) %>%
  ##   dplyr::filter(year == yy) %>%
  ##   group_by(lat, lon) %>%
  ##   summarize(val = mean(val, na.rm=TRUE))  %>%
  ##   ungroup() %>%
  ##   coarsen(4)
  ## this_month_dat %>% plot_dat(add_map = TRUE, hide_legend = FALSE)

  ## dplyr::mutate(val = pmin(val, 1))
  return(dat)
}

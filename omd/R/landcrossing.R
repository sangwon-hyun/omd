##' Takes a long matrix in \code{coords}, and indices of lands \code{land_inds}
##'
##' Calculate a binary matrix of size \code{nrow(coords)} x \code{nrow(coords)}
##' (the same size as a cost matrix) that marks where the problem
##'
##' @export
check_cross_land_outer <- function(coords, land_coords, maxdist){

  ## Basic checks
  stopifnot(all(c("lat", "lon") %in% names(coords)))

  ## Make the crossing matrix
  crossmat = matrix(NA, nrow = nrow(coords), ncol = nrow(coords))
  for(istart in 1:nrow(coords)){
    printprogress(istart, nrow(coords))
    for(iend in 1:nrow(coords)){
      if(istart < iend){
        start = coords[istart,] %>% select(lon, lat)##rename(lon = lon, lat = lat)
        end = coords[iend,] %>% select(lon, lat)##rename(lon = lon, lat = lat)
        if(start$lat < 25 & end$lat < 25 ){
          cross = 0
        } else if (start$lon < -160 & end$lon < -140 ){
          cross = 0
        } else {
          cross = check_cross_land(start$lon, start$lat,
                                   end$lon, end$lat,
                                   land_coords, maxdist)
        }
        crossmat[istart, iend] = crossmat[iend, istart] = cross
      }
    }
  }
  return(crossmat)
}


##' Checks if two (lon, lat) points -- (x0, y0) and (x1, y1) -- crosses any land
##' in \code{land_coords}.
##'
##' @param x0 x0
##' @param y0 y0
##' @param x1 x1
##' @param y1 y1
##' @param land_coords A data frame whose names are \code{lon} and \code{lat}
##' @param maxdist The maximum distance allowed between any land coordinates
##'
##' @export
check_cross_land <- function(x0, y0, x1, y1, land_coords, maxdist){

  ## Basic checks
  assertthat::assert_that(all(c("lon", "lat") %in% colnames(land_coords)))

  ## All distances between the line (x0, y0) and (x1, y1)
  land_coords = land_coords %>% select(lon, lat)
  nn = nrow(land_coords)
  dists = sapply(1:nn, function(ii){
    land_coord =  land_coords[ii,] %>% select(lon, lat) %>% unlist()
    dist2d(land_coord, c(x0, y0), c(x1, y1))
  })
  return(min(dists) < maxdist)
}

##' Distance of point a to line connected between b and c
dist2d <- function(a, b, c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
}

##' Obtain the indices of the set difference between two data frames of lat &
##' lons. If rows 1,3,6 of \code{complete_grid} are missing in
##' \code{incomplete_grid}, we find those indices.
##'
##' @param complete grid A complete data frame of lat & lon columns.
##' @param incomplete grid An incomplete data frame with lat & lon columns. All
##'   rows must be contained in \code{complete_grid}.
##'
##' @return Row indices of \code{complete_grid} of the rows that are missing in
##'   \code{incomplete_grid}.
##'
my_setdiff <- function(incomplete_grid, complete_grid){

  ## Basic checks
  stopifnot(is.data.frame(incomplete_grid))
  stopifnot(is.data.frame(complete_grid))
  stopifnot(sort(names(incomplete_grid)) == c("lat", "lon"))
  stopifnot(sort(names(complete_grid)) == c("lat", "lon"))


  ## Get indices of set difference of the rows
  diff_ind = my_basic_setdiff(a = incomplete_grid,
                              b = complete_grid,
                              check_if_a_in_b = TRUE)

  return(diff_ind)
}

##' More basic version of \code{my_setdiff()}.
my_basic_setdiff <- function(a, b, check_if_a_in_b = TRUE){
  ic = interaction(a)
  cc = interaction(b)
  if(check_if_a_in_b) stopifnot(all(ic %in% cc))
  diff_ind = match(setdiff(cc, ic), cc) %>% na.omit()
  return(diff_ind)
}


##' @import spData
##' @import sf
mark_land <- function(lon, lat, overreact=FALSE){

  ## Coordinates of land
  ## library(sf)
  ## library(spData) ## For `world`, an sf MULTIPOLYGON object
  pts <- st_as_sf(data.frame(lon = lon, lat = lat),
                  coords = 1:2, crs = 4326)
  if(!overreact){
    ## ## Find which points fall over land
    world <- ne_countries(scale = "medium", returnclass = "sf")
    ii <- !is.na(as.numeric(st_intersects(pts, world)))
    return(ii)

  } else {

    world <- ne_countries(scale = "medium", returnclass = "sf")
    ROI = ne_countries(scale = "medium", returnclass = 'sf') %>% st_combine()
    ##KM/earth circumference * degrees in circle
    buffer_in_km <- 250
    buffer_as_arc_degrees<- buffer_in_km/40075*360

    coastalWaters = ROI %>%
      st_buffer(buffer_as_arc_degrees)  %>%
      st_wrap_dateline()

    ii1 <- !is.na(as.numeric(st_intersects(pts, coastalWaters)))
    ii2 <- !is.na(as.numeric(st_intersects(pts, world)))
    ii1[which(lat<25)] = FALSE
    ii1[which(lat<30 & lon< -140)] = FALSE
    ii1[which(lat<40 & lon< -150)] = FALSE
    ii1[which(lat>50 & lon< -160)] = TRUE
    ii = (ii1 | ii2)

  ## tibble(lon, lat, val=ii1) %>% plot_dat(add_map = FALSE, hide_legend = FALSE)
  ## tibble(lon, lat, val=ii2) %>% plot_dat(add_map = FALSE, hide_legend = FALSE)
  ## tibble(lon, lat, val=ii) %>% plot_dat(add_map = FALSE, hide_legend = FALSE)

    return(ii)
  }


  ## plot(st_geometry(world), ylim = latrange, xlim = lonrange)

  ## plot(pts, col=1+ii, pch=16, add=TRUE)
  ## if(!is.null(window)){

  ##   lons = c(lon - window, lon, lon + window)
  ##   lats = c(lat - window, lat, lat + window)

  ##   approx_lonlats = expand.grid(lon=lons, lat=lats)%>%as_tibble()
  ##   inland = sapply(1:nrow(approx_lonlats), function(irow){
  ##     pts <- st_as_sf(data.frame(lon = approx_lonlats[irow,"lon"],
  ##                                lat = approx_lonlats[irow,"lat"]),
  ##                     coords = 1:2, crs = 4326)
  ##     return(!is.na(as.numeric(st_intersects(pts, world))))
  ##   })
  ##   ii = any(inland)
  ## }
}



##' Takes a data matrix with columns: month, lat, lon, val. Then, it makes the
##' grid a complete rectangle; more specifically, it adds new rows with lat/lon
##' values that were missing in \code{dat}, and adds rows with empty (NA) values
##' to the matrix.
##'
##' @param dat Data matrix.
##'
##' @return Completed data
complete_rectangle <- function(dat){

  ## Basic check
  mo = unique(dat$month)
  ty = unique(dat$dat_type)
  stopifnot(length(mo) == 1)
  stopifnot(length(ty) == 1)

  ## Add rows with NA values to make a complete rectanble of a grid.
  complete_grid = expand.grid(lat = dat$lat %>% unique(),
                              lon = dat$lon %>% unique()) %>% as_tibble()
  incomplete_grid = dat %>% select(lat, lon)

  rest_of_grid = dplyr::anti_join(as_tibble(complete_grid), incomplete_grid)
  if(nrow(rest_of_grid) > 0){
    rest_of_dat = data.frame(month = mo, lat = rest_of_grid$lat,
                             lon = rest_of_grid$lon, val = NA,
                             dat_type = ty)
    dat = dat %>% bind_rows(rest_of_dat)
  } else {
    print("Rectangle was complete already! Nothing was changed.")
  }

  return(dat)
}

##' Same as \code{complete_rectangle} but doesn't need \code{month} or
##' \code{dat_type},.
complete_rectangle_barebones <- function(dat){
  complete_grid = expand.grid(lat = dat$lat %>% unique(),
                              lon = dat$lon %>% unique()) %>% as_tibble()
  incomplete_grid = dat %>% dplyr::select(lat, lon)

  rest_of_grid = dplyr::anti_join(as_tibble(complete_grid), incomplete_grid)
  if(nrow(rest_of_grid) > 0){
    rest_of_dat = data.frame(lat = rest_of_grid$lat,
                             lon = rest_of_grid$lon, val = NA)
    dat = dat %>% dplyr::bind_rows(rest_of_dat)
  } else {
    print("Rectangle was complete already! Nothing was changed.")
  }
  return(dat)
}


##' Make the two data sources share the same lat/lon points; discard the
##' non-overlapping ones. Basically erases a-b from a.
##' @export
intersect_latlon <- function(a, b){
  a_minus_b = my_basic_setdiff(b %>% select(lat, lon),
                               a %>% select(lat, lon), FALSE)
  if(length(a_minus_b) > 0){
    a = a[-a_minus_b,]
  }
  return(a)
}



##' From a long data frame with lat, lon, obtain the land indices by (1)
##' completing the rectangle, and (2) marking the overlaps with land using
##' \code{mark_land()}.
##'
##' @param longdat Long data matrix containing \code{lat}, \code{lon}.
##'
##' @return Subset of the rows of \code{dat} that are land.
##' @export
get_land_coords <- function(dat){
  dat %>% complete_rectangle() %>% mutate(land = mark_land(lon, lat)) %>% filter(land == TRUE)
}

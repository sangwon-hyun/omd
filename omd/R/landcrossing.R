##' Takes a long matrix in \code{coords}, and indices of lands \code{land_inds}
##'
##' Calculate a binary matrix of size \code{nrow(coords)} x \code{nrow(coords)}
##' (the same size as a cost matrix) that marks where the problem.
##'
##' @param coords data frame with lon & lat.
##' @param land_coords coordinates with lon & lat.
##'
##' @export
check_cross_land_outer <- function(coords, land_coords, lonrange = NULL, latrange = NULL,
                                   verbose = FALSE,
                                   mc.cores = 1){

  ## Basic checks
  stopifnot(setequal(c("lon", "lat"), names(coords)))


  ## Make the crossing matrix
  crossmat = matrix(NA, nrow = nrow(coords), ncol = nrow(coords))
  coords = coords %>% dplyr::select(lon, lat)
  for(istart in 1:nrow(coords)){
    if(verbose) printprogress(istart, nrow(coords))

    ## do an MCLAPPLY instead here.
    onecol = parallel::mclapply(1:nrow(coords), function(iend){
      ## printprogress(iend, nrow(coords))
      if(istart >= iend){
        return(NA)
      } else {
        ## Get start and end indices
        start = coords[istart,]
        end = coords[iend,]
        assertthat::assert_that(!all(start == end))

        ## Check whether (start --- end) crosses any land
        cross = crosses_land(start, end, latrange, lonrange, land_coords)
        ## crossmat[istart, iend] = crossmat[iend, istart] =
        return(cross)
      }
    }, mc.cores = mc.cores)

    ## for(iend in 1:nrow(coords)){
    ##   if(istart < iend){

    ##     ## Get start and end indices
    ##     start = coords[istart,]
    ##     end = coords[iend,]
    ##     assertthat::assert_that(!all(start == end))

    ##     ## Check whether (start --- end) crosses any land
    ##     cross = crosses_land(start, end, latrange, lonrange, land_coords)
    ##     crossmat[istart, iend] = crossmat[iend, istart] = cross
    ##   }
    ## }
    onecol = unlist(onecol)
    crossmat[istart,] = onecol
  }

  ## Symmetrize this.
  M <- crossmat
  for(i in 1:nrow(M)) {for(j in 1:i) {M[i,j]=M[j,i] }}
  return(M)

  ## return(crossmat)
}


##' This does the same thing as \code{check_cross_land_outer()} without the
##' mclapply, which seems to be wacky.
check_cross_land_outer2 <- function(coords, land_coords, lonrange = NULL, latrange = NULL,
                                   verbose = FALSE){

  ## Basic checks
  stopifnot(setequal(c("lon", "lat"), names(coords)))


  ## Make the crossing matrix
  crossmat = matrix(NA, nrow = nrow(coords), ncol = nrow(coords))
  coords = coords %>% dplyr::select(lon, lat)
  for(istart in 1:nrow(coords)){
    if(verbose) printprogress(istart, nrow(coords))

    ## do an MCLAPPLY instead here.
    ## onecol = mclapply(1:nrow(coords), function(iend){
    for(iend in 1:nrow(coords)){
      if(istart >= iend){
        next
      } else {
        ## Get start and end indices
        start = coords[istart,]
        end = coords[iend,]
        assertthat::assert_that(!all(start == end))

        ## Check whether (start --- end) crosses any land
        cross = crosses_land(start, end, latrange, lonrange, land_coords)
        crossmat[istart, iend] = crossmat[iend, istart] =cross
      }
    }
  }

  ## ## Symmetrize this.
  ## M <- crossmat
  ## for(i in 1:nrow(M)) {for(j in 1:i) {M[i,j]=M[j,i] }}
  ## return(M)

  return(crossmat)
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



##' Takes a data matrix with columns: lat, lon, val. Then, it makes the
##' grid a complete rectangle; more specifically, it adds new rows with lat/lon
##' values that were missing in \code{dat}, and adds rows with empty (NA) values
##' to the matrix.
##'
##' @param dat Data matrix.
##'
##' @return Completed data
complete_rectangle <- function(dat){

  ## ## Basic check
  ## mo = unique(dat$month)
  ## ty = unique(dat$dat_type)
  ## stopifnot(length(mo) == 1)
  ## stopifnot(length(ty) == 1)

  ## Basic check
  stopifnot(setequal(c("lon", "lat", "val"), colnames(dat)))
  lonlat = dat %>% dplyr::select(lon, lat)
  stopifnot(nrow(unique(lonlat)) == nrow(lonlat))

  ## Add rows with NA values to make a complete rectanble of a grid.
  complete_grid = expand.grid(lat = dat$lat %>% unique(),
                              lon = dat$lon %>% unique()) %>% as_tibble()
  incomplete_grid = dat %>% dplyr::select(lat, lon)

  rest_of_grid = dplyr::anti_join(as_tibble(complete_grid), incomplete_grid)
  if(nrow(rest_of_grid) > 0){
    rest_of_dat = data.frame(lat = rest_of_grid$lat,
                             lon = rest_of_grid$lon,
                             val = NA)
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
  a_minus_b = my_basic_setdiff(b %>% dplyr::select(lat, lon),
                               a %>% dplyr::select(lat, lon), FALSE)
  if(length(a_minus_b) > 0){
    a = a[-a_minus_b,]
  }
  return(a)
}



##' From a long data frame with lat, lon, obtain the land indices by
##'
##' (1) completing the rectangle, and
##'
##' (2) marking the overlaps with land using \code{mark_land()}.
##'
##' @param longdat Long data matrix containing \code{lat}, \code{lon}.
##'
##' @return Subset of the rows of \code{dat} that are land.
##' @export
get_land_coords <- function(dat){
  dat %>% complete_rectangle() %>% mutate(land = mark_land(lon, lat)) %>% filter(land == TRUE)
}



##' Angle btw vector x and y
angle <- function(x,y){
  dot.prod <- x%*%y
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}

##' In a triangle with base a=(x0, y0) and b=(x1, y1), is c=|land_coord| in the
##' "middle" i.e. are both angles (c-b-a) and (c-a-b) smaller than 90 degrees?
check_angles <- function(a, b, land_coord){
  x0 = a[1]
  y0 = a[2]
  x1 = b[1]
  y1 = b[2]
  angle1 = angle(c(x0, y0) - land_coord, c(x1, y1) - c(x0, y0)) * 180 / pi
  angle2 = angle(c(x1, y1) - land_coord, c(x0, y0) - c(x1, y1)) * 180 / pi
  ## print(x0)
  ## print(y0)
  ## print(x1)
  ## print(y1)
  ## print(land_coord)
  ## print(angle1)## < 90) & (angle2 < 90)
  ## print(angle2)
  return((angle1 < 90) & (angle2 < 90))
}




##' Get all in between points (that touch the grid).
##'
##' @param x0 Starting lon.
##' @param y0 Starting lat.
##' @param x1 Ending lon.
##' @param y1 Ending lat.
##' @param one_unit Length of one unit, in lon/lat.
##'
##' @return All points in between that touch the grid (two-column matrix).
##'
get_all_in_btw_points <- function(x0, y0, x1, y1, one_unit){
  if(abs(x1 - x0) >= one_unit){
    ## xpoints = seq(from = x0 + one_unit, to = x1 - one_unit, by = one_unit * sign(x1-x0))
    xpoints = seq(from = x0, to = x1, by = one_unit * sign(x1-x0))
    points1 = lapply(xpoints, function(x){
      xinc = x - x0
      slope = (y1-y0) / (x1-x0)
      y_counterpart = xinc * slope + y0
      return(data.frame(lon=x, lat=y_counterpart))
    }) %>% bind_rows()
  } else {
    points1 = NULL
  }

  if(abs(y1 - y0) >= one_unit){
    ## ypoints = seq(from = y0 + one_unit, to = y1 - one_unit, by = one_unit * sign(y1-y0))
    ypoints = seq(from = y0, to = y1, by = one_unit * sign(y1-y0))
    points2 = lapply(ypoints, function(y){
      yinc = y - y0
      slope = (x1-x0)/(y1-y0)
      x_counterpart = yinc * slope + x0
      return(c(lon=as.numeric(x_counterpart), lat=y))
    }) %>% bind_rows()
  } else {
    points2 = NULL
  }
  all_points = rbind(points1, points2)
  return(all_points)
}



##' Helper function
##'
##' @examples
##' \dontrun{
##' ## Sample settings
##' land_coords = rbind(c(lon = -115.75, lat = 30.25),
##'                     c(lon = -115.75, lat = 32.25),
##'                     c(lon = -115.75, lat = 34.25),
##'                     c(lon = -115.75, lat = 36.25),
##'                     c(lon = -115.75, lat = 38.25),
##'                     c(lon = -115.75, lat = 40.25),
##'                     c(lon = -115.75, lat = 42.25),
##'                     c(lon = -117.75, lat = 34.25),
##'                     c(lon = -117.75, lat = 36.25),
##'                     c(lon = -117.75, lat = 38.25),
##'                     c(lon = -117.75, lat = 40.25),
##'                     c(lon = -117.75, lat = 42.25),
##'                     c(lon = -119.75, lat = 36.25),
##'                     c(lon = -119.75, lat = 38.25),
##'                     c(lon = -119.75, lat = 40.25),
##'                     c(lon = -119.75, lat = 42.25))
##'
##' one_unit = 2
##' orig = c(x0 = -115.75,
##'          y0 = 26.25,
##'          x1 = -115.75,
##'          y1 = 32.25)
##' list_of_ports = list(orig + one_unit * c(-1, 3, 1, 0),
##'                      orig + one_unit * c(-1, 3, 1, 3),
##'                      orig + one_unit * c(-2, 3, 1, 4),
##'                      orig + one_unit * c(-2, 3, 1, 10),
##'                      orig + one_unit * c(-2, 3, 1, -4),
##'                      orig + one_unit * c(-2, 3, 1, -3),
##'                      orig + one_unit * c(-2, 3, 1, -2),
##'                      c(x0 = -115.75, y0 = 42.25, x1 = -119.75, y1 = 40.25))
##'
##' par(mfrow = c(2,4))
##' for(ii in 1:8){
##'   ports = list_of_ports[[ii]]
##'   ylim = c(24.25, 42.25)
##'   xlim = c(-119.75, -111.75)
##'   plot(x = c(ports["x0"], ports["x1"]),
##'        y = c(ports["y0"], ports["y1"]),
##'        type = 'o',
##'        xlim = xlim, ylim = ylim,
##'        ylab = "lat",
##'        xlab = "lon", lwd=2)
##'   abline(v=seq(from=xlim[1], to=xlim[2], by=one_unit), col=rgb(0,0,0,0.5))
##'   abline(h=seq(from=ylim[1], to=ylim[2], by=one_unit), col=rgb(0,0,0,0.5))
##'   all_points = get_all_in_btw_points(ports["x0"],
##'                                      ports["y0"],
##'                                      ports["x1"],
##'                                      ports["y1"], one_unit)
##'   points(all_points, pch = 3)
##'   points(land_coords, pch = 16, col='green', cex = 2)
##'
##'   nn = nrow(land_coords)
##'   all_close_to_land = sapply(1:nn, function(ii){
##'     land_coord =  land_coords[ii,] %>% dplyr::select(lon, lat) %>% unlist()
##'     close_to_land = apply(all_points, 1, function(one_point){
##'       check_if_close_to_land(x = one_point["lon"],
##'                              y = one_point["lat"], land_coord, one_unit/2)
##'     })
##'     points(all_points[which(close_to_land),], col='red', pch=3, cex=3)
##'     return(any(close_to_land))
##'   })
##'   legend("bottomright", col = c("green", "red"),
##'          pch = c(16, 3),
##'          cex = c(1, 1),
##'          bg = "white",
##'          legend = c("Land", "Too close to land"))
##' }
##' }
check_if_close_to_land <- function(x, y, land_coord, dst){
  sum((land_coord - c(x,y))^2) <= dst
}



##' Helper function that checks if (start) and (end) are in box define by
##' (latrange) and (lonrange).
both_start_and_end_are_in_box <- function(start, end, latrange, lonrage){
  if(is.null(latrange) | is.null(lonrange)){
    return(TRUE)
  } else {
    return((latrange[1] < start$lat  &  start$lat < latrange[2]) &
           (latrange[1] <   end$lat  &    end$lat < latrange[2]) &
           (lonrange[1] < start$lon  &  start$lon < lonrange[2]) &
           (lonrange[1] <   end$lon  &    end$lon < lonrange[2]))
  }
}



##' Checks if the line between (start ---- end) crosses any land.
crosses_land <- function(start, end, latrange, lonrange, land_coords){
  if(!both_start_and_end_are_in_box(start, end, latrange, lonrange)){
    return(FALSE)
  } else {

    ## Get all in-between points touching the grid
    all_points = get_all_in_btw_points(start$lon, start$lat, end$lon,
                                       end$lat, one_unit)

    ## Check if close to all land coordinates
    cross = FALSE
    nn = nrow(land_coords)
    ## land_coords = land_coords[sample(1:nrow(land_coords), replace = FALSE), ]
    for(ii in 1:nn){
      land_coord =  land_coords[ii,] %>% dplyr::select(lon, lat) %>% unlist()
      close_to_land = apply(all_points, 1, function(one_point){
        check_if_close_to_land(x = one_point["lon"],
                               y = one_point["lat"], land_coord, one_unit/2)
      })
      if(any(close_to_land)){ break }
    }
    crosses_land_at_least_once = (ii < nn)
    return(crosses_land_at_least_once)
  }
}

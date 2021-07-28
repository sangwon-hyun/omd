##' Base R plotting of the optimal transport. Produces three figures: (1) The
##' "before" image, (2) the "after" image, and (3) the optimal transports.
##'
##' @param obj An omd object.
##' @param colfun A function that takes in a numeric vector of numbers between 0
##'   and 1, and produces a vector of colors, to use for plotting the grid.
##' @param main1 Title of first image.
##' @param main2 Title of second image.
##' @param cex Size of the squares to plot in a grid.
##' @param maxcol Maximum value to use for color evaluation
##'
##' @export
plot_omd <- function(obj, lwd_max = 10, colfun = NULL, cex = 3, maxcol = NULL,
                     main1 = NULL, main2 = NULL, ...){


  ## Get all the transfers
  mat = obj$transport_obj
  all_transfers = mat[,"mass"]
  ## cutoff = quantile(all_transfers, 0.1)

  ## Setup
  if(is.null(maxcol)) maxcol = max(obj$M1_long[,"val"], obj$M2_long[,"val"])
  nr = max(obj$M2_long$lat)
  if(is.null(main1))main1 = "From"
  if(is.null(main2))main2 = "To"


  ## Make the three plots.
  par(mfrow = c(1, 3))
  base_r_plot(obj$M1_long[,"lon", drop=TRUE],
              obj$M1_long[,"lat", drop=TRUE],
              max(obj$M1_long[,"lat", drop=TRUE]),
              val = obj$M1_long[,"val", drop=TRUE],
              main = main1,
              colfun = colfun,
              cex_dat = cex_dat,
              maxcol = maxcol,
              ...)

  base_r_plot(obj$M2_long[,"lon", drop=TRUE],
              obj$M2_long[,"lat", drop=TRUE],
              max(obj$M2_long[,"lat", drop=TRUE]),
              val = obj$M2_long[,"val", drop=TRUE],
              main = main2,
              colfun = colfun,
              cex_dat = cex_dat,
              maxcol = maxcol,
              ...)

  bw_colfun <- function(cols){ cols %>% sapply(., function(col)rgb(0, 0, 0, col))}
  base_r_plot(obj$M1_long[,"lon", drop=TRUE],
              obj$M1_long[,"lat", drop=TRUE],
              max(obj$M1_long[,"lat", drop=TRUE]),
              val = obj$M1_long[,"val", drop=TRUE],
              main = "Mass transports",
              colfun = bw_colfun,
              cex_dat = cex_dat,
              ...)

  ## Add the Arrows
  lwd = mat[,"mass", drop=TRUE]
  lwd = lwd/max(lwd)*lwd_max
  for(ii in 1:nrow(mat)){
    one_transfer = mat[ii,, ]

    ## Map the transfers
    coord_from = obj$M1_long[one_transfer[,"from"], c("lon", "lat")] %>% unlist()
    coord_to = obj$M1_long[one_transfer[,"to"], c("lon", "lat")] %>% unlist()

    ## Skip the transfers to the same coordinate
    if(all(coord_from == coord_to)) next

    arrows(x0 = coord_from['lon'], y0 = nr-coord_from['lat'],
           x1 = coord_to['lon'],   y1 = nr-coord_to['lat'],
           col = "red" %>% adjustcolor(alpha = 0.3),
           length = 0.1/2,
           lwd = lwd[ii])
  }
}



##' Base R plot for plotting latitude and longitude.
##'
##' @param lon Longitude.
##' @param lat Latitude.
##' @param maxlat Maximum latitude value.
##' @param val Values to plot.
##' @inheritParams plot_omd
##'
##' @export
base_r_plot <- function(lon, lat, maxlat, val, colfun = NULL, maxcol = NULL,
                        cex_dat = 3,
                        ...){


  ## Setup
  if(is.null(maxcol)){ maxcol = max(val) }

  ## Define color function.
  if(is.null(colfun)){
    ramp <- colorRamp(c("blue", "white", "red"))
    colfun <- function(cols){
      ramp(cols) %>% apply(1, function(r.g.b)rgb(r.g.b[1], r.g.b[2], r.g.b[3], max = 255))
    }
  }

  ## Define color function.
  ## val = obj$M1_long[,"val"]
  cols = val %>% pmax(0)
  cols = (cols/maxcol) %>% sapply(., colfun)

  ## Longitude
  plot(x = lon,
       y = maxlat - lat,
       xlab = "lon",
       ylab = "lat",
       col = cols,
       pch = 15,
       cex = cex_dat,
       ## xlab = "",
       ## ylab = "",
       ...)
}



##' ggplot-based plotting of the optimal transport. Produces three figures: (1)
##' The "before" image, (2) the "after" image, and (3) the optimal transports.
##'
##' @param obj An omd object.
##'
##' @import ggplot2
##' @import sf
##' @import rnaturalearth
plot_omd_ggplot <- function(obj, plot_type = c("one", "four"),
                            name_from = NULL, name_to = NULL,
                            add_map = FALSE,
                            classify_quantile = TRUE,
                            sample_arrows = FALSE){

  ## Setup
  ## stopifnot(class(obj) == "omd")
  plot_type = match.arg(plot_type)
  if(!classify_quantile & (plot_type=="four")) stop("If you don't want to classify quantiles, use plot_type='one'.")

  ## Get the transport table, get rid of same-cell transfers.
  tab = obj$transport_object %>% filter(from!=to)
  if(sample_arrows) tab = tab[seq(from = 1, to = nrow(tab), length=nrow(tab)/3), ]
  nr = nrow(tab)

  ## Transform to lat/lon coordinates on things
  all_rows = lapply(1:nr, function(irow){
    fr = tab[irow,"from"]
    to = tab[irow,"to"]
    mass =  tab[irow, "mass"]
    fr_coord = obj$M1_long %>% .[fr,] %>% dplyr::select(fr_lon = lon, fr_lat = lat)
    to_coord = obj$M1_long %>% .[to,] %>% dplyr::select(to_lon = lon, to_lat = lat)
    cbind(fr_coord, to_coord, mass = mass)
  })
  all_rows = do.call(rbind, all_rows) %>% as_tibble()

  if(classify_quantile){
    ## Categorizing the according to QUARtiles in the mass transfers.
    breaks = all_rows$mass %>% quantile(c(0, 0.25, 0.5, 0.75,1))
    breaks = breaks + c(-1E-20, 0, 0, 0,1E-20)
    all_rows = all_rows %>% mutate(rng = cut(mass, breaks)) %>%
      mutate(rng = factor(rng, labels=c("Smallest mass (Q1)", "Q2", "Q3", "Largest mass (Q4)")))

    ## ## Categorizing the according to top 90$ and the rest in the mass transfers.
    ## p + geom_sf(data = world) +
    ## labs(x = "Longitude", y = "Latitude") +
    ## coord_sf(xlim = lonrange, ylim = latrange, expand = FALSE)

    breaks = all_rows$mass %>% quantile(c(0, 0.95,1))
    all_rows = all_rows %>% mutate(rng2 = cut(mass, breaks)) %>%
    mutate(rng2 = factor(rng2, labels = c("rest", "Largest mass (90%+)")))
  } else {
    all_rows = all_rows %>% add_column(rng2 = "red")
  }

  ## Make the transfer plot
  before_dat = obj$M1_long
  after_dat = obj$M2_long ## Not used


  if(plot_type == "one"){
    p = all_rows %>% ggplot() +
      geom_raster(aes(x = lon, y = lat, fill = val),
                  dat = before_dat) +
      scale_fill_gradientn(colours = c("white", "black"), guide="colorbar") +
      ## facet_wrap(~rng) +
      geom_segment(aes(x = fr_lon,
                       y = fr_lat,
                       xend = to_lon,
                       yend = to_lat,
                       size = mass,
                       col = rng2),
                   arrow = arrow(length = unit(0.05, "inches"))) +
      theme_minimal() +
      theme(legend.position = "none") +
      scale_size_continuous(range = c(1E-3, 2)) +
      ylab("Latitude") +
      xlab("Longitude")  +
      coord_fixed(ratio = 1)

    if(classify_quantile){
      p = p + scale_color_manual(values=c("rest"=rgb(0,0,1,0.01), "Largest mass (90%+)"=rgb(1,0,0, 0.4)))
    } else {
      ## p = p + scale_color_manual("red"= rgb(1,0,0,0.2))
    }
  }

  if(plot_type == "four") {

    ## Four plots
    p = all_rows %>% ggplot() +
      geom_raster(aes(x = lon, y = lat, fill = val), dat = before_dat) +
      scale_fill_gradientn(colours = c("white", "black"), guide="colorbar") +
      facet_wrap(~rng) +
      geom_segment(aes(x = fr_lon, y = fr_lat, xend = to_lon, yend = to_lat, size = mass),
                   arrow = arrow(length = unit(0.1, "inches")),
                   col = rgb(1,0,0,0.2)) +
      theme_minimal() +
      theme(legend.position = "none") +
      scale_size_continuous(range = c(1E-3, 2)) +
      coord_fixed(ratio = 1) +
      ylab("Latitude") +
      xlab("Longitude")

  }

  mytitle = paste0("Optimal transport (", obj$p, "-Wasserstein)")
  if(!is.null(name_from)) mysubtitle = paste0(name_from, " to ", name_to) else mysubtitle = NULL
  p = p + ggtitle(label = mytitle, subtitle = mysubtitle)

  ## Other style elements
  p = p + theme(plot.title = element_text(size = rel(1.3), colour = "black", face = "bold"))
  p = p + theme(axis.text = element_text(size = rel(0.9)),
                axis.title = element_text(size = rel(1.3), face="bold"))
  p = p +  theme(plot.title = element_text(hjust = 0.5),
                 plot.subtitle = element_text(hjust = 0.5))

  p = p + theme(strip.text.x = element_text(size = rel(1.3)))

  ## Add map elements
  if(add_map){
  latrange = obj$M1_long %>% pull(lat) %>% range()
  lonrange = obj$M1_long %>% pull(lon) %>% range()
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  p = p + geom_sf(data = world, inherit.aes = F)+
    coord_sf(xlim = lonrange, ylim = latrange, expand = FALSE,
             label_axes = list(bottom = "E", right = "N"))
    ## theme(axis.title.y = element_text(vjust=-20))
  }

  return(p)
}


##' Make a multi-panel plot of the data.
##'
##' @param longdat A matrix with columns: \code{lat}, \code{lon},
##'   \code{val}. Optionally, \code{longdat} can have a column named
##'   \code{dat_type}. Also otionally, provide \code{mat}, which will be
##'   converted to long format, then plotted.
##' @param mat Defaults to NULL. Will be converted to long format, then plotted
##'   in a single panel
##' @param valname Optionally provide name of the measurement
##'   (e.g. \code{valname = "Chlorophyll"}.
##'
##' @return ggplot object.
##'
##' @export
plot_dat <- function(longdat = NULL, mat = NULL, valname = NULL, standardize = TRUE, hide_legend = TRUE,
                     add_map = FALSE, colours = NULL, lonrange = NULL, latrange = NULL){

  ## Basic checks
  assertthat::assert_that(!is.null(mat) | !is.null(longdat))
  assertthat::assert_that(xor(!is.null(mat), !is.null(longdat)))
  if(!is.null(mat)){
    longdat = image_to_long_format(mat)
  }
  assertthat::assert_that(all(c("lon", "lat", "val") %in% names(longdat)))

  ## If not dat_type exists, create one
  if(!("dat_type" %in% colnames(longdat))) longdat = longdat %>% add_column(dat_type = "")

  ## Optionally, standardize the values within each group
  if(standardize) longdat = longdat %>% group_by(dat_type) %>%
                    mutate(val = val/sum(val, na.rm=TRUE)) %>% ungroup()

  ## Change value name (this is a cosmetic change)
  if(is.null(valname)) valname = "val"
  if(!is.null(valname)) longdat = longdat %>% rename(!!valname := val)


  ## Colors setup
  if(is.null(colours)) colours = c("blue", "red", "yellow")

  ## Make a facet plot
  p = longdat %>%
    group_by(dat_type) %>%
    ggplot() +
    facet_grid(.~dat_type)  +
    geom_raster(aes(x = lon, y = lat, fill = !!sym(valname))) +
    theme_minimal() +
    scale_fill_gradientn(colours = colours, guide="colorbar") +
    ylab("Latitude") + xlab("Longitude") +
    coord_fixed(ratio = 1)

  if(hide_legend) p = p + theme(legend.position = "none")

  ## Other style elements
  p = p + theme(strip.text = element_text(size = rel(1.3), colour = "black", face = "bold"))
  p = p + theme(axis.text = element_text(size = rel(0.9)),
                axis.title = element_text(size = rel(1.3), face="bold"))
  p = p + theme(strip.text.x = element_text(size = rel(1.3)))

  ## Add map elements
  if(add_map){
    if(is.null(latrange)) latrange = longdat %>% pull(lat) %>% range()
    if(is.null(lonrange)) lonrange = longdat %>% pull(lon) %>% range()
    world <- ne_countries(scale = "medium", returnclass = "sf")
    p = p + geom_sf(data = world, inherit.aes = TRUE) +
      coord_sf(xlim = lonrange, ylim = latrange, expand = FALSE,
               label_axes = list(bottom = "E", right = "N")) +
      theme(panel.grid.major = element_line(color = "gray60", linetype = "dashed", size = 0.25),
            panel.background = element_rect(fill = "grey80"))
      ## theme(axis.title.y = element_text(vjust=-20))
  }


  return(p)
}

##' Add blob.
##'
##' @param longdat Data.
##' @param center lon, lat of center of the blob.
##'
##' @return abc
add_blob <- function(longdat, center = c(-150, 15),
                     sigma = NULL,
                     fac = 60,
                     rotate = 0){

  if(is.null(sigma)) sigma = c(50,0,0,10) %>% matrix(ncol=2)
  if(rotate == 0) rotatemat = diag(c(1,1))
  if(rotate !=0) rotatemat = matrix(c(cos(rotate), -sin(rotate),
                                      sin(rotate), cos(rotate)), ncol = 2,
                                    byrow = TRUE)
  sigma = rotatemat %*% sigma %*% t(rotatemat)

  for(ii in 1:nrow(longdat)){
      pt = longdat[ii,c("lon", "lat")]
      distance = sqrt(sum((pt - center)^2))
      to_add = mvtnorm::dmvnorm(pt, mean = center, sigma = sigma) * fac
      longdat[ii,"val"] = longdat[ii,"val"] + to_add
  }
  return(longdat)
}

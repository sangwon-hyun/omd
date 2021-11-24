##' ggplot-based plotting of the optimal transport.
##'
##' @param obj An "omd" class object.
##' @param plot_type One or four plots.
##' @param add_map If \code{TRUE}, add a map of the coastline.
##' @param classify_quantile If \code{TRUE}, classify into blue and red arrows;
##'   otherwise just use blue arrows whose thickness.
##' @param sample_arrows If \code{TRUE}, only plot half of the arrows.
##' @param map_projection If \code{TRUE}, project onto the "globe", so the plot
##'   appears curved.
##' @param map_orientation Orientation to input to \code{coord_map()} from
##'   \code{dplyr}. \url{https://ggplot2.tidyverse.org/reference/coord_map.html}
##' @param return_arrows If \code{TRUE}, return a table containing information
##'   about arrows.
##' @param arrow_range Range of arrow thicknesses.
##'
##' @import ggplot2
##' @import sf
##' @import rnaturalearth
##'
##' @export
plot_omd_ggplot <- function(obj, plot_type = c("one", "four"),
                            name_from = NULL, name_to = NULL,
                            add_map = FALSE,
                            classify_quantile = TRUE,
                            sample_arrows = FALSE,
                            map_projection = TRUE,
                            map_orientation = c(-36.92, 200,0),
                            return_arrows = FALSE,
                            arrow_range = c(1E-3, 2)
                            ){

  ## Setup
  ## stopifnot(class(obj) == "omd")
  plot_type = match.arg(plot_type)
  if(!classify_quantile & (plot_type=="four")) stop("If you don't want to classify quantiles, use plot_type='one'.")

  ## Get the transport table, get rid of same-cell transfers.
  tab = obj$transport_object %>% filter(from!=to)
  if(sample_arrows){
    arrow_inds = seq(from = 1, to = nrow(tab), length=nrow(tab)/2)
    tab = tab[arrow_inds, ]
  }
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
  if(return_arrows) return(all_rows)

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
    all_rows = all_rows %>%
      mutate(rng2 = cut(mass, breaks)) %>%
      mutate(rng2 = factor(rng2, labels = c("rest", "Largest mass (90%+)")))
  } else {
    all_rows = all_rows %>% add_column(rng2 = "all") ##adjustcolor("blue", alpha = .5))
  }

  ## Make the transfer plot
  before_dat = obj$M1_long
  after_dat = obj$M2_long ## Not used


  if(plot_type == "one"){
    p = all_rows %>% ggplot() +
      geom_tile(aes(x = lon, y = lat, fill = val),
                  dat = before_dat) +
      scale_fill_gradientn(colours = c("white", "black"), guide="colorbar") +
      ## facet_wrap(~rng) +
      geom_segment(aes(x = fr_lon,
                       y = fr_lat,
                       xend = to_lon,
                       yend = to_lat,
                       size = mass,
                       col = rng2),
                   ## arrow = arrow(length = unit(0.05, "inches"))) +
                   arrow = arrow(length = unit(0.025, "inches"))) +
      theme_minimal() +
      theme(legend.position = "none") +
      scale_size_continuous(range = arrow_range) +
      ylab("Latitude") +
      xlab("Longitude")  +
      coord_fixed(ratio = 1)

    if(classify_quantile){
      p = p + scale_color_manual(values=c("rest"=rgb(0,0,1,0.1), "Largest mass (90%+)"=rgb(1,0,0, 0.4)))
    } else {
      p = p + scale_color_manual(values = c("all" = rgb(1,0,0,0.2)))
    }
  }

  if(plot_type == "four") {

    ## Four plots
    p = all_rows %>% ggplot() +
      geom_tile(aes(x = lon, y = lat, fill = val), dat = before_dat) +
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
  ## p = p + theme(plot.title = element_text(size = rel(1.3), colour = "black", face = "bold"))
  ## p = p + theme(axis.text = element_text(size = rel(0.9)),
  ##               axis.title = element_text(size = rel(1.3), face="bold"))
  ## p = p +  theme(plot.title = element_text(hjust = 0.5),
  ##                plot.subtitle = element_text(hjust = 0.5))

  ## p = p + theme(strip.text.x = element_text(size = rel(1.3)))


  p = p + theme(strip.text = element_text(size = rel(1.2), colour = "black", face = "bold"))

  p = p + theme(axis.text = element_text(size = rel(0.9)),
                axis.title = element_text(size = rel(1.2)))##, face="bold"))

  p = p + theme(strip.text.x = element_text(size = rel(1.2)))

  ## Add map elements
  if(add_map){
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sp")
    p = p + geom_polygon(data = world,
                         aes(x = long, y = lat, group = group),
                         fill = "white", colour = "grey70")
  }

  if(map_projection){
    latrange = obj$M1_long %>% pull(lat) %>% range()
    lonrange = obj$M1_long %>% pull(lon) %>% range()
    p = p + coord_map("azequalarea",
                      orientation = map_orientation,
                      ylim = latrange, xlim = lonrange)
  }

  return(p)
}

plot_omd = plot_omd_ggplot




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
plot_dat <- function(longdat = NULL, mat = NULL, valname = NULL,
                     standardize = TRUE, hide_legend = TRUE, add_map = FALSE,
                     colours = NULL, lonrange = NULL, latrange = NULL,
                     map_projection = TRUE, map_orientation = c(-36.92, 200,0),
                     limits = NULL, breaks = waiver(),## for color scale
                     legend_name = waiver()){

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
    geom_tile(aes(x = lon, y = lat, fill = !!sym(valname))) +
    scale_fill_gradientn(colours = colours, limits = limits, breaks = breaks, name = legend_name)+
    theme_minimal() +
    ylab("Latitude") + xlab("Longitude") +
    coord_fixed(ratio = 1)

  if(hide_legend) p = p + theme(legend.position = "none")

  ## Other style elements
  p = p + theme(strip.text = element_text(size = rel(1.2), colour = "black"))##, face = "bold"))

  p = p + theme(axis.text = element_text(size = rel(0.9)),
                axis.title = element_text(size = rel(1.2)))##, face="bold"))

  p = p + theme(strip.text.x = element_text(size = rel(1.2)))

  ## Add map elements
  if(add_map){
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sp")
    p = p + geom_polygon(data = world,
                         aes(x = long, y = lat, group = group),
                         fill = "white", colour = "grey70")
  }

  if(map_projection){
    if(is.null(latrange)) latrange = longdat %>% pull(lat) %>% range()
    if(is.null(lonrange)) lonrange = longdat %>% pull(lon) %>% range()
    p = p + coord_map("azequalarea",
                      orientation = map_orientation,##c(-36.92, 200,0),
                      ylim = latrange, xlim = lonrange)
  }


  return(p)
}




##' helper to add patch (blob) to existing plot.
##'
##' @param longdat Data.
##' @param center lon, lat of center of the blob.
##' @param sigma Covariance of Gaussian blob.
##' @param fac How strong should this blob be.
##' @param rotate how much to rotate the blob.
##'
##' @return The resulting long dataset.
##'
##' @export
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

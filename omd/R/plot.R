##' Base R plotting of the optimal transport.
##'
##' @export
plot_omd <- function(obj, lwd_max = 10, cex_dat = 2, colfun = NULL){


  ## Get all the transfers
  mat = obj$transport_obj
  all_transfers = mat[,"mass"]
  ## cutoff = quantile(all_transfers, 0.1)

  ## Make a plot of the "from" matrix
  maxcol = max(obj$M1_long[,"val"], obj$M2_long[,"val"])
  nr = max(obj$M2_long$lat)

  par(mfrow=c(1,3))
  base_r_plot(obj$M1_long[,"lon"],
              obj$M1_long[,"lat"],
              max(obj$M1_long[,"lat"]),
              val = obj$M1_long[,"val"],
              main = "From")

  base_r_plot(obj$M2_long[,"lon"],
              obj$M2_long[,"lat"],
              max(obj$M2_long[,"lat"]),
              val = obj$M2_long[,"val"],
              main = "To")

  bw_colfun <- function(cols){ cols %>% sapply(., function(col)rgb(0, 0, 0, col))}
  base_r_plot(obj$M1_long[,"lon"],
              obj$M1_long[,"lat"],
              max(obj$M1_long[,"lat"]),
              val = obj$M1_long[,"val"],
              main = "Mass transports",
              colfun = bw_colfun)

  ## Add the Arrows
  lwd = mat[,"mass"]
  lwd = lwd/max(lwd)*lwd_max
  for(ii in 1:nrow(mat)){
    one_transfer = mat[ii,]

    ## Skip the transfers to the same coordinate
    if(one_transfer[,"from"] == one_transfer[,"to"]) next

    ## Map the transfers
    coord_from = obj$M1_long[one_transfer[,"from"], c("lon", "lat")]
    coord_to = obj$M1_long[one_transfer[,"to"], c("lon", "lat")]

    if(all(coord_from == coord_to)) browser()
    ## if(x0 == x1 & y0 == y1) browser()
    arrows(x0 = coord_from[,'lon'], y0 = nr-coord_from[,'lat'],
           x1 = coord_to[,'lon'],   y1 = nr-coord_to[,'lat'],
           col = "red" %>% adjustcolor(alpha = 0.5),
           length = 0.1/2,
           lwd = lwd[ii])
  }
}



##' Base R plot.
##' @export
base_r_plot <- function(lon, lat, maxlat, val, colfun = NULL, maxcol = NULL,
                        cex_dat = 3,
                        ...){

  ## Setup
  if(is.null(maxcol)){ maxcol = max(val) }

  ## Define color function.
  if(is.null(colfun)){
    ramp <- colorRamp(c("blue", "red"))
    colfun <- function(cols){
      ramp(cols) %>% apply(1, function(r.g.b)rgb(r.g.b[1], r.g.b[2], r.g.b[3], max = 255))
    }
  }

  ## Define color function.
  ## val = obj$M1_long[,"val"]
  cols = val %>% pmax(0)
  cols = (cols/maxcol) %>% colfun()

  ## Longitude
  plot(x = lon,
       y = maxlat - lat,
       col = cols,
       pch = 15,
       cex = cex_dat,
       xlab = "",
       ylab = "",
       ...)
}

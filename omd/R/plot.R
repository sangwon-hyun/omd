##' Base R plotting of the optimal transport
##'  d = 1
##' @export
plot_omd <- function(obj, lwd_max = 10, cex_dat = 2){

  ## Get all the transfers
  mat = obj$transport_obj
  all_transfers = mat[,"mass"]
  ## cutoff = quantile(all_transfers, 0.1)

  ## Make a plot of the "from" matrix
  maxcol = max(obj$M1_long[,"val"], obj$M2_long[,"val"])
  nr = max(obj$M2_long$lat)

  par(mfrow=c(1,3))
  cols = obj$M1_long[,"val"] %>% pmax(0)
  cols = cols/maxcol
  cols = sapply(cols, function(col) rgb(0,0,0,col))
  plot(x = obj$M1_long[,"lon"],
       y = nr-obj$M1_long[,"lat"],
       col = cols,
       pch = 15, cex=cex_dat,
       xlab = "",
       ylab = "")
  title(main = "From")

  cols = obj$M2_long[,"val"] %>% pmax(0)
  cols = cols/maxcol
  cols = sapply(cols, function(col) rgb(0,0,0,col))
  plot(x = obj$M2_long[,"lon"],
       y = nr-obj$M2_long[,"lat"],
       col = cols,
       pch = 15, cex=cex_dat,
       xlab = "",
       ylab = "")
  title(main = "To")

  cols = obj$M1_long[,"val"] %>% pmax(0)
  cols = cols/maxcol
  cols = sapply(cols, function(col) rgb(0,0,0,col))
  plot(x = obj$M1_long[,"lon"],
       y = nr-obj$M1_long[,"lat"],
       col = cols,
       pch = 15, cex=cex_dat,
       xlab = "",
       ylab = "")
  title(main = "Mass transports")

  ## Add the colors
  lwd = mat[,"mass"]
  lwd = lwd/max(lwd)*lwd_max
  for(ii in 1:nrow(mat)){
    one_transfer = mat[ii,]

    ## Skip the transfers to the same coordinate
    if(one_transfer[,"from"] == one_transfer[,"to"]) next

    ## Map the transfers
    coord_from = obj$M1_long[one_transfer[,"from"], c("lon", "lat")]
    coord_to = obj$M1_long[one_transfer[,"to"], c("lon", "lat")]
    ## if(one_transfer[,"mass"] > cutoff){
    ##   col = 'red'
    ## } else {
    ##   col = 'green'
    ## }
    ## lines(rbind(coord_from, coord_to), lwd=2, col="red" %>% adjustcolor(alpha = 0.5))
    if(all(coord_from == coord_to)) browser()
    ## if(x0 == x1 & y0 == y1) browser()
    arrows(x0 = coord_from[,'lon'], y0 = nr-coord_from[,'lat'],
           x1 = coord_to[,'lon'],   y1 = nr-coord_to[,'lat'],
           col = "red" %>% adjustcolor(alpha = 0.5),
           length = 0.1/2,
           lwd = lwd[ii])
  }
}

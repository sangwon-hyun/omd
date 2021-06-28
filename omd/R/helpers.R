##' Plots a vector field that characterizes the optimal transport.
##'
##' @param jan.dat one dataset.
##' @param feb.dat another dataset
##'
##' @return The transport result.
##' @export
plot_vf <- function(jan.dat, feb.dat, res){


  ## Make point sizes
  jan.cex = jan.dat[,"Chl"]
  jan.cex = jan.cex/sum(jan.cex) * length(jan.cex)
  feb.cex = feb.dat[,"Chl"]
  feb.cex = feb.cex/sum(feb.cex) * length(feb.cex)

  ## Make points
  plot(x=jan.dat[,"lon"], y=jan.dat[,"lat"], cex=jan.cex, col='blue', pch=16)
  points(x=feb.dat[,"lon"], y=feb.dat[,"lat"], cex=feb.cex, col='red', pch=16)

  ## Overlay the smaller January points once again (commented out for now)
  inds = which(feb.cex < jan.cex)
  jan.cex[inds] = 0
  points(x = jan.dat[,"lon"], y = jan.dat[,"lat"],
         cex = jan.cex, col = 'blue', pch = 16)

  ## make the line segments
  from_indices <- res$from##a$from[nonzero]
  to_indices <- res$to##a$to[nonzero]


  ## Get masses
  masses = res$mass
  nomove = which(from_indices == to_indices)
  masses[-nomove] = 0
  masses = masses / max(masses) * 10
  for (ii in 1:length(from_indices)){
    if(ii %in% nomove) next
    shift = rnorm(1,0,0.1)

    yco <- rev(gg[, 1])
    xco <- gg[, 2]

    wh <- which(res$from != res$to)
    for(whi in wh){
    arrows(xco[rs$from[whi]], yco[rs$from[whi]],
           xco[rs$to[whi]], yco[rs$to[whi]], angle = 5,
           length, col = arrcols[whi], lwd = lwd)
    }
  }

  return(res)
}


##' Compute the wasserstein distance using the algorithm by Aurenhammer,
##' Hoffmann and Aronov (1998) for finding optimal transference plans in terms
##' of the squared Euclidean distance in two dimensions.
##'
##' @param p1 The first \code{pgrid} object.
##' @param p2 The second \code{pgrid} object.
##'
##' @return The Wasserstein distance
##'
##' @export
compute_wasserstein <- function(p1, p2){

  res <- transport::transport(p1, p2, p = 2, method = "aha")
  transport::wasserstein(p1, p2, p = 2, tplan = res)

}


##' Make a vector of colors that are all blue, except for the top \code{prob}
##' quantile of mass transfers.
##'
##' @param res An object from \code{transform::transform()}.
##' @param prob Defaults to 0.1.
##'
##' @return A vector of colors, mostly blue, except for the top \code{prob}
##'   quantile.
##'
##' @export
make_colors <- function(res, prob = 0.1){
  cutoffs = quantile(res$mass, probs = c(1-prob, 1))
  ii = 1
  massrange = cutoffs[ii:(ii+1)]
  cols = rep(rgb(0,0,1,0.2), length(res$mass)) ## blue for all else
  cols[which(massrange[1] < res$mass & res$mass < massrange[2])] = rgb(1,0,0,0.7) ## red for the large ones
  cols
}

##' Fill in NA from surrounding values.
##'
##' @param mat A matrix with possibly missing values
##'
##' @return The same object.
##'
##' @export
fill_na <- function(mat){
  if(any(is.na(mat))){
    all_na_inds = which(is.na(mat), arr.ind=TRUE)
    for(irow in 1:nrow(all_na_inds)){
      inds = all_na_inds[irow,]
      surrounding_values = get_surround(mat, inds)
      mat[inds[1], inds[2]] = mean(surrounding_values)
    }
  }
  ## stopifnot(is.na(mat) %>% sum() == 0)
  return(mat)
}



##' Taken from stackoverflow..
##'
##' @export
get_surround <- function(data, index, type="all"){
  row_height = dim(data)[1]
  col_height = dim(data)[2]
  if(type!= "all" && type!="direct"){
    print(" type has to take values of either 'all' or 'direct' ")
  }else{
    if(type=="all"){
      index1     = c(index[1]-1,index[1]-1,index[1]-1,index[1]  ,index[1]+1,index[1]+1,index[1]+1,index[1]  )
      index2     = c(index[2]-1,index[2]  ,index[2]+1,index[2]+1,index[2]+1,index[2]  ,index[2]-1,index[2]-1)
    }
    if(type=="direct"){
      index1     = c(index[1]-1,index[1]+1,index[1]  ,index[1]  )
      index2     = c(index[2]  ,index[2]  ,index[2]-1,index[2]+1)
    }
  }
  adj_ind1   = index1[which(index1>=1 & index1<=row_height)[which(index1>=1 & index1<=row_height) %in% which(index2>=1 & index2<=col_height)]]
  adj_ind2   = index2[which(index1>=1 & index1<=row_height)[which(index1>=1 & index1<=row_height) %in% which(index2>=1 & index2<=col_height)]]
  empty_ind  = c()
  for(i in 1:length(adj_ind1)){
    empty_ind = c(data[adj_ind1[i],adj_ind2[i]], empty_ind)
    }
  return(empty_ind)
}


##' Taken from Justin's convenience functions..
##'
##' @export
drawmat_precise <- function(mat, contour = FALSE, colfun = NULL, num_color_ramp = 100, ...){

## Dummy data
## data <- matrix(runif(100, 0, 5) , 10 , 10)

  if(is.null(colnames(mat))){
    colnames(mat) <- paste(rep("col\n",ncol(mat)),
                            c(1:ncol(mat)) , sep=" ")
    rownames(mat) <- paste(rep("row",nrow(mat)),
                           c(1:nrow(mat)) , sep=" ")
  }

  ## Color function
  if(is.null(colfun)){
    colfun = colorRampPalette(c("blue", "red"))
  }

  # plot it flipping the axis
  lattice::levelplot(t(mat[c(nrow(mat):1) , ]),
                     col.regions = colfun(num_color_ramp),
                     contour = contour,
                     ## xaxt = 'n',
                     las = 2,
                     ...)
}


##' Makes a pgrid object, from a matrix whose row names and column names contain
##' *equally spaced* latitudes and longitudes.
##'
##' @export
make_pgrid <- function(dat){
  lat = rownames(dat) %>% as.numeric() ##+ c(1, rep(0,nrow(dat)-1))
  lon = colnames(dat) %>% as.numeric()
  generator = list(lat, lon)## %>% print()
  p1 = transport::pgrid(dat, generator = generator)
  p1
}


##' @export
matsplitter<-function(M, r, c) {
    rg <- (row(M)-1)%/%r+1
    cg <- (col(M)-1)%/%c+1
    rci <- (rg-1)*max(cg) + cg
    N <- prod(dim(M))/r/c
    cv <- unlist(lapply(1:N, function(x) M[rci==x]))
    dim(cv)<-c(r,c,N)
    cv
}


##' @export
avg <- function(num){
  nn = length(num)
  starts = seq(from=1, to=nn, by=2)
  sapply(starts, function(start){
    num[c(start, start+1)] %>% mean()
  })
}

##' @export
avg_coarsen <- function(d1){
  resol = nrow(d1)
  d2 = matrix(apply(matsplitter(d1, 2, 2), 3, mean), ncol=resol/2)
  rownames(d2) = rownames(d1) %>% as.numeric %>% avg()
  colnames(d2) = colnames(d1) %>% as.numeric %>% avg()
  return(d2)
}


##' @export
enlarge <- function(d){
  d = d[,rep(1:ncol(d), each=2)]
  d = d[rep(1:nrow(d), each=2), ]
  return(d)
}



##' A sliding window average, of window size (\code{size} x \code{size}).
##'
##' @param M Matrix.
##' @param size Size of sliding window.
##'
##' @return Smoothed matrix of the same size
##'
##' @export
smoothmat <- function(M, size){
  delta = (size - 1) / 2
  M_smoothed = ma.matrix(M, av = "mean", delta = delta, edgeNA = FALSE) ##%>% drawmat_precise()
  stopifnot(all(dim(M_smoothed) == dim(M)))
  return(M_smoothed)
  ##   smoothie::kernel2dsmooth(M, kernel.type="boxcar", n=bw)  %>% drawmat_precise()
}


##' Helper function taken directly from the OLIN package. Takes a sliding window
##' average i.e. (j,k)'th pixel's value is taken from a square window extending
##' +-delta out in either direction. When the window of size (delta*2+1) x
##' (delta*2+1) centered at (j,k) tries to extend to outside of the matrix, it's
##' not allowed to. This is equal to a local regression of zero'th order.
##'
##' @param X Data matrix
##' @param av What function to use for summarizing. Defaults to median.
##' @param delta Window size.
##' @param edgeNA if TRUE, then any pixel whose (delta*2+1) x (delta*2+1) window
##'   extends outside of the matrix is assigned zero.
##'
##' @export
ma.matrix <- function (X, av = "median", delta = 2, edgeNA = FALSE)
{
    Xav <- matrix(NA, nrow = dim(X)[[1]], ncol = dim(X)[[2]])
    if (av == "mean") {
        average <- mean
    }
    else {
        average <- median
    }

    delta1 = floor(delta)
    delta2 = ceiling(delta)

    #### SLIDING WINDOW
    for (j in 1:dim(X)[[1]]) {
        for (k in 1:dim(X)[[2]]) {
            a <- (j - delta1)
            c <- (j + delta1)
            b <- (k - delta2)
            d <- (k + delta2)
            if (a < 1) {
                a <- 1
                c <- 2 * delta + 1
            }
            if (b < 1) {
                b <- 1
                d <- 2 * delta + 1
            }
            if (c > dim(X)[[1]]) {
                a <- dim(X)[[1]] - 2 * delta
                c <- dim(X)[[1]]
            }
            if (d > dim(X)[[2]]) {
                b <- dim(X)[[2]] - 2 * delta
                d <- dim(X)[[2]]
            }
            Xav[j, k] <- average(X[a:c, b:d], na.rm = TRUE)
        }
    }

    ### TREATMENT OF EDGES
    if (edgeNA) {
        Xav[1:(delta), ] <- NA
        Xav[, 1:(delta)] <- NA
        Xav[, (dim(Xav)[[2]] - delta + 1):dim(Xav)[[2]]] <- NA
        Xav[(dim(Xav)[[1]] - delta + 1):dim(Xav)[[1]], ] <- NA
    }
    Xav
}


##' Nonoverlapping matrix average.
##'
##' @param window (window x window) sized non-overlapping bins.
##' @param m matrix
##'
##' @return Matrix of same size but with non-overlapping average.
##' @export
##'
##' @examples \dontrun{
##' m = matrix(rnorm(100), nrow=10, ncol=10)
##' stopifnot(all(m %>% bin_matrix(1) == m))
##' }
binmat <- function(m, window){
  nr = nrow(m)
  nc = ncol(m)
  ir = c(seq(from=1, to=nr, by=window), nr+1) %>% unique()
  ic = c(seq(from=1, to=nc, by=window), nc+1) %>%unique()
  new_m <- matrix(NA, nr, nc)
  for(iir in 1:(length(ir)-1)){
    for(iic in 1:(length(ic)-1)){
      jr = ir[iir]
      jc = ic[iic]
      jr_next = ir[iir+1]-1
      jc_next = ic[iic+1]-1
      new_m[jr:jr_next, jc:jc_next] = mean(m[jr:jr_next, jc:jc_next])
    }
  }
  return(new_m)
}



##' Drawing a 3d plot of a 2d matrix, where the Z values are the pixel values.
##'
##' @param M Matrix.
##' @param zfac Factor for z limits for the 3d plots.
##'
##' @export
drawmat_3d <- function(M, zfac = 30){

  col.facet = level.colors(M, at = do.breaks(range(M), 20),
                           col.regions = cm.colors,
                           colors = TRUE)
  col.facet = col.facet %>% adjustcolor(alpha.f = 0.9)
  cloud(M,
        panel.3d.cloud = panel.3dbars,
        col="white",                      # white borders for bars
        xbase = 1, ybase = 1,
        zlim = c(min(M), max(M)) * zfac,                              # No space around the bars
        scales = list(arrows = FALSE, just = "right"),
        ## col.facet = "grey",
        col.facet = col.facet,
        screen = list(z = 65, x = -65),
        xlim = c(-1, n+1),
        ylim = c(-1,n+1),
        ylab = "lat",
        xlab = "lon",
        zlab = NULL)
}


##' Helper to add dimension names to matrix \code{M}.
##' @param M Matrix that contains image.
##' @export
add_dimnames <- function(M){
  colnames(M) = sapply(1:ncol(M), toString)
  rownames(M) = sapply(nrow(M):1, toString)
  return(M)
}



##' Helper only used in \code{add_checkerboard}.
##' @export
cond <- function(input, div, remainders = 0, reverse = FALSE){
  a = (input %% div) %in% remainders
  ifelse(reverse, !a, a)
}


##' Helper to add global pattern.
##' @export
add_global <- function(M, sig, offset=0){
  n = ncol(M)
  M_add = matrix(-sig, nrow=n, ncol=n)
  M_add[seq(n/2 + 1 + offset, n),] = sig
  return(M + M_add)
}

##' Helper to add local pattern.
##' @export
add_local <- function(M, shift=0, sig = 1){

  n = ncol(M)
  M_shifted = M
  ii = n/2
  M_shifted[(ii - n/8):(ii + n/8) + shift, (ii - n/8):(ii + n/8)] = M_shifted[(ii - n/8):(ii + n/8) + shift, (ii - n/8):(ii + n/8)] + sig
  return(M_shifted)
}

##' Takes in a (2^k x 2^k) matrix, and adds a checkerboard pattern whose squares
##' are sized \code{div}, by adding and minusing a value \code{sig} on the
##' alternating checkboard colors.
##'
##' @param sig how much to add and subgract, in applying the checkerboard
##'   pattern.
##' @param div some power of 2 (i.e. 2^0, 2^1, ...), which is the size of the
##'   checkerboard square you want to create.
##' @param M original matrix to apply checkerboard pattern to.
##'
##' @return A matrix of the same size as the input \code{M}, but with a
##'   checkerboard pattern mean added to it.
##'
##' @export
add_checkerboard <- function(M, div, sig = 1){

  ## Setup (very crude, since brain is lazy today)
  div = div * 2
  stopifnot(ncol(M) == nrow(M))
  n = ncol(M)
  stopifnot(div > 1)
  stopifnot(round(log2(div), 1000) == log2(div))
  stopifnot(round(log2(n), 1000) == log2(n))
  remainder = 1: (div / 2)

  ## Make checkberboard
  ## for(reverse in c(TRUE, FALSE)){
  M_add = matrix(-sig, ncol=n, nrow=n)
  reverse=TRUE
  for(ii in 1:n){
    for(jj in 1:n){
      for(reverse in c(TRUE, FALSE)){
        if( ii %>% cond(div, remainder, reverse)  &  jj %>% cond(div, remainder, reverse)){
          M_add[ii, jj] = sig
        }
      }
    }
  }
  stopifnot(mean(M_add) == 0)
  return(M + M_add)
}

##' Helper to make some colors.
##'
##' @export
make_cols <- function(val){
  col = val
  col = col + min(col);
  col = col/max(col)
  col = sapply(col, function(a) rgb(0,0,1,a/2))
  return(col)
}


##' Fill in empty by taking the average of the 4 x 4 window around it.
##'
##' @param mat Matrix.
##' @param wind Window size, defaults to 4.
##'
##' @return Matrix of same size as \code{mat}.
##' @export
fill_in_empty <- function(mat, wind = 4){
  na.ind = which(is.na(mat), arr.ind=TRUE)
  new_mat = mat
  if(length(na.ind) > 0){
    for(ii in 1:nrow(na.ind)){
      ind = na.ind[ii,]
      rows = ind["row"] + (-wind):wind
      cols = ind["col"] + (-wind):wind
      rows = rows %>% pmin(nrow(mat))
      cols = cols %>% pmin(ncol(mat))
      rows = rows %>% pmax(1) %>% unique()
      cols = cols %>% pmin(1) %>% unique()
      ## print('before')
      ## new_mat[ind["row"], ind["col"]] %>% print()
      new_mat[ind["row"], ind["col"]] = mean(mat[rows, cols], na.rm=TRUE)
      ## print('after')
      ## new_mat[ind["row"], ind["col"]] %>% print()
    }
  }
  return(new_mat)
}



##' Takes data that contains "month", "lon" and "lat" as columns, and isolates
##' it to one month \code{mo} and to a latitude/longitutde range \code{lonrange}
##' and \code{latrange}.
##'
##' @param dat data matrix. Must include "lat" "lon" and "month"
##' @param mo Month.
##' @param latrange Latitude range.
##' @param lonrange Longitude range.
##'
##' @export
get_time_space_box <- function(dat, mo, latrange, lonrange, fill_na = FALSE){

  ##latmin, latmax, lonmin, lonmax){
  stopifnot(length(mo)==1)
  stopifnot(mo %in% 1:12)
  ## lonrange = c(lonmin, lonmax)
  ## latrange = c(latmin, latmax)
  dat = dat %>% subset(month %in% mo) %>% filter(lonrange[1] < lon,
                                                 lon < lonrange[2],
                                                 latrange[1] < lat,
                                                 lat < latrange[2]) %>%
    make_mat()
  if(fill_na) dat = dat %>% fill_in_empty(2)
  return(dat)
}






##' Scales \code{x} into between 0 and 1.
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}


##' Scales \code{x} into between \code{min} and \code{max}.
range_custom <- function(x, min, max){
  stopifnot(min <= x & x <= max)
  (x - min) / (max - min)
}


##' From a numeric vector between 0 and 1, make Red-Yellow-Blue colors.
colfun <- function(vec){
  colfun_val <- colorRamp(RColorBrewer::brewer.pal(11,"RdYlBu"))
  colfun_val(vec) %>% apply(1, function(r.g.b)rgb(r.g.b[1], r.g.b[2], r.g.b[3], max = 255))
}


## ##' Checks crossing
## check_crossing <- function(lat1, lat2, lon1, lon2){

##   ## Two coordinates
##   coord1 = c(lon1, lat1)
##   coord2 = c(lon2, lat2)

##   ## Find all coordinates to check.
##   get_all_land_coordinates <- funciton(){
##   }

##   lat, lon

## }


##' Combine two dataframes
combine <- function(d1, d2, name1="from", name2="to"){

  ## Setup
  cols = c(name1, name2)

  ## Merge (to retain only the shared lat-lons) then row-bind.
  twodat = merge(d1, d2, by = c("lat", "lon"),
                 suffixes = c(".from", ".to")) %>% as_tibble() %>%
    dplyr::rename(!!(name1) := val.from,
           !!(name2) := val.to) %>%
    tidyr::pivot_longer(cols = cols, names_to = "dat_type", values_to = "val") %>%
    dplyr::select(lat, lon, dat_type, val) %>%
    mutate(dat_type = factor(dat_type, levels = c(name1, name2)))

  return(twodat)
}

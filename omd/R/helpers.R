

##' Plots a vector field that characterizes the optimal transport.
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
    ## segments(jan.dat[from_indices[i], 'lon'],
    ##          feb.dat[from_indices[i], 'lat'],
    ##          jan.dat[to_indices[i], 'lon'],
    ##          feb.dat[to_indices[i], 'lat'],
    ##          lty = 1,
    ##          lwd = 0.5)

    shift = rnorm(1,0,0.1)
    ## shift = 0

    yco <- rev(gg[, 1])
    xco <- gg[, 2]

    wh <- which(res$from != res$to)
    for(whi in wh){
    arrows(xco[rs$from[whi]], yco[rs$from[whi]],
           xco[rs$to[whi]], yco[rs$to[whi]], angle = 5,
           length, col = arrcols[whi], lwd = lwd)
    }

    ## arrows(x0=jan.dat[from_indices[ii], 'lon'] + shift,
    ##        y0=jan.dat[from_indices[ii], 'lat'] + shift,
    ##        x1=feb.dat[to_indices[ii], 'lon']   + shift,
    ##        y1=feb.dat[to_indices[ii], 'lat']   + shift,
    ##        length = 0.1,
    ##        lwd = 1)
           ## lty = 1,
           ## lwd = masses[ii])
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
get_surround = function(data, index, type="all"){
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
drawmat_precise <- function(mat, contour = FALSE, ...){

## Dummy data
## data <- matrix(runif(100, 0, 5) , 10 , 10)

  if(is.null(colnames(mat))){
    colnames(mat) <- paste(rep("col\n",ncol(mat)),
                            c(1:ncol(mat)) , sep=" ")
    rownames(mat) <- paste(rep("row",nrow(mat)),
                           c(1:nrow(mat)) , sep=" ")
  }

  ## Color function
  colfun = colorRampPalette(c("blue", "red"))

  # plot it flipping the axis
  lattice::levelplot(t(mat[c(nrow(mat):1) , ]),
                     col.regions = colfun(100),
                     contour = contour,
                     ## xaxt = 'n',
                     las = 2,
                     ...)
}


##' Makes a pgrid object, from a matrix whose row names and column names contain
##' *equally spaced* latitudes and longitudes.
make_pgrid <- function(dat){
  lat = rownames(dat) %>% as.numeric() ##+ c(1, rep(0,nrow(dat)-1))
  lon = colnames(dat) %>% as.numeric()
  generator = list(lat, lon)## %>% print()
  p1 = transport::pgrid(dat, generator = generator)
  p1
}



matsplitter<-function(M, r, c) {
    rg <- (row(M)-1)%/%r+1
    cg <- (col(M)-1)%/%c+1
    rci <- (rg-1)*max(cg) + cg
    N <- prod(dim(M))/r/c
    cv <- unlist(lapply(1:N, function(x) M[rci==x]))
    dim(cv)<-c(r,c,N)
    cv
}


avg <- function(num){
  nn = length(num)
  starts = seq(from=1, to=nn, by=2)
  sapply(starts, function(start){
    num[c(start, start+1)] %>% mean()
  })
}


avg_coarsen <- function(d1){
  resol = nrow(d1)
  d2 = matrix(apply(matsplitter(d1, 2, 2), 3, mean), ncol=resol/2)
  rownames(d2) = rownames(d1) %>% as.numeric %>% avg()
  colnames(d2) = colnames(d1) %>% as.numeric %>% avg()
  return(d2)
}

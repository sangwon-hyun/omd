##' Given image |mat| whose number of entries is N, calculate cost matrix (of
##' size (M x M).
##'
##' @param mat Input matrix (image).
##' @param land_inds Indices of |mat| that correspond to land.
##'
##' @return Cost matrix of size (M x M).
get_costm <- function(mat, land_inds = NULL){

  matlong = mat %>% image_to_long_format() %>% as_tibble() %>% filter(!is.na(val))
  coords = matlong %>% select(lat, lon)

  ## Check which entries in cost matrix represent land crossings
  crossmat = matrix(NA, nrow = nrow(coords), ncol = nrow(coords))
  for(istart in 1:nrow(coords)){
    printprogress(istart, nrow(coords))
    for(iend in 1:nrow(coords)){
      if(istart > iend) next
      start = coords[istart,] %>% rename(x = lon, y = lat)
      end = coords[iend,] %>% rename(x = lon, y = lat)
      cross = check_cross_land(start$x, start$y, end$x, end$y, land_inds)
      crossmat[istart, iend] =  crossmat[iend, istart] = cross
    }
  }
  diag(crossmat) = 0

  ## Calculate cost matrix
  costm = matlong %>% form_cost_matrix()
  ## costm %>% drawmat_precise()

  ## Mark land crossings
  costm[crossmat] = Inf  ##Inf##1E10##Inf
  ## costm %>% drawmat_precise()

  ## Check that cost matrix is symmetric
  stopifnot(isSymmetric(costm))

  return(costm)
}


##' Given image |mat| whose number of entries is N, calculate cost matrix (of
##' size (M x M).
##'
##' @param mat Input matrix (image).
##' @param land_inds Indices of |mat| that correspond to land.
##'
##' @return Cost matrix of size (M x M).
get_costm_geodesic <- function(mat ){

  ## ## Cast to long matrix
  ## matlong = mat %>% image_to_long_format() %>% as_tibble() %>% filter(!is.na(val))
  ## coords = matlong %>% select(lat, lon)

  ## ## Obtain geodesic distances
  ## distmat = geodist::geodist(
  ##     x,
  ##     y,
  ##     paired = FALSE,
  ##     sequential = FALSE,
  ##     pad = FALSE,
  ##     measure = "cheap")
  ## ??geodist

  ## ## Temporary
  ## n <- 50
  ## x <- cbind (-10 + 20 * runif (n), -10 + 20 * runif (n))
  ## y <- cbind (-10 + 20 * runif (2 * n), -10 + 20 * runif (2 * n))
  ## colnames (x) <- colnames (y) <- c ("x", "y")
  ## install.packages("geodist")
  ## library(geodist)
  ## d0 <- geodist (x)
  ## d0 %>% drawmat_precise()
  ## ## End of temporary

  ## geodist::geodist

}

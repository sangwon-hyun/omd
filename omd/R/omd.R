##' Converts a matrix with `n` rows and with columns `lat`, `lon`, into a cost
##' matrix of size (n x n).
##'
##' @param Matrix whose rows are data points, and columns are named "lat" and
##'   "lon".
##'
##' @importFrom magrittr "%>%"
##'
##' @return p x p cost matrix.
##'
##' @export
form_cost_matrix <- function(dat){
  ## dat = image_to_long_format(img1)
  stopifnot("lat" %in% colnames(dat))
  stopifnot("lon" %in% colnames(dat))
  return(dat %>% select(lat, lon) %>% dist() %>% as.matrix())
}


##' Calculate Sinkhorn distance.
##'
##' @param costmat (m x n) Cost matrix
##' @param lambda Regularization parameter; large value is slow but gives
##'   solution closer to the exact Wasserstein distance.
##' @param invec Input vector (size m)
##' @param outvec Output vector (size n)
##' @param eps Tolerance for stopping the algorithm, in terms of row sums of P.
##' @param verbose Set to \code{TRUE} if you want loud code. Defaults to \code{TRUE}
##'
##' @importFrom magrittr "%>%"
##'
##' @return List containing three things; `dist`, the sinkhorn distance, `costm`
##'   the original cost matrix, `transports` is an (m x n) matrix of the optimal
##'   (regularized) transports.
##'
##' @export
sinkhorn <- function(costm, lambda, invec, outvec, eps = 1E-3, verbose=FALSE){

  ## Basic checks
  check_sinkhorn_inputs(costm, lambda, invec, outvec)

  ## Run Sinkhorn-type algorithm
  n = nrow(costm)
  P = exp(-lambda * costm)
  P = P / rowSums(P)
  u = rep(0, n)
  while(max(abs(u - rowSums(P)) > eps)){
    u = rowSums(P)

    ## Scale the rows
    P = P %>% sweep(MARGIN = 1, (invec / u), "*")
    stopifnot(all(abs(rowSums(P)- invec <  1E-8)))

    ## Scale the columns
    fac = outvec / colSums(P)
    P = P %>% sweep(MARGIN = 2, fac, "*")
    stopifnot(all(abs(colSums(P) - outvec < 1E-8)))

    if(verbose) abs(u - rowSums(P)) %>% round(3) %>% print()
  }

  dist = P * costm
  dist[is.nan(dist)] = 0
  dist = sum(dist)

  return(list(transports = P,
              costm = costm,
              dist = dist))
}

##' Basic checks for the sinkhorn inputs.
##' @inheritParams sinkhorn
check_sinkhorn_inputs <- function(costm, lambda, invec, outvec){
  assertthat::assert_that(nrow(costm) == length(invec))
  assertthat::assert_that(ncol(costm) == length(outvec))
  assertthat::assert_that(all.equal(sum(invec), 1)==TRUE)
  assertthat::assert_that(all.equal(sum(outvec), 1)==TRUE)
  assertthat::assert_that(lambda >= 0)
}

##' Helper function to transform matrix of pixel values (image) to a matrix of
##' long format, where one row is one data point.
##'
##' @param img matrix of pixel values; the row names and col names are assumed
##'   to contain lat and lon points.
##'
##' @importFrom magrittr "%>%"
##'
##' @return Matrix with three columns, named `lat`, `lon`, and `val`.
##'
##' @export
image_to_long_format <- function(img){
  lat = rep(rownames(img), ncol(img)) %>% as.numeric()
  lon = rep(colnames(img), each=nrow(img)) %>% as.numeric()
  return( data.frame(lat = lat, lon = lon, val = as.numeric(img)))
}

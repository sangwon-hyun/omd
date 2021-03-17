##' Calculate earthmovers distance between two images. They don't have to add up
##' to 1.
##'
##' @param M1 One matrix containing pixel values (m x m).
##' @param M2 Another matrix.
##' @param costm Cost matrix \code{m^2 x m^2}.
##' @param type Defaults to "transport".
##'
##' @return Not written yet.
##'
##' @export
omd <- function(M1, M2, costm = NULL, M1_long = NULL, M2_long = NULL,
                type = c("transport", "sinkhorn"),
                p = 1,
                sinkhorn_lambda = 30,
                sinkhorn_eps = 1E-15,
                sinkhorn_rel_eps = NULL
                ){

  ## Setup
  type = match.arg(type)

  ## Reformat and make into density
  M1 = M1 %>% add_dimnames()
  M2 = M2 %>% add_dimnames()
  M1 = M1 / sum(M1)
  M2 = M2 / sum(M2)

  ## Make sure they are the same size
  assertthat::assert_that(all(dim(M1) == dim(M2)))

  ## Image to long format
  if(is.null(M1_long)) M1_long = image_to_long_format(M1)
  if(is.null(M2_long)) M2_long = image_to_long_format(M2)

  ## Make cost matrix
  if(is.null(costm)) costm = form_cost_matrix(M1_long, pow = p)

  if(type == "sinkhorn"){

    ## Calculate sinkhorn
    sinkhorn_res = sinkhorn(costm = costm,
                             invec = M1_long %>% pull(val),
                             outvec = M2_long %>% pull(val),
                             lambda = sinkhorn_lambda,
                             eps = sinkhorn_eps,
                             rel_eps = sinkhorn_rel_eps)
    dist = (sinkhorn_res$dist)^(1/p)
    transport_res = NULL
  }
  if(type == "transport"){

    transport_res <- transport::transport(M1_long %>% pull(val),
                                M2_long %>% pull(val),
                                costm = costm)
    dist = transport::wasserstein(M1_long %>% pull(val),
                                  M2_long %>% pull(val),
                                  costm = costm,
                                  tplan = transport_res)
    sinkhorn_res = NULL
  }

  obj = list(dist = dist,
             transport_object = transport_res,
             sinkhorn_object = sinkhorn_res,
             costm = costm,
             M1_long = M1_long,
             M2_long = M2_long)
  return(obj)
}


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
form_cost_matrix <- function(dat, pow=2){
  ## dat = image_to_long_format(img1)
  stopifnot("lat" %in% colnames(dat))
  stopifnot("lon" %in% colnames(dat))
  return(dat %>% dplyr::select(lat, lon) %>% dist(p=pow) %>% as.matrix())
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

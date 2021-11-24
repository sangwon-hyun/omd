##' Calculate Wasserstein's distance between two images (either two matrices
##' representing 2d pixelized images, or two long data frames with columns
##' "lon", "lat", and "val"). The data does not need to be normalized. "omd"
##' stands for Ocean Mover's Distance.
##'
##' @param M1 One matrix containing pixel values (m x n). The row names and
##'   column names of this matrix are assumed to contain longitude and latitude
##'   information. Otherwise, lon/lat values of \code{1:m} and \code{1:n} are
##'   sassumed.
##' @param M2 Another matrix; see \code{M1}.
##' @param costm Cost matrix \code{mm x mn}. If not provided, one is calculated
##'   as Euclidean distance.
##' @param type Defaults to "transport".
##' @param M1_long Long matrix, with a column called "val". If provided, this
##'   supercedes \code{M1}.
##' @param M2_long Another long matrix, with a column called "val". If provided,
##'   this supercedes \code{M2}.
##' @param coordnames Defaults to \code{c("lon", "lat")}.
##' @param geodesic If \code{TRUE}, use geodesic distances instead of euclidean
##'   distances.
##' @param type Whether to use Sinkhorn approximation, or the \code{R}
##'   \code{transport} package.
##' @param p Power of Wasserstein's distance.
##' @param sinkhorn_lambda Regularization parameter for Sinkhorn.
##' @param sinkhorn_eps Tolerance for Sinkhorn objective value.
##' @param sinkhorn_rel_eps Relative tolerance for Sinkhorn objective value.
##' @param geodesic Whether to use geodesic distance instead of Euclidean
##'   distance.
##' @return A list object with the results and data, of S3 class "omd".
##'
##' @importFrom assertthat assert_that
##' @import dplyr
##'
##' @export
omd <- function(M1 = NULL,
                M2 = NULL,
                costm = NULL,
                M1_long = NULL,
                M2_long = NULL,
                coordnames = c("lon", "lat"),
                type = c("transport", "sinkhorn"),
                p = 1,
                sinkhorn_lambda = 30,
                sinkhorn_eps = 1E-15,
                sinkhorn_rel_eps = NULL,
                geodesic = FALSE ##
                ){

  ## Setup
  type = match.arg(type)

  ## Make sure they are the same size
  if(!is.null(M1) & !is.null(M2)){
    assertthat::assert_that(all(dim(M1) == dim(M2)))

    ## Image to long format; then, reformat and make into density.
    if(is.null(colnames(M1))) M1 = M1 %>% add_dimnames()
    M1_long = image_to_long_format(M1)
    coordnames = c("lon", "lat")
    if(is.null(colnames(M2)))M2 = M2 %>% add_dimnames()
    M2_long = image_to_long_format(M2)
    coordnames = c("lon", "lat")
  }

  ## Normalize
  M1_long_orig = M1_long
  M2_long_orig = M2_long
  ## M1_long = M1_long %>% as_tibble() %>% filter(!is.na(val))
  ## M2_long = M2_long %>% as_tibble() %>% filter(!is.na(val))

  ## Get rid of any NAs in the both datasets.
  M1_long_na = M1_long %>% as_tibble() %>% filter(is.na(val)) ##%>% dplyr::select(lon, lat)
  M2_long_na = M2_long %>% as_tibble() %>% filter(is.na(val)) ##%>% dplyr::select(lon, lat)
  if(nrow(M1_long_na) > 0 | nrow(M2_long_na) > 0){
    na_coords = bind_rows(M1_long_na, M2_long_na) %>% dplyr::select(!!!syms(coordnames)) %>% unique()
    M1_long = anti_join(M1_long, na_coords, by = coordnames)
    M2_long = anti_join(M2_long, na_coords, by = coordnames)

    ## ## If cost matrix is provided, take subsets of rows & columns
    M1_remove_row = left_join(M1_long_orig, add_column(na_coords, na = TRUE)) %>% pull(na)
    M2_remove_row = left_join(M2_long_orig, add_column(na_coords, na = TRUE)) %>% pull(na)
    costm = costm[-which(M1_remove_row), -which(M2_remove_row)]
  }

  ## Check if the coordinates match up exactly (this isn't actually necessary..)
  M1_long_coords = M1_long %>% dplyr::select(!!!syms(coordnames))
  M2_long_coords = M2_long %>% dplyr::select(!!!syms(coordnames))
  stopifnot(all( M1_long_coords == M2_long_coords ))

  ## Scale to sum to 1
  M1_long$val = M1_long$val / sum(M1_long$val)
  M2_long$val = M2_long$val / sum(M2_long$val)

  ## Make cost matrix
  if(is.null(costm)){
    if(geodesic){
      costm = form_geodesic_cost_matrix(M1_long)
    } else {
      costm = form_cost_matrix(M1_long)
    }
  }


  ## ## Temporary
  ## if(nrow(M1_long_na) > 0 | nrow(M2_long_na) > 0){
  ##   na_coords = bind_rows(M1_long_na, M2_long_na) %>% dplyr::select(!!!syms(coordnames))
  ##   ## M1_long = anti_join(M1_long, na_coords, by = coordnames)
  ##   ## M2_long = anti_join(M2_long, na_coords, by = coordnames)

  ##   ## ## If cost matrix is provided, take subsets of rows & columns
  ##   browser()
  ##   M1_remove_row = full_join(M1_long_orig, na_coords %>% add_column(na = TRUE)) %>% pull(na)
  ##   M2_remove_row = full_join(M2_long_orig, na_coords %>% add_column(na = TRUE)) %>% pull(na)
  ##   costm = costm[-which(M1_remove_row), -which(M2_remove_row)]
  ## }
  ## ## End of Temporary


  if(type == "sinkhorn"){

    ## Calculate sinkhorn
    sinkhorn_res = sinkhorn(costm = costm^p,
                            invec = M1_long %>% pull(val),
                            outvec = M2_long %>% pull(val),
                            lambda = sinkhorn_lambda,
                            eps = sinkhorn_eps,
                            rel_eps = sinkhorn_rel_eps)
    dist = (sinkhorn_res$dist)##^(1/p)
    transport_res = NULL
  }
  if(type == "transport"){
    transport_res <- transport::transport(M1_long %>% pull(val),
                                          M2_long %>% pull(val),
                                          costm = costm^p)
    dist = transport::wasserstein(M1_long %>% pull(val),
                                  M2_long %>% pull(val),
                                  costm = costm^p,
                                  tplan = transport_res)
    sinkhorn_res = NULL
  }

  obj = list(dist = dist^(1/p),
             transport_object = transport_res,
             sinkhorn_object = sinkhorn_res,
             costm = costm,
             M1_long = M1_long,
             M2_long = M2_long,
             M1_long_orig = M1_long_orig,
             M2_long_orig = M2_long_orig,
             p = p)
  class(obj) = "omd"
  return(obj)
}


##' Converts a matrix with `n` rows and with columns `lat`, `lon`, into a cost
##' matrix of size (n x n).
##'
##' @param dat Matrix whose rows are data points, and columns are named "lat"
##'   and "lon".
##'
##' @importFrom magrittr "%>%"
##'
##' @return p x p cost matrix.
##'
##' @export
form_cost_matrix <- function(dat, geodesic = FALSE){
  ## dat = image_to_long_format(img1)
  stopifnot("lat" %in% colnames(dat))
  stopifnot("lon" %in% colnames(dat))
  if(!geodesic){
    return(form_geodesic_cost_matrix(dat))
  } else {
    ## return(dat %>% dplyr::select(lat, lon) %>% dist() %>% as.matrix())
    return(dat %>% dplyr::select(lat, lon) %>% dist() %>% as.matrix())
  }
}

##' Converts a matrix with `n` rows and with columns `lat`, `lon`, into a cost
##' matrix of size (n x n).
##'
##' @param data Matrix whose rows are data points, and columns are named "lat"
##'   and "lon".
##'
##' @importFrom magrittr "%>%"
##'
##' @return p x p cost matrix in geodesic distance (unit: kilometers).
##'
##' @export
form_geodesic_cost_matrix <- function(dat){
 ## dat = image_to_long_format(img1)
  stopifnot("lat" %in% colnames(dat))
  stopifnot("lon" %in% colnames(dat))
  stopifnot(nrow(dat)>0)
  distmat = (dat %>% dplyr::select(lat, lon) %>% geodist::geodist(measure="geodesic") %>% as.matrix())
  distmat = distmat / 1000 ## converting to KM
  return(distmat)
}


##' Helper function to transform matrix of pixel values (image) to a matrix of
##' long format, where one row is one data point.
##'
##' @param img matrix of pixel values; the row names and col names are assumed
##'   to contain latitude and longitude points.
##'
##' @importFrom magrittr "%>%"
##'
##' @return Matrix with three columns, named `lat`, `lon`, and `val`.
##'
##' @export
image_to_long_format <- function(img){
  assert_that(!is.null(rownames(img)))
  assert_that(!is.null(colnames(img)))
  lat = rep(rownames(img), ncol(img)) %>% as.numeric()
  lon = rep(colnames(img), each=nrow(img)) %>% as.numeric()
  dat = data.frame(lat = lat, lon = lon, val = as.numeric(img))
  dat = dat %>% dplyr::arrange(lon, lat)

  ## ## This is another way to do it
  ## lats = as.numeric(rownames(mat))
  ## mat %>%  as_tibble() %>%
  ##   add_column(lat = as.character(lats)) %>%
  ##   pivot_longer(names_to = "lon", values_to = "val", -lat)

  return(dat)
}


##' The reverse of \code{image_to_long_format()}.
##'
##' @param dat A data frame (or matrix) whose columns are named: \code{lon},
##'   \code{lat}, \code{val}.
##'
##' @return An image, with values in pixels where values exist, and NA's
##'   otherwise.
##'
long_format_to_image <- function(dat){

  ## Size of gridded image
  nlat = length(unique(dat$lat))
  nlon = length(unique(dat$lon))
  lons_orig = sort(unique(dat$lon))
  lats_orig = sort(unique(dat$lat), decreasing = TRUE)
  counts = dat %>% group_by(lat, lon) %>% count() %>% ungroup() %>% pull(n)
  stopifnot(all(counts == 1))

  ## Ensure the column names are lat/lon/val
  stopifnot(all(colnames(dat) %in% c("lon", "lat", "val")))
  dat = dat %>% dplyr::select(lon, lat, val)

  ## Convert first two columns as lon-lat and third as value
  dfr <- raster::rasterFromXYZ(dat)
  lons_new = raster::xFromCol(dfr)
  lats_new = raster::yFromRow(dfr)
  img = raster::as.matrix(dfr)
  rownames(img) = lats_new
  colnames(img) = lons_new

  ## For sanity checks, compare the three plots:
  ## drawmat_precise(mat)
  ## raster::plot(dfr)
  ## plot_dat(dat)

  ## Checking the reverse operation
  dat2 = img %>%
    image_to_long_format() %>%
    as_tibble() %>%
    dplyr::select(lon, lat, val) %>%
    dplyr::arrange(lon, lat) %>%
    drop_na()

  stopifnot(identical((dat %>% drop_na() %>% dplyr::arrange(lon, lat)), (dat2 %>% dplyr::arrange(lon, lat))))

  ## Final checks
  stopifnot((nrow(img) == nlat) & (ncol(img) == nlon))
  stopifnot(all(colnames(img) == lons_orig))
  stopifnot(all(as.numeric(rownames(img)) == lats_orig))

  return(img)
}


##' Print function for omd object.
##'
##' @param obj .
print.omd <- function(obj){
  cat("Your OMD object came from two datasets of size",
      nrow(obj$M1_long), "x",
      ncol(obj$M1_long), fill=TRUE)
  cat("The", obj$p, "- Wasserstein's distance is ", obj$dist, fill=TRUE)}

is.omd <- function(x) inherits(x, "omd")

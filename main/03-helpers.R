##' @param shift How much to shift one of the means (top-right corner) to the
##'   left.
##' @param wiggle How much to wiggle the three means?
##' @param noisy If noisy, draw one from {1,2,3} from the relative probabilities
##'   of the three Gaussians with three spread-out centers. If not noisy, then
##'   just assign {1,2,3} to the cluster with the largest relative probability
##'   (density).
##'
##' @
get_mat <- function(shift=NULL, wiggle=NULL, random_means=FALSE, noisy = FALSE){
  y = matrix(NA, ncol = 30, nrow=30)
  mn1 = c(1,1)
  mn2 = c(1,30)
  mn3 = c(30,15)
  if(random_means){
    mn1 = c(sample(1:30,1), sample(1:30,1))
    mn2 = c(sample(1:30,1), sample(1:30,1))
    mn3 = c(sample(1:30,1), sample(1:30,1))
  }
  Sigma = c(10, 0, 0, 10) %>% matrix(ncol = 2)


  ## If necessary, shift or wiggle the means
  if(!is.null(shift)){
    mn2 = mn2 + c(0, -shift)
  }
  if(!is.null(wiggle)){
    mn1 = mn1 + rmvnorm(1, c(0,0), diag(c(wiggle^2, wiggle^2)))
    mn2 = mn2 + rmvnorm(1, c(0,0), diag(c(wiggle^2, wiggle^2)))
    mn3 = mn3 + rmvnorm(1, c(0,0), diag(c(wiggle^2, wiggle^2)))
  }

  for(ii in 1:30){
    for(jj in 1:30){
      x = c(ii,jj)
      prob = c(dmvnorm(x, mn1, Sigma),
               dmvnorm(x, mn2, Sigma),
               dmvnorm(x, mn3, Sigma))
      prob = prob/sum(prob)
      if(noisy){
        y[ii, jj] = sample(1:3, 1, prob = prob) + rnorm(1,0, .05)
      } else {
        y[ii, jj] = which.max(prob)
      }
    }
  }
  drawmat_precise(y)
  mat = get_cluster_boundary_mat(y)
  return(mat)
}


##' In a matrix y of size 30 x 30, do a K=3 Kmeans clustering (starting means
##' being 1,2, and 3). Then, do horizontal and vertical scan of the matrix to
##' find the breakpoints of the cluster labels. This is a crude way to get 1's
##' marked on the cluster boundaries.
##'
##' @param y 30 x 30 matrix, presumed to have three mean pixel levels of 1 and 2
##'   and 3.
##'
##' @return Matrix of the same size, with 1's on the cluster boundaries, and 0
##'   everywhere else.
get_cluster_boundary_mat <- function(y){

  stopifnot(nrow(y) == 30)
  stopifnot(ncol(y) == 30)

  cluster_res = y %>% as.numeric() %>% kmeans(3, centers=1:3)  ##%>% pull(means)
  three_mns = cluster_res$centers
  labels = cluster_res$cluster ##%>% matrix(ncol = 30)

  labelmat = labels %>% matrix(ncol = 30) ##%>% drawmat_precise()
  labelmat %>% apply(1, function(myrow){
    a = rep(0, length(myrow))
    ind = which(diff(myrow) !=0)
    a[ind] = 1
    a
  }) %>% t() %>% drawmat_precise()

  myfun = function(myrow){
    a = rep(0, length(myrow))
    ind = which(diff(myrow) !=0)
    a[ind] = 1
    a
  }

  mat1 = labelmat %>% apply(1, myfun) %>% t()
  mat2 = labelmat %>% apply(2, myfun)
  mat = (mat1 | mat2) %>% as.matrix()
  return(mat)
}



get_cluster_boundary_mat_barebones <- function(labelmat){

  ## cluster_res = y %>% as.numeric() %>% kmeans(3, centers=1:3)  ##%>% pull(means)
  ## three_mns = cluster_res$centers
  ## labels = cluster_res$cluster ##%>% matrix(ncol = 30)

  ## labelmat = labels %>% matrix(ncol = 30) ##%>% drawmat_precise()
  labelmat %>% apply(1, function(myrow){
    a = rep(0, length(myrow))
    ind = which(diff(myrow) !=0)
    a[ind] = 1
    a
  }) %>% t() ##%>% drawmat_precise()

  myfun = function(myrow){
    a = rep(0, length(myrow))
    ind = which(diff(myrow) !=0)
    a[ind] = 1
    a
  }

  mat1 = labelmat %>% apply(1, myfun) %>% t()
  mat2 = labelmat %>% apply(2, myfun)
  mat = (mat1 | mat2) %>% as.matrix()
  return(mat)
}



## ##' Get a (long-format) data frame with cluster boundaries marked as 1; 0
## ##' otherwise.
## ##'
## ##' @param dat A matrix with columns: \code{lat}, \code{lon},
## ##'   \code{val}.
## get_cluster_boundary_mat_new <- function(dat){

##   ## Make sure data is in right format
##   assertthat::assert_that(all(c("lon", "lat", "val") %in% names(longdat)))



##   ## ## Visualize
##   ## dat %>% ggplot() + geom_line(aes(x=lat, y=val)) + facet_wrap(~lon)
##   ## realdat %>% ggplot() + geom_line(aes(x=lat, y=val)) + facet_wrap(~lon)
##   pipeline = . %>% group_by(lon) %>%
##     mutate(clustered_val = as.factor(kmeans(val, centers=c(0, 0.3))$cluster) ) %>%
##     group_by(lon, clustered_val) %>%
##     mutate(meanval = mean(val)) %>%
##     ungroup() %>%
##     group_by(lon) %>%
##     mutate(changepoint = go_down_lat_and_mark_changepoint(clustered_val)) %>%
##     arrange(lon)

## }
## ```

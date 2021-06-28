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
sinkhorn <- function(costm, lambda, invec, outvec,
                     eps = 1E-3, rel_eps = 1E-8,
                     niter = 1E5){
  ## Basic checks
  check_sinkhorn_inputs(costm, lambda, invec, outvec)

  ## Run Sinkhorn-type algorithm
  n = nrow(costm)
  P = exp(-lambda * costm)
  P = P / rowSums(P)
  ## P might start as a diag(rep(1,m)) matrix.
  assertthat::assert_that(!any(is.nan(P)))
  u = rep(0, n)

  gap = 1E10 ## Just some large number
  gaps = dists = rep(NA, niter)
  iter = 1
  for(iter in 1:niter){
    ## cat("iter", iter, fill=TRUE)
    u = rowSums(P)

    ## Scale the rows
    P = P %>% sweep(MARGIN = 1, (invec / u), "*")
    stopifnot(all(abs(rowSums(P)- invec <  1E-8)))

    ## Scale the columns
    fac = outvec / colSums(P)
    P = P %>% sweep(MARGIN = 2, fac, "*")
    stopifnot(all(abs(colSums(P) - outvec < 1E-8)))

    dist = P * costm
    dist[is.nan(dist)] = 0
    dist = sum(dist)

    gap = max(abs(u - rowSums(P)))

    ## Temporary: Recording the gap and distance
    gaps[iter] = gap
    dists[iter] = dist
    if(iter > 1) dist_ratio = (dists[iter] - dists[iter-1]) / dists[iter-1]

    ## Stopping rule
    if(is.null(rel_eps)){
      if(gap <= eps) break
    } else {
      if(iter > 1){
        if(dist_ratio < rel_eps) break
      }
    }

    ## gap = (max(abs(u - rowSums(P))))
    ## cat("gap", signif(gap, 10), fill=TRUE)

    ## (P * costm) %>% round(2) %>% print()
  }

  return(list(transports = P,
              costm = costm,
              dist = dist,
              ## Things to track for convergence;
              gaps = na.omit(gaps),
              dists = na.omit(dists),
              max_iter = iter
              ))
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

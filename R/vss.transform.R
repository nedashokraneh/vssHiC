#' vss transformation
#'
#' This function takes measurements from two replicates and a trained
#' mean-variance relationship and outputs a two-column data.frame
#' object including variance-stabilized counts.
#'
#' @param cts A two-column data.frame including counts from two replicates.
#'
#' @param vss_fit A smooth.spline object.
#'
#' @param log_scale Whether spline is learned using log-scaled
#' means and variances.
#'
#' @param log_converge If TRUE, variance-stabilizing transformation
#' converges to the log function asymptotically.
#'
#' @return A variance-stabilized InteractionSet object.
#'
#' @import dplyr
#'
#' @export

vss.transform <- function(cts,
                          vss_fit,
                          log_scale = TRUE,
                          log_converge = FALSE){

  colnames(cts) = c('rep1', 'rep2')
  max_ct = max(cts)
  xg = sinh( seq(asinh(0), asinh(max_ct),
                 length.out = 1000) )[-1]
  if (log_scale){
    xg_stdevs = exp(stats:::predict.smooth.spline(vss_fit,
                                                  log(xg + 1))$y) - 1
  } else{
    xg_stdevs = stats:::predict.smooth.spline(vss_fit, xg)$y
  }
  xg_inv_stdevs = 1 / xg_stdevs
  splf <- splinefun(
    asinh( (xg[-1] + xg[-length(xg)]) / 2 ),
    cumsum(
      (xg[-1] - xg[-length(xg)]) *
        (xg_inv_stdevs[-1] + xg_inv_stdevs[-length(xg_inv_stdevs)]) / 2
    )
  )
  if (! log_converge){
    vss_cts = data.frame(rep1 = splf( asinh( cts$rep1 ) ),
                          rep2 = splf (asinh( cts$rep2)))
  } else{
    h1 <- quantile( pmax(cts$rep1, cts$rep2), .95 )
    h2 <- quantile( pmax(cts$rep1, cts$rep2), .99999 )
    log_asymp_params <- get_log_asymp_spline_params(splf, h1, h2)
    eta = log_asymp_params$eta
    xi = log_asymp_params$xi
    vss_cts = data.frame(rep1 = eta * splf( asinh( cts$rep1 ) ) + xi,
                          rep2 = eta * splf( asinh( cts$rep2 ) ) + xi)
  }

  return (vss_cts)
}


get_log_asymp_spline_params <- function(splineFunc,
                                        h1,
                                        h2){
  eta <- ( log2(h2) - log2(h1) ) / ( splineFunc(asinh(h2)) - splineFunc(asinh(h1)) )
  xi <- log2(h1) - eta * splineFunc(asinh(h1))
  #print(eta)
  #print(xi)
  return (list(eta = eta, xi = xi))
}


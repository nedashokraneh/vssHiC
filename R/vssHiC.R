#' vss transformation for Hi-C data
#'
#' This function uses Hi-C contact counts from two replicates
#' (provided as an InteractionSet object) and a trained
#' mean-variance relationship and outputs an InteractionSet
#' object including variance-stabilized counts.
#'
#' @param IntSet An InteractionSet object including two replicates of
#' Hi-C contact counts.
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

vssHiC.transform <- function(IntSet,
                             vss_fit,
                             log_scale = TRUE,
                             log_converge = FALSE){

  cts = data.frame(assay(IntSet))
  vss_cts = vss.transform(cts, vss_fit, log_scale, log_converge)
  vss_cts = as(vss_cts, 'matrix')
  transformed_IntSet = InteractionSet(assays=list(counts=vss_cts),
                                      interactions(IntSet))
  return (transformed_IntSet)
}




#' Learning the empirical mean-variance relationship from
#' two Hi-C replicates.
#'
#' This function uses Hi-C contact counts from two replicates
#' (provided as an InteractionSet object) to calculate
#' observed variances among different count intensities by VSS
#' algorithm (https://academic.oup.com/bioinformatics/article/37/23/4383/6308936).
#'
#'
#' @param IntSet An InteractionSet object including two replicates of Hi-C contact counts.
#'
#' @param aux_type The type of auxiliary counts in VSS algorithm.
#' 'One_rep' uses one replicate to construct base signals and
#' another one to construct auxiliary signals. 'mix' uses
#' both replicates to construct each of base and auxiliary
#' signals. 'One_rep' should be used when library sizes are different
#' and data is not normalized.
#'
#' @param group_type The strategy to group auxiliary counts to
#' calculate empirical means and variances. 'bin' divides sorted
#' auxiliary counts into equal size bins and calculate mean and
#' variance in each bin. 'value' groups auxiliary counts based
#' on their corresponding base counts and calculate mean and
#' variance in each group. For rare signals (large signal values)
#' where the size of groups decrease, bin strategy is used for
#' more robust estimations.
#'
#' @param zero_bin Due to excess of zero counts, there are
#' significant number of bins corresponding to base count = 0.
#' Including large amount of such mean and variance estimations
#' might cause a problem for fitting the mean-variance relationship.
#' zero_bin = TRUE creates one auxiliary bin corresponding to
#' zero base counts.
#'
#' @param bin_size The size of the bin in a grouping step.
#'
#' @param log_scale If TRUE, spline is learned using log-scaled
#' means and variances.

#' @return A list of two objects, a two-columns dataframe
#' including means and variances used for training and a
#' smooth.spline object.
#'
#' @import dplyr
#'
#' @export

vssHiC.train <- function(IntSet,
                         aux_type = 'rep1',
                         group_type = 'bin',
                         zero_bin = TRUE,
                         bin_size = 100,
                         log_scale = TRUE){

  cts = data.frame(assay(IntSet))
  fit.res = vss.train(cts, aux_type, group_type, zero_bin, bin_size, log_scale)
  return(fit.res)

}


#' The variance instability measurement
#'
#' This function calculates a variance instability (VI) measurement based on
#' observed differences between replicates.
#'
#' @param IntSet An InteractionSet object including two replicates of
#' Hi-C contact counts.
#'
#' @return VI value.
#'
#' @export


VI <- function(IntSet, bin){
  cts = assay(IntSet)
  # calculating the variance of each replicate
  var1 = var(cts[, 1])
  var2 = var(cts[, 2])
  scores = data.frame(cts)
  colnames(scores) = c('Base_score', 'Aux_score')
  # binning the pairs based on the value of replicate 1
  scores = scores[order(scores$Base_score),]
  scores$index <- 1:dim(scores)[1]
  scores$index <- ceiling(scores$index / bin)
  # calculating the variance of variance
  vars <- data.frame(scores %>% group_by(index) %>%
                       summarise(var = var(Aux_score)))

  VarOfVars = var(vars$var)
  # calculating the normalized variance of difference variances
  VI = VarOfVars/(var1*var2)
  return(VI)
}



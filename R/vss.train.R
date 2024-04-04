
#' Learning the empirical mean-variance relationship from
#' two replicates.
#'
#' This function takes measurements from two replicates (as a two-column
#' data.frame) to calculate observed variances across different signal intensities
#' by VSS algorithm
#' (https://academic.oup.com/bioinformatics/article/37/23/4383/6308936).
#'
#'
#' @param cts A two-column data.frame including counts from two replicates.
#'
#' @param aux_type The type of auxiliary signals in the VSS algorithm.
#' 'rep1' (or 'rep2') uses 'rep1' (or 'rep2') to construct base signals and
#' 'rep2' (or 'rep1') to construct auxiliary signals. 'mix' uses
#' both replicates to construct each of the base and auxiliary
#' signals. 'rep1' or 'rep2' should be used when library sizes differ
#' and data is not normalized.
#'
#' @param group_type The strategy to group auxiliary counts to
#' calculate empirical means and variances. 'bin' divides a sorted
#' auxiliary counts into equal size bins and calculate mean and
#' variance in each bin. 'value' groups auxiliary counts based
#' on their corresponding base counts and calculate mean and
#' variance in each group. For rare signals (large signal values)
#' where the size of groups decreases, bin strategy is used for
#' more robust estimations.
#'
#' @param zero_bin Due to excess of zero counts, there are
#' significant number of bins corresponding to base count = 0.
#' Including large amounts of such mean and variance estimations
#' might cause a problem in fitting the mean-variance relationship.
#' zero_bin = TRUE creates one auxiliary bin corresponding to
#' zero base counts.
#'
#' @param bin_size The size of the bin in a grouping step.
#'
#' @param log_scale If TRUE, spline is learned using log-scaled
#' means and variances.

#' @return A list of two objects, a three-column data.frame
#' including estimated means, variances, and predicted variances by smooth.spline
#' and a fitted smooth.spline object.
#'
#' @import dplyr
#'
#' @export

vss.train <- function(cts,
                      aux_type = 'rep1',
                      group_type = 'bin',
                      zero_bin = TRUE,
                      bin_size = 100,
                      log_scale = TRUE){

  cts = cts[rowSums(cts) > 0, ]
  row.names(cts) = seq(dim(cts)[1])
  colnames(cts) = c('rep1', 'rep2')
  base_aux_scores = get_base_aux(cts,
                                 aux_type)
  mean_sd = get_vss_mean_sd(base_aux_scores,
                            group_type,
                            zero_bin,
                            bin_size)
  fit.res = vss.fit(mean_sd,
                    log_scale)
  return(list(vss_mean_sd = fit.res$mean_sd, vss_fit = fit.res$vss_fit))
}


get_base_aux <- function(cts,
                         aux_type = 'rep1'){

  # if (aux_type == 'one_rep'){
  #   if (sum(cts$rep1) > sum(cts$rep2)){
  #     scores = data.frame(base = cts$rep2, aux = cts$rep1)
  #   }
  #   else{
  #     scores = data.frame(base = cts$rep1, aux = cts$rep2)
  #   }
  # }
  if (aux_type == 'rep1'){
    scores = data.frame(base = cts$rep2, aux = cts$rep1)
  }
  else if (aux_type == 'rep2'){
    scores = data.frame(base = cts$rep1, aux = cts$rep2)
  }
  else if (aux_type == 'mix'){
    scores = data.frame(base = c(cts$rep1, cts$rep2),
                        aux = c(cts$rep2, cts$rep1))
  }
  else{
    stop('Invalid aux type. Options: rep1, rep2, mix.')
  }

  scores = scores[sample(1:nrow(scores)),]
  scores = scores[order(scores$base),]
  return(scores)
}

get_vss_mean_sd <- function(scores,
                            group_type = 'bin',
                            zero_bin = TRUE,
                            bin_size = 100){

  if (group_type == 'bin'){
    if (zero_bin){
      zero_scores = scores[scores$base == 0, ]
      zero_mean = mean(zero_scores$aux)
      zero_sd = sd(zero_scores$aux)
      nonzero_scores = scores[scores$base != 0, ]
      nonzero_scores$index = 1:nrow(nonzero_scores)
      nonzero_scores$index <- ceiling(nonzero_scores$index / bin_size)
      mean_sd <- data.frame(nonzero_scores %>%
                              group_by(index) %>%
                              summarise(mean = mean(aux),
                                        sd = sd(aux)))
      mean_sd = rbind(data.frame(mean = zero_mean, sd = zero_sd),
                      mean_sd[,c('mean', 'sd')])
    }
    else{
      scores$index = 1:nrow(scores)
      scores$index <- ceiling(scores$index / bin_size)
      mean_sd <- data.frame(scores %>% group_by(index) %>%
                              summarise(mean = mean(aux),
                                        sd = sd(aux)))
    }
  }

  else if (group_type == 'value'){
    means = c()
    stdevs = c()
    curr_base_value = 0
    curr_bin = c()
    for (i in 1:nrow(scores)){
      base_value = scores[i,1]
      aux_value = scores[i,2]
      if (base_value != curr_base_value &
          length(curr_bin) >= bin_size){
        means = c(means, mean(curr_bin))
        stdevs = c(stdevs, sd(curr_bin))
        curr_bin = c(aux_value)
        curr_base_value = base_value
      }
      else{
        curr_bin = c(curr_bin, aux_value)
        curr_base_value = base_value
      }
    }
    mean_sd = data.frame(mean = means, sd = stdevs)
  }
  else{
    stop('Invalid group type. Options: bin, value.')
  }
  mean_sd = mean_sd[order(mean_sd$mean),]
  return (mean_sd)
}


vss.fit <- function(training_mean_sd,
                    log_scale = TRUE){
  training_mean_sd = training_mean_sd[!duplicated(round(training_mean_sd$mean, 2)),]
  training_mean_sd = na.omit(training_mean_sd)
  if (log_scale){
    vss_fit <- smooth.spline(log(training_mean_sd$mean + 1),
                             log(training_mean_sd$sd + 1),
                             lambda = 0.1)
    training_mean_sd$pred_sd = exp(vss_fit$y) - 1
  }
  else{
    vss_fit <- smooth.spline(training_mean_sd$mean,
                             training_mean_sd$sd,
                             cv = TRUE)
    training_mean_sd$pred_sd = vss_fit$y
  }

  return(list(mean_sd = training_mean_sd, vss_fit = vss_fit))
}




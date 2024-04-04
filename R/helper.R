
### downsample replicate (rep_num) from iset with a
### ratio of ds_ratio

downsample <- function(iset, rep_num, ds_ratio){
  as = assay(iset)
  rep = as[, rep_num]
  rep_ls = sum(rep)
  ds_rep = rmultinom(n = 1, size = as.integer(rep_ls * ds_ratio),
                     prob = rep / rep_ls)[,1]
  as[, rep_num] = ds_rep
  ds_iset = data <- InteractionSet(assays=list(counts=as),
                                   interactions(iset))
  return(ds_iset)
}

### downsample replicates 1 and 2 from iset to ss1 and ss2
### sizes respectively

downsample2 <- function(iset, ss1, ss2){
  as = assay(iset)
  ls = colSums(as)
  as[, 1] = rmultinom(n = 1, size = ss1,
                      prob = as[, 1] / ls[1])[,1]
  as[, 2] = rmultinom(n = 1, size = ss2,
                      prob = as[, 2] / ls[2])[,1]
  ds_iset = data <- InteractionSet(assays=list(counts=as),
                                   interactions(iset))
  return(ds_iset)
}


### Given a set of trained mean-var, create a dataframe
### all to compare

mul_mean_var <- function(train_list,
                         max_ct = 500,
                         suffix = ""){
  xg = sinh( seq(asinh(0), asinh(max_ct),
                 length.out = 1000) )[-1]
  mul_mean_var = data.frame(mean = xg)
  for (nm in names(train_list)){
    xg_stdevs = exp(stats:::predict.smooth.spline(train_list[[nm]]$vss_fit,
                                                  log(xg + 1))$y) - 1
    mul_mean_var[, paste0(nm, suffix)] = xg_stdevs
  }
  return(melt(mul_mean_var, id.var = 'mean'))
}



hft_transform <- function(IntSet, chromoR = FALSE, ex_diag = FALSE){
  rg = regions(interactions(IntSet))
  chroms = unique(seqnames(rg))
  if (length(chroms) == 1){
    rg1 = rg
    rg2 = rg
    type = "intra"
  }
  else if(length(chroms) == 2){
    rg1 = rg[seqnames(rg) == chroms[1]]
    rg2 = rg[seqnames(rg) == chroms[2]]
    type = "inter"
  }
  else{
    stop("Only supports InteractionSet including one pair of chromosomes")
  }
  transformed_isets = list()
  for (rep_num in c(1, 2)){
    m = inflate(IntSet, rg1, rg2, assay = 1, sample = rep_num)
    m = as.matrix(m)
    na_idx = is.na(m)
    m[na_idx] = 0
    if (type == "intra" & ex_diag){
      y = list(diag = diag(m), nondiag = m[upper.tri(m)])
    } else if(type == "intra" & !ex_diag){
      y = list(whole = m[upper.tri(m, diag = T)])
    }
    else{
      y = list(whole = as.vector(m))
    }

    transformed_y = lapply(y, function(y){
      orig_length = length(y)
      padding_length = 2 ** ceiling(log2(orig_length)) - orig_length
      yy = c(rep(0, padding_length), y)
      yHft = hft(yy)
      if (chromoR){
        # Denoise - shrink coefficients for denoising
        yHftWd = wd(yHft , filter.number=1, family="DaubExPhase")
        yHftWd = threshold(yHftWd, policy = "cv", dev=madmad)
        yHftWr = wr(yHftWd) # inverse shranked coeff.
        # Apply the inverse Fisz Haar transform to recontrcit the corrected sequence
        final_y = yHftWr[(padding_length + 1) :
                           (padding_length + orig_length)]
      }
      else{
        final_y = yHft[(padding_length + 1) :
                         (padding_length + orig_length)]
      }
      return(final_y)
    })

    if (type == "intra" & ex_diag){
      transformed_m = matrix(data = NA, nrow = nrow(m), ncol = ncol(m))
      diag(transformed_m) = transformed_y$diag
      transformed_m[upper.tri(transformed_m)] = transformed_y$nondiag
      transformed_m[na_idx] = NA
    } else if (type == "intra" & !ex_diag){
      transformed_m = matrix(data = NA, nrow = nrow(m), ncol = ncol(m))
      transformed_m[upper.tri(transformed_m, diag = T)] = transformed_y$whole
      transformed_m[na_idx] = NA
    }
    else{
      transformed_m = matrix(transformed_y$whole, nrow = nrow(m), ncol = ncol(m))
    }
    transformed_CM = ContactMatrix(transformed_m, rg1, rg2)
    transformed_isets[[rep_num]] = deflate(transformed_CM)
  }
  transformed_iset = do.call(cbind, transformed_isets)
  transformed_iset = sort(transformed_iset)

}

# hft for multiple isets
hft_transform2 <- function(IntSet_list, chromoR = FALSE,
                           type = "inter"){
  transformed_isets = list()
  for (name in names(IntSet_list)){
    transformed_isets[[name]] = list()
  }
  for (rep_num in c(1, 2)){
    y = list(whole = c())
    for (iset in IntSet_list){
      rg = regions(interactions(iset))
      chroms = unique(seqnames(rg))
      if (type == "intra"){
        rg1 = rg
        rg2 = rg
      }
      else{
        rg1 = rg[seqnames(rg) == chroms[1]]
        rg2 = rg[seqnames(rg) == chroms[2]]
      }
      m = inflate(iset, rg1, rg2, assay = 1, sample = rep_num)
      m = as.matrix(m)
      na_idx = is.na(m)
      m[na_idx] = 0
      if(type == "intra"){
        y$whole = c(y$whole, m[upper.tri(m, diag = T)])
      }
      else{
        y$whole = c(y$whole, as.vector(m))
      }
    }
    transformed_y = lapply(y, function(y){
      orig_length = length(y)
      padding_length = 2 ** ceiling(log2(orig_length)) - orig_length
      yy = c(rep(0, padding_length), y)
      yHft = hft(yy)
      if (chromoR){
        sample_idx = sample((length(yHft) - 10000):length(yHft),10)
        # Denoise - shrink coefficients for denoising
        yHftWd = wd(yHft , filter.number=1, family="DaubExPhase")
        yHftWd = threshold(yHftWd, policy = "cv", dev=madmad)
        yHftWr = wr(yHftWd) # inverse shranked coeff.
        # Apply the inverse Fisz Haar transform to recontrcit the corrected sequence
        final_y = yHftWr[(padding_length + 1) :
                           (padding_length + orig_length)]

      }
      else{
        final_y = yHft[(padding_length + 1) :
                         (padding_length + orig_length)]
      }
      return(final_y)
    })

    offset = 0

    for (name in names(IntSet_list)){
      iset = IntSet_list[[name]]
      rg = regions(interactions(iset))
      if (type == "intra"){
        y_length = length(rg) * (length(rg) + 1) / 2
        transformed_m = matrix(data = NA, nrow = length(rg), ncol = length(rg))
        transformed_m[upper.tri(transformed_m, diag = T)] =
          transformed_y$whole[(offset + 1) : (offset + y_length)]
        offset = offset + y_length
        transformed_m[na_idx] = NA
        transformed_CM = ContactMatrix(transformed_m, rg, rg)
      }
      else{
        chroms = unique(seqnames(rg))
        rg1 = rg[seqnames(rg) == chroms[1]]
        rg2 = rg[seqnames(rg) == chroms[2]]
        y_length = length(rg1) * length(rg2)
        transformed_m = matrix(transformed_y$whole[(offset + 1) : (offset + y_length)],
                               nrow = length(rg1), ncol = length(rg2))
        offset = offset + y_length
        transformed_CM = ContactMatrix(transformed_m, rg1, rg2)
      }

      transformed_isets[[name]][[rep_num]] = deflate(transformed_CM)

    }
  }

  output_transformed_isets = list()
  for (name in names(transformed_isets)){
    output_transformed_isets[[name]] = do.call(cbind, transformed_isets[[name]])
    output_transformed_isets[[name]] = sort(output_transformed_isets[[name]])
  }
  return(output_transformed_isets)
}

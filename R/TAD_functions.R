### wrapper functions for TAD callers


#' Apply HiCseg on InteractionSet
#'
#' This is a wrapper function for HiCseg to take InteractionSet
#' object as input and return TADs as GRanges object.
#'
#' @param iset Input InteractionSet.
#' @param rep_num The replicate's number in InteractionSet to
#' be used for TAD annotation.
#' @param max_cp HiCseg's parameter: maximum number of
#' change points.
#' @param model HiCSeg's parameter: "D" for block-diagonal
#' and "Dplus" for the extended block-diagonal model.
#' @param dist HiCseg's parameter: Distribution of the data:
#' "B" is for Negative Binomial distribution,
#' "P" is for the Poisson distribution and
#' "G" is for the Gaussian distribution.
#'
#' @return A GRanges object including identified TADs' ranges.
#'
#' @importFrom HiCseg HiCseg_linkC_R
#' @import InteractionSet
#' @importFrom InteractionSet inflate
#' @export

iset_hicseg <- function(iset,
                        rep_num,
                        max_cp = 40,
                        model = "D",
                        dist = "G"){
  rg = regions(interactions(iset))
  CM = inflate(iset, rg, rg, assay = 1, sample = rep_num)
  CM = as.matrix.Vector(CM)
  #CM = as(CM, "Matrix")
  CM[is.na(CM)] = 0
  CM_size = dim(CM)[1]
  result = HiCseg_linkC_R(CM_size, max_cp, dist, CM, model)
  #return(list(matrix = CM, result = result))
  t_hat = c(result$t_hat[result$t_hat != 0])
  starts = start(ranges(rg)) - 1
  boundaries = starts[t_hat]
  chrom = unique(seqnames(regions(interactions(iset))))
  g_range = GRanges(seqnames = rep(chrom, length(boundaries) - 1),
                    IRanges(start = boundaries[1: (length(boundaries) - 1)],
                            end = boundaries[2: length(boundaries)]))
  return(g_range)
}

#' Apply TopDom on InteractionSet
#'
#' This is a wrapper function for TopDom to take InteractionSet
#' object as input and return TADs as GRanges object.
#'
#' @param iset Input InteractionSet.
#' @param rep_num The replicate's number in InteractionSet to
#' be used for TAD annotation.
#' @param res resolution of iset.
#' @param window_size TopDom's parameter: The number of bins
#' to extend (as a non-negative integer).
#'
#' @return A GRanges object including identified TADs' ranges.
#'
#' @import InteractionSet
#' @importFrom InteractionSet inflate
#' @export

iset_topdom <- function(iset,
                        rep_num,
                        res,
                        window_size){

  rg = regions(interactions(iset))
  starts = start(ranges(regions(interactions(iset)))) - 1
  CM = inflate(iset, rg, rg, assay = 1, sample = rep_num)
  CM = as.matrix.Vector(CM)
  CM[is.na(CM)] = 0
  window_size = window_size / res
  tads <- TopDom(CM, window_size)
  tads$start = starts[tads$bin]
  #return(tads)
  boundaries = tads$start[tads$pvalue_boundaries]
  chrom = unique(seqnames(regions(interactions(iset))))
  g_range = GRanges(seqnames = rep(chrom, length(boundaries) - 1),
                    IRanges(start = boundaries[1: (length(boundaries) - 1)],
                            end = boundaries[2: length(boundaries)]))
  return(g_range)
}

#' Apply SpectralTAD on InteractionSet
#'
#' This is a wrapper function for SpectralTAD to take
#' InteractionSet object as input and return TADs as
#' GRanges object.
#'
#' @param iset Input InteractionSet.
#' @param rep_num The replicate's number in InteractionSet to
#' be used for TAD annotation.
#' @param res resolution of iset.
#' @param window_size SpectralTAD's parameter: The size of
#' the sliding window for calculating TADs.
#'
#' @return A GRanges object including identified TADs' ranges.
#'
#' @import InteractionSet
#' @importFrom InteractionSet inflate
#' @export

iset_spectralTAD <- function(iset,
                             rep_num,
                             res,
                             window_size){

  rg = regions(interactions(iset))
  CM = inflate(iset, rg, rg, assay = 1, sample = rep_num)
  CM = as.matrix.Vector(CM)
  CM[is.na(CM)] = 0
  colnames(CM) = start(ranges(rg)) - 1 # seq(0:(dim(CM)[1]-1))*res
  window_size = window_size / res
  chrom = unique(seqnames(regions(interactions(iset))))
  tads <- SpectralTAD(CM, chrom, resolution = res,
                      min_size = 1, window_size = window_size)
  #return(tads$Level_1)
  boundaries <- unique(c(tads$Level_1$start, tads$Level_1$end))
  boundaries = sort(boundaries)
  #return(boundaries)
  chrom = unique(seqnames(regions(interactions(iset))))
  g_range = GRanges(seqnames = rep(chrom, length(boundaries) - 1),
                    IRanges(start = boundaries[1: (length(boundaries) - 1)],
                            end = boundaries[2: length(boundaries)]))
  return(g_range)
}

### TAD evaluation functions

around_boundary2 <- function(boundaries,
                             peaks_vec,
                             res_step,
                             window_size,
                             chr_size){

  chr_size <- ceiling(chr_size/(res_step*1000))
  num_step = window_size / (res_step*1000)
  around_bound = rep(0,(2*num_step)+1)
  for (bound in boundaries){
    bound_id = bound / (res_step*1000)
    if (bound_id > num_step & bound_id <= chr_size - num_step){
      around_bound = around_bound + peaks_vec[(bound_id-num_step):(bound_id+num_step)]
    }
  }
  around_bound = around_bound/length(boundaries)
  return(around_bound)

}

enrichment_df2 <- function(boundaries,
                           list_of_peaks,
                           TF_names,
                           res_step,
                           window_size,
                           chr_size){

  num_step = window_size / (res_step * 1000)
  enrichment_df = data.frame(pos=c(-num_step:num_step))
  for (i in seq(1:length(TF_names))){
    enrichment_df[TF_names[i]] = around_boundary2(boundaries, list_of_peaks[[i]],
                                                  res_step, window_size, chr_size)
  }
  enrichment_df = melt(enrichment_df, id.vars='pos')
  colnames(enrichment_df) = c('pos', 'TF', 'avg_peak')
  return(enrichment_df)
}

read_peak_file <- function(peak_filepath, res_step, chroms, chr_sizes){
  peaks <- read.table(peak_filepath,
                      quote = '', stringsAsFactors=FALSE,
                      sep = "\t", header = FALSE)[,c(1,2,3)]
  peaks$V2 <- round(peaks$V2/(res_step*1000))
  peaks$V3 <- round(peaks$V3/(res_step*1000))
  peaks_vecs = hash()
  for (chrom in chroms){
    chr_size <- chr_sizes[chr_sizes[,1]==chrom,2]
    chr_size <- ceiling(chr_size/(res_step*1000))
    peaks_vecs[[chrom]] <- rep(0,chr_size)
    chr_peaks = peaks[peaks$V1==chrom,]
    for (i in c(1:nrow(chr_peaks))){
      peaks_vecs[[chrom]][chr_peaks[i,'V2']:chr_peaks[i,'V3']] =
        peaks_vecs[[chrom]][chr_peaks[i,'V2']:chr_peaks[i,'V3']] + 1
    }
  }
  return(peaks_vecs)
}

around_boundary <- function(boundaries_df,
                            peaks_vec,
                            res_step,
                            window,
                            chr_sizes){

  num_step = window / (res_step * 1000)
  around_bound = rep(0, (2 * num_step) + 1)
  num_bounds = 0
  for (i in c(1:nrow(boundaries_df))){
    bound_chr = boundaries_df[i, 1]
    bound_pos = boundaries_df[i, 2] / (res_step * 1000)
    chr_size = chr_sizes[chr_sizes$V1 == bound_chr, 2]
    chr_size <- ceiling(chr_size / (res_step * 1000))

    if (bound_pos > num_step & bound_pos <= chr_size - num_step){

      around_bound = around_bound + peaks_vec[[bound_chr]][(bound_pos - num_step) : (bound_pos + num_step)]
      num_bounds = num_bounds + 1
    }
  }
  around_bound = around_bound / num_bounds
  return(around_bound)
}

# input: boundaries is a 2-column (chr, boundary) dataframe

enrichment_df <- function(boundaries_df,
                          list_of_peaks,
                          TF_names,
                          res_step,
                          window,
                          chr_sizes){

  window_size = as.integer(window / (res_step * 1000))
  enrichment_df = data.frame(pos=c(-window_size : window_size))
  for (i in seq(1 : length(TF_names))){
    enrichment_df[TF_names[i]] = around_boundary(boundaries_df,
                                                 list_of_peaks[[i]],
                                                 res_step,
                                                 window,
                                                 chr_sizes)
  }
  enrichment_df = melt(enrichment_df, id.vars='pos')
  colnames(enrichment_df) = c('pos', 'TF', 'avg_peak')
  return(enrichment_df)
}


all_boundary_fc <- function(boundaries_df, TF_peaks, res,
                            res_step, chr_sizes){
  around_bins = ceiling(100000 / (res_step*1000))
  control_start = 400000 / (res_step*1000)
  control_end = 500000 / (res_step*1000)
  boundary_fc <- function(boundary_chr, boundary_pos){
    chr_size = chr_sizes[chr_sizes$V1 == boundary_chr, 2]
    chr_size = ceiling(chr_size / (res_step * 1000))
    bound_id = boundary_pos / (res_step * 1000)
    if ((bound_id > control_end) & (bound_id <= chr_size - control_end)){
      around_bound = mean(TF_peaks[[boundary_chr]][(bound_id - around_bins):
                                                  (bound_id + around_bins)])
      us_avg = mean(TF_peaks[[boundary_chr]][(bound_id - control_end):
                                            (bound_id - control_start)])
      ds_avg = mean(TF_peaks[[boundary_chr]][(bound_id + control_start):
                                            (bound_id + control_end)])
      distant_avg = max(mean(us_avg, ds_avg), 1/40)

      fc = (around_bound / distant_avg) - 1
      return(fc)
      #return(c(around_bound, distant_avg))
    }
    else{
      return(NA)
    }
  }

  enrichments = mapply(function(chr, pos)
    boundary_fc(chr, pos),
    boundaries_df$chr, boundaries_df$boundary,
    SIMPLIFY = T)
  return(unlist(enrichments))
  #return(sum(enrichments[1,], na.rm = T) / sum(enrichments[2,], na.rm = T) - 1)
}

# boundary_fc <- function(boundaries, peaks_vec, res_step, chr_size){
#   chr_size <- ceiling(chr_size/(res_step*1000))
#   control_start = 400000 / (res_step*1000)
#   control_end = 500000 / (res_step*1000)
#   fcs = c()
#   for (bound in boundaries){
#     bound_id = bound / (res_step*1000)
#     if (bound_id > control_end & bound_id <= chr_size - control_end){
#       around_bound = mean(peaks_vec[(bound_id-2):(bound_id+2)])
#       us_avg = mean(peaks_vec[(bound_id-control_end):(bound_id-control_start)])
#       ds_avg = mean(peaks_vec[(bound_id+control_start):(bound_id+control_end)])
#       distant_avg = max(mean(us_avg,ds_avg),1/40)
#
#       fcs = c(fcs,around_bound/distant_avg)
#     }
#   }
#   return(fcs)
# }

### helper functions corresponding to TAD callers

# TopDom

get_TopDom_boundaries <- function(mat, window_size, res){
  tads <- TopDom(mat, window_size)
  boundaries = which(tads$pvalue_boundaries==TRUE)*res
  return (boundaries)
}

TopDom <- function(mat, window_size){
  n_bins = dim(mat)[1]
  binSignals = rep(0,n_bins)
  local.ext = rep(-0.5,n_bins)
  pvalue <- rep(1, times = n_bins)
  gaps <- rep(0,n_bins)
  proc_regions = data.frame()
  curr_start = NULL
  curr_end = NULL
  # calculating binSignals
  for (i in 1:n_bins-1){
    lb <- max( 1, i-window_size+1 )
    ub <- min( i+window_size, n_bins)
    binSignals[i] <- mean(mat[lb:i, (i+1):ub])
    is_gap <- (sum(mat[i, max(1, i-window_size):ub]) == 0)
    if(is_gap){
      gaps[i] = 1
    }
    if (is.null(curr_start) && !is_gap){
      curr_start = i
    }
    if (!is.null(curr_start) && is_gap){
      curr_end <- i-1
      if ((curr_end - curr_start) > 2){
        proc_regions <- rbind(proc_regions, c(curr_start, curr_end))
      }
      curr_start = NULL
      curr_end = NULL
    }
  }
  is_gap <- (sum(mat[n_bins, max(1, n_bins-window_size):n_bins]) == 0)
  if (!is.null(curr_start) && is_gap){
    curr_end = n_bins - 1
    proc_regions <- rbind(proc_regions, c(curr_start, curr_end))
  }
  if (!is.null(curr_start) && !is_gap){
    curr_end = n_bins
    proc_regions <- rbind(proc_regions, c(curr_start, curr_end))
  }
  colnames(proc_regions) <- c('start', 'end')
  # detecting boundaries
  for(i in 1:nrow(proc_regions))
  {
    start = proc_regions[i, "start"]
    end = proc_regions[i, "end"]
    local.ext[start:end] = Detect.Local.Extreme(x=binSignals[start:end])
  }
  scaled_mat = mat
  for( i in 1:(2*window_size) )
  {
    scaled_mat[ seq(1+(n_bins*i), n_bins*n_bins, 1+n_bins) ] =
      scale( mat[ seq(1+(n_bins*i), n_bins*n_bins, 1+n_bins) ] )
  }
  for( i in 1:nrow(proc_regions))
  {
    start = proc_regions[i, "start"]
    end = proc_regions[i, "end"]
    pvalue[start:end] <- Get.Pvalue(matrix.data=scaled_mat[start:end, start:end],
                                    size=window_size, scale=1)
  }
  ext_boundaries = rep(0,n_bins)
  ext_boundaries[local.ext==-1] = 1
  pvalue_boundaries = rep(0,n_bins)
  pvalue_boundaries[pvalue<0.05] = 1
  domains <- get_domains(proc_regions, ext_boundaries)
  return(data.frame(bin = seq(n_bins), gap = gaps, binSignals = binSignals,
                    ext_boundaries = ext_boundaries,
                    pvalue_boundaries = (ext_boundaries & pvalue_boundaries)
  ))
  #return(list(data.frame(gap = gaps, binSignals = binSignals, ext_boundaries = ext_boundaries,
  #                  pvalue_boundaries = pvalue_boundaries),domains))
}

get_domains <- function(proc_regions, boundaries){
  domains = data.frame()
  for (i in seq(nrow(proc_regions))){
    start = proc_regions[i,1]
    end = proc_regions[i,2]
    curr_boundaries = start - 1 + which(boundaries[start:end]==1)
    if (length(curr_boundaries)>1){
      for (b in seq(1:(length(curr_boundaries)-1))){
        domain = data.frame(start = curr_boundaries[b], end = curr_boundaries[b+1])
        domains = rbind(domains, domain)
      }
    }
  }
  return (domains)
}
# This function is from https://github.com/HenrikBengtsson/TopDom/blob/develop/R/TopDom_0.0.2.R
Detect.Local.Extreme <- function(x)
{
  n_bins = length(x)
  ret = rep(0, n_bins)
  x[is.na(x)]=0

  if(n_bins <= 3)
  {
    ret[which.min(x)]=-1
    ret[which.max(x)]=1

    return(ret)
  }
  # Norm##################################################3
  new.point = Data.Norm(x=1:n_bins, y=x)
  x=new.point$y
  ##################################################
  cp = Change.Point(x=1:n_bins, y=x)

  if( length(cp$cp) <= 2 ) return(ret)
  if( length(cp$cp) == n_bins) return(ret)
  for(i in 2:(length(cp$cp)-1))
  {
    if( x[cp$cp[i]] >= x[cp$cp[i]-1] && x[cp$cp[i]] >= x[cp$cp[i]+1] ) ret[cp$cp[i]] = 1
    else if(x[cp$cp[i]] < x[cp$cp[i]-1] && x[cp$cp[i]] < x[cp$cp[i]+1]) ret[cp$cp[i]] = -1

    min.val = min( x[ cp$cp[i-1] ], x[ cp$cp[i] ] )
    max.val = max( x[ cp$cp[i-1] ], x[ cp$cp[i] ] )

    if( min( x[cp$cp[i-1]:cp$cp[i]] ) < min.val ) ret[ cp$cp[i-1] - 1 + which.min( x[cp$cp[i-1]:cp$cp[i]] ) ] = -1
    if( max( x[cp$cp[i-1]:cp$cp[i]] ) > max.val ) ret[ cp$cp[i-1] - 1 + which.max( x[cp$cp[i-1]:cp$cp[i]] ) ] = 1
  }

  return(ret)
}

Data.Norm <- function(x, y)
{
  ret.x = rep(0, length(x))
  ret.y = rep(0, length(y))

  ret.x[1] = x[1]
  ret.y[1] = y[1]

  diff.x = diff(x)
  diff.y = diff(y)

  scale.x = 1 / mean( abs(diff(x) ) )
  scale.y = 1 / mean( abs( diff(y) ) )

  #print(scale.x)
  #print(scale.y)

  for(i in 2:length(x))
  {
    ret.x[i] = ret.x[i-1] + (diff.x[i-1]*scale.x)
    ret.y[i] = ret.y[i-1] + (diff.y[i-1]*scale.y)
  }

  return(list(x=ret.x, y=ret.y))
}

Change.Point <- function( x, y )
{
  if( length(x) != length(y))
  {
    print("ERROR : The length of x and y should be the same")
    return(0)
  }

  n_bins <- length(x)
  Fv <- rep(NA, n_bins)
  Ev <- rep(NA, n_bins)
  cp <- 1

  i=1
  Fv[1]=0
  while( i < n_bins )
  {
    j=i+1
    Fv[j] = sqrt( (x[j]-x[i])^2 + (y[j] - y[i] )^2 )

    while(j<n_bins)
    {
      j=j+1
      k=(i+1):(j-1)
      Ev[j] = ( sum( abs( (y[j]-y[i] )*x[k] - (x[j] -x[i])*y[k] - (x[i]*y[j]) + (x[j]*y[i]) ) ) / sqrt( (x[j]-x[i])^2 + (y[j] - y[i] )^2 ) )
      Fv[j] = sqrt( (x[j]-x[i])^2 + (y[j] - y[i])^2 ) - ( sum( abs( (y[j]-y[i] )*x[k] - (x[j] -x[i])*y[k] - (x[i]*y[j]) + (x[j]*y[i]) ) ) / sqrt( (x[j]-x[i])^2 + (y[j] - y[i] )^2 ) )

      #################################################
      #Not Original Code
      if( is.na(Fv[j]) || is.na(Fv[j-1]) ) {
        j = j-1
        cp <- c(cp, j)
        break
      }
      ####################################################3
      if(Fv[j] < Fv[j-1] ) {
        j = j - 1
        cp <- c(cp, j )
        break
      }
    }
    i=j
  }

  cp <- c(cp, n_bins)

  return(list(cp=cp, objF=Fv, errF=Ev))
}

Get.Pvalue <- function( matrix.data, size, scale=1 )
{
  n_bins = nrow(matrix.data)
  pvalue <- rep(1, n_bins)

  for( i in 1:(n_bins-1) )
  {
    dia = as.vector( Get.Diamond.Matrix2(matrix.data, i, size=size) )
    ups = as.vector( Get.Upstream.Triangle(matrix.data, i, size=size) )
    downs = as.vector( Get.Downstream.Triangle(matrix.data, i, size=size) )
    wil.test =  wilcox.test(x=dia*scale, y=c(ups, downs), alternative="less", exact=F)
    pvalue[i] = wil.test$p.value

    #print(paste(i, "=", wil.test$p.value) )
  }

  pvalue[ is.na(pvalue) ] = 1
  return(pvalue)
}

Get.Diamond.Matrix2 <- function(mat.data, i, size)
{
  n_bins = nrow(mat.data)
  new.mat = matrix(rep(NA, size*size), nrow=size, ncol=size)

  for(k in 1:size)
  {
    if(i-(k-1) >= 1 && i < n_bins)
    {
      lower = min(i+1, n_bins)
      upper = min(i+size, n_bins)

      new.mat[size-(k-1), 1:(upper-lower+1)] = mat.data[i-(k-1), lower:upper]
    }
  }

  return(new.mat)
}

Get.Upstream.Triangle <- function(mat.data, i, size)
{
  n_bins = nrow(mat.data)

  lower = max(1, i-size)
  tmp.mat = mat.data[lower:i, lower:i]
  return( tmp.mat[ upper.tri( tmp.mat, diag=F ) ] )
}

Get.Downstream.Triangle <- function(mat.data, i, size)
{
  n_bins = nrow(mat.data)
  if(i==n_bins) return(NA)

  upperbound = min(i+size, n_bins)
  tmp.mat = mat.data[(i+1):upperbound, (i+1):upperbound]
  return( tmp.mat[ upper.tri( tmp.mat, diag=F ) ] )
}

Convert.Bin.To.Domain.TMP <- function(bins, signal.idx, gap.idx, pvalues = NULL, pvalue.cut = NULL) {
  n_bins <- nrow(bins)
  ret <- data.frame(chr = character(0), from.id = numeric(0), from.coord = numeric(0), to.id = numeric(0), to.coord = numeric(0), tag = character(0), size = numeric(0), stringsAsFactors = FALSE)
  levels(x = ret[, "tag"]) <- c("domain", "gap", "boundary")

  rmv.idx <- setdiff(seq_len(n_bins), gap.idx)
  proc.region <- Which.process.region(rmv.idx, n_bins = n_bins, min.size = 0L)
  from.coord <- bins[proc.region[, "start"], "from.coord"]
  n_procs <- nrow(proc.region)
  zeros <- double(length = n_procs)
  gap <- data.frame(chr = rep(bins[1, "chr"], times = n_procs), from.id = zeros, from.coord = from.coord, to.id = zeros, to.coord = zeros, tag = rep("gap", times = n_procs), size = zeros, stringsAsFactors = FALSE)

  rmv.idx <- union(signal.idx, gap.idx)
  proc.region <- Which.process.region(rmv.idx, n_bins = n_bins, min.size = 0L)
  n_procs <- nrow(proc.region)
  from.coord <- bins[proc.region[, "start"], "from.coord"]
  zeros <- double(length = n_procs)
  domain <- data.frame(chr = rep(bins[1, "chr"], times = n_procs), from.id = zeros, from.coord = from.coord, to.id = zeros, to.coord = zeros, tag = rep("domain", times = n_procs), size = zeros, stringsAsFactors = FALSE)

  rmv.idx <- setdiff(seq_len(n_bins), signal.idx)
  proc.region <- as.data.frame(Which.process.region(rmv.idx, n_bins = n_bins, min.size = 1L))
  n_procs <- nrow(proc.region)
  if (n_procs > 0) {
    from.coord <- bins[proc.region[, "start"] + 1, "from.coord"]
    zeros <- double(length = n_procs)
    boundary <- data.frame(chr = rep(bins[1, "chr"], times = n_procs), from.id = zeros, from.coord = from.coord, to.id = zeros, to.coord = zeros, tag = rep("boundary", times = n_procs), size = zeros, stringsAsFactors = FALSE)
    ret <- rbind(ret, boundary)
  }

  if (nrow(domain) == 0L) {
    ret <- gap
  } else {
    ret <- rbind(gap, domain)
    ret <- ret[order(ret[, 3]), ]

    ## FIXME: Below code assumes nrow(ret) >= 2
    ret[, "to.coord"] <- c(ret[2:nrow(ret), "from.coord"], bins[n_bins, "to.coord"])
    ret[, "from.id"] <- match(ret[, "from.coord"], table = bins[, "from.coord"])
    ret[, "to.id"] <- match(ret[, "to.coord"], table = bins[, "to.coord"])
    ret[, "size"] <- ret[, "to.coord"] - ret[, "from.coord"]

    if (!is.null(pvalues) && !is.null(pvalue.cut)) {
      for (i in seq_len(nrow(ret))) {
        if (ret[i, "tag"] == "domain") {
          domain.bins.idx <- ret[i, "from.id"]:ret[i, "to.id"]
          p.value.constr <- which(pvalues[domain.bins.idx] < pvalue.cut)

          if (length(domain.bins.idx) == length(p.value.constr)) ret[i, "tag"] <- "boundary"
        }
      }
    }
  }

  ret
}

# SpectralTAD

SpectralTAD = function(cont_mat, chr, levels = 1, qual_filter = FALSE,
                       z_clust = FALSE, eigenvalues = 2, min_size = 5,
                       window_size = 25,
                       resolution = "auto", gap_threshold = 1,
                       grange = FALSE, out_format = "none", out_path = chr) {

  #Disable scientific notation
  options(scipen = 999)

  #Calculate the number of rows and columns of the contact matrix


  if (missing("chr")) {
    stop("Must specify chromosome")
  }

  row_test = dim(cont_mat)[1]
  col_test = dim(cont_mat)[2]

  if (row_test == col_test) {
    if (all(is.finite(cont_mat)) == FALSE) {
      stop("Contact matrix must only contain real numbers")
    }
  }

  if (col_test == 3) {

    if (!is.matrix(cont_mat)) {
      cont_mat = as.matrix(cont_mat)
    }

    #Convert sparse matrix to n x n matrix

    message("Converting to n x n matrix")

    if (nrow(cont_mat) == 1) {
      stop("Matrix is too small to convert to full")
    }
    cont_mat = HiCcompare::sparse2full(cont_mat)

    if (all(is.finite(cont_mat)) == FALSE) {
      stop("Contact matrix must only contain real numbers")
    }

    if (resolution == "auto") {
      message("Estimating resolution")
      resolution = as.numeric(names(table(as.numeric(colnames(cont_mat))-dplyr::lag(as.numeric(colnames(cont_mat)))))[1])
    }

  } else if (col_test-row_test == 3) {

    message("Converting to n x n matrix")

    #Find the start coordinates based on the second column of the bed file portion of matrix

    start_coords = cont_mat[,2]

    #Calculate resolution based on given bin size in bed file

    resolution = as.numeric(cont_mat[1,3])-as.numeric(cont_mat[1,2])

    #Remove bed file portion

    cont_mat = as.matrix(cont_mat[,-c(seq_len(3))])

    if (all(is.finite(cont_mat)) == FALSE) {
      stop("Contact matrix must only contain real numbers")
    }

    #Make column names correspond to bin start

    colnames(cont_mat) = start_coords

  } else if (col_test!=3 & (row_test != col_test) & (col_test-row_test != 3)) {

    #Throw error if matrix does not correspond to known matrix type

    stop("Contact matrix must be sparse or n x n or n x (n+3)!")

  } else if ( (resolution == "auto") & (col_test-row_test == 0) ) {
    message("Estimating resolution")

    #Estimating resolution based on most common distance between loci

    resolution = as.numeric(names(table(as.numeric(colnames(cont_mat))-dplyr::lag(as.numeric(colnames(cont_mat)))))[1])
  }

  if (resolution>200000) {
    stop("Resolution must be less than (or equal to) 200kb")
  }

  if (nrow(cont_mat) < 2000000/resolution) {
    stop("Matrix must be larger than 2 megabases divided by resolution")
  }

  #Performed window spectral clustering

  bed = .windowedSpec(cont_mat, chr = chr, resolution = resolution,
                      z_clust = z_clust, eigenvalues = eigenvalues,
                      min_size = min_size, window_size = window_size,
                      qual_filter = qual_filter, gap_threshold = gap_threshold) %>%
    mutate(Level = 1)

  #Calculate the end point of TADs based on bin instead of genomic coordinate

  coords = cbind(match(bed$start, as.numeric(colnames(cont_mat))), match(bed$end-resolution, as.numeric(colnames(cont_mat))))

  #Create a list of tad start and end points

  tads = apply(coords, 1, function(x) cont_mat[x[1]:x[2], x[1]:x[2]])


  called_tads = list(bed)

  #Initialize TAD level

  curr_lev = 2

  while (curr_lev != (levels + 1) ) {

    #Get a list of TAD coordinates at the previous level

    coords = cbind(match(called_tads[[curr_lev-1]]$start, as.numeric(colnames(cont_mat))), match(called_tads[[curr_lev-1]]$end-resolution, as.numeric(colnames(cont_mat))))

    #Get tads that are less than the twice the minmium length and thus not seperable

    less_5 = which( (coords[,2]-coords[,1])<min_size*2  )

    if (length(less_5)>0) {

      #Remove TADs which cannot be seperate

      pres_tads = called_tads[[curr_lev-1]][less_5,]

      coords = coords[-less_5, ]

    } else {
      pres_tads = c()
    }

    #Account for the situation where there is only 1 potential sub-tad

    if (is.null(nrow(coords))) {
      coords = t(as.matrix(coords))
    }

    tads = apply(coords, 1, function(x) cont_mat[x[1]:x[2], x[1]:x[2]])



    #Remove sub-tads with too many zeros

    zeros = which(unlist(lapply(tads, function(x) nrow(x)-sum(rowSums(x)==0)))<min_size*2)

    if (length(zeros)>0) {
      pres_tads = rbind(pres_tads, called_tads[[curr_lev-1]][zeros,])
      tads[zeros] = NULL
    }

    #Calculate sub-TADs for each seperable TAD

    sub_tads = lapply(tads, function(x) {
      .windowedSpec(x, chr =chr, resolution = resolution, qual_filter = qual_filter, z_clust = TRUE, min_size = min_size)
    })

    #Convert sub-TADs to BED format

    called_tads[[curr_lev]] = bind_rows(sub_tads, pres_tads) %>% mutate(Level = curr_lev) %>% arrange(start)

    curr_lev = curr_lev+1

  }

  #Assign names based on levels

  names(called_tads) = paste0("Level_", seq_len(levels))

  if ( !(out_format == "none")) {
    if (out_format %in% c("bedpe", "juicebox")) {
      #Get just coordinates
      bed_out = bind_rows(called_tads) %>%
        dplyr::select(chr,start,end)
      #Same object, different colnames
      bed_out1 <- bed_out
      colnames(bed_out1) <- c("chr1", "start1", "end1")
      #Combine into first six columns of bedpe and add extra columns
      bed_out = bind_cols(bed_out, bed_out1) %>%
        mutate(name = ".", score = ".", strand1 =".", strand2 = ".")
      colnames(bed_out)
      #Binding tads for color assignment
      bound_tads = bind_rows(called_tads)
      #Create vector of colors
      colors = c("0,0,0", "255,0,0", "0,255,0", "0,0,255")
      #Assign colors
      bed_out = bed_out %>%
        mutate(color =colors[bound_tads$Level])

      # bed_out = bed_out %>% mutate(start = format(start, scientific = FALSE), end = format(end, scientific = FALSE),
      #                              start1 = format(start1, scientific = FALSE), end1 = format(end1, scientific = FALSE))
      write.table(bed_out, out_path, quote = FALSE,
                  row.names = FALSE, sep = "\t", col.names = FALSE)
    } else if (out_format %in% c("bed", "hicexplorer")) {
      bed_out = bind_rows(called_tads) %>%
        dplyr::select(chr,start,end)
      # mutate(start = format(start, scientific = FALSE),
      #        end = format(end, scientific = FALSE))
      write.table(bed_out, out_path, quote = FALSE,
                  row.names = FALSE, sep = "\t", col.names = FALSE)
    } else {
      warning("No file output, unsupported output format chosen")
    }

  }

  if (grange == TRUE) {
    called_tads = lapply(called_tads, function(x) {
      GenomicRanges::GRanges(x)
    })
    called_tads = GenomicRanges::GRangesList(called_tads)
  }

  return(called_tads)
}




#Function to perform the actual sliding window Spectral clustering
#Used within SpectralTAD

.windowedSpec = function(cont_mat, resolution, chr,
                         gap_filter = TRUE,z_clust = FALSE,  qual_filter = TRUE,
                         eigenvalues = 2, min_size = 5,
                         window_size = ceiling(2000000/resolution),
                         gap_threshold = 1

) {

  #Set window sized based on biologically maximum TAD size of 2000000
  window_size = ceiling(window_size)

  #Find all regions which aren't completely zero and remove those that are


  #Get end point of the first window

  Group_over = dplyr::bind_rows()

  #Initialize first window

  start = 1
  end = window_size

  #Set parameter for determining end of loop

  end_loop = 0

  #Test if start+window is larger than the contact matrix and correct end point

  if (end+window_size>nrow(cont_mat)) {
    end = nrow(cont_mat)
  }

  #Begin sliding window clustering

  while (end_loop == 0) {

    #Subset matrix based on window size

    sub_filt = cont_mat[seq(start,end, 1), seq(start,end, 1)]
    #Remove columns and rows with % zeros higher than threshold
    zero_thresh = round(nrow(sub_filt)*(gap_threshold))
    non_gaps_within = which((colSums(sub_filt == 0))<zero_thresh)

    #Subset based on gaps
    sub_filt = sub_filt[non_gaps_within, non_gaps_within]

    #If matrix is empty then move window
    if (length(nrow(sub_filt)) == 0) {
      start = end
      end = start+window_size

      #If the new end is same as start end while loop
      if (start == nrow(cont_mat)) {
        end_loop = 1
        next
      }

      #If window overlaps with end of matrix make it move to last column
      if ( (end + (2000000/resolution)) > nrow(cont_mat) ) {
        end = nrow(cont_mat)

      }
      next
    }

    #Ignore if sub matrix if too small
    if (nrow(sub_filt) < min_size*2) {
      start = end
      end = start+window_size

      #If we reach the end of the matrix then end
      if (start == nrow(cont_mat)) {
        end_loop = 1
        next
      }

      if ( (end + (2000000/resolution)) > nrow(cont_mat) ) {
        end = nrow(cont_mat)

      }
      next
    }

    # sub_gaps = colSums(sub_filt)>0
    # sub_filt = sub_filt[sub_gaps, sub_gaps]

    #Calculate distance matrix for silhouette score

    dist_sub = 1/(1+sub_filt)

    #Get degree matrix

    dr = rowSums(abs(sub_filt))

    #Creating the normalized laplacian

    Dinvsqrt = diag((1/sqrt(dr)))

    P_Part1 = Matrix::crossprod(as.matrix(sub_filt), Dinvsqrt)
    sub_mat = Matrix::crossprod(Dinvsqrt, P_Part1)

    colnames(sub_mat) = colnames(cont_mat)[non_gaps_within]

    sub_mat[is.nan(sub_mat)] = 0

    #Get first k eigenvectors

    Eigen = PRIMME::eigs_sym(sub_mat, NEig = eigenvalues)

    eig_vals = Eigen$values
    eig_vecs = Eigen$vectors

    #Get order of eigenvalues from largest to smallest

    large_small = order(-eig_vals)

    eig_vals = eig_vals[large_small]
    eig_vecs = eig_vecs[,large_small]

    index = 1
    Group_mem = list()

    #Calculate the range of possible clusters

    clusters = seq_len(ceiling( (end-start+1)/min_size))

    #Normalize the eigenvectors from 0-1

    norm_ones = sqrt(dim(sub_mat)[2])

    for (i in seq_len(dim(eig_vecs)[2])) {
      eig_vecs[,i] = (eig_vecs[,i]/sqrt(sum(eig_vecs[,i]^2)))  * norm_ones
      if (eig_vecs[1,i] !=0) {
        eig_vecs[,i] = -1*eig_vecs[,i] * sign(eig_vecs[1,i])
      }
    }

    n = dim(eig_vecs)[1]
    k = dim(eig_vecs)[2]

    #Project eigenvectors onto a unit circle

    eig_vecs = crossprod(diag(diag(tcrossprod(eig_vecs))^(-1/2)), eig_vecs)

    #Get distance between points on circle

    point_dist = sqrt(rowSums( (eig_vecs-rbind(NA,eig_vecs[-nrow(eig_vecs),]))^2  ))

    #Use z-score to select significant gaps

    if (z_clust) {

      #Get statisticaly significant boundaries

      sig_bounds = which(scale(point_dist[-length(point_dist)])>2)

      #Remove boundaries within the minimum size

      sig_bounds = subset(sig_bounds, sig_bounds>min_size)

      #2*min_size is to offset and remove the second occurence

      dist_bounds = which(c(min_size*2,diff(sig_bounds))<min_size)

      #Remove bounds within the mininum size if they exist

      if (length(dist_bounds) > 0) {
        sig_bounds = sig_bounds[-dist_bounds]
      }

      #Create TADs using significant boundaries

      TAD_start = c(1, sig_bounds+1)

      TAD_end = c(sig_bounds, nrow(sub_filt))

      widths = (TAD_end-TAD_start)+1

      memberships = unlist(lapply(seq_len(length(TAD_start)), function(x) rep(x,widths[x])))

      #Create groups

      if (length(sig_bounds) == 0) {

        #Create empty set if non-significant

        end_group = dplyr::bind_rows()
      } else {

        sig_bounds = which(scale(point_dist[-length(point_dist)])>2)

        #Remove boundaries within the minimum size

        sig_bounds = subset(sig_bounds, sig_bounds>min_size)

        #2*min_size is to offset and remove the second occurence

        dist_bounds = which(c(min_size*2,diff(sig_bounds))<min_size)

        #Assign IDs based on coordinate and groups based on significant boundaries

        end_group = data.frame(ID = as.numeric(colnames(sub_filt)), Group = memberships)

        #Compile into bed file

        end_group = end_group %>% dplyr::mutate(group_place = Group) %>% dplyr::group_by(group_place) %>% dplyr::mutate(Group = last(ID)) %>% dplyr::ungroup() %>% dplyr::select(ID, Group)

      }


    } else {


      #Find largest gaps

      gap_order = order(-point_dist)

      #Remove boundaries occuring before minimum size at the very beginning of window

      #gap_order = gap_order[-which(gap_order<min_size)]

      #Initialize silhouette score

      sil_score = c()


      for (cluster in clusters) {

        #Loop through first k gaps and remove repeating boundaries

        #Set intial cutpoints to the number of clusters

        k = 1
        partition_found = 0
        first_run = TRUE
        cutpoints = c()

        #Loop through cluster numbers by iteratively adding new candidate boundaries and testing

        while(partition_found == 0) {

          #Get candidate gaps

          new_gap = gap_order[k]

          cutpoints = c(cutpoints, new_gap)

          #Identify gaps which are closer together than the minimum TAD size

          diff_points = which( abs(new_gap-cutpoints[-length(cutpoints)]) <= min_size)

          #If a point exists that is too close to another, remove it

          if (length(diff_points)>0) {
            cutpoints = cutpoints[-length(cutpoints)]
          }

          #If not these are final clusters

          if (length(cutpoints) == cluster) {
            partition_found = 1
          } else {
            k = k+1
          }
        }

        #If the new candidate cluster is an NA value, ignore

        if (any(is.na(cutpoints))) {
          next
        }

        #Order

        cutpoints = cutpoints[order(cutpoints)]

        #Combine cutpoints with start and end of window

        cutpoints = c(1, cutpoints, length(non_gaps_within)+1)

        #Find size of each cluster (TAD)

        group_size = diff(cutpoints)

        #Assign locations of the window memberships based on cutpoints

        memberships = c()
        for (i in seq_len(length(group_size))) {
          memberships = c(memberships, rep(i,times = group_size[i]))
        }

        #Get silhouette score for current number of clusters (TADs)

        sil = summary(cluster::silhouette(memberships,dist_sub))

        #Save silhouette scores for each configuration in vector

        sil_score = c(sil_score, sil$si.summary[4])

        #Save memberships in list

        Group_mem[[cluster]] = memberships

      }



      #Pull out the cutpoints which maximize silhouette score

      end_group = Group_mem[[which(diff(sil_score)<0)[1]]]

      #Put coordinates and group IDs into data frame

      if (length(end_group) == 0) {
        end_group = dplyr::bind_rows()
      } else {
        end_group = data.frame(ID = as.numeric(colnames(sub_filt)), Group = end_group)

        #Convert IDs to coordinates of endpoint to avoid overlaps

        end_group = end_group %>%dplyr::mutate(group_place = Group) %>%dplyr::group_by(group_place) %>%dplyr::mutate(Group = max(ID)) %>% ungroup() %>% dplyr::select(ID, Group)
      }
    }

    #End while loop if window reaches end of contact matrix

    if (end == nrow(cont_mat)) {
      Group_over = dplyr::bind_rows(Group_over, end_group)
      end_loop = 1
    } else {

      #Remove the last group (To account for overlap across windows) and set new start to start of group

      if (nrow(end_group)!=0) {
        end_IDs = which(end_group$Group == last(end_group$Group))
      } else {
        end_IDs = 1:window_size
      }

      start = end-length(end_IDs)+1

      #Account for cases when final TAD can't be removed

      if (length(start) == 0 ) {
        start = end
      }

      #Set new window end

      if (nrow(end_group != 0)) {
        end = start+window_size
      } else {
        end=start+window_size*2
      }
      #Remove final group to avoid repeating

      end_group = end_group[-end_IDs, ]

      #Combine TAD coordinates into single bed file
      Group_over = dplyr::bind_rows(Group_over, end_group)

      #Set end point to end of contact matrix if window is larger than end of matrix

      if ( (end + (2000000/resolution)) > nrow(cont_mat) ) {
        end = nrow(cont_mat)

      }
    }
  }


  #Organize final results based on options selected

  if (z_clust) {

    if (nrow(Group_over) > 0) {
      bed = Group_over %>% dplyr::group_by(Group) %>% dplyr::summarise(start = min(ID), end = max(ID) + resolution) %>%dplyr::mutate(chr = chr) %>% dplyr::select(chr, start, end) %>%
        dplyr::filter((end-start)/resolution >= min_size) %>%dplyr::arrange(start)
    } else {
      bed = Group_over
    }
  } else {

    if (qual_filter) {

      #Calculate an overall distance matrix for calculating silhouette score for filtering

      #Get range of values in the contact matrix

      fin_range = match(Group_over$ID,colnames(cont_mat))

      over_dist_mat = 1/(1+cont_mat[fin_range, fin_range])


      #Calculate group-wise silhouette

      sil = cluster::silhouette(Group_over$Group, over_dist_mat)

      ave_sil = summary(sil)$clus.avg.widths

      #Subset results based on silhouette score depending on qual_filter option

      bed = Group_over %>% dplyr::group_by(Group) %>% dplyr::summarise(start = min(ID), end = max(ID) + resolution) %>% dplyr::mutate(chr = chr) %>% dplyr::select(chr, start, end) %>%
        dplyr::mutate(Sil_Score = ave_sil) %>% dplyr::filter( ((end-start)/resolution >= min_size) & Sil_Score > .15)  %>%dplyr::arrange(start)
    } else {
      bed = Group_over %>% dplyr::group_by(Group) %>% dplyr::summarise(start = min(ID), end = max(ID) + resolution) %>% dplyr::mutate(chr = chr) %>% dplyr::select(chr, start, end) %>%dplyr::filter((end-start)/resolution >= min_size) %>% dplyr::arrange(start)
    }
  }

  return(bed)
}

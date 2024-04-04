#' Haar-Fisz variance-stabilizing transformation for InteractionSet
#'
#' This function applies Haar-Fisz variance-stabilizing transformation
#' to Hi-C contact counts from two replicates (provided as an
#' InteractionSet object)
#'
#' @param IntSet An InteractionSet object including two replicates of
#' Hi-C contact counts.
#'
#' @param chromoR If TRUE, it does denoising of coefficients as in chromoR.
#'
#' @return A variance-stabilized InteractionSet object.
#'
#' @import haarfisz
#' @import InteractionSet
#'
#' @export

# hft.transform <- function(IntSet){
#
#   cts = data.frame(assay(IntSet))
#   orig_length = dim(cts)[1]
#   padding_length = 2 ** ceiling(log2(orig_length)) - orig_length
#   hft_cts = matrix(nrow = orig_length, ncol = 2)
#   for (rep_num in c(1,2)){
#     rep = cts[,rep_num]
#     rep_order = order(rep)
#     rep = rep[rep_order]
#     rep = c(rep(0, padding_length), rep)
#     hft_rep = hft(rep)
#     hft_rep = hft_rep[(padding_length + 1) :
#                         (padding_length + orig_length)]
#     hft_rep = hft_rep[order(rep_order)]
#     hft_cts[, rep_num] = hft_rep
#   }
#   transformed_IntSet = InteractionSet(assays=list(counts = hft_cts),
#                                       interactions(IntSet))
#   return (transformed_IntSet)
# }

hft.transform <- function(IntSet, chromoR = FALSE){
  rg = regions(interactions(IntSet))
  transformed_isets = list()
  for (rep_num in c(1, 2)){
    m = inflate(IntSet, rg, rg, assay = 1, sample = rep_num)
    #print(class(m))
    #show(m)
    m = as.matrix(m)
    na_idx = is.na(m)
    m[na_idx] = 0
    y = list(diag = diag(m), nondiag = m[upper.tri(m)])
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
    transformed_m = matrix(data = NA, nrow = nrow(m), ncol = ncol(m))
    diag(transformed_m) = transformed_y$diag
    transformed_m[upper.tri(transformed_m)] = transformed_y$nondiag
    transformed_m[na_idx] = NA
    transformed_CM = ContactMatrix(transformed_m, rg, rg)
    transformed_isets[[rep_num]] = deflate(transformed_CM)
  }
  transformed_iset = do.call(cbind, transformed_isets)
  transformed_iset = sort(transformed_iset)
  return(transformed_iset)
}



#' Learning mean-dispersion relationship from
#' two Hi-C replicates.
#'
#' This function learns mean-dispersion function from DESeq2
#' given Hi-C contact counts from two replicates (provided as an
#' InteractionSet object)
#'
#' @param IntSet An InteractionSet object including two replicates of
#' Hi-C contact counts.
#'
#' @return A dispersionFunction object.
#'
#' @import DESeq2
#'
#' @export

vst.train <- function(IntSet){
  cts = assay(IntSet)
  colnames(cts) = c('rep1', 'rep2')
  dispF = get_dispersionFunction(cts)
  return(dispF)
  cell_types = c('A','A')
  coldata = data.frame(cell_type = cell_types,
                       row.names = c('rep1','rep2'))
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ 1)
  dds <- estimateSizeFactors(dds)
  suppressMessages({dispersionFunction(dds) <- dispF})
  dds_vst = varianceStabilizingTransformation(dds, blind=FALSE)
  transformed_IntSet = InteractionSet(assays=list(counts = assay(dds_vst)),
                                      interactions(IntSet))
  return(transformed_IntSet)
}

get_dispersionFunction <- function(cts,
                                   nsub=10000,
                                   fitType='local'){
  cell_types = c('A','A')
  coldata = data.frame(cell_type = cell_types,
                       row.names = c('rep1','rep2'))
  dds_object <- DESeqDataSetFromMatrix(countData = cts,
                                       colData = coldata,
                                       design = ~ 1)
  dds_object <- estimateSizeFactors(dds_object)
  baseMean <- rowMeans( counts(dds_object, normalized=TRUE) )
  dds.sub <- dds_object[baseMean > 1,]
  baseMean <- baseMean[baseMean > 1]
  o <- order(baseMean)
  idx <- o[round( seq(from = 1, to = length(o), length = nsub) )]
  dds.sub <- dds.sub[idx,]
  dds.sub <- estimateDispersionsGeneEst(dds.sub, quiet=TRUE)
  dds.sub <- estimateDispersionsFit(dds.sub, fitType=fitType, quiet=TRUE)
  return(dispersionFunction(dds.sub))
}

#' A variance-stabilizing transformation for InteractionSet based on 'vst'
#' function from 'DESeq2'
#'
#' This function takes a learned dispersionFunction and stabilize the variance of
#' Hi-C contact counts from two replicates (provided as an
#' InteractionSet object)
#'
#' @param IntSet An InteractionSet object including two replicates of
#' Hi-C contact counts.
#'
#' @param dispF A dispersionFunction object.
#'
#' @return variance-stabilized InteractionSet object.
#'
#' @import DESeq2
#'
#' @export

vst.transform <- function(IntSet, dispF){
  cts = assay(IntSet)
  colnames(cts) = c('rep1', 'rep2')
  cell_types = c('A','A')
  coldata = data.frame(cell_type = cell_types,
                       row.names = c('rep1','rep2'))
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ 1)
  dds <- estimateSizeFactors(dds)
  suppressMessages({dispersionFunction(dds) <- dispF})
  dds_vst = varianceStabilizingTransformation(dds, blind=FALSE)
  transformed_IntSet = InteractionSet(assays=list(counts = assay(dds_vst)),
                                      interactions(IntSet))
  return(transformed_IntSet)
}

#' Log transformation for InteractionSet
#'
#' This function applies log transformation
#' to Hi-C contact counts from two replicates (provided as an
#' InteractionSet object)
#'
#' @param IntSet An InteractionSet object including two replicates of
#' Hi-C contact counts.
#'
#' @return A variance-stabilized InteractionSet object.
#'
#' @export

log.transform <- function(IntSet){
  cts = assay(IntSet)
  log_cts = log2(cts + 1)
  transformed_IntSet = InteractionSet(assays=list(counts = log_cts),
                                      interactions(IntSet))
  return(transformed_IntSet)
}

#' Asinh transformation for InteractionSet
#'
#' This function applies asinh transformation
#' to Hi-C contact counts from two replicates (provided as an
#' InteractionSet object)
#'
#' @param IntSet An InteractionSet object including two replicates of
#' Hi-C contact counts.
#'
#' @return A variance-stabilized InteractionSet object.
#'
#' @export

asinh.transform <- function(IntSet){
  cts = assay(IntSet)
  asinh_cts = asinh(cts)
  transformed_IntSet = InteractionSet(assays=list(counts = asinh_cts),
                                      interactions(IntSet))
  return(transformed_IntSet)
}

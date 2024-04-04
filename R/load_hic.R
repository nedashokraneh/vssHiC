# adapted from https://github.com/dovetail-genomics/coolR/blob/master/R/coolR.R
# to load InteractionSet corresponding to only one pair of chromosomes and
# inter-chromosomal InteractionSet

#' Read .cool/.mcool sparse matrix into a InteractionSet object
#'
#' This function reads HiC contact matrix file(s) created by the
#' cooler application (https://github.com/mirnylab/cooler) stored in a
#' HDF5 data storage and converts them to an InteractionSet object for genomic bins
#' within a pair of chromosomes with the counts of reads in these bin pairs.
#'
#' @param files The list of paths to a HDF5 stored cools (uni-dimension)
#' or mcools (multi-dimension sparse matrix). If using a multi-dimesion,
#' the resolution of one of the dimension need to be passed to res.
#' @param chr1 The chromosome name including genomic bins of the first anchor.
#' @param chr2 The chromosome name including genomic bins of the second anchor.
#' @param res 'NULL' if using a uni-dimensional cool file or the resolution of one of layer in the mcool file
#'
#' @return An InteractionSet of relevant bin pairs with the number of interaction for each pair
#'
#' @import rhdf5
#' @import coolR
#' @importFrom coolR getBins
#'
#' @export


cool2IntSet <- function(files, chr1, chr2, res=NULL) {

  GIs = lapply(files, function(file){
    read.cool(file,
              res = res,
              chr1 = chr1,
              chr2 = chr2)
    })
  IntSet = gi2is(GIs, paste0("IF", seq(length(files))))
  IntSet = sort(IntSet)
  return(IntSet)
}

read.cool <- function(file,
                      res=NULL,
                      chr1=NULL,start1=NULL,end1=NULL,
                      chr2=NULL,start2=NULL,end2=NULL) {

  anchors <- getBins(file,res)

  slice <- getSlice(anchors,file,res,chr1,start1,end1,chr2,start2,end2)

  g.i <- GInteractions(anchors[slice$bin1_id],
                       anchors[slice$bin2_id],
                       count=slice$count)

  return(g.i)
}


gi2is <- function(gi.counts,col.names) {
  gi <- unique(do.call(c,gi.counts))
  S4Vectors::mcols(gi) <- c()

  ## Get a dataframe of counts for each interaction pairs
  counts <- do.call(cbind,lapply(gi.counts,function(x){
    ## Initialize a vector of counts 0 of the lenght of all interaction
    counts <- as.integer(rep(0,length(gi)))
    ## Find the interactions in the common set intersecting with the current GI
    mm <- findMatches(gi,x)
    ## Replace in the vector of counts for the universe, replace the intersection with the current counts
    counts[queryHits(mm)] <- x$count[subjectHits(mm)]
    ## Return
    return(counts)
  }))

  ## Asign the name of the Column to the sample
  colnames(counts) <- col.names

  ## Get the library size
  lib.data <- DataFrame(totals=colSums(counts))

  ## Create the interaction set
  data <- InteractionSet(assays=list(counts=counts),
                         gi,
                         colData=lib.data)

  return(data)
}

getSlice <- function(anchors,file,res,chr1,start1,end1,chr2,start2,end2){
  ################################################################################
  ## Get the chromosme bins as a GRanges
  ################################################################################
  if(is.null(start1)) start1 <- 1
  if(is.null(end1)) end1 <- seqlengths(anchors)[chr1]

  indexes.group <- ifelse(is.null(res),"/indexes",paste('resolutions',res,'indexes',sep='/'))
  chr.group <- ifelse(is.null(res),"/chroms",paste('resolutions',res,'chroms',sep='/'))
  pixels.group <- ifelse(is.null(res),"/pixels",paste('resolutions',res,'pixels',sep='/'))

  indexes <- list(chr=as.vector(h5read(file,name=paste(chr.group,'name',sep="/"))),
                  chr_idx=as.vector(h5read(file,name=paste(indexes.group,'chrom_offset',sep="/"))),
                  bin1_idx=as.vector(h5read(file,name=paste(indexes.group,'bin1_offset',sep="/")))
  )

  ################################################################################
  ## Reading chromosome chunk If chr1 is null, return the full cool file
  ################################################################################
  if (is.null(chr1)){
    chunk <- NULL
  } else {
    chr1.start <- GRanges(chr1,IRanges(start1,width=1))
    chr1.end <- GRanges(chr1,IRanges(end1,width=1))

    chr1.start.idx <- subjectHits(findOverlaps(chr1.start,anchors))
    chr1.end.idx <- subjectHits(findOverlaps(chr1.end,anchors))

    idx.chunk <- seq(chr1.start.idx,chr1.end.idx)

    bin1.idx <- as.vector(h5read(file,
                                 name=paste(indexes.group,'bin1_offset',sep="/"),
                                 index=list(idx.chunk)))

    slice <- sum(bin1.idx[-1] - bin1.idx[-length(bin1.idx)])-1

    chunk  <- seq(bin1.idx[1]+1,bin1.idx[1]+1+slice)

  }

  ################################################################################
  ## Reading the chunks from the cool file
  ################################################################################
  d.f <- data.frame(bin1_id = as.vector(h5read(file,
                                               name=paste(pixels.group,'bin1_id',sep="/"),
                                               index=list(chunk)))+1,
                    bin2_id = as.vector(h5read(file,
                                               name=paste(pixels.group,'bin2_id',sep="/"),
                                               index=list(chunk)))+1,
                    count = as.vector(h5read(file,
                                             name=paste(pixels.group,'count',sep="/"),
                                             index=list(chunk)))
  )

  ################################################################################
  ## If Chr2 is set, only return corresponding ranges
  ################################################################################
  if (!is.null(chr2)){
    if(is.null(start2)) start2 <- 1
    if(is.null(end2)) end2 <- seqlengths(anchors)[chr2]

    chr2.start <- GRanges(chr2,IRanges(start2,width=1))
    chr2.end <- GRanges(chr2,IRanges(end2,width=1))

    chr2.start.idx <- subjectHits(findOverlaps(chr2.start,anchors))
    chr2.end.idx <- subjectHits(findOverlaps(chr2.end,anchors))

    filter.bin2 <- d.f$bin2_id %in% chr2.start.idx:chr2.end.idx

    d.f <- d.f[filter.bin2,]
  }

  return(d.f)

}


# scipen is set to a large value to read resolutions as a whole number instead of
# scientific notation
.onLoad <- function(...) {
  options(scipen = 9)
}

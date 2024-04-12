
subset_GI <- function(GI, res, start, end){
  seqname = unfactor(seqnames(anchors(GI)$first)@values)
  query = GRanges(seqnames = seqname,
                  IRanges(seq(start + 1, end, res),
                          seq(start + res, end, res)))
  subset_idx = (anchors(GI)$first %within% query) &
    (anchors(GI)$second %within% query)
  subset_GI = GI[subset_idx, ]
  return(subset_GI)
}

iset2GI <- function(iset, rep_num){
  GI = interactions(iset)
  GI$count = assay(iset)[, rep_num]
  return(GI)
}

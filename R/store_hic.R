#' Storing IntSet assays
#'
#' This function writes interaction frequencies of assay from IntSet into 7-columns
#' tab-sep format (chr1, start1, end1, chr2, start2, end2, count)
#'
#' @param IntSet The InteractionSet object to store.
#' @param rep_num Number of the assay to store.
#' @param filepath The path to write a tab-sep file.
#'
#' @import data.table
#'
#' @export
#'

store_IntSet <- function(IntSet, rep_num, filepath){
  pair_info_cols = c("seqnames1", "start1", "end1", "seqnames2", "start2", "end2")
  interactions_df = interactions(IntSet) %>% data.frame()
  interactions_df = interactions_df[, pair_info_cols]
  interactions_df[, c("start1", "start2")] = interactions_df[, c("start1", "start2")] - 1
  interactions_df$count = assay(IntSet)[, rep_num]
  interactions_df = interactions_df[interactions_df$count != 0, ]
  data.table::fwrite(interactions_df, filepath, sep = "\t",
                     row.names = F, col.names = F)
}


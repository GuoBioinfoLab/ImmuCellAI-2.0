#' Merge reference profiles
#'
#' @param reference Gene x state matrix.
#' @param merge.groups Named list. Names are new merged states; values are original states.
#' @return Merged reference matrix.
#' @export
merge_reference <- function(reference, merge.groups) {
  reference <- as.matrix(reference)
  reference[is.na(reference)] <- 0
  reference[reference < 0] <- 0
  original.cells <- colnames(reference)
  used.cells <- unique(unlist(merge.groups))
  used.cells <- intersect(used.cells, original.cells)
  unmerged.cells <- setdiff(original.cells, used.cells)

  ref.list <- list()
  for (ct in unmerged.cells) ref.list[[ct]] <- reference[, ct]

  for (new.name in names(merge.groups)) {
    cells <- intersect(merge.groups[[new.name]], original.cells)
    if (length(cells) == 0) next
    ref.list[[new.name]] <- if (length(cells) == 1) reference[, cells] else rowMeans(reference[, cells, drop = FALSE])
  }
  ref.merged <- do.call(cbind, ref.list)
  ref.merged <- as.matrix(ref.merged)
  ref.merged[is.na(ref.merged)] <- 0
  ref.merged[ref.merged < 0] <- 0
  ref.merged
}

#' Build merge map
#'
#' @param original.cells Original cell-type names.
#' @param merge.groups Named merge groups.
#' @return Named character vector mapping original cells to merged cells.
#' @export
build_merge_map <- function(original.cells, merge.groups) {
  merge.map <- stats::setNames(original.cells, original.cells)
  if (length(merge.groups) == 0) return(merge.map)
  for (new.name in names(merge.groups)) {
    cells <- intersect(merge.groups[[new.name]], original.cells)
    merge.map[cells] <- new.name
  }
  merge.map
}

#' @title Circ.filter
#'
#' @description
#' Filter circRNA candidates based on expression level. User specify how many samples (filter.sample) need to have above
#' which minimum number (filter.count) of read count support. In addition, circRNA need to account for which percentage (percentage) of the total transcripts in at least one group.
#'
#' @param circ CircRNACount file. A file of circRNA read count table. First three columns are circRNA coordinates, and followed by columns for circRNA read counts, each sample per column.
#' @param linear LinearCount file. A file of circRNA host gene expression count table. Same configuration as CircRNACount file.
#' @param Nreplicates Number of replicates in your data. Expect each group have the same number of replicates.
#' @param filter.count The minimum read count used for filtering.
#' @param filter.sample The minimum number of samples need to have above filter.count number of circRNA supporting reads.
#' @param percentage The minimum percentage of circRNAs account for the total transcripts in at least one group.
#' @param circle_description Column indices which do not carry circle/linear read counts.
#' @export Circ.filter
#'

Circ.filter <- function(circ = circ, linear = linear, Nreplicates = 3, filter.sample = 4, filter.count = 5, percentage = 1, circle_description = c(1:3)) {

  del_row <- c()

  for (i in 1:nrow(circ)) {

    if (!is.null(Nreplicates)) {

      if (sum(circ[i, -circle_description] >= filter.count) < filter.sample | circ_max_perc(circ[i, -circle_description], linear[i, -circle_description], Nreplicates = Nreplicates) < percentage)
        del_row <- c(del_row, i)

    } else {
      if (sum(circ[i, -circle_description] >= filter.count) < filter.sample)
        del_row <- c(del_row, i)
    }


  }

  if (length(del_row) > 0) {
    new_dat = circ[-del_row,]
    return(new_dat)
  } else {
    return(circ)
  }
}

circ_max_perc <- function(circ = circ, linear = linear, Nreplicates = 3) {
  # convert to vector
  circ <- as.numeric(circ)
  linear <- as.numeric(linear)
  if (length(circ) != length(linear)) {
    stop('Number of samples in circRNA is not equal to Hostgene.')
  }
  Ngroups <- length(circ) / Nreplicates
  # calculate percentage
  circ_sum <- unname(tapply(circ, (seq_along(1:length(circ)) - 1) %/% Nreplicates, sum))
  linear_sum <- unname(tapply(linear, (seq_along(1:length(linear)) - 1) %/% Nreplicates, sum))
  perc <- max(circ_sum / (circ_sum + linear_sum), na.rm = T)
  return(perc)
}

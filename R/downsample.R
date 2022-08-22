#' Down sample a vector of type \code{numeric}; nonsensical for negative values.
#'
#' @param x A \code{numeric} vector.
#'
#' @param thres If \code{thres} >= 3,000, then in
#' all input values > \code{thres} are
#' downsampled.
#' If \code{thres} < 3,000,
#'  some input values > \code{thres} are
#'  downsampled.
#'
#' @return A vector of integers (type \code{numeric}) of the same
#'   length as \code{x}, downsampled as described in the
#'   documentation for the \code{thres} argument.
#'
#' @export
#'
downsample <- function(x, thres = 3000) {
  return(
    ifelse(x <= thres,
           x,
           ceiling(pmin(x, thres + 3000 * log10(x/thres)))))
}

#' Down sample a set of mutational spectra
#'
#' @param spec Input spectra as a numerical matrix or similar \code{data.frame};
#'   each column is a spectrum, each row is a mutation type (e.g. CAG -> CTG).
#'
#' @param thres See \code{\link{downsample}}.
#'
#' @return A numeric matrix with same shape as \code{spec},
#'   with the entries each column reduced based on
#'   the ratio of
#'   \code{\link{downsample}(colSums(spec))} to
#'   \code{colSums(spec)}.
#'
#' @export
#'
downsample_spectra <- function(spec, thres) {

  mut_sum <- colSums(spec)
  down_sum <- downsample(mut_sum, thres = thres)
  down_factor <- down_sum / mut_sum

  down_spec <- round(t(t(spec) * down_factor))

  # Return down-sampled spectra
  return(list(down_spec = down_spec, down_factor = down_factor))
}



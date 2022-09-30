#' Downsample a vector of type \code{numeric}; nonsensical for negative values.
#'
#' @param x A \code{numeric} vector.
#'
#' @param downsample_threshold
#' If \code{downsample_threshold} >= 3,000, then
#' all input values > \code{downsample_threshold} are
#' downsampled.
#' If \code{downsample_threshold} < 3,000,
#'  only some input values > \code{downsample_threshold} are
#'  downsampled but all input values > 3000
#'  are downsampled. See \code{\link{show_downsample_curves}}.
#'
#' @param thres Deprecated alternative to \code{downsample_threshold}.
#'
#' @return A vector of integers (type \code{numeric}) of the same
#'   length as \code{x}, downsampled as described in the
#'   documentation for the \code{downsample_threshold} argument.
#'
#' @export
#'
downsample <- function(x, thres = NULL, downsample_threshold = 3000) {
  if (!is.null(thres)) {
    downsample_threshold = thres
    warning("Use downsample_threshold insead of thres in downwample()")
  }
  return(
    ifelse(
      x <= downsample_threshold,
      x,
      ceiling(pmin(x, downsample_threshold + 3000 * log10(x/downsample_threshold)))))
}

#' Downsample a set of mutational spectra
#'
#' @param spec Input spectra as a numerical matrix or similar \code{data.frame};
#'   each column is a spectrum, each row is a mutation type (e.g. CAG -> CTG).
#'
#' @param thres Deprecated alternative to \code{downsample_threshold}.
#'
#' @param downsample_threshold See \code{\link{downsample}}
#'   and \code{\link{show_downsample_curves}}.
#'
#' @md
#'
#' @return A list with elements:
#'
#'  * `down_spec`: A numeric matrix with same shape as \code{spec},
#'   with the mutation counts in each column reduced based on
#'   the corresponding ratio in the second element of this list, `down_factor`.
#'
#'  * `down_factor` Numeric vector of the ratios of
#'   \code{\link{downsample}(colSums(spec))} to
#'   \code{colSums(spec)}.
#'
#' @export
#'
downsample_spectra <- function(spec, downsample_threshold = NULL, thres = NULL) {
  if (is.null(downsample_threshold))  {
    downsample_threshold = thres
    warning("Use downsample_threshold insead of thres in downwample_spectra()")
  }
  mut_sum <- colSums(spec)
  down_sum <- downsample(mut_sum, downsample_threshold = downsample_threshold)
  down_factor <- down_sum / mut_sum

  down_spec <- round(t(t(spec) * down_factor))

  # Return down-sampled spectra
  return(list(down_spec = down_spec, down_factor = down_factor))
}

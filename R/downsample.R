#' Down sample a vector of integers; nonsensical for negative values
#'
#' @param x Vector of integers
#'
#' @param thres Values $le$ \code{thres} are unmodified; for values of
#'   \code{thres} > 3,000 some of the return values will also not be
#'   reduced.
#'
#' @return A vector of integers (type \code{numeric}) of the same
#'   length as \code{x}, with elements that were $le$ \code{thres} in
#'   \code{x} reduced.
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
#'   with each column down-sampled by
#'   \code{\link{downsample}} based on its \code{colSums}.
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



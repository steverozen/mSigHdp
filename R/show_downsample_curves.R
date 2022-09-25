#' Show the effects of the downsample_threshold argument to \code{\link{RunHdpxParallel}}.
#'
#' @param ... One or more postive numbers to use as
#'            possible values for \code{downsample_threshold}.
#'
#' @export
#'
#' @examples show_downsample_curves(3000, 5000, 10000)
#'
#' @value Called for side effects.
#'
#' TODO -- update parameter names in downsample, downsample_spectra

show_downsample_curves <- function(...) {
  thresholds <- unlist(list(...))
  last <- length(thresholds)
  stopifnot(last >= 1)
  maxx <- max(thresholds) * 2
  xvec <- round(seq(0, maxx, by = maxx / 1000))
  pch = "."
  thresholds <- sort(thresholds)
  plot(xvec, downsample(xvec, thres = thresholds[last]),
       pch = pch,
       xlab = "Original number of mutations",
       ylab = "Number of mutatons after downsampling")
  if (length(thresholds) > 1) {
    thresholds <- thresholds[-last]
    for (tt in thresholds) {
      # browser()
      points(xvec, downsample(xvec, thres = tt), pch = pch, new = FALSE)
    }
  }
}


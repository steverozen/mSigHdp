#' Show the effects of the downsample_threshold argument to \code{\link{RunHdpxParallel}}.
#'
#' @param ... One or more postive numbers to use as
#'            possible values for \code{downsample_threshold}.
#'
#' @export
#'
#' @importFrom graphics lines text points
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
  maxt <- max(thresholds)
  maxx <- maxt * 1.5
  mint <- min(thresholds)
  xvec <- round(seq(0, maxx, by = maxx / 1000))
  pch <- "."
  lty <- 3
  start.line.fac <- 0.5
  thresholds <- sort(thresholds)
  plot(xvec, downsample(xvec, thres = maxt),
       pch = pch,
       xlab = "Original number of mutations",
       ylab = "Number of mutations after downsampling")
  lines(x = c(mint * start.line.fac, maxx), y = rep(maxt, 2), type="l", lty = lty)
  text(x = 0, y = maxt, maxt, adj = c(0, 0.5))
  text(x = 0, y = maxt * 1.07, "downsampling_threshold", adj = c(0, 0.5))
  # abline(h = thresholds[last], lty = lty)
  if (length(thresholds) > 1) {
    thresholds <- thresholds[-last]
    for (tt in thresholds) {
      # browser()
      points(xvec, downsample(xvec, thres = tt), pch = pch, new = FALSE)
      lines(x = c(mint * start.line.fac, maxx), y = rep(tt, 2), type="l", lty = lty)
      text(x = 0, y = tt, tt, adj = c(0, 0.5))
      # abline(h = tt, lty = lty)
    }
  }
}


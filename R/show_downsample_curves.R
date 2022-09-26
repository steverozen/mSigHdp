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
#' @return Called for side effects.

show_downsample_curves <- function(...) {
  th <- unlist(list(...))
  last <- length(th)
  stopifnot(last >= 1)
  maxt <- max(th)
  maxx <- maxt * 1.5
  mint <- min(th)
  xvec <- round(seq(0, maxx, by = maxx / 1000))
  pch <- "."
  lty <- 3
  line.start.x <- maxx * 0.08
  th <- sort(th)
  plot(xvec, downsample(xvec, downsample_threshold = maxt),
       pch = pch,
       xlab = "Original number of mutations",
       ylab = "Number of mutations after downsampling",
       ylim = c(0, 1.1 * maxt))
  lines(x = c(line.start.x, maxx),
        y = rep(maxt, 2), type="l", lty = lty)
  text(x = 0, y = maxt, maxt, adj = c(0, 0.5))
  text(x = 0, y = maxt * 1.07, "downsampling_threshold", adj = c(0, 0.5))
  if (length(th) > 1) {
    th <- th[-last]
    for (tt in th) {
      points(xvec, downsample(xvec, downsample_threshold = tt), pch = pch, new = FALSE)
      lines(x = c(line.start.x, maxx),
            y = rep(tt, 2), type="l", lty = lty)
      text(x = 0, y = tt, tt, adj = c(0, 0.5))
    }
  }
}


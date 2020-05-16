#' Return a function to calculate the unsigned Stirling numbers of the first kind
#'
#' @export
#'
#' @return A function to calculate a vector of unsigned Stirling numbers,
#' \eqn{s(n ,k), k = 1...n},
#' each divided by the maximum Stirling number in the series.
#' The returned function is a closure with state that includes
#' a list of all the unsigned Stirling number series \eqn{<=} the argument, \eqn{n},
#'
#' i.e. \eqn{[s(1, 1)], [s(2, 1), s(2, 2)], ..., [s(n, 1), ..., s(n, n)]}.
#' Memory usage could be substantial, but the stored state
#' does not include the many trailing zeros in the vectors.
#' For this to work within the hdp (https://github.com/nicolaroberts/hdp)
#' package the
#' function returned *must* be called \code{stir.closure}.

xmake.s <- function() {

  allss <- list(1)

  stir <- function(nn) {
    # cat(".\n")

    if (nn == 0) return(1)

    len.all <- length(allss)

    if (nn > len.all) {
      for (mm in (len.all + 1):nn) { # For each new Stirling number series, if any

        ss <- allss[[mm - 1]]

        newss <- c(ss * (mm - 1), 0) + c(0, ss)
        newss <- newss / max(newss)

        last.non.0 <- max(which(newss > 0))
        newss <- newss[1:last.non.0]

        allss[[mm]] <<- newss
      }
    }

    retval <- allss[[nn]]
    retval <- c(retval, rep(0, nn - length(retval)))
    stopifnot(length(retval) == nn)
    attr(retval, "closure") <- TRUE
    return(retval)
  }

  return(stir)

}


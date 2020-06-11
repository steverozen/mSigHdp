#' Run hdp extraction and attribution on a spectra catalog file
#' A function to do burn-in iteration only. This returns a list of hdp object.
#' This needs to be converted to a hdpState object before hdp_posterior
#' (hdpx:::as.hdpState(hdplist))
#'
#'
#' @inheritParams PrepInit
#'
#' @param seedNumber An integer that is used to generate separate
#'   random seeds for each call to \code{\link[hdpx]{dp_activate}},
#'   and each call of \code{\link[hdpx]{hdp_posterior}}; please see the code
#'   on how this is done. But repeated calls with same value of
#'   \code{seedNumber} and other inputs should produce the same results.
#'
#' @param post.burnin Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{burnin}.
#'
#' @param post.cpiter Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{cpiter}.
#'
#' @param post.verbosity Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{verbosity}.
#'
#' @param hdplist If a hdplist is provided, SetupAndActivate will be skipped.
#'
#' @return A list with hdp object after burn-in iteration and likelihood of iteration
#'
#' @export
#'
BurninIteration <-
  function(input.catalog,
           seedNumber     = 1,
           K.guess,
           multi.types    = FALSE,
           verbose        = TRUE,
           post.burnin    = 4000,
           post.cpiter    = 3,
           post.verbosity = 0,
           gamma.alpha    = 1,
           gamma.beta     = 1,
           hdplist        = NULL

  ) { # 8 arguments

    if(is.null(hdplist)){
      hdp.state <- SetupAndActivate(input.catalog = input.catalog,
                                    seedNumber    = seedNumber,
                                    K.guess       = K.guess,
                                    multi.types   = multi.types,
                                    verbose       = verbose,
                                    gamma.alpha   = gamma.alpha,
                                    gamma.beta    = gamma.beta)

      hdplist <- hdpx::as.list(hdp.state)
    }
    iterate <- utils::getFromNamespace(x = "iterate", ns = "hdpx")
    output <- iterate(hdplist, post.burnin, post.cpiter, post.verbosity)##burn-in first, then return the hdplist after burnt in.
    return(invisible(list(hdplist    = output[[1]],
                          likelihood = output[[2]])))
  }

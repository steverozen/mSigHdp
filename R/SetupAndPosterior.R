#' Generate an HDP Gibbs sampling chain from a spectra catalog.
#'
#' @inheritParams PrepInit
#'
#' @inheritParams SetupAndActivate
#'
#' @param post.burnin Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{burnin}.
#'
#' @param post.n Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{n}.
#'
#' @param post.space Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{space}.
#'
#' @param post.cpiter Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{cpiter}.
#'
#' @param post.verbosity Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{verbosity}.
#'
#' @return Invisibly, an \code{\link[hdpx]{hdpSampleChain-class}} object
#'  as returned from \code{\link[hdpx]{hdp_posterior}}.
#'
#' @export

SetupAndPosterior <-
  function(input.catalog,
           seedNumber          = 1,
           K.guess,
           multi.types         = FALSE,
           verbose             = TRUE,
           post.burnin         = 4000,
           post.n              = 50,
           post.space          = 50,
           post.cpiter         = 3,
           post.verbosity      = 0,
           gamma.alpha         = 1,
           gamma.beta          = 1,
           one.parent.hack     = FALSE)
{ # 12 arguments

    hdp.state <- SetupAndActivate(input.catalog = input.catalog,
                                  seedNumber    = seedNumber,
                                  K.guess       = K.guess,
                                  multi.types   = multi.types,
                                  verbose       = verbose,
                                  gamma.alpha   = gamma.alpha,
                                  gamma.beta    = gamma.beta,
                                  one.parent.hack = one.parent.hack)

    if (verbose) message("calling hdp_posterior, seed = ",
                         seedNumber, " ", Sys.time())
    posterior.time <- system.time(
      sample.chain <- hdpx::hdp_posterior(
        hdp       = hdp.state,
        verbosity = post.verbosity,
        burnin    = post.burnin,
        n         = post.n,
        space     = post.space,
        cpiter    = post.cpiter,
        seed      = seedNumber)
    )

    if (verbose) {
      message("compute sample.chain time: ")
      for (xn in names(posterior.time)) {
        message(" ", xn, " ", posterior.time[[xn]])
      }
    }

    return(invisible(sample.chain))

  }

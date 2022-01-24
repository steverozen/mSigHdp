#' Generate an HDP Gibbs sampling chain from a spectra catalog
#'
#' @inheritParams PrepInit
#'
#' @inheritParams SetupAndActivate
#'
#' @inheritParams Burnin
#'
#' @param post.n The number of posterior samples to collect.
# Pass to \code{\link[hdpx]{hdp_posterior_sample}} \code{n}.
#'
#' @param post.space Pass to \code{\link[hdpx]{hdp_posterior_sample}}
#'      \code{space}. The number of iterations between collected samples.
#'
#' @param post.cpiter The number of iterations of concentration
#'        parameter samplings to perform after each iteration.
#        Pass to \code{\link[hdpx]{hdp_posterior_sample}} and
#         \code{\link[hdpx]{hdp_burnin}} \code{cpiter}.
#'
#' @param post.verbosity Verbosity of debugging statements.
#'       No need to change except for development purposes.
#        Pass to \code{\link[hdpx]{hdp_posterior_sample}}
#        \code{verbosity}.
#'
#' @param checkpoint If \code{TRUE}, then \itemize{
#'      \item Checkpoint each final Gibbs sample
#'        chain to the current working directory, in a file called
#'        mSigHdp.sample.checkpoint.*x*.Rdata, where
#'        *x* depends on \code{seedNumber}.
#'      \item Periodically checkpoint the burnin state
#'        to the current working directory, in files called
#'        mSigHdp.burnin.checkpoint.*x*.Rdata,
#'        where *x* depends on the \code{seedNumber}.
#'  }
#'
#' @return Invisibly, an \code{\link[hdpx]{hdpSampleChain-class}} object
#'  as returned from \code{\link[hdpx]{hdp_posterior}}.
#'
#' @keywords internal

SetupAndPosterior <-
  function(input.catalog,
           seedNumber          = 1,
           K.guess,
           multi.types         = FALSE,
           verbose             = TRUE,
           burnin              = 5000,
           post.n              = 50,
           post.space          = 50,
           post.cpiter         = 3,
           post.verbosity      = 0,
           gamma.alpha         = 1,
           gamma.beta          = 20,
           burnin.multiplier   = 2,
           checkpoint          = TRUE)
  {

    hdp.state <- SetupAndActivate(input.catalog = input.catalog,
                                  seedNumber    = seedNumber,
                                  K.guess       = K.guess,
                                  multi.types   = multi.types,
                                  verbose       = verbose,
                                  gamma.alpha   = gamma.alpha,
                                  gamma.beta    = gamma.beta)

    if (verbose) message("calling hdp_posterior, seed = ",
                         seedNumber, " ", Sys.time())

    burnin.output <- Burnin(
      hdp.state         = hdp.state,
      burnin.verbosity  = post.verbosity,
      burnin            = burnin,
      cpiter            = post.cpiter,
      seedNumber        = seedNumber,
      burnin.multiplier = burnin.multiplier,
      checkpoint        = checkpoint)

    posterior.time <- system.time(
      sample.chain <- hdpx::hdp_posterior_sample(post.input     = burnin.output,
                                                 post.n         = post.n,
                                                 post.space     = post.space,
                                                 post.cpiter    = post.cpiter,
                                                 seed           = seedNumber,
                                                 post.verbosity = post.verbosity,
                                                 checkpoint     = checkpoint)
    )

    if (checkpoint) {
      save(sample.chain,
           file = paste0("mSigHdp.sample.checkpoint.", seedNumber, ".Rdata"))
    }

    if (verbose) {
      message("compute sample.chain time: ")
      for (xn in names(posterior.time)) {
        message(" ", xn, " ", posterior.time[[xn]])
      }
    }

    return(invisible(sample.chain))

  }



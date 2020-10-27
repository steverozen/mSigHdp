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
#' @param checkpoint.1.chain If \code{TRUE} checkpoint the sample
#'      chain to current working directory, in a file called
#'      sample.chain.*seed_number*.Rdata.
#'
#' @return Invisibly, an \code{\link[hdpx]{hdpSampleChain-class}} object
#'  as returned from \code{\link[hdpx]{hdp_posterior}}.
#'
#' @export

SetupAndPosterior2 <-
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
           gamma0.alpha        = gamma.alpha,
           gamma0.beta         = gamma.beta,
           checkpoint.1.chain  = TRUE,
           burnin.multiplier   = 1,
           burnin.checkpoint   = NULL,
           prior.signatures    = NULL,
           prior.pseudocounts  = NULL)
{
    hdp.state <- SetupAndActivate(input.catalog = input.catalog,
                                  seedNumber    = seedNumber,
                                  K.guess       = K.guess,
                                  multi.types   = multi.types,
                                  verbose       = verbose,
                                  gamma.alpha   = gamma.alpha,
                                  gamma.beta    = gamma.beta,
                                  gamma0.alpha  = gamma0.alpha,
                                  gamma0.beta   = gamma0.beta)

    if (verbose) message("calling hdp_posterior, seed = ",
                         seedNumber, " ", Sys.time())

    burnin.output <- ActivateAndBurnin2( # Change to burnin
      hdp       = hdp.state,
      verbosity = post.verbosity,
      burnin    = post.burnin,
      n         = post.n,
      space     = post.space,
      cpiter    = post.cpiter,
      seed      = seedNumber)

    sample.chain <- hdp_posterior_sample(burnin.output, XXXXXXXXXXX) # hdpx::hdp_posterior no longer needed?

    if (checkpoint.1.chain) {
      save(sample.chain, file = paste0("sample.chain.", seedNumber, ".Rdata"))
    }

    if (verbose) {
      message("compute sample.chain time: ")
      for (xn in names(posterior.time)) {
        message(" ", xn, " ", posterior.time[[xn]])
      }
    }

    return(invisible(sample.chain))

  }

#' Generate an HDP Gibbs sampling chain from a spectra catalog.
#'
#' @inheritParams PrepInit
#'
#' @inheritParams SetupAndActivate
#'
#' @inheritParams PriorSetupAndActivate
#'
#' @inheritParams ChainBurnin
#'
#' @param post.n Pass to \code{\link[hdpx]{hdp_posterior_sample}}
#'      \code{n}.The number of posterior samples to collect.
#'
#' @param post.space Pass to \code{\link[hdpx]{hdp_posterior_sample}}
#'      \code{space}. The number of iterations between collected samples.
#'
#' @param post.cpiter Pass to \code{\link[hdpx]{hdp_posterior_sample}} and
#'        \code{\link[hdpx]{hdp_burnin}} \code{cpiter}.The number of iterations of concentration
#'        parameter sampling to perform after each iteration
#'
#' @param post.verbosity Pass to \code{\link[hdpx]{hdp_posterior_sample}}
#'      \code{verbosity}. Verbosity of debugging statements.
#'       No need to change unless for development purpose
#'
#' @param checkpoint If \code{TRUE}, then \itemize{
#'      \item Checkpoint each final Gibbs sample
#'        chain to the current working directory, in a file called
#'        mSigHdp.sample.checkpoint.*seed_number*.Rdata.
#'      \item Periodically checkpoint the burnin state
#'        to the current working directory, in files called
#'        mSigHdp.burnin.checkpoint.*seed_number*.Rdata.
#'  }
#'
#'  All checkpoint files go in current working directory.
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
           gamma0.alpha        = gamma.alpha,
           gamma0.beta         = gamma.beta,
           burnin.multiplier   = 2,
           checkpoint          = TRUE,
           prior.sigs          = NULL,
           prior.pseudoc       = NULL)
  {

    if(!is.null(prior.sigs)){
      if(verbose) message('Prior signatures are found')
      stop("Computation with prior.sigs not currently supported")
      if(is.null(prior.pseudoc)){
        stop('Prior signatures pseudo counts are not set')
      } else {
        hdp.state <- PriorSetupAndActivate(input.catalog = input.catalog,
                                           seedNumber    = seedNumber,
                                           K.guess       = K.guess,
                                           multi.types   = FALSE, ##multi.types=T doesn't work for now
                                           verbose       = verbose,
                                           gamma.alpha   = gamma.alpha,
                                           gamma.beta    = gamma.beta,
                                           gamma0.alpha  = gamma0.alpha,
                                           gamma0.beta   = gamma0.beta,
                                           prior.sigs    = prior.sigs,
                                           prior.pseudoc = prior.pseudoc)
      }

    } else { # No prior signatures (is.null(prior.sigs))
      hdp.state <- SetupAndActivate(input.catalog = input.catalog,
                                    seedNumber    = seedNumber,
                                    K.guess       = K.guess,
                                    multi.types   = multi.types,
                                    verbose       = verbose,
                                    gamma.alpha   = gamma.alpha,
                                    gamma.beta    = gamma.beta,
                                    gamma0.alpha  = gamma0.alpha,
                                    gamma0.beta   = gamma0.beta)
    }

    if (verbose) message("calling hdp_posterior, seed = ",
                         seedNumber, " ", Sys.time())


    burnin.output <- ChainBurnin(
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



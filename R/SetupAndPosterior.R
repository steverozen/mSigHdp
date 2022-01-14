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
#' @inheritParams GibbsSamplingAfterBurnin
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
           checkpoint          = TRUE,
           prior.sigs          = NULL,
           prior.pseudoc       = NULL)
  {

    if(!is.null(prior.sigs)) {
      stop("Computation with prior.sigs not currently supported")

      if(verbose) message('Prior signatures are found')

      if(is.null(prior.pseudoc)){
        stop('Prior signatures pseudo counts are not set')
      } else {
        hdp.state <- PriorSetupAndActivate(input.catalog = input.catalog,
                                           seedNumber    = seedNumber,
                                           K.guess       = K.guess,
                                           multi.types   = FALSE, ## multi.types=T doesn't work for now
                                           verbose       = verbose,
                                           gamma.alpha   = gamma.alpha,
                                           gamma.beta    = gamma.beta,
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
                                    gamma.beta    = gamma.beta)
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

      sample.chain <- GibbsSamplingAfterBurnin(burnin.output     = burnin.output,
                                                     post.n            = post.n,
                                                     post.space        = post.space,
                                                     post.cpiter       = post.cpiter,
                                                     post.verbosity    = post.verbosity)
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



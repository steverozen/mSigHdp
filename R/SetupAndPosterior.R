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
#' @param checkpoint.1.chain If \code{TRUE} checkpoint the sample
#'      chain to current working directory, in a file called
#'      sample.chain.*seed_number*.Rdata.
#'
#' @param posterior.checkpoint If \code{TRUE} checkpoint the posterior sampling after every 10
#'                             posterior samples collected
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
           burnin              = 5000,
           post.n              = 50,
           post.space          = 50,
           post.cpiter         = 3,
           post.verbosity      = 0,
           gamma.alpha         = 1,
           gamma.beta          = 20,
           gamma0.alpha        = gamma.alpha,
           gamma0.beta         = gamma.beta,
           checkpoint.1.chain  = TRUE,
           burnin.multiplier   = 2,
           burnin.checkpoint   = TRUE,
           prior.sigs          = NULL,
           prior.pseudoc       = NULL,
           posterior.checkpoint= F)
  {

    if(!is.null(prior.sigs)){
      if(verbose) message('Prior signatures are found')
      if(is.null(prior.pseudoc)){
        stop('Prior signatures pseudo counts are not set')
      }else{
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

    }else{
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
      burnin.checkpoint = burnin.checkpoint)

    posterior.time <- system.time(

      sample.chain <- hdpx::hdp_posterior_sample(post.input     = burnin.output,
                                                 post.n         = post.n,
                                                 post.space     = post.space,
                                                 post.cpiter    = post.cpiter,
                                                 seed           = seedNumber,
                                                 post.verbosity = post.verbosity,
                                                 checkpoint     = posterior.checkpoint)

    )



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



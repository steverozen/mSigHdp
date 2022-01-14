#' Start Gibbs sampling after burnin
#'
#'
#' @param burnin.output A path to burnin checkpoint Rdata or An S4 object from \code{ChainBurnin}
#'
#'@param post.n Pass to \code{\link[hdpx]{hdp_posterior_sample}}
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
#'       No need to change unless for development purpose. Default is 0.
#'
#' @return  A hdpSampleChain object with the salient information from each
#'  posterior sample. See \code{\link{hdpSampleChain-class}}
#' @keywords internal

GibbsSamplingAfterBurnin<- function(burnin.output     = burnin.output,
                                    post.n         = post.n,
                                    post.space     = post.space,
                                    post.cpiter    = 3,
                                    post.verbosity = 0,
                                    seedNumber     = NULL){
  if(is.character(burnin.output)){
    if(grepl("mSigHdp.burnin.checkpoint", burnin.output)){
      load(burnin.output)
    }else{
      stop("The input is not a burnin checkpoint")
    }

  }

  #We allow the user to use a different seed when running Gibbs sampling, but we don't recommend.
  if(is.null(seedNumber)){
    seedNumber <- burnin.output$hdplist$seed_activate
  }
  ##call hdpx::hdp_posterior_sample function to run Gibbs sampling
  sample.chain <- hdpx::hdp_posterior_sample(post.input     = burnin.output,
                                             post.n         = post.n,
                                             post.space     = post.space,
                                             post.cpiter    = post.cpiter,
                                             seed           = seedNumber,
                                             post.verbosity = post.verbosity)
  return(sample.chain)

}

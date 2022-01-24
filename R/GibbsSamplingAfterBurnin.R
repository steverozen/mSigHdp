#' Start Gibbs sampling on one chain after burnin
#'
#' This function might be used to start Gibbs sampling after
#' \code{\link{ExtendBurnin}}.
#'
#' @param burnin.output A path to burnin checkpoint Rdata
#'  or to an S4 object from \code{\link{Burnin}}.
#'
#' @param post.n The number of posterior samples to collect.
#   (Passed to argument \code{n} in
#   \code{\link[hdpx]{hdp_posterior_sample}}.)
#'
#' @param post.space The number of iterations between collected samples.
#       (Passed to argument \code{space} in
#        \code{\link[hdpx]{hdp_posterior_sample}}.)
#'
#' @param post.cpiter The number of iterations of concentration
#'  parameter sampling
#'  to perform after each main Gibbs-sample iteration. (See Teh et al.,
#'  "Hierarchical Dirichlet Processes", Journal of the American Statistical
#'  Association 2006;101(476):1566-1581
#'  (https://doi.org/10.1198/016214506000000302).)
#        (Passed to argument \code{cpiter} in
#        \code{\link[hdpx]{hdp_posterior_sample}} and
#        \code{\link[hdpx]{hdp_burnin}}.)
#'
#' @param post.verbosity Verbosity of debugging statements.
#'       No need to change unless for testing or debugging.
#       (Passed to argument \code{verbosity} in
#        \code{\link[hdpx]{hdp_posterior_sample}}.)
#'
#' @return  An \code{hdpSampleChain}
#'  S4 object with the salient information from each
#'  posterior sample. See \code{\link{hdpSampleChain-class}}
#'  in package hdpx.
#'
#' @param seedNumber A random seed that ensures reproducible
#'   results.
#'
#'
#' @export

GibbsSamplingAfterBurnin<- function(burnin.output,
                                    post.n,
                                    post.space,
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

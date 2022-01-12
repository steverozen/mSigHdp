#' Generate an HDP Gibbs sampling chain from a spectra catalog.
#' @inheritParams GeneratePriorppindex
#'
#' @param K.guess Suggested initial value of the number of
#' signatures, passed to \code{\link[hdpx]{dp_activate}} as
#' \code{initcc}.
#'
#' @param verbose If \code{TRUE} then \code{message} progress information.
#' @param seedNumber A random seeds passed to \code{\link[hdpx]{dp_activate}}.
#' @param gamma.alpha Shape parameter of
#'   the gamma distribution prior for the Dirichlet process concentration
#'   parameters; in this
#'   function the gamma distributions for all Dirichlet processes,
#'   except possibly the top level process, are the same.
#'
#' @param gamma.beta Inverse scale parameter (rate parameter) of
#'   the gamma distribution prior for the Dirichlet process concentration
#'   parameters; in this
#'   function the gamma distributions for all Dirichlet processes, except
#'   possibly the top level process, are the same.
#'
#' @param gamma0.alpha See figure B.1 from Nicola Robert's thesis.
#'   The shape parameter (\eqn{\alpha_0}) of the gamma
#'   distribution priors for
#'   the Dirichlet process concentration parameters  (\eqn{\gamma_0})
#'   for \eqn{G_0}.
#'
#' @param gamma0.beta See figure B.1 from Nicola Robert's thesis.
#'   Inverse scale parameter (rate parameter, \eqn{\beta_0}) of the gamma
#'   distribution priors for
#'   the Dirichlet process concentration parameters (\eqn{\gamma_0})
#'   for \eqn{G_0}.
#'
#' @param prior.sigs DELETE ME LATER, NOT SUPPORTED. A matrix containing prior signatures.
#'
#' @param prior.pseudoc DELETE ME LATER, NOT SUPPORTED. A numeric list. Pseudo counts of each prior signature. Recommended is 1000. In practice,
#'                       it may be advisable to put lower weights on prior signatures that you do not expect
#'                       to be present in your dataset, or even exclude some priors entirely.
#'
#' @return Invisibly, an \code{\link[hdpx]{hdpState-class}} object
#'  as returned from \code{\link[hdpx]{dp_activate}}.
#'
#' @keywords internal

PriorSetupAndActivate <- function(prior.sigs,
                                  prior.pseudoc,
                                  gamma.alpha       = 1,
                                  gamma.beta        = 1,
                                  K.guess,
                                  gamma0.alpha      = gamma.alpha,
                                  gamma0.beta       = gamma.beta,
                                  multi.types       = F,
                                  input.catalog,
                                  verbose           = TRUE,
                                  seedNumber        = 1){


  input.catalog <- GetInputCatalogAsMatrix(input.catalog)
  nps <- ncol(prior.sigs)
  number.channels <- nrow(prior.sigs)
  convSpectra <- t(input.catalog)
  if(verbose) message(nps," prior signatures are defined")

  hdp.prior <- hdpx::hdp_prior_init(prior_distn = prior.sigs,
                                    prior_pseudoc = prior.pseudoc,
                                    hh=rep(1, number.channels),
                                    alphaa=c(gamma.alpha, gamma.alpha),
                                    alphab=c(gamma.beta, gamma.beta))

  ppindex <- GeneratePriorppindex(multi.types = multi.types,input.catalog = input.catalog,nps = nps)


  dp.levels <- length(unique(ppindex$ppindex))

  alphaa <- c(gamma0.alpha, rep(gamma.alpha, dp.levels - 1))
  alphab <- c(gamma0.beta,  rep(gamma.beta,  dp.levels - 1))

  hdp.prior <- hdpx::hdp_addconparam(hdp.prior,
                                     alphaa = c(gamma0.alpha, rep(gamma.alpha, dp.levels - 1)),
                                     alphab = c(gamma0.beta,  rep(gamma.beta,  dp.levels - 1)))



  hdp.prior <- hdpx::hdp_adddp(hdp.prior,
                               numdp = length(ppindex$ppindex), ##number of samples + 1
                               ppindex = ppindex$ppindex,
                               cpindex = ppindex$cpindex)


  num.process <- hdpx::numdp(hdp.prior)


  hdp.prior <- hdpx::hdp_setdata(hdp.prior,
                                 dpindex = (1+nps+1 + ppindex$num.tumor.types + 1):num.process,
                                 convSpectra)

  hdp.state <- hdpx::dp_activate(hdp.prior,
                                 dpindex =  (1+nps+1):num.process,#don’t activate the frozen pseudo-count nodes for the prior signatures
                                 initcc = nps + K.guess,
                                 seed = seedNumber)

  return(invisible(hdp.state))

}

#' Generate an HDP Gibbs sampling chain from a spectra catalog.
#' @inheritParams GeneratePriorppindex
#'
#' @param K.guess Suggested initial value of the number of
#' signatures, passed to \code{\link[hdpx]{dp_activate}} as
#' \code{initcc}.
#'
#' @param verbose If \code{TRUE} then \code{message} progress information.
#' @param seedNumber A random seeds passed to \code{\link[hdpx]{dp_activate}}.
#'
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
#' @return Invisibly, an \code{\link[hdpx]{hdpState-class}} object
#'  as returned from \code{\link[hdpx]{dp_activate}}.
#'
#' @keywords internal

PriorSetupAndActivate <- function(prior.sigs,
                                  prior.pseudoc,
                                  gamma.alpha       = 1,
                                  gamma.beta        = 20,
                                  K.guess,
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

  ppindex <- GeneratePriorppindex(
    multi.types = multi.types,
    input.catalog = input.catalog,
    nps = nps)

  dp.levels <- length(unique(ppindex$ppindex))

  alphaa <- rep(gamma.alpha, dp.levels)
  alphab <- rep(gamma.beta,  dp.levels)

  hdp.prior <- hdpx::hdp_addconparam(hdp.prior,
                                     alphaa = alphaa,
                                     alphab = alphab)

  hdp.prior <- hdpx::hdp_adddp(hdp.prior,
                               numdp   = length(ppindex$ppindex),
                               ppindex = ppindex$ppindex,
                               cpindex = ppindex$cpindex)

  num.process <- hdpx::numdp(hdp.prior)

  hdp.prior <- hdpx::hdp_setdata(hdp.prior,
                                 dpindex = (1+nps+1 + ppindex$num.tumor.types + 1):num.process,
                                 convSpectra)

  hdp.state <- hdpx::dp_activate(hdp.prior,
                                 dpindex =  (1+nps+1):num.process, # don’t activate the frozen pseudo-count nodes for the prior signatures
                                 initcc = nps + K.guess,
                                 seed = seedNumber)

  return(invisible(hdp.state))
}

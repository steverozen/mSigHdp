#' Prepare an \code{\link[hdpx]{hdpState-class}} object and run the Gibbs sampling burnin.
#'
#' @inheritParams PrepInit
#'
#' @param seedNumber An integer that is used to generate separate
#'   random seeds for the call to \code{\link[hdpx]{dp_activate}},
#'   and before the call of \code{\link[hdpx]{hdp_burnin}}.
#'
#' @param burnin Pass to \code{\link[hdpx]{hdp_burnin}}
#'      \code{burnin}.
#'
#' @param cpiter Pass to \code{\link[hdpx]{hdp_burnin}}
#'      \code{cpiter}.
#'
#' @param burnin.verbosity Pass to \code{\link[hdpx]{hdp_burnin}}
#'      \code{verbosity}.
#'
#' @return A list with 2 elements: \describe{
#' \item{\code{hdplist}}{A list representation of
#'    an \code{\link[hdpx]{hdpState-class}} object.}
#' \item{likelihood}{A numeric vector with the likelihood at each iteration.}
#' }
#'
#' @export
#'
ActivateAndBurnin2 <-
  function(input.catalog,
           seedNumber        = 1,
           K.guess,
           multi.types       = FALSE,
           verbose           = TRUE,
           burnin            = 4000,
           cpiter            = 3,
           burnin.verbosity  = 0,
           gamma.alpha       = 1,
           gamma.beta        = 1,
           gamma0.alpha      = 1,
           gamma0.beta       = 1,
           burnin.multiplier = 1,
           burnin.checkpoint = FALSE
  ) { # 10 arguments

    ### DELETE
    hdp.state <- SetupAndActivate(input.catalog = input.catalog,
                                  seedNumber    = seedNumber,
                                  K.guess       = K.guess,
                                  multi.types   = multi.types,
                                  verbose       = verbose,
                                  gamma.alpha   = gamma.alpha,
                                  gamma.beta    = gamma.beta,
                                  gamma0.alpha  = gamma0.alpha,
                                  gamma0.beta   = gamma0.beta)


    set.seed(seedNumber) ###?

    for (ii in 1:burnin.multiplier) {
      burnin.checkpoint <- hdpx::hdp_burnin(hdp         = hdp.state,
                                 burnin      = burnin,
                                 cpiter      = cpiter,
                                 verbosity   = burnin.verbosity)
      if (burnin.checkpoint) {
         save(burnin.checkpoint,
              file = paste0("burnin.checkpoint.", seedNumber, ".Rdata"))
      }
    }

    return(invisible(output))
  }

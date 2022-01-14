
#'
#' Prepare an \code{\link[hdpx]{hdpState-class}} object and run the Gibbs sampling burnin.
#'
#' @param hdp.state An \code{\link[hdpx]{hdpState-class}} object or a list
#'   representation of an \code{\link[hdpx]{hdpState-class}} object.
#'
#' @param seedNumber An integer that is used to generate separate
#'   random seeds for the call to \code{\link[hdpx]{dp_activate}},
#'   and before the call of \code{\link[hdpx]{hdp_burnin}}.
#'
#' @param burnin Pass to \code{\link[hdpx]{hdp_burnin}}
#'      \code{burnin}. The number of burn-in iterations
#'
#' @param burnin.multiplier Run \code{burnin.multiplier} rounds of \code{burnin}
#'   iterations.
#'   If \code{checkpoint} is \code{TRUE}, save the burnin chain (see
#'   parameter \code{checkpoint}.)
#'   The diagnostic plot \code{diagnostics.likelihood.pdf} can help
#'   determine if the chain is stationary.
#'   The burnin can be continued from a checkpoint file
#'   with \code{\link{ExtendBurnin}} (see argument \code{checkpoint}).
#'
#' @param cpiter Pass to \code{\link[hdpx]{hdp_burnin}}
#'  \code{cpiter}. The number of iterations of concentration
#'  parameter sampling
#'  to perform after each iteration. (See Teh et al.,
#" Hierarchical Dirichlet Processes", Journal of the American Statistical
#' Association 2006;101(476):1566-1581
#' (https://doi.org/10.1198/016214506000000302).)
#'
#' @param burnin.verbosity Pass to \code{\link[hdpx]{hdp_burnin}}
#'      \code{verbosity}.Verbosity of debugging statements.
#'#'
#' @param checkpoint If \code{TRUE}, a checkpoint
#' file called mSigHdp.burnin.checkpoint.*seed_number*.Rdata.
#'  will be created
#'  in the current working directory.
#'
#' @return A list with 2 elements: \describe{
#' \item{\code{hdplist}}{A list representation of
#'    an \code{\link[hdpx]{hdpState-class}} object.}
#' \item{likelihood}{A numeric vector with the likelihood at each iteration.}
#' }
#'
#' @keywords internal
#'
ChainBurnin <-
  function(hdp.state,
           seedNumber        = 1,
           burnin            = 5000,
           cpiter            = 3,
           burnin.verbosity  = 0,
           burnin.multiplier = 2,
           checkpoint        = TRUE
  ) {

    set.seed(seedNumber)

    my.checkpoint <- function() {
      if (checkpoint) {
        save(burnin.output,
             file = paste0("mSigHdp.burnin.checkpoint.", seedNumber, ".Rdata"))
      }
    }

    burnin.output <- hdpx::hdp_burnin(hdp         = hdp.state,
                                      burnin      = burnin,
                                      cpiter      = cpiter,
                                      verbosity   = burnin.verbosity)
    my.checkpoint()

    if(burnin.multiplier>1){
      for (ii in 2:burnin.multiplier) {
        old.likelihood <- burnin.output$likelihood
        burnin.output <- hdpx::hdp_burnin(hdp         = burnin.output$hdplist,
                                          burnin      = burnin,
                                          cpiter      = cpiter,
                                          verbosity   = burnin.verbosity)
        burnin.output$likelihood <- c(old.likelihood,burnin.output$likelihood)
        my.checkpoint
      }
    }

    return(invisible(burnin.output))
  }

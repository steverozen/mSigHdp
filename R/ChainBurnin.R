
#'
#' Prepare an \code{\link[hdpx]{hdpState-class}} object and run the Gibbs sampling burnin.
#'
#' @param hdp.state An \code{\link[hdpx]{hdpState-class}} object or a list representation of an \code{\link[hdpx]{hdpState-class}} object.
#'
#' @param seedNumber An integer that is used to generate separate
#'   random seeds for the call to \code{\link[hdpx]{dp_activate}},
#'   and before the call of \code{\link[hdpx]{hdp_burnin}}.
#'
#' @param burnin Pass to \code{\link[hdpx]{hdp_burnin}}
#'      \code{burnin}. The number of burn-in iterations
#'
#' @param burnin.multiplier A checkpoint setting. \code{burnin.multiplier} rounds of \code{burnin} iterations will be run.
#'        After each round, a burn-in chain will be save for checkpoint.
#'        A total number of 10,000 iterations is recommended for most analysis. However, number of iterations
#'        can be adjusted based on the size of dataset. According to our experience, 20,000 iterations are
#'        needed when analyzing all PCAWG7 genomes (2,780 samples).
#'        The burnin can be continued from a checkpoint file with \code{\link{ExtendBurnin}}.
#'
#' @param cpiter Pass to \code{\link[hdpx]{hdp_burnin}}
#'      \code{cpiter}. The number of iterations of concentration parameter sampling
#'  to perform after each iteration.
#'
#' @param burnin.verbosity Pass to \code{\link[hdpx]{hdp_burnin}}
#'      \code{verbosity}.Verbosity of debugging statements.
#'
#'
#' @param burnin.checkpoint Default is False. If True, a checkpoint for burnin will be created.
#'
#' @return A list with 2 elements: \describe{
#' \item{\code{hdplist}}{A list representation of
#'    an \code{\link[hdpx]{hdpState-class}} object.}
#' \item{likelihood}{A numeric vector with the likelihood at each iteration.}
#' }
#'
#' @export
#'
ChainBurnin <-
  function(hdp.state,
           seedNumber        = 1,
           burnin            = 5000,
           cpiter            = 3,
           burnin.verbosity  = 0,
           burnin.multiplier = 2,
           burnin.checkpoint = TRUE
  ) { # 10 arguments

    set.seed(seedNumber)

    burnin.output <- hdpx::hdp_burnin(hdp         = hdp.state,
                                      burnin      = burnin,
                                      cpiter      = cpiter,
                                      verbosity   = burnin.verbosity)
    if (burnin.checkpoint) {
      save(burnin.checkpoint,
           file = paste0("burnin.checkpoint.", seedNumber, ".Rdata"))
    }


    if(burnin.multiplier>1){
      for (ii in 2:burnin.multiplier) {
        old.likelihood <- burnin.output$likelihood
        burnin.output <- hdpx::hdp_burnin(hdp         = burnin.output$hdplist,
                                          burnin      = burnin,
                                          cpiter      = cpiter,
                                          verbosity   = burnin.verbosity)
        burnin.output$likelihood <- c(old.likelihood,burnin.output$likelihood)

        if (burnin.checkpoint) {
            save(burnin.output,
               file = paste0("burnin.checkpoint.", seedNumber, ".Rdata"))
        }
      }
    }



    return(invisible(burnin.output))
  }

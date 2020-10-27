
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
#'      \code{burnin}.
#'
#' @param cpiter Pass to \code{\link[hdpx]{hdp_burnin}}
#'      \code{cpiter}.
#'
#' @param burnin.verbosity Pass to \code{\link[hdpx]{hdp_burnin}}
#'      \code{verbosity}.
#'
#' @param burnin.multiplier A checkpoint setting. \code{burnin.multiplier} rounds of \code{burnin} iterations will be run.
#'        After each round, a burn-in chain will be save for checkpoint.
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
           burnin            = 4000,
           cpiter            = 3,
           burnin.verbosity  = 0,
           burnin.multiplier = 1,
           burnin.checkpoint = FALSE
  ) { # 10 arguments

    set.seed(seedNumber) ###?

    burnin.output <- hdpx::hdp_burnin(hdp         = hdp.state,
                                      burnin      = burnin,
                                      cpiter      = cpiter,
                                      verbosity   = burnin.verbosity)
    if (burnin.checkpoint) {
      save(burnin.checkpoint,
           file = paste0(out.dir,"/checkpoint.Rdatas/","burnin.checkpoint.", seedNumber, ".Rdata"))
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
               file = paste0(out.dir,"/checkpoint.Rdatas/","burnin.checkpoint.", seedNumber, ".Rdata"))
        }
      }
    }



    return(invisible(burnin.output))
  }

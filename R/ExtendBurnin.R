#' Extend Burn in iteration for a list representation of
#'    an \code{\link[hdpx]{hdpState-class}} object. This list
#'    is an output from \code{\link[hdpx]{hdp_burnin}} or
#'    \code{ActivateandBurnin}.
#'
#' @param seedNumber A random seed for setting the environment of
#'   \code{\link[hdpx]{hdp_burnin}}.
#'
#' @param burnin Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{burnin}.
#'
#' @param cpiter Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{cpiter}.
#'
#' @param verbosity Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{verbosity}.
#'
#' @param hdplist A list representation of
#'    an \code{\link[hdpx]{hdpState-class}} object
#'
#' @return A list with hdp object after burn-in iteration and likelihood of iteration
#'
#' @export
#'
ExtendBurnin <-
  function(hdplist,
           seedNumber     = 1,
           burnin    = 4000,
           cpiter    = 3,
           verbosity = 0
  ) { # 11 arguments

    as.hdpState <- utils::getFromNamespace(x = "as.hdpState", ns = "hdpx")

    hdp.state <- as.hdpState(hdplist)

    set.seed(seedNumber)
    output <- hdpx::hdp_burnin (hdp         = hdp.state,
                                burnin      = burnin,
                                cpiter      = cpiter,
                                verbosity   = verbosity)
    return(invisible(output))
  }

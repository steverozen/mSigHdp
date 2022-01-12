#' Extend burnin iterations generated from \code{\link{ChainBurnin}}.
#'
#' @param seedNumber A random seed for setting the environment of
#'   \code{\link[hdpx]{hdp_burnin}}.
#'
#' @param burnin Passed to \code{\link[hdpx]{hdp_posterior}}
#'      \code{burnin}.
#'
#' @param cpiter Passed to \code{\link[hdpx]{hdp_posterior}}
#'      \code{cpiter}. Please see that documentation.
#'
#' @param verbosity Passed to \code{\link[hdpx]{hdp_posterior}}
#'      \code{verbosity}.
#'
#' @param hdplist A list representation of
#'    an \code{\link[hdpx]{hdpState-class}} object;
#'    a return value \code{\link{ChainBurnin}}.
#'
#' @return The same type of object as returned from
#'  \code{\link{ChainBurnin}.}
#'
#' The envisioned application is extending burnins from burnin checkpoints.
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

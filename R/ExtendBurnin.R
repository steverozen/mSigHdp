#' Extend burnin iterations generated from \code{\link{Burnin}}
#'
#' @param seedNumber A random seed for reproducible results.
#'
#' @param burnin The number of burnin iterations to perform.
# Passed \code{\link[hdpx]{hdp_posterior}} argument
#      \code{burnin} (package hdpx).
#'
#' @param cpiter The number of iterations of concentration
#'  parameter sampling
#'  to perform after each main Gibbs-sample iteration. (See Teh et al.,
#' "Hierarchical Dirichlet Processes", Journal of the American Statistical
#' Association 2006;101(476):1566-1581
#' (https://doi.org/10.1198/016214506000000302).)
#  Passed to argument \code{cpiter} in \code{\link[hdpx]{hdp_burnin}}
#  in package hdpx.
#'
#' @param verbosity Passed to \code{\link[hdpx]{hdp_posterior}}
#'      \code{verbosity}.
#'
#' @param hdplist A list representation of
#'    an \code{\link[hdpx]{hdpState-class}} object;
#'    a return value \code{\link{Burnin}}.
#'
#' @return The same type of object as returned from
#'  \code{\link{Burnin}.}
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
  ) {

    as.hdpState <- utils::getFromNamespace(x = "as.hdpState", ns = "hdpx")

    hdp.state <- as.hdpState(hdplist)

    set.seed(seedNumber)
    output <- hdpx::hdp_burnin (hdp         = hdp.state,
                                burnin      = burnin,
                                cpiter      = cpiter,
                                verbosity   = verbosity)
    return(invisible(output))
  }

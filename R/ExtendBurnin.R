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
#' @param burnin.verbosity Number that controls whether progress
#'    messages are printed.
#'
#' @param previous.burnin.output Output from
#'    \code{\link{Burnin}} or the file path of
#'    a checkpoint file written by \code{\link{Burnin}}.
#'
#' @return The same type of object as returned from
#'  \code{\link{Burnin}.}
#'
#' The envisioned application is extending burnins from burnin checkpoints.
#'
#' @export
#'
ExtendBurnin <-
  function(previous.burnin.output,
           burnin    = 5000,
           cpiter    = 3,
           burnin.verbosity = 0,
           seedNumber = NULL
  ) {

    if(is.character(previous.burnin.output)){
      if(grepl("mSigHdp.burnin.checkpoint", previous.burnin.output)) {
        tmp.env <- new.env()
        load(previous.burnin.output, envir = tmp.env)
        vars <- ls(tmp.env)
        if (length(vars) > 1) {
          stop("The input file does not seem to be a burnin checkpoint")
        }
        previous.burnin.output <- tmp.env[[vars]]

      }else{
        stop("The input is not a burnin checkpoint")
      }
    }

    if(is.null(seedNumber)){
      seedNumber <- previous.burnin.output$hdplist$seed_activate
    }

    set.seed(seedNumber)

    output <- hdpx::hdp_burnin (hdp         = previous.burnin.output$hdplist,
                                burnin      = burnin, # This is the number of burnin iterations
                                cpiter      = cpiter,
                                verbosity   = burnin.verbosity)
    output$likelihood <- c(previous.burnin.output$likelihood, output$likelihood)
    return(invisible(output))
  }

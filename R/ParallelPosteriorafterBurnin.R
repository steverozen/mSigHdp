#' Generate an HDP Gibbs sampling chain from a spectra catalog.
#'
#' @param retval A list object containing hdplist after burn-in iteration
#'               and likelihood from \code{BurninIteration}.
#'
#' @param seedNumber Pass to \code{\link[hdpx]{hdp_posterior}}
#'
#' @param post.burnin Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{burnin}. This can be set to a small number
#'
#' @param post.n Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{n}.
#'
#' @param post.space Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{space}.
#'
#' @param post.cpiter Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{cpiter}.
#'
#' @param post.verbosity Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{verbosity}.
#'
#' @return Invisibly, an \code{\link[hdpx]{hdpSampleChain-class}} object
#'  as returned from \code{\link[hdpx]{hdp_posterior}}.
#'
#' @export

ParallelPosteriorafterBurnin <-
  function(retval,
           seedNumber          = 1,
           verbose             = TRUE,
           post.burnin         = 4000,
           post.n              = 50,
           post.space          = 50,
           post.cpiter         = 3,
           post.verbosity      = 0,
           num.child.process   = 2,
           CPU.cores           = 2
          )
{ # 12 arguments

    hdp.state <- hdpx:::as.hdpState(retval$hdplist)

    run.posterior <- function(seedNumber) {

      if (verbose) message("Runing run.posterior on ", seedNumber)

      if (verbose) message("calling hdp_posterior, seed = ",
                           seedNumber, " ", Sys.time())
      posterior.time <- system.time(
        sample.chain <- hdpx::hdp_posterior(
          hdp       = hdp.state,
          verbosity = post.verbosity,
          burnin    = post.burnin,
          n         = post.n,
          space     = post.space,
          cpiter    = post.cpiter,
          seed      = seedNumber)
      )

      if (verbose) {
        message("compute sample.chain time: ")
        for (xn in names(posterior.time)) {
          message(" ", xn, " ", posterior.time[[xn]])
        }
      }
      return(invisible(sample.chain))

    }

    chlist <- parallel::mclapply(
      # Must choose a different seed for each of the chains
      X = (seedNumber + 1:num.child.process * 10^6) ,
      FUN = run.posterior,
      mc.cores = CPU.cores)

    clean.chlist <- CleanChlist(chlist)

    return(invisible(clean.chlist))

  }

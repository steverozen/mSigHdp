#' Initialize hdp object
#' Allocate process index for hdp initialization.
#' Prepare for \code{\link[hdpx]{hdp_init}}
#'
#' @inheritParams Generateppindex
#'
#'
#' @param K.guess Suggested initial value of the number of
#' signatures, passed to \code{\link[hdpx]{dp_activate}} as
#' \code{initcc}.
#'
#' @param verbose If \code{TRUE} then \code{message} progress information.
#'
#' @param gamma.alpha Shape parameter of gamma distribution from which
#'   the Dirichlet process concentration parameters are drawn; in this
#'   function the gamma distributions for all Dirichlet processes are the same.
#'
#' @param gamma.beta Inverse scale parameter (rate parameter) of gamma distribution
#'   from which the Dirichlet process concentration parameters are drawn; in this
#'   function the gamma distributions for all Dirichlet processes are the same.
#'
#' @export

PrepInit <- function(multi.types,
                     input.catalog,
                     verbose,
                     K.guess,
                     gamma.alpha=1,
                     gamma.beta=1,
                     one.parent.hack = FALSE){

  if (mode(input.catalog) == "character") {
    if (verbose) message("Reading input catalog file ", input.catalog)
    input.catalog <- ICAMS::ReadCatalog(input.catalog, strict = FALSE)
  } else {
    input.catalog <- input.catalog
  }

  convSpectra <- t(input.catalog)

  number.channels <- nrow(input.catalog)
  number.samples  <- ncol(input.catalog)

  if (verbose) {
    message("Guessed number of signatures ",
            "(= Dirichlet process data clusters) = ", K.guess)
  }

  ppindex <- Generateppindex(multi.types = multi.types,
                             one.parent.hack = one.parent.hack,
                             input.catalog = input.catalog)

  # cpindex (concentration parameter)
  cpindex <- 1 + ppindex$ppindex

  ## Calculate the number of levels in the DP node tree.
  dp.levels <- length(unique(ppindex$ppindex))

  if (verbose) {
    message("Gamma distribution was set to shape = ", gamma.alpha,
            " rate (inverse scale) = ", gamma.beta)
  }

  alphaa <- rep(gamma.alpha,dp.levels)
  alphab <- rep(gamma.beta,dp.levels)

  invisible(list(num.tumor.types = ppindex$num.tumor.types,
                 number.channels = number.channels,
                 convSpectra     = convSpectra,
                 ppindex         = ppindex$ppindex,
                 cpindex         = cpindex,
                 alphaa          = alphaa,
                 alphab          = alphab))
}

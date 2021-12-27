#' Initialize hdp object
#' Allocate process index for hdp initialization.
#' Prepare for \code{\link[hdpx]{hdp_init}}
#'
#' @inheritParams Generateppindex
#'
#'
#' @param K.guess Suggested initial value of the number of
#'                clusters. Usually, the number of clusters is two times of the number
#'                of extracted signatures. Passed to \code{\link[hdpx]{dp_activate}} as
#'                \code{initcc}.
#'
#' @param verbose If \code{TRUE} then \code{message} progress information.
#'
#' @param gamma.alpha Shape parameter of
#'   the gamma distribution prior for the Dirichlet process concentration
#'   parameters; in this
#'   function the gamma distributions for all Dirichlet processes,
#'   except possibly the top level process, are the same.
#'
#' @param gamma.beta Inverse scale parameter (rate parameter) of
#'   the gamma distribution prior for the Dirichlet process concentration
#'   parameters; in this
#'   function the gamma distributions for all Dirichlet processes, except
#'   possibly the top level process, are the same.
#'
#'   We recommend gamma.alpha = 1 and gamma.beta = 20 for single-base-substitution signatures extraction;
#'   gamma.alpha = 1 and gamma.beta = 50  for doublet-base-substitution/INDEL signature extraction
#'
#' @param gamma0.alpha See figure B.1 from Nicola Robert's thesis.
#'   The shape parameter (\eqn{\alpha_0}) of the gamma
#'   distribution priors for
#'   the Dirichlet process concentration parameters  (\eqn{\gamma_0})
#'   for \eqn{G_0}.
#'
#' @param gamma0.beta See figure B.1 from Nicola Robert's thesis.
#'   Inverse scale parameter (rate parameter, \eqn{\beta_0}) of the gamma
#'   distribution priors for
#'   the Dirichlet process concentration parameters (\eqn{\gamma_0})
#'   for \eqn{G_0}.
#'
#' @param IS.ICAMS If \code{TRUE} using ICAMS functions to plot, read and write signatures. Set to \code{FALSE}
#'   if your input cannot be taken by ICAMS.
#'
#' @export

PrepInit <- function(multi.types,
                     input.catalog,
                     IS.ICAMS     = TRUE,
                     verbose      = TRUE,
                     K.guess,
                     gamma.alpha  = 1,
                     gamma.beta   = 1,
                     gamma0.alpha = gamma.alpha,
                     gamma0.beta  = gamma.beta){

  if (mode(input.catalog) == "character") {
    if (verbose) message("Reading input catalog file ", input.catalog)
    if(IS.ICAMS==T){
      input.catalog <- ICAMS::ReadCatalog(input.catalog)
    }else{
      input.catalog <- data.table::fread(input.catalog)
    }
  } else {
    input.catalog <- input.catalog
  }
  input.catalog <- input.catalog[,colSums(input.catalog)>0] ##0 catalogs will cause problems in later steps

  convSpectra <- t(input.catalog)

  number.channels <- nrow(input.catalog)
  number.samples  <- ncol(input.catalog)

  if (verbose) {
    message("Guessed number of signatures ",
            "(= Dirichlet process data clusters) = ", K.guess)
  }

  ppindex <- Generateppindex(multi.types = multi.types,
                             input.catalog = input.catalog)

  # cpindex (concentration parameter)
  cpindex <- 1 + ppindex$ppindex

  ## Calculate the number of levels in the DP node tree.
  dp.levels <- length(unique(ppindex$ppindex))

  if (verbose) {
    message("Gamma distribution was set to shape = ", gamma.alpha,
            " rate (inverse scale) = ", gamma.beta)
    message("gamma0.alpha=", gamma0.alpha, "; gamma0.beta=", gamma0.beta)
  }

  alphaa <- c(gamma0.alpha, rep(gamma.alpha, dp.levels - 1))
  alphab <- c(gamma0.beta,  rep(gamma.beta,  dp.levels - 1))

  invisible(list(num.tumor.types = ppindex$num.tumor.types,
                 number.channels = number.channels,
                 convSpectra     = convSpectra,
                 ppindex         = ppindex$ppindex,
                 cpindex         = cpindex,
                 alphaa          = alphaa,
                 alphab          = alphab))
}

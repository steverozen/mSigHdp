#' Initialize hdp object
#' Allocate process index for hdp initialization.
#' Prepare for \code{\link[hdpx]{hdp_init}}
#'
#' @inheritParams Generateppindex
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
#'   parameters (\eqn{\alpha_0} and all \eqn{\alpha_j} in
#'   Figure B.1 of
#'   https://www.repository.cam.ac.uk/bitstream/handle/1810/275454/Roberts-2018-PhD.pdf
#'
#' @param gamma.beta Inverse scale parameter (rate parameter) of
#'   the gamma distribution prior for the Dirichlet process concentration
#'   parameters(\eqn{\beta_0} and all \eqn{\beta_j} in
#'   Figure B.1 of
#'   https://www.repository.cam.ac.uk/bitstream/handle/1810/275454/Roberts-2018-PhD.pdf
#'
#'   We recommend gamma.alpha = 1 and gamma.beta = 20 for single-base-substitution signatures extraction;
#'   gamma.alpha = 1 and gamma.beta = 50  for doublet-base-substitution/INDEL signature extraction
#'
#' @keywords internal

PrepInit <- function(multi.types,
                     input.catalog,
                     verbose      = TRUE,
                     K.guess,
                     gamma.alpha  = 1,
                     gamma.beta   = 1) {

  input.catalog <- GetInputCatalogAsMatrix(input.catalog)

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
  }

  alphaa <- rep(gamma.alpha, dp.levels)
  alphab <- rep(gamma.beta,  dp.levels)

  invisible(list(num.tumor.types = ppindex$num.tumor.types,
                 number.channels = number.channels,
                 convSpectra     = convSpectra,
                 ppindex         = ppindex$ppindex,
                 cpindex         = cpindex,
                 alphaa          = alphaa,
                 alphab          = alphab))
}

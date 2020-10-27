#' test.spectra
#'
#' Synthetic SBS96 spectra for testing.
#'
#' @format An \code{\link[ICAMS]{ICAMS}} catalog (each column is a sample, each row is a mutation type,
#'   e.g. ACT > AGT).
#'
#' @name test.spectra
"test.spectra"


#' test.ground.truth.sig
#'
#' Ground truth signature for testing.
#'
#' @format An \code{\link[ICAMS]{ICAMS}} catalog (each column is a COSMICv3 signature, each row is a mutation type,
#'   e.g. ACT > AGT).
#'
#' @name test.ground.truth.sig
"test.ground.truth.sig"


#' test.ground.truth.exposure
#'
#' Ground truth exposure for testing.
#'
#' @format An exposure matrix (each row is a COSMICv3 signature used to generate \code{test.spectra} , each column is a sample).
#'
#' @name test.ground.truth.sig
"test.ground.truth.exposure"

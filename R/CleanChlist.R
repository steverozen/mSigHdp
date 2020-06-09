#' If the job of Gibbs sampling from \code{MultipleSetupAndPosterior}
#' has an error caught by R, the corresponding element
#' of chlist has class try-error.
#' If the job is stopped with, e.g. a segfault,
#' the \code{chlist} element is NULL.
#'
#'
#'
#'
#' @param chlist A list of \code{\link[hdpx]{hdpSampleChain-class}} objects.
#'
#' @return Invisibly, the clean, non-error \code{chlist}
#'    This is a list of \code{\link[hdpx]{hdpSampleChain-class}} objects.
#'


CleanChlist <- function(chlist, verbose = FALSE) {

  clean.chlist <- list()

  if (!is.list(chlist)) chlist <- list(chlist)

  ii <- 1
  for (i in 1:length(chlist)) {
    cclass <- class(chlist[[i]])
    if (verbose) message("chlist element", i, "has class ", cclass)
    if ("try-error" %in% cclass) {
      warning("class of element", i, "is try-error\n")
    } else {
      if (!("hdpSampleChain" %in% cclass)) {
        warning("An incorrect class for i = ", i, " ", cclass)
      } else{
        clean.chlist[[ii]] <- chlist[[i]]
        ii <- ii + 1
      }
    }
  }
  if (length(clean.chlist) == 0) stop("No usable result in chlist")
  return(clean.chlist)

}

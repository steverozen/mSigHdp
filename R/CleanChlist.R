#' If the job of Gibbs sampling from \code{MultipleSetupAndPosterior}
#' has an error caught by R, the corresponding element
#' of chlist has class try-error.
#' If the job is stopped with, e.g. a segfault,
#' the \code{chlist} element is NULL.
#'
#' @param chlist A list of \code{\link[hdpx]{hdpSampleChain-class}} objects.
#'
#' @param verbose If \code{TRUE} then \code{message} progress information.
#'
#' @return Invisibly, the clean, non-error \code{chlist}
#'    This is a list of \code{\link[hdpx]{hdpSampleChain-class}} objects.
#'
#' @export
#'
#'
CleanChlist <- function(chlist, verbose = FALSE) {

  clean.chlist <- list()

  if (!is.list(chlist)) chlist <- list(chlist)

  ii <- 1
  for (i in 1:length(chlist)) {
    my.error <- chlist[[i]]
    if (is.null(my.error)) {
      warning("Child process for element ", i, " likely crashed R")
    } else {
      cclass <- class(my.error)
      # if (verbose) message("chlist element ", i, " has class ", cclass)
      if ("try-error" %in% cclass) {
        warning("Child process for element ", i,
                " generated this error:\n       ", my.error[1])
      } else {
        if (!("hdpSampleChain" %in% cclass)) {
          warning("Got incorrect class from child process for element ", i, " (", cclass, ")")
        } else{
          clean.chlist[[ii]] <- chlist[[i]]
          ii <- ii + 1
        }
      }
    }
  }
  if (length(clean.chlist) == 0)
    stop("No usable results in from parallel child processes",
         " (from using mclapply)\n",
         "See warnings() for details")
  return(clean.chlist)

}

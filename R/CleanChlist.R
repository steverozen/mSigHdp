# If a child dies with an error caught by R,
# the corresponding element of chlist has class try-error,
# if the child dies with e.g. a segfault the chlist element
# is NULL.
#
# We filter these out and generate a warning. This is a bit
# tricky and I am not sure I have anticipated all possible
# returns, so I do this in a loop.

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

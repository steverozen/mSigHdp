# Generate the original multi_chain for the sample
# if (verbose) message("calling hdp_multi_chain ", Sys.time())
# If a child dies the corresponding element of chlist has
# class try-error.
#
# We filter these and generate a warning. This is a bit
# tricky and I am not sure I have anticipated all possible
# returns, so I do this in a loop.
# @import

CleanChlist <- function(chlist){

  clean.chlist <- list()

  if(length(chlist)>1){
    ii <- 1
    for (i in 1:length(chlist)) {
      cclass <- class(chlist[[i]])
      cat("chlist element", i, "has class ", cclass, "\n")
      if ("try-error" %in% cclass) {
        warning("class of element", i, "is try-error\n")
      } else {
        if (!("hdpSampleChain" %in% cclass)) {
          warning("A different incorrect class for i =", i, cclass)
        } else{
          clean.chlist[[ii]] <- chlist[[i]]
          ii <- ii + 1
        }
      }
    }
    return(clean.chlist)

  }else if(length(chlist)==1){
    cclass <- class(chlist)
    cat("chlist element", 1, "has class ", cclass, "\n")
    if ("try-error" %in% cclass) {
      warning("class of element", 1, "is try-error\n")
    } else {
      if (!("hdpSampleChain" %in% cclass)) {
        warning("A different incorrect class for i =", 1, cclass)
      } else{
        clean.chlist <- chlist
      }
    }
    return(clean.chlist)

  }else if(length(clean.chlist) == 0) {
    stop("No usable result in chlist")
  }


}

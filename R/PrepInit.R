#' Initialize hdp object
#' Allocate process index for hdp initialization.
#' Prepare for hdp_init
#' @param multi.types TODO
#' @param input.catalog TODO
#' @param verbose TODO
#' @param K.guess TODO
#  @import
PrepInit <- function(multi.types,
                      input.catalog,
                      verbose,
                      K.guess){
  if (!exists("stir.closure", envir = .GlobalEnv)) {
    assign("stir.closure", hdpx::xmake.s(), envir = .GlobalEnv)
  }

  # hdp gets confused if the class of its input is not matrix.
  convSpectra <- t(input.catalog)
  # class(convSpectra) <- "matrix"
  # convSpectra <- t(convSpectra)

  number.channels <- nrow(input.catalog)
  number.samples  <- ncol(input.catalog)

  if (verbose) {
    message("Guessed number of signatures ",
            "(= Dirichlet process data clusters) = ", K.guess)
  }


  if (multi.types == FALSE) { # All tumors belong to one tumor type
    num.tumor.types <- 1
    process.index <- c(0,1,rep(2,number.samples))
  } else {
    if (multi.types == TRUE) {
      sample.names <- colnames(input.catalog)
      if (!all(grepl("::", sample.names)))
        stop("Every sample name needs to be of",
             " the form <sample_type>::<sample_id>")

      tumor.types <- sapply(
        sample.names,
        function(x) {strsplit(x, split = "::", fixed = T)[[1]][1]})

      num.tumor.types <- length(unique(tumor.types))
    } else if (is.character(multi.types)) {
      num.tumor.types <- length(unique(multi.types))
      tumor.types <- multi.types
    } else {
      stop("multi.types should be TRUE, FALSE, or a character vector of tumor types")
    }
    # 0 refers to the grandparent Dirichelet process node.
    # There is a level-one node for each tumor type, indicated by a 1.
    process.index <- c(0, rep(1, num.tumor.types))

    # Each tumor type gets its own number.
    process.index <- c(process.index, 1 + as.numeric(as.factor(tumor.types))) # To do, update this with the more transparent code
    cat(process.index, "\n")
    # process.index is now something like
    # c(0, 1, 1, 2, 2, 2, 3, 3)
    # 0 is grandparent
    # 1 is a parent of one type (there are 2 types)
    # 2 indcates tumors of the first type
    # 3 indicates tumors of second type
  }

  ## Specify ppindex as process.index, TODO, why introduce a new variable here?
  ## and cpindex (concentration parameter) as 1 + process.index
  ppindex <- process.index
  cpindex <- 1 + process.index

  ## Calculate the number of levels in the DP node tree.
  dp.levels <- length(unique(ppindex))

  al <- rep(1,dp.levels)
  invisible(list(num.tumor.types       = num.tumor.types,
                 number.channels        = number.channels,
                 convSpectra      = convSpectra,
                 ppindex    = ppindex,
                 cpindex = cpindex,
                 al = al))
}

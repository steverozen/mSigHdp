#' Generate index for a HDP structure and num.tumor.types for other functions
#'
#' @param input.catalog Input spectra catalog as a matrix or
#' in \code{\link[ICAMS]{ICAMS}} format.
#'
#' @param multi.types A logical scalar or
#' a character vector.
#'
#' If \code{FALSE}, The HDP analysis
#'   will regard all input spectra as one tumor type, and the
#'   HDP structure will have one parent node for all tumors.
#'
#' If \code{TRUE}, Sample IDs in \code{input.catalog} must have the
#'   form \emph{sample_type}::\emph{sample_id}.
#'
#' If a character vector, then its length must be \code{ncol(input.catalog}),
#' and each value is the sample type of the corresponding column in \code{input.catalog},
#' e.g. \code{c(rep("Type-A", 23), rep("Type-B", 10))} for 23 Type-A samples and 10 Type-B
#' samples.
#'
#' If not \code{FALSE}, HDP will have one parent node for each sample type
#' and one grandparent node.
#'
#' @export


Generateppindex <- function(multi.types, input.catalog) {

  number.samples <- ncol(input.catalog)

  if (FALSE == multi.types) {
    num.tumor.types <- 0 # For input to hdpx::hdp_setdata
    ppindex <- c(0, rep(1, number.samples))
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

    # 0 refers to the grandparent Dirichlet process node.
    # There is a level-one node for each tumor type, indicated by a 1.
    ppindex <- c(0, rep(1, num.tumor.types))

    # Each tumor type gets its own number.
    tumor.index <- 1

    for(TumorType in unique(tumor.types)){
      tumor.index <- tumor.index+1
      num.sample.each <- sum(tumor.types==TumorType)
      ppindex <- c(ppindex, rep(tumor.index,num.sample.each))
      message(paste0(num.sample.each," samples in ",TumorType))
    }
    # cat(ppindex, "\n")
    # ppindex is now something like
    # c(0, 1, 1, 2, 2, 2, 3, 3)
    # 0 is grandparent
    # 1 is a parent of one type (there are 2 types)
    # 2 indicates tumors of the first type
    # 3 indicates tumors of second type
  }

  return(invisible(list(ppindex          = ppindex,
                        num.tumor.types  = num.tumor.types)))
}

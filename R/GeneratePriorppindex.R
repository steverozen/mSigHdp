#' Generate index for a HDP structure and num.tumor.types for other functions for hdp_prior_init
#'
#' @param input.catalog Input spectra catalog as a matrix or
#' in \code{\link[ICAMS]{ICAMS}} format.
#'
#' @param multi.types A logical scalar or
#' a character vector.
#' If \code{FALSE}, The HDP analysis
#'   will regard all input spectra as one tumor type.
#'
#' If \code{TRUE}, the HDP analysis
#'   will infer tumor types based on the string before "::" in their names.
#' e.g. tumor type for "SA.Syn.Ovary-AdenoCA::S.500" would be "SA.Syn.Ovary-AdenoCA"
#'
#' If \code{multi.types} is a character vector, then it should be of the same length
#' as the number of columns in \code{input.catalog}, and each value is the
#' name of the tumor type of the corresponding column in \code{input.catalog}.
#'
#' e.g. \code{c("SA.Syn.Ovary-AdenoCA", "SA.Syn.Kidney-RCC")}.
#'
#' @param nps Number of prior signatures
#'
#' @export

GeneratePriorppindex <- function(multi.types, input.catalog, nps){
  number.samples <- ncol(input.catalog)

  if (FALSE == multi.types) {
    num.tumor.types <- 0 # For input to hdpx::hdp_setdata
    ppindex <- c(1, rep(1+nps+1, number.samples))
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
    ppindex <- c(1, rep(1+nps+1, num.tumor.types))
    cpindex <- c(3, rep(4,num.tumor.types))

    tumor.index <- 1+nps+1
    cp.index <- 4

    for(TumorType in unique(tumor.types)){
      tumor.index <- tumor.index+1
      cp.index    <- cp.index+1
      num.sample.each <- sum(tumor.types==TumorType)
      ppindex <- c(ppindex, rep(tumor.index,num.sample.each))
      cpindex <- c(cpindex, rep(cp.index,num.sample.each))
      message(paste0(num.sample.each," samples in ",TumorType))
    }
    # Each tumor type gets its own number.

    # cat(ppindex, "\n")
    # ppindex is now something like
    # c(0, 1, 1, 2, 2, 2, 3, 3)
    # 0 is grandparent
    # 1 is a parent of one type (there are 2 types)
    # 2 indicates tumors of the first type
    # 3 indicates tumors of second type
  }

  return(invisible(list(ppindex          = ppindex,
                        cpindex          = cpindex,
                        num.tumor.types  = num.tumor.types)))
}

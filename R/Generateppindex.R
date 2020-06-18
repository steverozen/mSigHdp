#' Generate index for a HDP structure and num.tumor.types for other functions
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
#' @param one.parent.hack IF TRUE, the vector of parents for the Dirichlet
#'   processes looks like c(0, 1, 1, 1, 1, ...), not c(0, 1, 2, 2, 2, ....).
#'
#' @export


Generateppindex <- function(multi.types,
                             one.parent.hack,
                             input.catalog){

  number.samples <- ncol(input.catalog)

  if (one.parent.hack){ # if one.parent.hack is TRUE, ignore multi.types

    num.tumor.types <- 0 #to match with hdp_setdata
    ppindex <- c(0, rep(1, number.samples))

  }else{

    if (multi.types == FALSE) { # All tumors belong to one tumor type. This structure is
      # different from one.parent.hack
      num.tumor.types <- 1
      ppindex <- c(0,1,rep(2,number.samples)) # Parent Dirichlet process index

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
  }
  return(invisible(list(ppindex          = ppindex,
                        num.tumor.types  = num.tumor.types)))
}

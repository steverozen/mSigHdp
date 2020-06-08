#' Run hdp extraction and attribution on a spectra catalog file
#' A function to do burn-in iteration only. This returns a list of hdp object.
#' This needs to be converted to a hdpState object before hdp_posterior
#' (hdpx:::as.hdpState(hdplist))
#'
#' @param input.catalog Input spectra catalog as a matrix or
#' in \code{\link[ICAMS]{ICAMS}} format.
#'
#'
#' @param seedNumber An integer that is used to generate separate
#'   random seeds for each call to \code{\link[hdpx]{dp_activate}},
#'   and each call of \code{\link[hdpx]{hdp_posterior}}; please see the code
#'   on how this is done. But repeated calls with same value of
#'   \code{seedNumber} and other inputs should produce the same results.
#'
#' @param K.guess Suggested initial value of the number of
#' signatures, passed to \code{\link[hdpx]{dp_activate}} as
#' \code{initcc}.
#'
#' @param multi.types A logical scalar or
#' a character vector.
#' If \code{FALSE}, hdp will regard all input spectra as one tumor type.
#'
#' If \code{TRUE}, hdp will infer tumor types based on the string before "::" in their names.
#' e.g. tumor type for "SA.Syn.Ovary-AdenoCA::S.500" would be "SA.Syn.Ovary-AdenoCA"
#'
#' If \code{multi.types} is a character vector, then it should be of the same length
#' as the number of columns in \code{input.catalog}, and each value is the
#' name of the tumor type of the corresponding column in \code{input.catalog},
#' e.g. \code{c("SA.Syn.Ovary-AdenoCA", "SA.Syn.Ovary-AdenoCA", "SA.Syn.Kidney-RCC")}.
#'
#' @param verbose If \code{TRUE} then \code{message} progress information.
#'
#'
#' @param post.burnin Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{burnin}.
#'
#' @param post.cpiter Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{cpiter}.
#'
#' @param post.verbosity Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{verbosity}.
#'
#'
#' @return A list with hdp object after burn-in iteration
#'
#' @export
#'
#'

BurninIteration <-
  function(input.catalog,
           seedNumber          = 1,
           K.guess,
           multi.types         = FALSE,
           verbose             = TRUE,
           post.burnin         = 4000,
           post.cpiter         = 3,
           post.verbosity      = 0

  ) { # 8 arguments

    prep_val <- PrepInit(multi.types = multi.types,
                         input.catalog = input.catalog,
                         verbose       = verbose,
                         K.guess       = K.guess)

    if (verbose) message("calling hdp_init ", Sys.time())
    hdpObject <- hdpx::hdp_init(ppindex = prep_val$ppindex,
                                cpindex = prep_val$cpindex,
                                hh      = rep(1,prep_val$number.channels),
                                alphaa  = prep_val$al,
                                alphab  = prep_val$al)

    # num.process is the number of samples plus number of cancer types plus 1 (grandparent)
    num.process <- hdpx::numdp(hdpObject)

    if (verbose) message("calling hdp_setdata ", Sys.time())

    # (hdp/hdpx)::hdp_setdata generates the warning:
    # In if (!class(data) %in% c("matrix", "data.frame")) { :
    #     the condition has length > 1 and only the first element will be used
    # We circumvent this here
    tmp.cs <- prep_val$convSpectra
    attr(tmp.cs, "class") <- "matrix"
    hdpObject <-
      hdpx::hdp_setdata(hdpObject,
                        (1 + prep_val$num.tumor.types + 1):num.process,
                        tmp.cs)
    rm(tmp.cs)

    seed <- seedNumber
    if (verbose) message("calling dp_activate ", Sys.time())
    # dp_activate requires that stir.closure exists in .GlobalEnv;
    # see above in this function.
    hdp.state <- hdpx::dp_activate(hdpObject,
                                   1:num.process,
                                   initcc = K.guess,
                                   seed = seed + 3e6)

    hdplist <- hdpx::as.list(hdp.state)
    iterate <- utils::getFromNamespace(x = "iterate", ns = "hdpx")
    output <- hdpx:::iterate(hdplist, post.burnin, post.cpiter, post.verbosity)##burn-in first, then return the hdplist after burnt in.
    hdplist <- output[[1]]
    return(hdplist)

  }

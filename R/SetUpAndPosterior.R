#' Generate HDP gibs sampling chains from a spectra catalog.
#'
#' @param input.catalog Input spectra catalog as a matrix or
#' in \code{\link[ICAMS]{ICAMS}} format.
#'
#' @param CPU.cores Number of CPUs to use in running
#'    \code{\link[hdpx]{hdp_posterior}}; this is used to parallelize
#'    running the posterior sampling chains, so there is no
#'    point in making this larger than \code{num.posterior}.
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
#' @param num.posterior Number of posterior sampling chains; can set to
#'   1 for testing.
#'
#' @param post.burnin Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{burnin}.
#'
#' @param post.n Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{n}.
#'
#' @param post.space Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{space}.
#'
#' @param post.cpiter Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{cpiter}.
#'
#' @param post.verbosity Pass to \code{\link[hdpx]{hdp_posterior}}
#'      \code{verbosity}.
#'
#' @param checkpoint.aft.post If non-\code{NULL}, a file path to checkpoint
#'      the list of values returned from the calls to \code{\link[hdpx]{hdp_posterior}}
#'      as a .Rdata file.
#'
#' @param stop.after.hdp.posterior If non-\code{NULL}, then
#'      a file path to checkpoint
#'      the return value from the call to \code{\link[hdpx]{hdp_posterior}}
#'      as a .Rdata file. The function will then also invisibly return
#'      this value.
#'
#' @return Invisibly,
#'    the clean
#'    \code{chlist} (output of the hdp_posterior calls).
#'
#' @export

SetUpAndPosterior <-
  function(input.catalog,
           CPU.cores           = 1,
           seedNumber          = 1,
           K.guess,
           multi.types         = FALSE,
           verbose             = TRUE,
           num.posterior       = 4,
           post.burnin         = 4000,
           post.n              = 50,
           post.space          = 50,
           post.cpiter         = 3,
           post.verbosity      = 0,
           checkpoint.aft.post = NULL,
           stop.after.hdp.posterior = NULL
) { # 14 arguments

    prep_val <- prep_init(multi.types,input.catalog)

    if (verbose) message("calling hdp_init ", Sys.time())
    hdpObject <- hdpx::hdp_init(ppindex = prep_val$ppindex,
                               cpindex = prep_val$cpindex,
                               hh = rep(1,prep_val$number.channels),
                               alphaa = prep_val$al,
                               alphab = prep_val$al)

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

    # Run num.posterior independent sampling chains
    activate.and.sample <- function(my.seed) {

      if (verbose) message("calling dp_activate ", Sys.time())
      # dp_activate requires that stir.closure exists in .GlobalEnv;
      # see above in this function.
      hdp.state <- hdpx::dp_activate(hdpObject,
                                    1:num.process,
                                    initcc = K.guess,
                                    seed = my.seed + 3e6)

      if (verbose) message("calling hdp_posterior ", Sys.time())
      sample.chain <- hdpx::hdp_posterior(
        hdp       = hdp.state,
        verbosity = post.verbosity,
        burnin    = post.burnin,
        n         = post.n,
        space     = post.space,
        cpiter    = post.cpiter,
        seed      = my.seed)
      return(sample.chain)
    }

    parallel.time <- system.time(
      chlist <- activate.and.sample(seedNumber)
    )
    if (verbose) {
      message("compute chlist time: ")
      for (xn in names(parallel.time)) {
        message(" ", xn, " ", parallel.time[[xn]])
      }
    }



    return(invisible(chlist))

  }

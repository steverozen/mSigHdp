#' Generate an HDP Gibbs sampling chain from a spectra catalog.
#'
#' @inheritParams PrepInit
#'
#' @param seedNumber A random seeds passed to \code{\link[hdpx]{dp_activate}}.
#'
#' @return Invisibly, an \code{\link[hdpx]{hdpState-class}} object
#'  as returned from \code{\link[hdpx]{dp_activate}}.
#'
#' @export

SetupAndActivate <- function(input.catalog,
                             seedNumber          = 1,
                             K.guess,
                             multi.types         = FALSE,
                             verbose             = TRUE,
                             gamma.alpha         = 1,
                             gamma.beta          = 1,
                             gamma0.alpha        = gamma.alpha,
                             gamma0.beta         = gamma.beta)
{

    prep_val <- PrepInit(multi.types = multi.types,
                         input.catalog = input.catalog,
                         verbose       = verbose,
                         K.guess       = K.guess,
                         gamma.alpha   = gamma.alpha,
                         gamma.beta    = gamma.beta,
                         gamma0.alpha  = gamma0.alpha,
                         gamma0.beta   = gamma0.beta)

    if (verbose) {
      message("SetupAndActivate gamma0.alpha = ",
              gamma0.alpha, " gamma0.beta = ",
              gamma0.beta)
    }
    if (verbose) message("calling hdp_init ", Sys.time())
    hdpObject <- hdpx::hdp_init(ppindex = prep_val$ppindex,
                                cpindex = prep_val$cpindex,
                                hh      = rep(1,prep_val$number.channels),
                                alphaa  = prep_val$alphaa,
                                alphab  = prep_val$alphab)

    # num.process is the number of samples plus number of cancer types plus 1 (grandparent)
    num.process <- hdpx::numdp(hdpObject)

    if (verbose) message("calling hdp_setdata ", Sys.time())

    hdpObject <-
      hdpx::hdp_setdata(hdpObject,
                        (1 + prep_val$num.tumor.types + 1):num.process,
                        prep_val$convSpectra)

    if (verbose) message("calling dp_activate ", Sys.time())

    hdp.state <- hdpx::dp_activate(hdp     = hdpObject,
                                   dpindex = 1:num.process,
                                   initcc  = K.guess,
                                   seed    = seedNumber)

    return(invisible(hdp.state))

  }

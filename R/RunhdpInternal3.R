#' Run hdp extraction and attribution on a spectra catalog file
#'
#' @param input.catalog A catalog of spectra catalog
#' in \code{\link[ICAMS]{ICAMS}} format.
#'
#' @param CPU.cores Number of CPUs to use in running
#'    \code{\link[hdp]{hdp_posterior}}; this is used to parallize
#'    running the posterior sampling chains, so there is no
#'    point in making this larger than \code{num.posterior}.
#'
#' @param seedNumber An integer that is used to generate separate
#'   random seeds for each call to \code{\link[hdp]{dp_activate}},
#'   and each call of \code{\link[hdp]{hdp_posterior}}; please see the code
#'   on how this is done. But repeated calls with same value of
#'   \code{seedNumber} and other inputs should produce the same results.
#'
#' @param K.guess Suggested initial value of the number of
#' signatures, passed to \code{\link[hdp]{dp_activate}} as
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
#' @param post.burnin Pass to \code{\link[hdp]{hdp_posterior}}
#'      \code{burnin}.
#'
#' @param post.n Pass to \code{\link[hdp]{hdp_posterior}}
#'      \code{n}.
#'
#' @param post.space Pass to \code{\link[hdp]{hdp_posterior}}
#'      \code{space}.
#'
#' @param post.cpiter Pass to \code{\link[hdp]{hdp_posterior}}
#'      \code{cpiter}.
#'
#' @param post.verbosity Pass to \code{\link[hdp]{hdp_posterior}}
#'      \code{verbosity}.
#'
#' @param cos.merge The cosine similarity threshhold for merging raw clusters
#'      from the posterior sampling chains into "components" i.e. signatures;
#'      passed to \code{\link[hdp]{hdp_extract_components}}.
#'
#' @param min.sample A "component" (i.e. signature) must have at least
#'      this many samples; passed to \code{\link[hdp]{hdp_extract_components}}.
#'
#' @return A list with the following elements:\describe{
#' \item{signature}{The extracted signature profiles as a matrix;
#'             rows are mutation types, columns are
#'             samples (e.g. tumors).}
#' \item{exposure}{The inferred exposures as a matrix of mutation counts;
#'            rows are signatures, columns are samples (e.g. tumors).}
#' \item{exposure.p}{\code{exposure} converted to proportions.}
#'
#' \item{multi.chains}{A \code{\link[hdp]{hdpSampleMulti-class}} object.
#'     This object has the method \code{\link[hdp]{chains}} which returns
#'     a list of \code{\link[hdp]{hdpSampleChain-class}} objects. Each of these
#'     sample chains objects has a method \code{\link[hdp]{final_hdpState}}
#'     (actually the methods seems to be just \code{hdp})
#'     that returns the \code{hdpState} from which it was generated.}
#' }
#'
#' @export

RunhdpInternal3 <-
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
           cos.merge           = 0.9,
           min.sample          = 1
) { # 14 arguments

    if (!exists("stir.closure", envir = .GlobalEnv)) {
      assign("stir.closure", hdp::make.stirling(), envir = .GlobalEnv)
    }

    # hdp gets confused if the class of its input is not matrix.
    convSpectra <- input.catalog
    class(convSpectra) <- "matrix"
    convSpectra <- t(convSpectra)

    number.channels <- nrow(input.catalog)
    number.samples  <- ncol(input.catalog)

    if (verbose) {
      message("Guessed number of signatures ",
              "(= Dirichlet process data clusters) = ", K.guess)
    }

    # Initialize hdp object
    # Allocate process index for hdp initialization.

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
      process.index <- c(process.index, 1 + as.numeric(as.factor(tumor.types)))
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

    ## initialise hdp
    if (verbose) message("calling hdp_init")
    hdpObject <- hdp::hdp_init(ppindex = ppindex,
                               cpindex = cpindex,
                               hh = rep(1,number.channels),
                               alphaa = al,
                               alphab = al)

    # num.process is the number of samples plus number of cancer types plus 1 (grandparent)
    num.process <- hdp::numdp(hdpObject)

    if (verbose) message("calling hdp_setdata")
    hdpObject <- hdp::hdp_setdata(
      hdpObject,
      (1 + num.tumor.types + 1):num.process,
      convSpectra)

    # Run num.posterior independent sampling chains
    activate.and.sample <- function(my.seed) {

      if (verbose) message("calling dp_activate")
      # dp_activate requires that stir.closure exists in .GlobalEnv;
      # see above in this function.
      hdp.state <- hdp::dp_activate(hdpObject,
                                    1:num.process,
                                    initcc = K.guess,
                                    seed = my.seed + 3e6)

      if (verbose) message("calling hdp_posterior")
      sample.chain <- hdp::hdp_posterior(
        hdp       = hdp.state,
        verbosity = post.verbosity,
        burnin    = post.burnin,
        n         = post.n,
        space     = post.space,
        cpiter    = post.cpiter,
        seed      = my.seed)
      return(sample.chain)
    }

    chlist <- parallel::mclapply(
      # Must choose a different seed for each of the chains
      X = (seedNumber + 1:num.posterior * 10^6) ,
      FUN = activate.and.sample,
      mc.cores = CPU.cores)

    # Generate the original multi_chain for the sample
    if (verbose) message("calling hdp_multi_chain")
    multi.chains <- hdp::hdp_multi_chain(chlist)

    if (verbose) message("calling hdp_extract_components")
    # Group raw "clusters" into "components" (i.e. signatures).
    multi.chains <-
      hdp::hdp_extract_components(multi.chains,
                                  cos.merge  = cos.merge,
                                  min.sample = min.sample)

    if (verbose) message("calling hdp::comp_categ_distn")
    extractedSignatures <- t(hdp::comp_categ_distn(multi.chains)$mean)

    rownames(extractedSignatures) <- rownames(input.catalog)
    # Set signature names to "hdp.0","hdp.1","hdp.2", ...
    colnames(extractedSignatures) <-
      paste("hdp", colnames(extractedSignatures), sep = ".")

    ## Calculate the exposure probability of each signature (component) for each
    ## tumor sample (posterior sample corresponding to a dirichlet process node).
    ## This is the probability distribution of signatures (components) for all
    ## tumor samples (DP nodes); exposureProbs is the normalized
    ## signature exposure all tumor samples # TODO Wuyang, what do you mean
    # by normalize?

    if (verbose) message("Calling hdp::comp_dp_distn")
    exposureProbs <- hdp::comp_dp_distn(multi.chains)$mean

    # Remove columns corresponding to parent or grandparent nodes
    # (leaving only columns corresponding to samples.
    # Transpose so it conforms to SynSigEval format
    exposureProbs <- t(exposureProbs[-(1:(num.tumor.types + 1)), ])
    # Now rows are signatures, columns are samples

    # Calculate exposure counts from exposure probabilies and total mutation
    # counts
    exposureCounts <- exposureProbs %*% diag(rowSums(convSpectra))
    colnames(exposureCounts) <- colnames(input.catalog)
    rownames(exposureCounts) <- colnames(extractedSignatures)

    invisible(list(signature       = extractedSignatures,
                   exposure        = exposureCounts,
                   exposure.p      = exposureProbs,
                   multi.chains    = multi.chains))
  }

#' Run hdp extraction and attribution on a spectra catalog file
#' This repeats what Nicola do in her thesis. It starts four independent initial chains with post.burnin iterations,
#' then pick up from the end of each initial chain,
#' started another num.posterior MCMC chains for another post.burnin iterations
#' and then collected post.n posterior samples at intervals of post.space iterations.
#' In total, this collects 4 times num.posterior times post.n posterior samples from 4 times y separate chains.
#'
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
#' @param cos.merge The cosine similarity threshold for merging raw clusters
#'      from the posterior sampling chains into "components" i.e. signatures;
#'      passed to \code{\link[hdpx]{hdp_extract_components}}.
#'
#' @param min.sample A "component" (i.e. signature) must have at least
#'      this many samples; passed to \code{\link[hdpx]{hdp_extract_components}}.
#'
#' @param checkpoint.aft.post If non-\code{NULL}, a file path to checkpoint
#'      the list of values returned from the calls to \code{\link[hdpx]{hdp_posterior}}
#'      as a .Rdata file.
#'
#' @return A list with the following elements:\describe{
#' \item{signature}{The extracted signature profiles as a matrix;
#'             rows are mutation types, columns are
#'             samples (e.g. tumors).}
#' \item{exposure}{The inferred exposures as a matrix of mutation counts;
#'            rows are signatures, columns are samples (e.g. tumors).}
#' \item{exposure.p}{\code{exposure} converted to proportions.}
#'
#' \item{multi.chains}{A \code{\link[hdpx]{hdpSampleMulti-class}} object.
#'     This object has the method \code{\link[hdpx]{chains}} which returns
#'     a list of \code{\link[hdpx]{hdpSampleChain-class}} objects. Each of these
#'     sample chains objects has a method \code{\link[hdpx]{final_hdpState}}
#'     (actually the methods seems to be just \code{hdp})
#'     that returns the \code{hdpState} from which it was generated.}
#' }
#'
#' @export
#'
#'

ExtendIterationAndPosterior <-
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
           min.sample          = 1,
           checkpoint.aft.post = NULL
  ) { # 15 arguments

    # if (!exists("stir.closure", envir = .GlobalEnv)) {
    #  assign("stir.closure", hdpx::xmake.s(), envir = .GlobalEnv)
    # }

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



    if (verbose) message("calling hdp_init ", Sys.time())
    hdpObject <- hdpx::hdp_init(ppindex = ppindex,
                                cpindex = cpindex,
                                hh = rep(1,number.channels),
                                alphaa = al,
                                alphab = al)

    # num.process is the number of samples plus number of cancer types plus 1 (grandparent)
    num.process <- hdpx::numdp(hdpObject)

    if (verbose) message("calling hdp_setdata ", Sys.time())

    # (hdp/hdpx)::hdp_setdata generates the warning:
    # In if (!class(data) %in% c("matrix", "data.frame")) { :
    #     the condition has length > 1 and only the first element will be used
    # We circumvent this here

    tmp.cs <- convSpectra
    attr(tmp.cs, "class") <- "matrix"
    hdpObject <-
      hdpx::hdp_setdata(hdpObject,
                        (1 + num.tumor.types + 1):num.process,
                        tmp.cs)
    rm(tmp.cs)


    # Run num.posterior independent sampling chains for burned-in hdp
    hdp_posterior_sample <- function(my.seed) {

      if (verbose) message("calling hdp_posterior ", Sys.time())
      print(my.seed)
      sample.chain <- hdpx::hdp_posterior(
        hdp       = hdp.state.burned,
        verbosity = post.verbosity,
        burnin    = post.burnin,
        n         = post.n,
        space     = post.space,
        cpiter    = post.cpiter,
        seed      = my.seed)

      return(sample.chain)
    }

    chlist <- {}
    for(i in 1:4){
      print(paste0("init_chain",i))
      seed <- seedNumber + i
      if (verbose) message("calling dp_activate ", Sys.time())
      # dp_activate requires that stir.closure exists in .GlobalEnv;
      # see above in this function.
      hdp.state <- hdpx::dp_activate(hdpObject,
                                     1:num.process,
                                     initcc = K.guess,
                                     seed = seed + 3e6)

      hdplist <- hdpx::as.list(hdp.state)
      iterate <- utils::getFromNamespace(x = "iterate", ns = "hdpx")
      output <- iterate(hdplist, post.burnin, post.cpiter, post.verbosity)##burn-in first, then return the hdplist after burnt in.
      hdplist <- output[[1]]
      as.hdpState <- utils::getFromNamespace(x = "as.hdpState", ns = "hdpx")
      hdp.state.burned <- as.hdpState(hdplist)

      parallel.time <- system.time(
        chlist <- c(chlist,parallel::mclapply(
          # Must choose a different seed for each of the chains
          X = (seed + 1:num.posterior * 10^6) ,
          FUN = hdp_posterior_sample,
          mc.cores = CPU.cores))

      )
    }
    if (verbose) {
      message("compute chlist time: ")
      for (xn in names(parallel.time)) {
        message(" ", xn, " ", parallel.time[[xn]])
      }
    }

    if (!is.null(checkpoint.aft.post)) {
      save(chlist, file = checkpoint.aft.post)
    }

    # Generate the original multi_chain for the sample
    if (verbose) message("calling hdp_multi_chain ", Sys.time())
    # If a child dies the corresponding element of chlist has
    # class try-error.
    #
    # We filter these and generate a warning. This is a bit
    # tricky and I am not sure I have anticipated all possible
    # returns, so I do this in a loop.
    ii <- 1
    clean.chlist <- list()
    for (i in 1:length(chlist)) {
      cclass <- class(chlist[[i]])
      cat("chlist element", i, "has class", cclass, "\n")
      if ("try-error" %in% cclass) {
        warning("class of element", i, "\n")
      } else {
        if (!("hdpSampleChain" %in% cclass)) {
          warning("A different incorrect class for i =", i, cclass)
        } else{
          clean.chlist[[ii]] <- chlist[[i]]
          ii <- ii + 1
        }
      }
    }

    if (length(clean.chlist) == 0) {
      fname <- "chlist.from.aborted.run.of.RunhdpInternal4.Rdata"
      save(chlist, file = fname)
      stop("No usable result in chlist, look in ", fname)
    }

    multi.chains <- hdpx::hdp_multi_chain(clean.chlist)
    rm(chlist)
    rm(clean.chlist)

    if (verbose) message("calling hdp_extract_components ", Sys.time())
    # Group raw "clusters" into "components" (i.e. signatures).
    extract.time <- system.time(
      multi.chains <-
        hdpx::hdp_extract_components(multi.chains,
                                     cos.merge  = cos.merge,
                                     min.sample = min.sample)
    )
    if (verbose) {
      message("hdp_extract_components time: ")
      for (xn in names(extract.time)) {
        message(" ", xn, " ", extract.time[[xn]])
      }
    }

    if (verbose) message("calling hdpx::comp_categ_distn ", Sys.time())
    extractedSignatures <- t(hdpx::comp_categ_distn(multi.chains)$mean)

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

    if (verbose) message("Calling hdpx::comp_dp_distn ", Sys.time())
    exposureProbs <- hdpx::comp_dp_distn(multi.chains)$mean

    # Remove columns corresponding to parent or grandparent nodes
    # (leaving only columns corresponding to samples.
    # Transpose so it conforms to SynSigEval format
    exposureProbs <- t(exposureProbs[-(1:(num.tumor.types + 1)), ])
    # Now rows are signatures, columns are samples

    # Calculate exposure counts from exposure probabilities and total mutation
    # counts
    exposureCounts <- exposureProbs %*% diag(rowSums(convSpectra))
    colnames(exposureCounts) <- colnames(input.catalog)
    rownames(exposureCounts) <- colnames(extractedSignatures)

    invisible(list(signature       = extractedSignatures,
                   exposure        = exposureCounts,
                   exposure.p      = exposureProbs,
                   multi.chains    = multi.chains))
  }

# Test status
# devtools::test(filter = "Runhdp3-fast96")


#' Run hdp extraction and attribution on a spectra catalog file
#'
#' @inheritParams RunhdpInternal4
#'
#' @param input.catalog.file File containing a spectra catalog
#' in \code{\link[ICAMS]{ICAMS}} format.
#'
#' @param out.dir Directory that will be created for the output;
#'   if \code{overwrite} is \code{FALSE} then
#'   abort if \code{out.dir} already exits.
#'
#' @param remove.noise Deprecated; ignored
#'
#' @param test.only If > 0, only analyze the first \code{test.only} columns
#'  in \code{input.catalog.file}.
#'
#' @return The same list as returned by \code{\link{RunhdpInternal}}.
#'
#' @details Creates several files in \code{out.dir}. These are:
#'  TODO(Steve): list the files
#'
#' @export

Runhdp4 <-
  function(input.catalog.file,
           out.dir,
           CPU.cores           = 1,
           seedNumber          = 1,
           K.guess,
           multi.types         = FALSE,
           remove.noise        = FALSE,
           test.only           = 0,
           overwrite           = FALSE,
           verbose             = TRUE,
           num.posterior       = 4,
           post.burnin         = 4000,
           post.n              = 50,
           post.space          = 50,
           post.cpiter         = 3,
           post.verbosity      = 0,
           cos.merge           = 0.9,
           min.sample          = 1) {

    if (verbose) message("Reading input catalog file ", input.catalog.file)
    spectra <- ICAMS::ReadCatalog(input.catalog.file, strict = FALSE)
    if (test.only > 0) spectra <- spectra[ , 1:test.only]

    ## Create output directory
    if (dir.exists(out.dir)) {
      if (!overwrite) stop(out.dir, " already exits")
      if (verbose) message("Using existing out.dir ", out.dir)
    } else {
      dir.create(out.dir, recursive = T)
      if (verbose) message("Created new out.dir", out.dir)
    }

    utils::capture.output(date(), cat("\n\n"),
                          match.call(), cat("\n\n"),
                          sessionInfo(),
                          file = file.path(out.dir, "call.and.session.info.txt"))

    retval <- mSigHdp::RunhdpInternal4(
      input.catalog   = spectra,
      CPU.cores       = CPU.cores,
      seedNumber      = seedNumber,
      K.guess         = K.guess,
      multi.types     = multi.types,
      num.posterior   = num.posterior,
      verbose         = verbose,
      post.burnin     = post.burnin,
      post.n          = post.n,
      post.space      = post.space,
      post.cpiter     = post.cpiter,
      post.verbosity  = post.verbosity,
      cos.merge       = cos.merge,
      min.sample      = min.sample
    ) # 14 Arguments

    save(retval, file = file.path(out.dir, "Runhdp3.retval.Rdata"))

    multi <- retval[["multi.chains"]] # class hdpSampleMulti
    chains <- hdpx::chains(multi)      # list of hdpSampleChain

    # Plot the diagnostics of sampling chains.
    if (verbose) message("Writing hdp.diagnostics.pdf")
    grDevices::pdf(file = paste0(out.dir,"/hdp.diagnostics.pdf"))
    graphics::par(mfrow=c(2,2), mar=c(4, 4, 2, 1))

    # This is the likelihood plot along each chain
    lapply(chains, hdpx::plot_lik, bty = "L")

    # This is the number of raw clusters sampled along each chain
    lapply(chains, hdpx::plot_numcluster, bty = "L")

    # This is the number of mutations assigned as a function of
    # the number of raw clusters
    lapply(chains, hdpx::plot_data_assigned, bty = "L")

    graphics::par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
    # Components were already extracted, so this call will work
    hdpx::plot_comp_size(multi, bty="L")

    graphics::par(mfrow=c(8, 1), mar = c(1, 1, 1, 1))
    # This plots the component (signature) profiles with
    # 95% credibility intervals
    hdpx::plot_comp_distn(multi)

    # TODO, need argument dpinices and col_comp;
    # Need to return the hdp object (perhaps) from RunhdpInternal
    # to get the required values.
    if (FALSE) { # Not finished
    num.dpindices <- length(chains[[1]]@hdp@ppindex)
    hdpx::plot_dp_comp_exposure(
      multi, dpindices = 3:num.dpindices,
      col_comp = myCol[1:ncol(retval$signature)],
      dpnames = colnames(retval$exposure))
    }

    grDevices::dev.off()

    if (verbose) message("Writing signatures")
    extractedSignatures <- ICAMS::as.catalog(retval$signature,
                                             region       = "unknown",
                                             catalog.type = "counts.signature")
    ICAMS::WriteCatalog(extractedSignatures,
                        paste0(out.dir,"/extracted.signatures.csv"))

    if (verbose) message("Writing exposures")
    SynSigGen::WriteExposure(retval$exposure.p,
                  paste0(out.dir,"/exposure.probs.csv"))
    SynSigGen::WriteExposure(retval$exposure,
                  paste0(out.dir,"/inferred.exposures.csv"))

    invisible(retval)
  }

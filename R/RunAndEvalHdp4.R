#' Run and evaluate hdp
#'
#' @inheritParams Runhdp4
#'
#' @param ground.truth.exposure.file Path to file with ground truth exposures.
#'
#' @param ground.truth.sig.file Path to file with ground truth signatures.
#'
#' @param ground.truth.sig.catalog \code{\link[ICAMS]{ICAMS}} catalog with signatures
#'   used to construct the ground truth spectra.  Specify only one of
#'   \code{ground.truth.sig.file.path} or \code{ground.truth.sig.catalog}.
#'
#' @export

RunAndEvalHdp4 <- function(
  input.catalog.file,
  ground.truth.exposure.file,
  ground.truth.sig.file      = NULL,
  ground.truth.sig.catalog   = NULL,
  out.dir,
  CPU.cores                  = 1,
  seedNumber                 = 1,
  K.guess,
  multi.types                = FALSE,
  remove.noise               = FALSE,
  test.only                  = 0,
  overwrite                  = FALSE,
  verbose                    = TRUE,
  num.posterior              = 4,
  post.burnin                = 4000,
  post.n                     = 50,
  post.space                 = 50,
  post.cpiter                = 3,
  post.verbosity             = 0,
  cos.merge                  = 0.9,
  min.sample                 = 1) {

  if (!is.null(ground.truth.sig.catalog)
      && !is.null(ground.truth.sig.file)) {
    stop("Specify only one of ground.truth.sig.catalog or ",
         "ground.truth.sig.file.path")
  }

  if (is.null(ground.truth.sig.catalog)) {
    ground.truth.sig.catalog <- ICAMS::ReadCatalog(ground.truth.sig.file)
  }
  stopifnot(is.matrix(ground.truth.sig.catalog))

  # Do this early to catch any possible error before we do a lot
  # of computation
  ground.truth.exp <- SynSigGen::ReadExposure(ground.truth.exposure.file)

  retval <- Runhdp4(
    input.catalog.file         = input.catalog.file,
    out.dir                    = out.dir,
    CPU.cores                  = CPU.cores,
    seedNumber                 = seedNumber,
    K.guess                    = K.guess,
    multi.types                = multi.types,
    test.only                  = test.only,
    overwrite                  = overwrite,
    verbose                    = verbose,
    num.posterior              = num.posterior,
    post.burnin                = post.burnin,
    post.n                     = post.n,
    post.space                 = post.space,
    post.cpiter                = post.cpiter,
    post.verbosity             = post.verbosity,
    cos.merge                  = cos.merge,
    min.sample                 = min.sample)

  sig.without.0 <- retval$signature
  hdp.0.col <- which(colnames(sig.without.0) == "hdp.0")
  stopifnot(length(hdp.0.col) == 1)
  sig.without.0 <- sig.without.0[ , -hdp.0.col]
  ICAMS::WriteCatalog(ICAMS::as.catalog(sig.without.0),
                      file.path(out.dir,"extracted.signatures.no.0.csv"))

  sigAnalysis <- SynSigEval::MatchSigsAndRelabel(
    ex.sigs  = sig.without.0,
    gt.sigs  = ground.truth.sig.catalog,
    exposure = ground.truth.exp)

  ICAMS::PlotCatalogToPdf(ICAMS::as.catalog(sigAnalysis$gt.sigs,
                                            catalog.type = "counts.signature"), # Need to fix this
                          file.path(out.dir, "ground.truth.sigs.pdf"))

  ICAMS::PlotCatalogToPdf(ICAMS::as.catalog(sigAnalysis$ex.sigs,
                                            catalog.type = "counts.signature"),
                          file.path(out.dir, "extracted.sigs.pdf"))

  # Writes bi-directional matching and cos.sim calculation
  utils::write.csv(sigAnalysis$match1, file = file.path(out.dir, "match1.csv"))
  utils::write.csv(sigAnalysis$match2, file = file.path(out.dir, "match2.csv"))

  utils::capture.output(
    cat("Call\n"),
    match.call(),
    cat("\nAverage cosine similarity\n"),
    sigAnalysis$averCosSim,
    cat("\nAverage cosine similarity to each ground-truth signature\n"),
    sigAnalysis$cosSim,
    cat("\nNumber of ground-truth signatures\n"),
    ncol(sigAnalysis$gt.sigs),
    cat("\nNumber of extracted signatures\n"),
    ncol(sigAnalysis$ex.sigs),
    cat("\nsigAnalysis$extracted.with.no.best.match\n"),
    sigAnalysis$extracted.with.no.best.match,
    cat("\nsigAnalysis$ground.truth.with.no.best.match\n"),
    sigAnalysis$ground.truth.with.no.best.match,
    file = file.path(out.dir,"other.results.txt"))

  if (FALSE) {
    # See file AssessExtracted.R
    exposureDiff <- ReadAndAnalyzeExposures(
      extracted.sigs = NULL,
      ground.truth.sigs = NULL ,
      inferred.exp.path = NULL,
      ground.truth.exposures = NULL)
  }

  invisible(retval)

}

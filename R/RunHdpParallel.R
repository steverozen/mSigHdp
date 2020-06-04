#' A caller function combines MultipleSetupandPosterior, CombinePosteriorChains and AnalyzeAndPlotretval
#'
#' @inheritParams MultipleSetupAndPosterior
#' @inheritParams AnalyzeAndPlotretval
#' @param out.dir Directory that will be created for the output
#' @return  If \code{out.dir} is not NULL, output including data and plots
#'          will be saved in \code{out.dir}.
#'          Else,invisibly, a list with the following elements:\describe{
#' \item{signature}{The extracted signature profiles as a matrix;
#'             rows are mutation types, columns are
#'             samples (e.g. tumors).}
#' \item{exposure}{The inferred exposures as a matrix of mutation counts;
#'            rows are signatures, columns are samples (e.g. tumors).}

#'
#' \item{multi.chains}{A \code{\link[hdpx]{hdpSampleMulti-class}} object.
#'     This object has the method \code{\link[hdpx]{chains}} which returns
#'     a list of \code{\link[hdpx]{hdpSampleChain-class}} objects. Each of these
#'     sample chains objects has a method \code{\link[hdpx]{final_hdpState}}
#'     (actually the methods seems to be just \code{hdp})
#'     that returns the \code{hdpState} from which it was generated.}
#' }
#'
#'
#' @export

RunHdpParallel <- function(input.catalog,
                           seedNumber          = 1,
                           K.guess,
                           multi.types         = FALSE,
                           verbose             = TRUE,
                           post.burnin         = 4000,
                           post.n              = 50,
                           post.space          = 50,
                           post.cpiter         = 3,
                           post.verbosity      = 0,
                           CPU.cores           = 1,
                           num.child.process   = 4,
                           ground.truth.sig,
                           ground.truth.exp,
                           out.dir             =NULL
){

  ##Step 1: Activate hierarchical Dirichlet processes and run posterior sampling in parallel.
  ##This returns a list of multiple posterior chains

  multi.chlist <- MultipleSetupAndPosterior(input.catalog,
                                            seedNumber          = 1,
                                            K.guess,
                                            multi.types         = FALSE,
                                            verbose             = TRUE,
                                            post.burnin         = 4000,
                                            post.n              = 50,
                                            post.space          = 50,
                                            post.cpiter         = 3,
                                            post.verbosity      = 0,
                                            CPU.cores           = 1,
                                            num.child.process   = 4)

  ##Step 2: Combine the multiple posterior chains and extract signatures and exposures
  ##This returns signatures and exposures

  retval <- CombinePosteriorChains(multi.chlist,
                                   input.catalog=input.catalog,
                                   multi.types=multi.types)
  ##Step 3: Plot diagnostic plots, signatures, exposurs and compare with ground truth signature and exposures
  if(!is.null(out.dir)){
    AnalyzeAndPlotretval(retval,
                         out.dir,
                         ground.truth.sig,
                         ground.truth.exp,
                         verbose,
                         overwrite)
  }else{
    return(invisible(retval))
  }

}

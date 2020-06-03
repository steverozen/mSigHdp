#' Run hdp extraction and attribution on a spectra catalog file
#' This assign one dp_activate and one hdp_posterior to one node, and return a list,
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
# @param cos.merge The cosine similarity threshold for merging raw clusters
#      from the posterior sampling chains into "components" i.e. signatures;
#      passed to \code{\link[hdpx]{hdp_extract_components}}.
#'
# @param min.sample A "component" (i.e. signature) must have at least
#      this many samples; passed to \code{\link[hdpx]{hdp_extract_components}}.
#
# @param checkpoint.aft.post If non-\code{NULL}, a file path to checkpoint
#      the list of values returned from the calls to \code{\link[hdpx]{hdp_posterior}}
#      as a .Rdata file.
#'
#' @param stop.after.hdp.posterior If non-\code{NULL}, then
#'      a file path to checkpoint
#'      the return value from the call to \code{\link[hdpx]{hdp_posterior}}
#'      as a .Rdata file. The function will then also invisibly return
#'      this value.
#'
#' @return If \code{stop.after.hdp.posterior} is not NULL, then invisibly,
#'    the clean
#'    \code{chlist} (output of the hdp_posterior calls).
#'    Otherwise, invisibly, a list with the following elements:\describe{
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


CallSetupAndPosterior <- function(input.catalog,
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
                                  stop.after.hdp.posterior = NULL){

  run.setup.and.posterior <- function(my.seed) {
    print(my.seed)
    temp.chlist <-SetUpAndPosterior(
      input.catalog,
      CPU.cores           = 1,
      seedNumber          = my.seed,
      K.guess = K.guess,
      multi.types         = FALSE,
      verbose             = TRUE,
      num.posterior       = 4,
      post.burnin         = 4000,
      post.n              = 50,
      post.space          = 50,
      post.cpiter         = 3,
      post.verbosity      = 0,
      checkpoint.aft.post = NULL,
      stop.after.hdp.posterior = NULL)
    return(temp.chlist)
  }


  chlist <- parallel::mclapply(
    # Must choose a different seed for each of the chains
    X = (seedNumber + 1:num.posterior * 10^6) ,
    FUN = run.setup.and.posterior,
    mc.cores = CPU.cores)

  clean.chlist <- CleanChlist(chlist)

  if (!is.null(stop.after.hdp.posterior)) { ##if the stop.after.hdp.posterior is set, then save the chlist without proceeding
    save(clean.chlist, file = stop.after.hdp.posterior)
  }else{
    retval <- CombinePosteriorChains(clean.chlist=clean.chlist,
                                     input.catalog=input.catalog,
                                     multi.types=multi.types)
    return(invisible(retval))
  }

}

#' Generate average pattern of clusters of each posterior chain
#'  from combined list of multiple posterior sample chains
#'
#'
#' @param clean.chlist A list of multiple (or one) posterior sample chains.
#'
#' @return A list of matrices containing the average pattern of clusters within each posterior chain
#'         and a list of matrices containing the sum of each cluster in each posterior chain
#'

GenerateAverageCluster <- function(clean.chlist){


  x <- hdpx::hdp_multi_chain(clean.chlist)
  # input checks
  if (class(x)=="hdpSampleChain") {
    warning('Extracting components on single posterior sampling chain. Recommend switching to multiple independent chains in a hdpSampleMulti object, see ?hdp_multi_chain')
    is_multi <- FALSE
  } else if (class(x)=="hdpSampleMulti") {
    is_multi <- TRUE
  } else {
    stop("x must have class hdpSampleChain or hdpSampleMulti")
  }

  if (is_multi) {
    # list of hdpSampleChain objects
    chlist <- x@chains
    nch <- length(chlist)

    # get final state and number of posterior samples
    finalstate <- hdpx::final_hdpState(chlist[[1]])

  } else {
    #get final state and number of posterior samples
    finalstate <- hdpx::final_hdpState(x)

  }

  # number of categories, DPs,data items at each DP, and frozen priors
  ncat <- hdpx::numcateg(finalstate) ##number of channel
  ndp <- hdpx::numdp(finalstate) ##number of dp

  rm(finalstate)


  # Step (1)
  # Make each ccc (clust_categ_counts) and
  # cdc (clust_dp_counts) matrix have the
  # same number of columns

  if(is_multi){

    maxclust <- max(sapply(chlist, function(x) max(hdpx::numcluster(x))))
    clust_label <- 1:maxclust

    ccc_0 <- lapply(chlist, function(ch){
      lapply(hdpx::clust_categ_counts(ch), function(x){
        ans <- cbind(x, matrix(0, nrow=ncat, ncol=(maxclust-ncol(x)+1)))
        return(ans[, -ncol(ans)])
      })
    })
  }else{

    maxclust <- max(hdpx::numcluster(x))
    clust_label <- 1:maxclust

    ccc_0 <- lapply(hdpx::clust_categ_counts(x), function(x){
      ans <- cbind(x, matrix(0, nrow=ncat, ncol=(maxclust-ncol(x)+1)))
      return(ans[, -ncol(ans)])
    })

  }



  ccc_raw_avg_per_ch <- lapply(ccc_0, function(matlist){ Reduce('+', matlist)/length(matlist) })

  ccc_raw_sum_per_ch <- lapply(ccc_0, function(matlist){ Reduce('+', matlist)})

  ##remove columns with only 0s
  for(i in 1:length(ccc_raw_avg_per_ch)){
    non.zero.columns <- which(colSums(ccc_raw_avg_per_ch[[i]])>0)
    ccc_raw_avg_per_ch[[i]] <- ccc_raw_avg_per_ch[[i]][,non.zero.columns]
    ccc_raw_sum_per_ch[[i]] <- ccc_raw_sum_per_ch[[i]][,non.zero.columns]
  }

  return(list(average.pattern = ccc_raw_avg_per_ch,
              sum.patterm     = ccc_raw_sum_per_ch))
}

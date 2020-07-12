function (x, cluster.method = "kmedians", categ.CI = 0.95, exposure.CI = 0.9,
          cos.merge = 0.9, min.sample = 1)
{
  if (class(x) == "hdpSampleChain") {
    warning("Extracting components on single posterior sampling chain. Recommend switching to multiple independent chains in a hdpSampleMulti object, see ?hdp_multi_chain")
    is_multi <- FALSE
  }
  else if (class(x) == "hdpSampleMulti") {
    is_multi <- TRUE
  }
  else {
    stop("x must have class hdpSampleChain or hdpSampleMulti")
  }
  if (!validObject(x))
    stop("x not a valid object")
  if (class(cos.merge) != "numeric" | cos.merge >= 1 | cos.merge <=
      0) {
    stop("cos.merge must be between 0 and 1")
  }
  if (min.sample%%1 != 0 | min.sample < 1) {
    stop("min.sample must be a positive integer")
  }
  if (is_multi) {
    chlist <- x@chains
    nch <- length(chlist)
    set.seed(sampling_seed(chlist[[1]]), kind = "Mersenne-Twister",
             normal.kind = "Inversion")
    finalstate <- final_hdpState(chlist[[1]])
    nsamp <- sum(sapply(chlist, function(x) hdp_settings(x)$n))
  }
  ncat <- numcateg(finalstate)
  ndp <- numdp(finalstate)
  numdata <- sapply(dp(finalstate), numdata)
  pseudo <- pseudoDP(finalstate)
  rm(finalstate)
  is_prior <- length(pseudo) > 0
  if (is_prior) {
    priorcc <- 1:length(pseudo)
  }
  ccc_0 <- lapply(chlist, function(ch) {
    lapply(clust_categ_counts(ch), function(x) {
      ans <- cbind(x)
      return(ans[, -ncol(ans)])
    })
  })

  ### ccc_0 a list, 1 element per posterior chain
  ## each element is a list of matrices, one matrix
  ## per Gibb's sample, e.g. View(ccc_0[[1]][[1]]).
  ## Each matrix is "aggregated partial spectrum" a.k.a
  ## raw cluster

  cdc_0 <- lapply(chlist, function(ch) {
    lapply(clust_dp_counts(ch), function(x) {
      ans <- cbind(x)
      return(ans[, -ncol(ans)])
    })
  })

  ## Analogous, but each matrix is an exposure, columns
  ## raw clusters, e.g. View(as.matrix(cdc_0[[1]][[1]]))
  ## (actually a matrix like class dgCMatrix)

  clust.number <- {
  }
  for (i in 1:length(ccc_0)) {
    ccc_temp <- ccc_0[[i]]
    cdc_temp <- cdc_0[[i]]
    for (j in 1:length(ccc_temp)) {
      clust_label <- 1:ncol(ccc_temp[[j]])
      clust_cos <- lsa::cosine(ccc_temp[[j]])
      clust_same <- (clust_cos > 0.95 & lower.tri(clust_cos))
      same <- which(clust_same, arr.ind = TRUE)
      if (length(same) > 0) {
        for (index in 1:nrow(same)) {
          clust_label[same[index, 1]] <- clust_label[same[index,
                                                          2]]
        }
      }
      ccc_temp[[j]] <- merge_cols(ccc_temp[[j]], clust_label)
      cdc_temp[[j]] <- merge_cols(cdc_temp[[j]], clust_label)
      clust.number <- c(clust.number, ncol(ccc_temp[[j]]))
    }
    ccc_0[[i]] <- ccc_temp
    cdc_0[[i]] <- cdc_temp
  }
  maxclust <- max(clust.number)
  for (i in 1:length(ccc_0)) {
    ccc_temp <- ccc_0[[i]]
    cdc_temp <- cdc_0[[i]]
    for (j in 1:length(ccc_temp)) {
      ccc_temp[[j]] <- cbind(ccc_temp[[j]], matrix(0,
                                                   nrow = ncat, ncol = (maxclust - ncol(ccc_temp[[j]]) +
                                                                          1)))
      cdc_temp[[j]] <- cbind(cdc_temp[[j]], matrix(0,
                                                   nrow = ndp, ncol = (maxclust - ncol(cdc_temp[[j]]) +
                                                                         1)))
    }
    ccc_0[[i]] <- ccc_temp
    cdc_0[[i]] <- cdc_temp
  }
  ccc_raw_avg_per_ch <- lapply(ccc_0, function(matlist) {
    Reduce("+", matlist)/length(matlist)
  })
  mclust <- ncol(ccc_raw_avg_per_ch[[1]])
  rapch_unlist <- t(do.call(cbind, ccc_raw_avg_per_ch))
  rapch_gf <- rep(1:nch, each = mclust)
  rapch_ic <- rep(1:mclust, times = nch)
  rapch_clust <- flexclust::kcca(rapch_unlist, k = rapch_ic,
                                 group = rapch_gf, family = flexclust::kccaFamily(cluster.method,
                                                                                  groupFun = "differentClusters"))
  rapch_label <- split(flexclust::clusters(rapch_clust), rapch_gf)
  ccc_1 <- Reduce("c", mapply(function(matlist, rank) {
    lapply(matlist, function(mat) {
      ans <- mat[, order(rank)]
      return(ans)
    })
  }, ccc_0, rapch_label, SIMPLIFY = FALSE))
  cdc_1 <- Reduce("c", mapply(function(matlist, rank) {
    lapply(matlist, function(mat) {
      ans <- mat[, order(rank)]
      return(ans)
    })
  }, cdc_0, rapch_label, SIMPLIFY = FALSE))
  remove(ccc_0, cdc_0, ccc_raw_avg_per_ch, rapch_unlist, rapch_gf,
         rapch_ic, rapch_clust, rapch_label, mclust)
  mclust <- ncol(ccc_1[[1]])
  if (mclust == 1) {
    ccc_label <- rep(1, length(ccc_1))
  }
  else {
    ccc_unlist <- t(do.call(cbind, ccc_1))
    for (i in 1:nrow(ccc_unlist)) {
      if (sum(ccc_unlist[i, ]) > 0) {
        ccc_unlist[i, ] <- ccc_unlist[i, ]/sum(ccc_unlist[i,
        ])
      }
    }
    groupfactor <- rep(1:(nsamp), each = mclust)
    initial_clust <- rep(1:mclust, times = nsamp)
    ccc_clust <- flexclust::kcca(ccc_unlist, k = initial_clust,
                                 group = groupfactor, family = flexclust::kccaFamily(cluster.method,
                                                                                     groupFun = "differentClusters"))
    ccc_label <- split(flexclust::clusters(ccc_clust), groupfactor)
    remove(ccc_unlist, groupfactor, initial_clust, ccc_clust)
  }
  ccc_2 <- mapply(function(ccc, label) {
    colnames(ccc) <- label
    ccc[, order(as.numeric(colnames(ccc)))]
  }, ccc_1, ccc_label, SIMPLIFY = FALSE)
  cdc_2 <- mapply(function(cdc, label) {
    colnames(cdc) <- label
    cdc[, order(as.numeric(colnames(cdc)))]
  }, cdc_1, ccc_label, SIMPLIFY = FALSE)
  maxclust <- mclust
  clust_label <- 1:maxclust
  remove(ccc_1, cdc_1, ccc_label)
  avgdistn <- matrix(0, nrow = ncat, ncol = maxclust)
  for (i in 1:maxclust) {
    distns <- sapply(ccc_2, function(x) x[, i]/sum(x[, i]))
    avgdistn[, i] <- rowMeans(distns, na.rm = T)
  }
  clust_cos <- lsa::cosine(avgdistn)
  clust_same <- (clust_cos > cos.merge & lower.tri(clust_cos))
  same <- which(clust_same, arr.ind = TRUE)
  if (length(same) > 0) {
    for (i in 1:nrow(same)) {
      clust_label[same[i, 1]] <- clust_label[same[i, 2]]
    }
  }
  avgdistn_ccc3 <- merge_cols(avgdistn, clust_label)
  ccc_3 <- lapply(ccc_2, merge_cols, clust_label)
  cdc_3 <- lapply(cdc_2, merge_cols, clust_label)
  clust_label <- colnames(ccc_3[[1]])
  if (any(clust_label != colnames(cdc_3)))
    stop("problem in step 3!")
  remove(avgdistn, distns, clust_cos, clust_same, same, ccc_2,
         cdc_2)
  clust_hdp0_ccc4 <- data.frame(matrix(ncol = 0, nrow = ncat))
  use_clust <- c()
  for (ii in 1:ncol(ccc_3[[1]])) {
    compii <- sapply(ccc_3, function(x) x[, ii])
    lowerb <- apply(compii, 1, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in%
          c(0, 1)) {
        NaN
      }
      else {
        round(coda::HPDinterval(samp, categ.CI)[1],
              3)
      }
    })
    if (any(lowerb > 0)) {
      use_clust <- c(use_clust, colnames(ccc_3[[1]])[ii])
    }
    else {
      clust_hdp0_ccc4 <- cbind(clust_hdp0_ccc4, rowSums(compii))
      colnames(clust_hdp0_ccc4)[ncol(clust_hdp0_ccc4)] <- paste("ccc_3_",
                                                                colnames(ccc_3[[1]])[ii], sep = "")
    }
  }
  clust_label[which(!clust_label %in% use_clust)] <- "0"
  ccc_4 <- lapply(ccc_3, merge_cols, clust_label)
  cdc_4 <- lapply(cdc_3, merge_cols, clust_label)
  if (!"0" %in% clust_label) {
    ccc_4 <- lapply(ccc_4, function(x) {
      ans <- cbind(0, x)
      colnames(ans) <- c(0, colnames(x))
      return(ans)
    })
    cdc_4 <- lapply(cdc_4, function(x) {
      ans <- cbind(0, x)
      colnames(ans) <- c(0, colnames(x))
      return(ans)
    })
  }
  avgdistn_ccc4 <- matrix(0, nrow = ncat, ncol = ncol(ccc_4[[1]]))
  for (i in 1:ncol(ccc_4[[1]])) {
    distns <- sapply(ccc_4, function(x) x[, i]/sum(x[, i]))
    avgdistn_ccc4[, i] <- rowMeans(distns, na.rm = T)
  }
  clust_label <- colnames(cdc_4[[1]])
  if (any(clust_label != colnames(cdc_4)))
    stop("problem in step 4!")
  remove(compii, ccc_3, cdc_3, ii, lowerb, use_clust)
  use_clust <- c()
  disregard <- if (is_prior)
    union(which(numdata == 0), pseudo)
  else which(numdata == 0)
  clust_hdp0_ccc5 <- data.frame(matrix(ncol = 0, nrow = ncat))
  for (ii in 1:ncol(cdc_4[[1]])) {
    compii <- sapply(cdc_4, function(x) x[, ii])
    lowerb <- apply(compii[-disregard, ], 1, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in%
          c(0, 1)) {
        NaN
      }
      else {
        round(coda::HPDinterval(samp, exposure.CI)[1],
              3)
      }
    })
    if (sum(lowerb > 0) >= min.sample) {
      use_clust <- c(use_clust, colnames(cdc_4[[1]])[ii])
    }
    else {
      ccc_compii <- sapply(ccc_4, function(x) x[, ii])
      clust_hdp0_ccc5 <- cbind(clust_hdp0_ccc5, rowSums(ccc_compii))
      colnames(clust_hdp0_ccc5)[ncol(clust_hdp0_ccc5)] <- paste("ccc_4_",
                                                                colnames(ccc_4[[1]])[ii], sep = "")
    }
  }
  clust_label[which(!clust_label %in% use_clust)] <- 0
  ccc_5 <- lapply(ccc_4, merge_cols, clust_label)
  cdc_5 <- lapply(cdc_4, merge_cols, clust_label)
  clust_label <- colnames(ccc_5[[1]])
  if (any(clust_label != colnames(cdc_5)))
    stop("problem in step 5!")
  remove(compii, ccc_4, cdc_4, ii, lowerb, use_clust, disregard)
  avg_ndi <- rowMeans(sapply(ccc_5, colSums))
  colorder <- c(1, setdiff(order(avg_ndi, decreasing = T),
                           1))
  ccc_6 <- lapply(ccc_5, function(x) {
    x <- x[, colorder]
    if (is_prior) {
      update <- setdiff(which(!grepl("P", colnames(x))),
                        1)
      if (length(update) > 0) {
        colnames(x)[update] <- paste0("N", 1:length(update))
      }
    }
    else {
      colnames(x) <- 0:(ncol(x) - 1)
    }
    return(x)
  })
  cdc_6 <- lapply(cdc_5, function(x) {
    x <- x[, colorder]
    if (is_prior) {
      update <- setdiff(which(!grepl("P", colnames(x))),
                        1)
      if (length(update) > 0) {
        colnames(x)[update] <- paste0("N", 1:length(update))
      }
    }
    else {
      colnames(x) <- 0:(ncol(x) - 1)
    }
    return(x)
  })
  ncomp <- length(colorder)
  remove(ccc_5, cdc_5, avg_ndi, colorder)
  ccc_ans <- rep(list(matrix(0, nrow = nsamp, ncol = ncat)),
                 ncomp)
  for (i in 1:ncomp) {
    ccc_ans[[i]] <- t(sapply(ccc_6, function(x) x[, i]))
  }
  names(ccc_ans) <- colnames(ccc_6[[1]])
  cdc_ans <- rep(list(matrix(0, nrow = nsamp, ncol = ncomp)),
                 ndp)
  for (i in 1:ndp) {
    cdc_ans[[i]] <- t(sapply(cdc_6, function(x) x[i, ]))
  }
  remove(ccc_6, cdc_6)
  ccc_norm <- lapply(ccc_ans, function(x) x/rowSums(x, na.rm = TRUE))
  ccc_mean <- t(sapply(ccc_norm, colMeans, na.rm = TRUE))
  ccc_credint <- lapply(ccc_norm, function(x) {
    apply(x, 2, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in%
          c(0, 1)) {
        c(NaN, NaN)
      }
      else {
        round(coda::HPDinterval(samp, 0.95), 4)
      }
    })
  })
  cdc_norm <- lapply(cdc_ans, function(x) x/rowSums(x, na.rm = TRUE))
  cdc_mean <- t(sapply(cdc_norm, colMeans, na.rm = TRUE))
  cdc_credint <- lapply(cdc_norm, function(x) {
    apply(x, 2, function(y) {
      samp <- coda::as.mcmc(y)
      if (min(sum(!is.na(samp)), sum(!is.nan(samp))) %in%
          c(0, 1)) {
        c(NaN, NaN)
      }
      else {
        round(coda::HPDinterval(samp, 0.95), 4)
      }
    })
  })
  x@numcomp <- as.integer(ncomp - 1)
  avcount <- colMeans(sapply(ccc_ans, rowSums, na.rm = TRUE),
                      na.rm = TRUE)
  x@prop.ex <- round(1 - avcount[1]/sum(avcount), 3)
  x@comp_cos_merge <- cos.merge
  x@comp_categ_counts <- ccc_ans
  x@comp_dp_counts <- lapply(cdc_ans, as, "dgCMatrix")
  x@comp_categ_distn <- list(mean = ccc_mean, cred.int = ccc_credint,
                             aggregated_raw_clusters_after_cos_merge = avgdistn_ccc3,
                             aggregated_raw_clusters_after_nonzero_categ = avgdistn_ccc4,
                             clust_hdp0_ccc4 = clust_hdp0_ccc4, clust_hdp0_ccc5 = clust_hdp0_ccc5)
  x@comp_dp_distn <- list(mean = cdc_mean, cred.int = cdc_credint)
  if (!validObject(x))
    warning("Not a valid hdpSampleChain/Multi object.")
  return(x)
}

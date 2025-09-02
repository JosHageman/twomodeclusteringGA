#' Validate a two-mode clustering partition by global and per-block significance
#'
#' @description
#' Given a numeric matrix and a full two-mode partition (exclusive row and column clusters),
#' this function tests whether the fitted block-means model explains more structure
#' than expected under a no-structure null. The global test uses an F-statistic based
#' on SS_fit and SSE derived from your fitness definition. Optionally, it also reports
#' per-block chi-square tests and a fast Monte Carlo p-value using random partitions
#' (no GA reruns).
#'
#' @param myMatrix Numeric matrix or coercible data.frame.
#' @param rowClusters Integer vector of length nrow(myMatrix) with cluster labels (1..kR, arbitrary labels allowed).
#' @param colClusters Integer vector of length ncol(myMatrix) with cluster labels (1..kC, arbitrary labels allowed).
#' @param center Logical, center the matrix by its global mean before testing (default TRUE).
#'   Centering aligns the null with zero-mean noise and generally stabilizes inference.
#' @param perBlock Logical, compute per-block tests (default TRUE).
#' @param monteCarlo Integer, number of random partitions to draw for a MC p-value (default 0 disables).
#' @param fixBlockSizes Logical, if TRUE keep row and column cluster sizes equal to the observed sizes
#'   when generating random partitions (default TRUE). If FALSE, only kR and kC are fixed.
#' @param storeNull Logical, store the vector of null F statistics from random partitions (default FALSE).
#'   If FALSE, only quantiles are stored.
#' @param seed Optional integer seed for reproducibility (default NULL).
#'
#' @return A list of class "twomodeValidation" with elements:
#' \itemize{
#'   \item \code{nR}, \code{nC}, \code{kR}, \code{kC}
#'   \item \code{dfModel}, \code{dfResid}
#'   \item \code{ssTot}, \code{ssFit}, \code{sse}, \code{sigma2Hat}, \code{r2}
#'   \item \code{fStat}, \code{pValue} (global F test)
#'   \item \code{perBlock} (data.frame with per-block stats) if \code{perBlock=TRUE}
#'   \item \code{mc} (list with \code{nSim}, \code{pMonteCarlo}, \code{fNull} or \code{fNullQuantiles}) if \code{monteCarlo>0}
#' }
#'
#' @importFrom stats pf pchisq p.adjust
#' 
#' @export
#'
#' @examples
#' data("twomodeToy")
#' myMatrix_s <- scale(twomodeToy)
#'
#' #Run the GA-based two-mode clustering
#' result <- twomodeClusteringGA(
#'   myMatrix = myMatrix_s,
#'   nRowClusters = 2,
#'   nColClusters = 3,
#'   seeds = 1,
#'   maxiter = 200,
#'   popSize = 30,
#'   elitism = 1,
#'   validate = FALSE,
#'   verbose = TRUE
#' )
#' 
#' result$validation <- validateTwomodePartition(myMatrix_s, 
#'                                 rowClusters=result$rowClusters, 
#'                                 colClusters=result$colClusters)
#'
#' #Inspect the result
#' print(result)
#' summary(result)
#' myTwomodeResult <- as.data.frame(result)
#' head(myTwomodeResult)
#' 
#' #Plot the clustered heatmap
#' plotTwomodeClustering(
#'   myMatrix = myMatrix_s,
#'   result   = result,
#'   title    = "Two-mode clustering Toy example",
#'   fixAspect = FALSE
#' )
validateTwomodePartition <- function(myMatrix,
                                     rowClusters,
                                     colClusters,
                                     center = TRUE,
                                     perBlock = TRUE,
                                     monteCarlo = 0,
                                     fixBlockSizes = TRUE,
                                     storeNull = FALSE,
                                     seed = NULL) {
  dat <- as.matrix(myMatrix)
  if (!is.numeric(dat)) stop("myMatrix must be numeric or coercible to numeric.")
  storage.mode(dat) <- "double"
  
  nR <- nrow(dat)
  nC <- ncol(dat)
  n <- nR * nC
  
  if (length(rowClusters) != nR) stop("rowClusters length must equal nrow(myMatrix).")
  if (length(colClusters) != nC) stop("colClusters length must equal ncol(myMatrix).")
  
  # Relabel clusters to 1..k (robust to arbitrary labels)
  a <- match(rowClusters, sort(unique(rowClusters)))
  b <- match(colClusters, sort(unique(colClusters)))
  kR <- max(a)
  kC <- max(b)
  
  if (center) {
    grandMean <- mean(dat)
    dat <- dat - grandMean
  } else {
    grandMean <- 0
  }
  
  # Helper: compute SS_tot, SS_fit, SSE, S (block sums), M (block counts)
  computeBlockFit <- function(x, aLab, bLab, kRval, kCval) {
    ssTot <- sum(x * x)
    
    # Block sums via two rowsum() calls (fast and robust)
    sRow <- rowsum(x, group = factor(aLab, levels = seq_len(kRval)), reorder = TRUE)        # kR x nC
    sCol <- rowsum(t(sRow), group = factor(bLab, levels = seq_len(kCval)), reorder = TRUE)  # kC x kR
    S <- t(sCol)  # kR x kC
    
    nPerRowCluster <- tabulate(aLab, nbins = kRval)
    nPerColCluster <- tabulate(bLab, nbins = kCval)
    M <- outer(nPerRowCluster, nPerColCluster, "*")
    
    ssFit <- sum((S * S) / M)
    sse <- ssTot - ssFit
    
    list(ssTot = ssTot, ssFit = ssFit, sse = sse, S = S, M = M,
         nPerRowCluster = nPerRowCluster, nPerColCluster = nPerColCluster)
  }
  
  fit <- computeBlockFit(dat, a, b, kR, kC)
  
  dfModel <- kR * kC
  dfResid <- n - dfModel
  if (dfResid <= 0) stop("Not enough degrees of freedom: nR*nC must exceed kR*kC.")
  
  sigma2Hat <- fit$sse / dfResid
  r2 <- if (fit$ssTot > 0) fit$ssFit / fit$ssTot else NA_real_
  
  fNum <- fit$ssFit / dfModel
  fDen <- fit$sse / dfResid
  fStat <- fNum / fDen
  pValue <- 1 - pf(fStat, dfModel, dfResid)
  
  result <- list(
    nR = nR, nC = nC, kR = kR, kC = kC,
    dfModel = dfModel, dfResid = dfResid,
    centered = center, grandMean = grandMean,
    ssTot = fit$ssTot, ssFit = fit$ssFit, sse = fit$sse,
    sigma2Hat = sigma2Hat, r2 = r2,
    fStat = fStat, pValue = pValue
  )
  
  # Optional: per-block tests
  if (isTRUE(perBlock)) {
    S <- fit$S
    M <- fit$M
    # For each block (u,v), test (S_uv^2)/(sigma2Hat*M_uv) ~ chi-square_1 under H0
    z2 <- (S * S) / (sigma2Hat * M)
    pBlock <- 1 - pchisq(z2, df = 1)
    
    blockMeans <- S / M
    perBlockDf <- data.frame(
      rowCluster = rep.int(seq_len(kR), times = kC),
      colCluster = rep(seq_len(kC), each = kR),
      nCells = as.integer(c(M)),
      sumValues = c(S),
      meanValue = c(blockMeans),
      effectSS = c((S * S) / M),
      chiSq1 = c(z2),
      pValue = c(pBlock),
      stringsAsFactors = FALSE
    )
    # Adjust for multiple testing (optional but useful)
    perBlockDf$pAdjBH <- p.adjust(perBlockDf$pValue, method = "BH")
    result$perBlock <- perBlockDf
  }
  
  # Optional: Monte Carlo using random partitions (no GA)
  if (is.numeric(monteCarlo) && monteCarlo > 0) {
    if (!is.null(seed)) {
      oldSeed <- .Random.seed
      set.seed(seed)
      on.exit({
        # Restore RNG state if it existed
        if (exists("oldSeed", inherits = FALSE)) {
          .Random.seed <<- oldSeed
        }
      }, add = TRUE)
    }
    
    fNull <- numeric(monteCarlo)
    
    if (isTRUE(fixBlockSizes)) {
      # Preserve observed cluster sizes
      rowSizes <- fit$nPerRowCluster
      colSizes <- fit$nPerColCluster
      # Precompute boundaries for assignment
      assignBySizes <- function(nTot, sizes) {
        idx <- sample.int(nTot, nTot, replace = FALSE)
        labs <- integer(nTot)
        start <- 1L
        for (k in seq_along(sizes)) {
          if (sizes[k] > 0) {
            labs[idx[start:(start + sizes[k] - 1L)]] <- k
            start <- start + sizes[k]
          }
        }
        labs
      }
      for (bIter in seq_len(monteCarlo)) {
        aRand <- assignBySizes(nR, rowSizes)
        bRand <- assignBySizes(nC, colSizes)
        fitRand <- computeBlockFit(dat, aRand, bRand, kR, kC)
        fNull[bIter] <- (fitRand$ssFit / dfModel) / (fitRand$sse / dfResid)
      }
    } else {
      # Only fix kR and kC; allow varying sizes
      for (bIter in seq_len(monteCarlo)) {
        aRand <- sample.int(kR, nR, replace = TRUE)
        bRand <- sample.int(kC, nC, replace = TRUE)
        # Ensure all clusters are nonempty; resample if needed (rare for moderate sizes)
        if (min(tabulate(aRand, nbins = kR)) == 0 || min(tabulate(bRand, nbins = kC)) == 0) {
          # simple retry loop
          tries <- 0L
          repeat {
            aRand <- sample.int(kR, nR, replace = TRUE)
            bRand <- sample.int(kC, nC, replace = TRUE)
            if (min(tabulate(aRand, nbins = kR)) > 0 && min(tabulate(bRand, nbins = kC)) > 0) break
            tries <- tries + 1L
            if (tries > 100L) stop("Failed to sample nonempty random partitions after 100 tries.")
          }
        }
        fitRand <- computeBlockFit(dat, aRand, bRand, kR, kC)
        fNull[bIter] <- (fitRand$ssFit / dfModel) / (fitRand$sse / dfResid)
      }
    }
    
    # Monte Carlo p-value (upper-tail)
    # add-one smoothing to avoid zero estimates
    pMc <- (sum(fNull >= fStat) + 1) / (length(fNull) + 1)
    
    if (isTRUE(storeNull)) {
      mcOut <- list(nSim = length(fNull), pMonteCarlo = pMc, fNull = fNull)
    } else {
      q <- stats::quantile(fNull, probs = c(0.5, 0.9, 0.95, 0.99), names = TRUE, type = 7)
      mcOut <- list(nSim = length(fNull), pMonteCarlo = pMc, fNullQuantiles = q)
    }
    result$mc <- mcOut
  }
  
  class(result) <- "twomodeValidation"
  return(result)
}

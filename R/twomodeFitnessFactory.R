#' Two-mode clustering genetic algorithm evaluation function (fast, robust)
#'
#' @description Fast evaluation of a two-mode clustering solution.
#' 
#' @param myMatrix Numeric matrix or coercible data.frame.
#' 
#' @return Function(string, ...) -> numeric fitness value = negative SSE (higher is better).
#' 
#' @export
twomodeFitnessFactory <- function(myMatrix) {
  dat <- as.matrix(myMatrix)
  if (!is.numeric(dat)) stop("myMatrix must be numeric or coercible to numeric.")
  storage.mode(dat) <- "double"
  nR <- nrow(dat)
  nC <- ncol(dat)
  constSumX2 <- sum(dat * dat)

  function(string, ...) { 

    # Relabel clusters to 1..k (robust tegen 0/NA/gaten)
    aRaw <- string[seq_len(nR)]
    bRaw <- string[nR + seq_len(nC)]
    if (any(!is.finite(aRaw)) || any(!is.finite(bRaw))) {
      stop("Cluster assignments must be finite.")
    }
    a <- match(aRaw, sort(unique(aRaw)))  # 1..kR
    b <- match(bRaw, sort(unique(bRaw)))  # 1..kC
    kR <- max(a)
    kC <- max(b)

    # Block sums S_{uv} via twee rowsum() calls
    sRow <- rowsum(dat, group = factor(a, levels = seq_len(kR)), reorder = TRUE)     # kR x nC
    sCol <- rowsum(t(sRow), group = factor(b, levels = seq_len(kC)), reorder = TRUE) # kC x kR
    S <- t(sCol)  # kR x kC

    # Block counts
    nPerRowCluster <- tabulate(a, nbins = kR)
    nPerColCluster <- tabulate(b, nbins = kC)
    M <- outer(nPerRowCluster, nPerColCluster, "*")  # > 0

    # SSE = sum(X^2) - sum(S^2 / M)
    sumS2OverM <- sum((S * S) / M)
    -(constSumX2 - sumS2OverM)
  }
}

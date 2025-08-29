#' Summary method for twomodeClustering objects
#'
#' Creates a summary of a twomodeClustering object, including matrix dimensions,
#' cluster sizes, fitness, optional bicluster summaries (if matrix available),
#' and optional validation highlights (if validation is present).
#'
#' @param object An object of class 'twomodeClustering'.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An object of class summary.twomodeClustering with components:
#' \describe{
#'   \item{matrixDim}{Named integer vector: rows, cols}
#'   \item{nRowClusters}{Number of row clusters}
#'   \item{nColClusters}{Number of column clusters}
#'   \item{rowClusterSizes}{Table of row cluster sizes}
#'   \item{colClusterSizes}{Table of column cluster sizes}
#'   \item{biclusters}{Data frame with bicluster summaries (if myMatrix present), possibly merged with validation per-block stats}
#'   \item{fitness}{Best fitness value if available, else NA}
#'   \item{validationGlobal}{List with r2, fStat, pValue, dfModel, dfResid, pMonteCarlo (if present), or NULL}
#'   \item{nSigBlocks}{Number of BH-significant blocks at 0.05 if available, else NULL}
#'   \item{rowContribution}{Data frame with total effectSS per row cluster (if available), else NULL}
#'   \item{colContribution}{Data frame with total effectSS per column cluster (if available), else NULL}
#' }
#' 
#' @export
summary.twomodeClustering <- function(object, ...) {
  if (!inherits(object, "twomodeClustering")) {
    stop("Object must be of class 'twomodeClustering'")
  }
  # basic structure
  rowCl <- object$rowClusters
  colCl <- object$colClusters
  nR <- length(rowCl); nC <- length(colCl)

  rowSizes <- table(rowCl)
  colSizes <- table(colCl)

  # matrix if available
  mat <- if (!is.null(object$myMatrix)) as.matrix(object$myMatrix) else NULL

  # bicluster summaries based on matrix (mean/sd)
  biclusters <- NULL
  if (!is.null(mat)) {
    biclusters <- data.frame(
      rowCluster = integer(0),
      colCluster = integer(0),
      nRows = integer(0),
      nCols = integer(0),
      nCells = integer(0),
      meanValue = numeric(0),
      sdValue = numeric(0),
      stringsAsFactors = FALSE
    )
    uR <- sort(unique(rowCl))
    uC <- sort(unique(colCl))
    for (r in uR) {
      rMask <- rowCl == r
      nRows <- sum(rMask)
      for (c in uC) {
        cMask <- colCl == c
        nCols <- sum(cMask)
        subMat <- mat[rMask, cMask, drop = FALSE]
        biclusters <- rbind(biclusters, data.frame(
          rowCluster = r,
          colCluster = c,
          nRows = nRows,
          nCols = nCols,
          nCells = as.integer(nRows * nCols),
          meanValue = mean(subMat),
          sdValue = stats::sd(subMat),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # fitness if present
  fitness <- if (!is.null(object$bestFitness)) object$bestFitness else NA_real_

  # validation highlights if present
  validationGlobal <- NULL
  nSigBlocks <- NULL
  rowContribution <- NULL
  colContribution <- NULL

  if (!is.null(object$validation) && is.list(object$validation)) {
    v <- object$validation
    # global
    if (!is.null(v$r2) && !is.null(v$fStat) && !is.null(v$pValue) &&
        !is.null(v$dfModel) && !is.null(v$dfResid)) {
      validationGlobal <- list(
        r2 = v$r2,
        fStat = v$fStat,
        pValue = v$pValue,
        dfModel = v$dfModel,
        dfResid = v$dfResid,
        pMonteCarlo = if (!is.null(v$mc) && !is.null(v$mc$pMonteCarlo)) v$mc$pMonteCarlo else NULL
      )
    }
    # per-block merge into biclusters if both exist
    if (!is.null(v$perBlock) && is.data.frame(v$perBlock)) {
      if (!is.null(biclusters)) {
        # merge by rowCluster, colCluster
        keepCols <- intersect(names(v$perBlock),
                              c("rowCluster","colCluster","nCells","sumValues","meanValue","effectSS","chiSq1","pValue","pAdjBH"))
        vb <- v$perBlock[keepCols]
        # avoid duplicate meanValue from both sides; rename validation mean as meanValue_fit if clash
        if ("meanValue" %in% names(vb) && "meanValue" %in% names(biclusters)) {
          names(vb)[names(vb) == "meanValue"] <- "meanValue_fit"
        }
        biclusters <- merge(biclusters, vb, by = c("rowCluster","colCluster"), all.x = TRUE, sort = FALSE)
      }
      # contributions and significance counts
      if ("effectSS" %in% names(v$perBlock)) {
        rowContribution <- stats::aggregate(effectSS ~ rowCluster, data = v$perBlock, sum)
        colContribution <- stats::aggregate(effectSS ~ colCluster, data = v$perBlock, sum)
      }
      if ("pAdjBH" %in% names(v$perBlock)) {
        nSigBlocks <- sum(is.finite(v$perBlock$pAdjBH) & v$perBlock$pAdjBH <= 0.05)
      } else if ("pValue" %in% names(v$perBlock)) {
        nSigBlocks <- sum(is.finite(v$perBlock$pValue) & v$perBlock$pValue <= 0.05)
      }
    }
  }

  out <- list(
    matrixDim = c(rows = nR, cols = nC),
    nRowClusters = length(rowSizes),
    nColClusters = length(colSizes),
    rowClusterSizes = rowSizes,
    colClusterSizes = colSizes,
    biclusters = biclusters,
    fitness = fitness,
    validationGlobal = validationGlobal,
    nSigBlocks = nSigBlocks,
    rowContribution = rowContribution,
    colContribution = colContribution
  )
  class(out) <- "summary.twomodeClustering"
  return(out)
}

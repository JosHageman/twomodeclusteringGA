#' Print method for summary.twomodeClustering objects
#'
#' Prints key information about a two-mode clustering result, including matrix dimensions,
#' cluster sizes, fitness, and (if available) validation highlights.
#'
#' @param x An object of class 'summary.twomodeClustering'.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns x.
#' 
#' @export
print.summary.twomodeClustering <- function(x, ...) {
  cat("Summary of two-mode clustering\n")
  cat(sprintf("Matrix: %d rows x %d cols\n", x$matrixDim["rows"], x$matrixDim["cols"]))
  cat(sprintf("Row clusters: %d (%s)\n",
              x$nRowClusters,
              paste(x$rowClusterSizes, collapse = ", ")))
  cat(sprintf("Col clusters: %d (%s)\n",
              x$nColClusters,
              paste(x$colClusterSizes, collapse = ", ")))
  if (!is.na(x$fitness)) cat(sprintf("Fitness (-SSE): %.6f\n", x$fitness))

  # validation highlights
  if (!is.null(x$validationGlobal)) {
    vg <- x$validationGlobal
    cat("Validation (global):\n")
    cat(sprintf("  R^2: %.4f\n", vg$r2))
    cat(sprintf("  F(%d,%d) = %.4f, p = %.6g\n",
                as.integer(vg$dfModel), as.integer(vg$dfResid),
                as.numeric(vg$fStat), as.numeric(vg$pValue)))
    if (!is.null(vg$pMonteCarlo)) {
      cat(sprintf("  Monte Carlo: p_MC = %.6g\n", as.numeric(vg$pMonteCarlo)))
    }
    if (!is.null(x$nSigBlocks)) {
      cat(sprintf("  Significant blocks (BH 0.05): %d\n", as.integer(x$nSigBlocks)))
    }
    # compact contributors
    if (!is.null(x$rowContribution) && nrow(x$rowContribution) > 0) {
      ordR <- order(-x$rowContribution$effectSS)
      topR <- x$rowContribution$rowCluster[ordR]
      cat(sprintf("  Top row contributions: %s\n",
                  paste(utils::head(topR, 3), collapse = ", ")))
    }
    if (!is.null(x$colContribution) && nrow(x$colContribution) > 0) {
      ordC <- order(-x$colContribution$effectSS)
      topC <- x$colContribution$colCluster[ordC]
      cat(sprintf("  Top col contributions: %s\n",
                  paste(utils::head(topC, 3), collapse = ", ")))
    }
  } else {
    # gentle hint if no validation present
    cat("Validation: not available. Consider validate=TRUE or run validateTwomodePartition().\n")
  }

  # biclusters: avoid dumping huge tables by default
  if (!is.null(x$biclusters)) {
    nB <- nrow(x$biclusters)
    cat(sprintf("\nBicluster summaries: %d blocks\n", nB))
    if (nB > 0) {
      toShow <- utils::head(x$biclusters, 10)
      print(toShow, row.names = FALSE, digits = 3)
      if (nB > 10) cat(sprintf("... %d more rows\n", nB - 10))
    }
  }
  invisible(x)
}

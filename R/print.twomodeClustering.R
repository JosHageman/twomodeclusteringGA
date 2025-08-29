#' Print method for twomodeClustering objects
#'
#' Prints a concise summary of a twomodeClustering object, including matrix dimensions,
#' cluster counts, fitness, and (if available) validation results.
#'
#' @param x An object of class 'twomodeClustering'.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns \code{x}.
#' @export
print.twomodeClustering <- function(x, ...) {
  if (!inherits(x, "twomodeClustering")) {
    stop("Object must be of class 'twomodeClustering'")
  }
  nR <- length(x$rowClusters)
  nC <- length(x$colClusters)
  nRowCl <- length(unique(x$rowClusters))
  nColCl <- length(unique(x$colClusters))
  fitness <- if (!is.null(x$bestFitness)) x$bestFitness else NA_real_

  cat("Two-mode clustering result\n")
  cat(sprintf("Matrix: %d rows x %d cols\n", nR, nC))
  cat(sprintf("Row clusters: %d\n", nRowCl))
  cat(sprintf("Column clusters: %d\n", nColCl))
  if (is.finite(fitness)) cat(sprintf("Fitness (-SSE): %.6f\n", fitness))

  # Optional validation block
  if (!is.null(x$validation) && is.list(x$validation)) {
    v <- x$validation
    hasGlobal <- !is.null(v$fStat) && !is.null(v$pValue) &&
                 !is.null(v$r2) && !is.null(v$dfModel) && !is.null(v$dfResid)
    if (isTRUE(hasGlobal)) {
      cat("Validation (global):\n")
      cat(sprintf("  R^2: %.4f\n", v$r2))
      cat(sprintf("  F(%d,%d) = %.4f, p = %.6g\n",
                  as.integer(v$dfModel), as.integer(v$dfResid),
                  as.numeric(v$fStat), as.numeric(v$pValue)))
    }
    # Monte Carlo, if present
    if (!is.null(v$mc) && is.list(v$mc) && !is.null(v$mc$pMonteCarlo)) {
      cat(sprintf("  Monte Carlo: p_MC = %.6g (n = %d)\n",
                  as.numeric(v$mc$pMonteCarlo),
                  if (!is.null(v$mc$nSim)) as.integer(v$mc$nSim) else NA_integer_))
    }
    # Per-block summary, if present
    if (!is.null(v$perBlock) && is.data.frame(v$perBlock) && nrow(v$perBlock) > 0) {
      if ("pAdjBH" %in% names(v$perBlock)) {
        nSig <- sum(is.finite(v$perBlock$pAdjBH) & v$perBlock$pAdjBH <= 0.05)
        cat(sprintf("  Significant blocks (BH FDR 0.05): %d of %d\n", nSig, nrow(v$perBlock)))
      } else if ("pValue" %in% names(v$perBlock)) {
        nSig <- sum(is.finite(v$perBlock$pValue) & v$perBlock$pValue <= 0.05)
        cat(sprintf("  Significant blocks (p <= 0.05): %d of %d\n", nSig, nrow(v$perBlock)))
      }
    }
  }

  invisible(x)
}

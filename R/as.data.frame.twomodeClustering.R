#' Convert a twomodeClustering object to a data.frame
#'
#' This function creates a data.frame representation of a twomodeClustering object,
#' listing the cluster assignments for both rows and columns.
#'
#' @param x An object of class 'twomodeClustering'.
#' @param row.names Optional vector of row names for the resulting data.frame.
#' @param optional Logical. If TRUE, allows optional parameters for data.frame.
#' @param myMatrix Optional matrix to provide row and column names.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A data.frame with columns: name, type (row/col), and cluster assignment.
#' 
#' @export
as.data.frame.twomodeClustering <- function(x, row.names = NULL, optional = FALSE, myMatrix = NULL, ...) {
  if (!inherits(x, "twomodeClustering")) {
    stop("Object must be of class 'twomodeClustering'")
  }
  
  nR <- length(x$rowClusters)
  nC <- length(x$colClusters)
  
  rowNames <- if (!is.null(myMatrix) && !is.null(rownames(myMatrix))) {
    rownames(myMatrix)
  } else {
    paste0("Row", seq_len(nR))
  }
  colNames <- if (!is.null(myMatrix) && !is.null(colnames(myMatrix))) {
    colnames(myMatrix)
  } else {
    paste0("Col", seq_len(nC))
  }
  
  dfRow <- data.frame(
    name = rowNames,
    type = "row",
    cluster = x$rowClusters,
    stringsAsFactors = FALSE
  )
  dfCol <- data.frame(
    name = colNames,
    type = "col",
    cluster = x$colClusters,
    stringsAsFactors = FALSE
  )
  
  df <- rbind(dfRow, dfCol)
  rownames(df) <- NULL
  df
}


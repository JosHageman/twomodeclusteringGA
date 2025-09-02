#' Plot two-mode clustering results (validation-aware, compact labels)
#'
#' @description Heatmap of the clustered matrix with clear cluster boundaries.
#' If \code{result$validation} is present, each block shows one label with
#' the chosen value plus significance stars.
#'
#' @param myMatrix Numeric matrix or coercible data.frame with the data.
#' @param result Result from \code{twomodeClusteringGA()}, with \code{rowClusters}, \code{colClusters},
#'   and optionally \code{validation}.
#' @param title Text for title.
#' @param xlabel Text for x-axis label.
#' @param ylabel Text for y-axis label.
#' @param varOrder Order of column clusters (0 = automatic).
#' @param objOrder Order of row clusters (0 = automatic).
#' @param palette Color scale: "diverging", "viridis", or "grey".
#' @param showBoundaries Logical; show cluster boundaries.
#' @param boundaryColor Color of the boundaries.
#' @param boundarySize Width of the boundaries.
#' @param showMeans Logical; show block labels (value + stars if validation).
#' @param fixAspect Logical; square cells.
#' @param showValidation Logical; use validation information if available.
#' @param value Which block statistic to label: "mean", "standardized", or "effectSS".
#'   For "standardized", sign(mean) * sqrt(chi^2_1) is shown if validation is available.
#' @param digits Number of decimals in the label.
#' @param sigLevels Thresholds for stars: c(0.001, 0.01, 0.05, 0.1).
#' @param showMarginal Logical; show "." for p < 0.1.
#' @param labelColor Color of the block labels.
#' @param showGlobal Logical; add global validation (R2, F, p, p_MC) to subtitle.
#'
#' @return A ggplot object.
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
#'   validate = TRUE,
#'   verbose = TRUE
#' )
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
#' 
plotTwomodeClustering <- function(myMatrix,
                                  result,
                                  title = "",
                                  xlabel = "",
                                  ylabel = "",
                                  varOrder = 0,
                                  objOrder = 0,
                                  palette = c("diverging","viridis","grey"),
                                  showBoundaries = TRUE,
                                  boundaryColor = "white",
                                  boundarySize = 1,
                                  showMeans = TRUE,
                                  fixAspect = TRUE,
                                  showValidation = TRUE,
                                  value = c("mean","standardized","effectSS"),
                                  digits = 2,
                                  sigLevels = c(0.001, 0.01, 0.05, 0.1),
                                  showMarginal = TRUE,
                                  labelColor = "white",
                                  showGlobal = TRUE) {
  
  palette <- match.arg(palette)
  value <- match.arg(value)
  
  mat <- as.matrix(myMatrix)
  if (!is.numeric(mat)) stop("myMatrix must be numeric or coercible to a numeric matrix.")
  nR <- nrow(mat); nC <- ncol(mat)
  
  if (!is.list(result) || is.null(result$rowClusters) || is.null(result$colClusters)) {
    stop("Argument 'result' must be the list returned by twomodeClusteringGA() with rowClusters and colClusters.")
  }
  rowCl <- as.integer(result$rowClusters)
  colCl <- as.integer(result$colClusters)
  if (length(rowCl) != nR || length(colCl) != nC) {
    stop("Length of rowClusters/colClusters does not match matrix dimensions.")
  }
  
  if (length(varOrder) == 1 && identical(varOrder[1], 0)) varOrder <- sort(unique(colCl))
  if (length(objOrder) == 1 && identical(objOrder[1], 0)) objOrder <- sort(unique(rowCl))
  
  colOrder <- order(factor(colCl, levels = varOrder))
  rowOrder <- order(factor(rowCl, levels = objOrder))
  matRe <- mat[rowOrder, colOrder, drop = FALSE]
  
  rowLabels <- rownames(mat); if (is.null(rowLabels)) rowLabels <- as.character(seq_len(nR))
  rowLabels <- rowLabels[rowOrder]
  colLabels <- colnames(mat); if (is.null(colLabels)) colLabels <- as.character(seq_len(nC))
  colLabels <- colLabels[colOrder]
  
  rr <- rep(seq_len(nR), times = nC)
  cc <- rep(seq_len(nC), each = nR)
  val <- as.vector(matRe)
  rowClOrd <- rep(rowCl[rowOrder], times = nC)
  colClOrd <- rep(colCl[colOrder], each = nR)
  df <- data.frame(row = rr, col = cc, value = val,
                   rowCluster = rowClOrd, colCluster = colClOrd)
  
  # Counts and centers per block
  colCounts <- as.numeric(table(factor(colCl[colOrder], levels = varOrder)))
  rowCounts <- as.numeric(table(factor(rowCl[rowOrder], levels = objOrder)))
  vBound <- if (length(colCounts) > 1) cumsum(colCounts)[seq_len(length(colCounts) - 1)] + 0.5 else numeric(0)
  hBound <- if (length(rowCounts) > 1) cumsum(rowCounts)[seq_len(length(rowCounts) - 1)] + 0.5 else numeric(0)
  rowCenters <- cumsum(rowCounts) - rowCounts/2 + 0.5
  colCenters <- cumsum(colCounts) - colCounts/2 + 0.5
  
  
  # Helper: p -> stars
  pToStars <- function(p, lv = sigLevels, showMarg = showMarginal) {
    if (!is.finite(p)) return("")
    if (p < lv[1]) return("***")
    if (p < lv[2]) return("**")
    if (p < lv[3]) return("*")
    if (showMarg && p < lv[4]) return(".")
    ""
  }
  
  # Validation present?
  hasValidation <- isTRUE(showValidation) &&
    !is.null(result$validation) && is.list(result$validation) &&
    !is.null(result$validation$perBlock) && is.data.frame(result$validation$perBlock)
  
  # Labels per block (value + stars if validation)
  meanDf <- NULL
  subtitleTxt <- NULL
  
  if (showMeans) {
    if (hasValidation) {
      vb <- result$validation$perBlock
      # choose p-column: BH > pValue if available
      pcol <- if ("pAdjBH" %in% names(vb)) "pAdjBH" else if ("pValue" %in% names(vb)) "pValue" else NULL
      pvec <- if (!is.null(pcol)) vb[[pcol]] else rep(NA_real_, nrow(vb))
      
      # choose label value
      labVal <- vb$meanValue
      if (identical(value, "standardized") && all(c("chiSq1","meanValue") %in% names(vb))) {
        sgn <- ifelse(vb$meanValue >= 0, 1, -1)
        labVal <- sgn * sqrt(pmax(vb$chiSq1, 0))
      } else if (identical(value, "effectSS") && "effectSS" %in% names(vb)) {
        labVal <- vb$effectSS
      }
      
      valStr <- sprintf(paste0("%.", digits, "f"), labVal)
      stars <- vapply(pvec, pToStars, character(1L))
      label <- ifelse(stars == "", valStr, paste0(valStr, stars))
      
      meanDf <- data.frame(
        rowCluster = vb$rowCluster,
        colCluster = vb$colCluster,
        label = label,
        stringsAsFactors = FALSE
      )
    } else {
      # no validation: show mean per block from the data
      agg <- stats::aggregate(value ~ rowCluster + colCluster, df, mean)
      showVal <- agg$value
      # for "standardized"/"effectSS" without validation we have no statistic; fallback to mean
      valStr <- sprintf(paste0("%.", digits, "f"), showVal)
      meanDf <- data.frame(
        rowCluster = agg$rowCluster,
        colCluster = agg$colCluster,
        label = valStr,
        stringsAsFactors = FALSE
      )
    }
    # positions
    meanDf$row <- rowCenters[match(meanDf$rowCluster, objOrder)]
    meanDf$col <- colCenters[match(meanDf$colCluster, varOrder)]
  }
  
  # Global validation text (R², F, p, p_MC) in subtitle
  if (hasValidation && isTRUE(showGlobal)) {
    r2 <- result$validation$r2
    fstat <- result$validation$fStat
    df1 <- result$validation$dfModel
    df2 <- result$validation$dfResid
    pvalG <- result$validation$pValue
    mcTxt <- if (!is.null(result$validation$mc) && !is.null(result$validation$mc$pMonteCarlo)) {
      sprintf("; p_MC = %.4g", result$validation$mc$pMonteCarlo)
    } else ""
    # pretty R²
    subtitleTxt <- sprintf("R\u00B2 = %.3f; F(%d,%d) = %.3f; p = %.4g%s",
                           r2, df1, df2, fstat, pvalG, mcTxt)
  }
  
  # Base heatmap
  p <- ggplot2::ggplot(df, ggplot2::aes(x = col, y = row, fill = value)) +
    ggplot2::geom_raster()
  
  if (palette == "diverging") {
    p <- p + ggplot2::scale_fill_gradient2(low = "#008000", mid = "#000000", high = "#ff0000",
                                           name = "Value")
  } else if (palette == "viridis") {
    p <- p + ggplot2::scale_fill_viridis_c(option = "C", name = "Value")
  } else if (palette == "grey") {
    p <- p + ggplot2::scale_fill_grey(name = "Value")
  }
  
  # Labels (value + stars or just value)
  if (showMeans && !is.null(meanDf)) {
    p <- p + ggplot2::geom_text(data = meanDf,
                                ggplot2::aes(x = col, y = row, label = label),
                                inherit.aes = FALSE,
                                color = labelColor, size = 3)
  }
  
  # Boundaries
  if (showBoundaries) {
    p <- p +
      ggplot2::geom_vline(xintercept = vBound, linewidth = boundarySize, color = boundaryColor) +
      ggplot2::geom_hline(yintercept = hBound, linewidth = boundarySize, color = boundaryColor)
  }
  
  if (fixAspect) p <- p + ggplot2::coord_fixed()
  
  p <- p +
    ggplot2::scale_y_reverse(expand = c(0, 0),
                             breaks = seq_len(nR),
                             labels = rowLabels) +
    ggplot2::scale_x_continuous(expand = c(0, 0),
                                breaks = seq_len(nC),
                                labels = colLabels) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
  
  if (!is.null(subtitleTxt)) {
    p <- p + ggplot2::labs(title = title, subtitle = subtitleTxt, x = xlabel, y = ylabel)
  } else {
    p <- p + ggplot2::labs(title = title, x = xlabel, y = ylabel)
  }
  
  return(p)
}

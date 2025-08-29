#' Two-mode clustering using genetic algorithm (with optional validation)
#'
#' @description Performs two-mode clustering on a numeric matrix using a genetic algorithm.
#' The algorithm simultaneously clusters rows and columns to minimize within-cluster
#' sum of squared errors (SSE). Optionally, a validation step is executed that tests
#' the statistical significance of the found partition using \code{validateTwomodePartition()}.
#'
#' @param myMatrix Numeric matrix or data.frame to be clustered. Must be coercible to numeric.
#' @param nColClusters Integer. Number of column clusters to form.
#' @param nRowClusters Integer. Number of row clusters to form.
#' @param seeds Integer vector. Random seeds for multiple GA runs. Default is 1:5.
#' @param verbose Logical. If TRUE, prints progress information. Default is FALSE.
#' @param maxiter Integer. Maximum number of GA iterations. Default is 1000.
#' @param popSize Integer. Population size for the GA. Default is 300.
#' @param pmutation Numeric. Probability of mutation (0-1). Default is 0.05.
#' @param pcrossover Numeric. Probability of crossover (0-1). Default is 0.7.
#' @param elitism Integer. Number of best individuals to preserve. If NULL, uses 5% of popSize.
#' @param interval Integer. Interval for progress monitoring when verbose=TRUE. Default is 100.
#' @param parallel Logical. Whether to use parallel processing. Default is FALSE.
#' @param run Integer. Number of consecutive generations without improvement before stopping.
#'   If NULL, runs for full maxiter iterations.
#'
#' @param validate Logical. If TRUE, run validation on the best partition and attach results
#'   under \code{$validation}. Default FALSE.
#' @param validateCenter Logical. Passed to \code{validateTwomodePartition(center=...)}. Default TRUE.
#' @param validatePerBlock Logical. Passed to \code{validateTwomodePartition(perBlock=...)}. Default TRUE.
#' @param validateMonteCarlo Integer. Number of random partitions for MC p-value.
#'   Passed to \code{validateTwomodePartition(monteCarlo=...)}. Default 0 (disabled).
#' @param validateFixBlockSizes Logical. Keep observed cluster sizes in MC. Default TRUE.
#' @param validateStoreNull Logical. Store full null vector from MC. Default FALSE.
#' @param validateSeed Optional integer seed for the validation step. Default NULL.
#'
#' @return A list of class \code{"twomodeClustering"} containing:
#' \describe{
#'   \item{bestGa}{The best GA object from all runs}
#'   \item{bestFitness}{Best fitness value achieved (negative SSE)}
#'   \item{bestSeed}{Seed that produced the best result}
#'   \item{rowClusters}{Integer vector of row cluster assignments}
#'   \item{colClusters}{Integer vector of column cluster assignments}
#'   \item{control}{List of control parameters used}
#'   \item{validation}{List returned by \code{validateTwomodePartition()} if \code{validate=TRUE}; otherwise NULL}
#' }
#'
#' @details The function runs multiple GA instances with different random seeds and
#' returns the best solution. The fitness function minimizes the sum of squared errors
#' within clusters. Row and column clusters are optimized simultaneously.
#'
#' @references
#' Hageman, J. A., van den Berg, R. A., Westerhuis, J. A., van der Werf, M. J., & Smilde, A. K. (2008).
#' Genetic algorithm based two-mode clustering of metabolomics data. *Metabolomics*, 4, 141–149. \doi{10.1007/s11306-008-0105-7}
#'
#' @seealso \code{\link[GA]{ga}} for the underlying genetic algorithm implementation
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
twomodeClusteringGA <- function(myMatrix,
                                nColClusters,
                                nRowClusters,
                                seeds = 1:5,
                                verbose = FALSE,
                                maxiter = 2000,
                                popSize = 100,
                                pmutation = 0.05,
                                pcrossover = 0.8,
                                elitism = 10,
                                interval = 100,
                                parallel = FALSE,
                                run = NULL,
                                # validation controls
                                validate = FALSE,
                                validateCenter = TRUE,
                                validatePerBlock = TRUE,
                                validateMonteCarlo = 0L,
                                validateFixBlockSizes = TRUE,
                                validateStoreNull = FALSE,
                                validateSeed = NULL) {

  X <- as.matrix(myMatrix)
  if (!is.numeric(X)) stop("myMatrix must be numeric or coercible to a numeric matrix.")
  nR <- nrow(X); nC <- ncol(X)

  if (length(nColClusters) != 1L || length(nRowClusters) != 1L) {
    stop("nColClusters and nRowClusters must be single integers.")
  }

  if (is.null(elitism)) elitism <- max(1L, round(0.05 * popSize))

  lower <- c(rep(1L, nR), rep(1L, nC))
  upper <- c(rep(as.integer(nRowClusters), nR), rep(as.integer(nColClusters), nC))

  fitnessFun <- twomodeFitnessFactory(X)
  monitorFun <- if (isTRUE(verbose)) monitorFactory(interval) else FALSE

  bestGa <- NULL; bestFitness <- -Inf; bestSeed <- NA_integer_

  for (i in seq_along(seeds)) {
    if (isTRUE(verbose)) cat(sprintf("Run %d (seed=%s):\n", i, as.character(seeds[i])))

    set.seed(seeds[i])

    gaArgs <- list(
      type       = "real-valued",
      fitness    = fitnessFun,
      lower      = lower,
      upper      = upper,
      population = gaintegerPopulation,
      mutation   = gaintegerMutation,
      crossover  = gaintegerTwoPointCrossover,
      selection  = GA::gareal_lrSelection,
      monitor    = monitorFun,
      parallel   = parallel,
      seed       = seeds[i],
      maxiter    = maxiter,
      popSize    = popSize,
      pmutation  = pmutation,
      pcrossover = pcrossover,
      elitism    = elitism
    )
    if (!is.null(run) && is.finite(run) && run > 0) gaArgs$run <- as.integer(run)

    gaObj <- do.call(GA::ga, gaArgs)

    fit <- gaObj@fitnessValue
    if (fit > bestFitness) { 
      bestFitness <- fit
      bestGa <- gaObj
      bestSeed <- seeds[i] 
      }
  }

  if (is.null(bestGa)) stop("GA failed to produce a solution.")

  bestSol <- as.numeric(bestGa@solution[1, ])
  bestSol <- pmin(upper, pmax(lower, round(bestSol)))
  rowClusters <- as.integer(bestSol[seq_len(nR)])
  colClusters <- as.integer(bestSol[nR + seq_len(nC)])

  if (isTRUE(verbose)) cat(sprintf("Best run: seed=%s, fitness=%.6f\n", as.character(bestSeed), bestFitness))

  # Optional validation hook
  validation <- NULL
  if (isTRUE(validate)) {
    validation <- validateTwomodePartition(
      myMatrix    = X,
      rowClusters = rowClusters,
      colClusters = colClusters,
      center      = validateCenter,
      perBlock    = validatePerBlock,
      monteCarlo  = as.integer(validateMonteCarlo),
      fixBlockSizes = validateFixBlockSizes,
      storeNull   = validateStoreNull,
      seed        = validateSeed
    )
    if (isTRUE(verbose)) {
      cat(sprintf("Validation: R^2=%.4f | F(%d,%d)=%.4f | p=%.6g\n",
                  validation$r2, validation$dfModel, validation$dfResid,
                  validation$fStat, validation$pValue))
      if (!is.null(validation$mc)) {
        cat(sprintf("Validation MC: p_MC=%.6g (n=%d)\n",
                    validation$mc$pMonteCarlo, validation$mc$nSim))
      }
    }
  }

  results <- list(
    bestGa      = bestGa,
    bestFitness = bestFitness,
    bestSeed    = bestSeed,
    rowClusters = rowClusters,
    colClusters = colClusters,
    control     = list(
      nRowClusters = as.integer(nRowClusters),
      nColClusters = as.integer(nColClusters),
      maxiter = maxiter, popSize = popSize, pmutation = pmutation,
      pcrossover = pcrossover, elitism = elitism, parallel = parallel,
      run = if (is.null(run)) 0L else as.integer(run)
    ),
    validation  = validation
  )

  results$myMatrix <- X
  class(results) <- "twomodeClustering"
  return(results)
}

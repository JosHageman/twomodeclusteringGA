#' One-point crossover for genetic algorithm with integer encoding
#'
#' @description Performs one-point crossover between two parent individuals in the
#' genetic algorithm, exchanging genetic material at a single randomly selected point.
#'
#' @param object GA object containing algorithm parameters.
#' @param parents Integer vector of length 2 containing indices of parent individuals.
#' @param ... Additional arguments (not used).
#'
#' @return List containing:
#'   \describe{
#'     \item{children}{Matrix with two rows representing the offspring}
#'     \item{fitness}{Vector of NA values (fitness will be calculated later)}
#'   }
#'
#' @export
gaintegerOnePointCrossover <- function(object, parents, ...) {  
  parents <- object@population[parents, , drop = FALSE]
  numCols <- ncol(parents)
  children <- parents
  
  crossOverPoint <- sample.int(numCols, size = 1)
  
  swap <- 1:crossOverPoint
  children[1, swap] <- parents[2, swap]
  children[2, swap] <- parents[1, swap]
  
  out <- list(children = children, fitness = rep(NA, 2))
  return(out)
}
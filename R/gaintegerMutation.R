#' Integer mutation for genetic algorithm
#'
#' @description Performs mutation on a genetic algorithm individual by randomly
#' changing cluster assignments with a specified probability.
#'
#' @param object GA object containing algorithm parameters.
#' @param parent Integer index of the parent individual to mutate.
#' @param ... Additional arguments (not used).
#'
#' @return Numeric vector representing the mutated individual.
#'
#' @importFrom stats runif
#' @export
#'
gaintegerMutation <- function(object, parent, ...) {
  maxVal <- object@upper
  mutate <- as.vector(object@population[parent, ])
  numVars <- length(mutate)
  w <- runif(numVars) < object@pmutation
  
  for (i in 1:numVars) {
    if (w[i] == TRUE) {
      mutate[i] <- sample.int(maxVal[i], size = 1, replace = TRUE)
    }
  }  
  return(mutate)
}
#' Integer population initialization for genetic algorithm
#'
#' @description Generates an initial population for the genetic algorithm where each
#' individual represents a clustering solution with integer cluster assignments.
#'
#' @param object GA object containing algorithm parameters.
#' @param ... Additional arguments (not used).
#'
#' @return Matrix where each row represents an individual in the population and
#'   each column represents a cluster assignment.
#'
#' @export
#'
gaintegerPopulation <- function(object, ...) {
  maxVal <- object@upper
  population <- matrix(as.double(NA), nrow = object@popSize, ncol = length(maxVal))
  
  for (i in 1:object@popSize) {
    for (j in 1:length(maxVal)) {
      population[i, j] <- sample.int(maxVal[j], size = 1, replace = TRUE)
    }
  }
  return(population)
}
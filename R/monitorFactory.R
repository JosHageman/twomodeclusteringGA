#' Two mode clustering monitoring function factory for GA progress
#'
#' Creates a monitoring function that prints the current generation and the best fitness score 
#' to the console at specified intervals. Intended for use as a \code{monitor} function in GA runs.
#'
#' @param interval An integer specifying the interval for printing progress updates. 
#' Default is 100 (prints every 100 generations).
#'
#' @return A monitoring function that can be used with GA. The returned function takes a GA object
#' and prints progress information at the specified interval.
#'
#' @examples
#' # Create monitor that prints every 100 generations (default)
#' monitor <- monitorFactory()
#' # ga(..., monitor = monitor)
#' 
#' # Create monitor that prints every 50 generations
#' monitor <- monitorFactory(50)
#' # ga(..., monitor = monitor)
#'
#' @export
monitorFactory <- function(interval = 100) { # nolint: object_name_linter.
  
  force(interval)

  return(function(object) {
    if (object@iter %% interval == 0) {
      cat("Generation:", object@iter,
          "Best fitness:", round(max(object@fitness), 2), "\n")
    }
    
    return(invisible(NULL))
  })
}

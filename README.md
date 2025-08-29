# twomodeclusteringGA

![R-CMD-check](https://github.com/joshageman/twomodeclusteringGA/actions/workflows/R-CMD-check.yaml/badge.svg)
![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
![CRAN status](https://www.r-pkg.org/badges/version/twomodeclusteringGA)
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/twomodeclusteringGA)


**twomodeclusteringGA** is an R package for performing **two-mode clustering (biclustering)** using a Genetic Algorithm (GA).
The package provides a flexible and efficient framework to detect row and column clusters simultaneously in numeric data matrices.

Key features:
- GA-based optimization tailored for two-mode clustering
- Built-in visualization tools for clustered heatmaps
- Functions for summarizing, comparing, and exporting results
- Easy-to-use workflow suitable for small and medium-sized datasets

## Installation

You can install the development version directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("joshageman/twomodeclusteringGA")
```

## Example

A toy dataset is included to demonstrate the functionality:

```r
data("twomodeToy")
myMatrix_s <- scale(twomodeToy)

# Run the GA-based two-mode clustering
result <- twomodeClusteringGA(
  myMatrix     = myMatrix_s,
  nRowClusters = 2,
  nColClusters = 3,
  seeds        = 1,
  maxiter      = 200,
  popSize      = 30,
  validate     = TRUE,
  verbose      = TRUE
)

# Inspect the result
print(result)
summary(result)
myTwomodeResult <- as.data.frame(result)
head(myTwomodeResult)

# Plot the clustered heatmap
plotTwomodeClustering(
  myMatrix  = myMatrix_s,
  result    = result,
  title     = "Two-mode clustering Toy example",
  fixAspect = FALSE
)
```

## References

This package builds on the literature on biclustering and genetic algorithms. Key references:

- Hageman, J. A., Wehrens, R., & Buydens, L. M. C. (2008). *Simplivariate Models: Ideas and First Examples*. PLoS ONE, 3(9), e3259. doi:10.1371/journal.pone.0003259
- Madeira, S. C., & Oliveira, A. L. (2004). *Biclustering algorithms for biological data analysis: a survey*. IEEE/ACM Transactions on Computational Biology and Bioinformatics, 1(1), 24–45. doi:10.1109/TCBB.2004.2

## Contributing

Contributions, suggestions, and bug reports are welcome. Please open an issue or submit a pull request.

## Author

Maintainer: **Jos Hageman** (jos.hageman@wur.nl)
Wageningen University & Research

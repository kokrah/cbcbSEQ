\name{batchSEQ}
\alias{batchSEQ}
\title{a pipeline for batch correction and mean-variance function calculation for count data}
\usage{
  batchSEQ(counts, design, batch, condition,
    lib.size = NULL, verbose = FALSE, plot = FALSE,
    filter = FALSE, ...)
}
\arguments{
  \item{counts}{Count matrix}

  \item{design}{model.matrix for differential expression
  testing}

  \item{batch}{factor indicating batch}

  \item{lib.size}{library sizes, if NULL or missing, uses
  the sum of quantile-normalized counts (default=NULL)}

  \item{verbose}{print extra information}

  \item{plot}{plot the mean-variance fit}

  \item{condition}{factor indicating biological condition}

  \item{filter}{filter genes with low expression}

  \item{...}{pass arguments to internal functions}
}
\value{
  list with components elist (result of calling voom) and
  combatEstimates (batch effect estimates from combat)
}
\description{
  The pipeline includes 1) quantile normalization 2)
  log-transformation of counts 3) combat batch correction
  4) voom calculation of weights for testing from
  mean-variance relationship
}


#' a pipeline for batch correction and mean-variance function calculation for count data
#' 
#' The pipeline includes
#'  1) quantile normalization
#'  2) log-transformation of counts
#'  3) combat batch correction
#'  4) voom calculation of weights for testing from mean-variance relationship
#'  
#' @param counts Count matrix
#' @param design model.matrix for differential expression testing
#' @param batch factor indicating batch
#' @param lib.size library sizes, if NULL or missing, uses the sum of quantile-normalized counts (default=NULL)
#' @param verbose print extra information
#' @param plot plot the mean-variance fit
#' @param condition factor indicating biological condition
#' @param filter filter genes with low expression
#' @param ... pass arguments to internal functions
#' @return list with components elist (result of calling voom) and combatEstimates (batch effect estimates from combat)
#' @export
batchSEQ <- function(counts, design, batch, condition, lib.size=NULL, verbose=FALSE, plot=FALSE, filter=FALSE, ...)
{
  if (filter) {
    # filter samples
    counts <- filterCounts(counts, lib.size, ...)
  }
  
  # quantile normalize counts
  qcounts = qNorm(counts, verbose=verbose)
  
  # convert to log(counts per million)
  res = log2CPM(qcounts, lib.size=lib.size)
  y = res$y
  lib.size = res$lib.size
  
  # run combat on the log(cpm) data
  res = combatMod(y, batch, condition, ...)
  
  y = res$bayesdata
  res = res[c("gamma.star", "delta.star")]
  
  # compute weights from voom
  voomRes = voomMod(y, design, lib.size, plot=plot)

  return(list(elist=voomRes, combatEstimates=res))
}
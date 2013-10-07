#' Filter low expression genes, using independent filtering
#' recommended in Nature Protocols paper
#' 
#' @param counts Count matrix
#' @param lib.size library sizes, if missing or sum uses the total counts per samples
#' @param thresh Minimum counts per million to determine expression
#' @param minSamples Minimum number of samples where gene is required to be expressed. This should be set to the numer of samples in the smallest group of interest.
#' @return Filtered count matrix
#' @export
filterCounts <- function(counts, lib.size = NULL,
                         thresh = 1, 
                         minSamples = 2)
{
  cpms <- 2^logCPM(counts, lib.size=lib.size)
  keep <- rowSums(cpms > thresh) >= minSamples
  counts <- counts[keep,]
  counts
}


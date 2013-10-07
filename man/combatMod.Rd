\name{combatMod}
\alias{combatMod}
\title{Call ComBat function with minor adjustments}
\usage{
  combatMod(dat, batch, mod, noScale = TRUE,
    prior.plots = FALSE)
}
\arguments{
  \item{dat}{Genomic measure matrix (dimensions probe x
  sample)}

  \item{batch}{Batch covariate (multiple batches allowed)}

  \item{mod}{bio factor of interest}

  \item{noScale}{do not adjust scale for batch (use
  default(TRUE) for now)}

  \item{prior.plots}{(Optional) TRUE gives prior plots.
  Used only when noScale is FALSE}
}
\value{
  list containing adjusted data(bayesdata), shift
  adjustments(gamma.star), and scale
  adjustments(delta.star)
}
\description{
  Call ComBat function with minor adjustments
}


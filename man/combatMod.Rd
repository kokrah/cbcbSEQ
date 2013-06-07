\name{combatMod}
\alias{combatMod}
\title{Call ComBat function with minor adjustments}
\usage{
  combatMod(dat, batch, mod, numCovs = NULL,
    par.prior = TRUE, prior.plots = FALSE, scal = FALSE)
}
\arguments{
  \item{dat}{Genomic measure matrix (dimensions probe x
  sample)}

  \item{batch}{Batch covariate (multiple batches allowed)}

  \item{mod}{Model matrix for outcome of interest and other
  covariates besides batch}

  \item{numCovs}{The column numbers of the variables in mod
  to be treated as continuous variables}

  \item{par.prior}{OPTIONAL (use default(TRUE) for now)}

  \item{prior.plots}{(Optional) TRUE gives prior plots}

  \item{scal}{sets batch scale adjustments equal to 1 when
  FALSE(default)}
}
\value{
  list containing adjusted data(bayesdata), shift
  adjustments(gamma.star), and scale
  adjustments(delta.star)
}
\description{
  Call ComBat function with minor adjustments
}


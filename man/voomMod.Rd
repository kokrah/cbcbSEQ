\name{voomMod}
\alias{voomMod}
\title{Call modified voom on ComBatLoc adjusted data}
\usage{
  voomMod(x, design, lib.size, plot = FALSE)
}
\arguments{
  \item{x}{ComBatLoc adjusted data}

  \item{design}{reparameterized design matrix}

  \item{lib.size}{library sizes of quantile counts obtained
  from calling logPCM()}

  \item{plot}{logical indicating whether or not to plot
  mean var trend}
}
\value{
  Elist object
}
\description{
  Call modified voom on ComBatLoc adjusted data
}


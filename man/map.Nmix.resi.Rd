\name{map.Nmix.resi}
\alias{map.Nmix.resi}
\encoding{UTF-8}
\title{
Mapping of residuals
}
\description{
Produces a map of the mean residuals from an N-mixture model fit by function \code{pcount} in \pkg{unmarked}. Used in Section 6.9.3 to produce maps of Switzerland with the residuals for each site.
}
\usage{
map.Nmix.resi(fm, x, y)
}
\arguments{
  \item{fm}{
the fitted model object
}
  \item{x}{
x coordinates of each site
}
  \item{y}{
y coordinates of each site
}
}

\value{
None. Used for its plotting side effects.
}
\references{
Kéry & Royle (2016) \emph{Applied Hierarachical Modeling in Ecology} Section 6.9.3
}
\author{
Marc Kéry
}
\note{
In previous versions, the defaults were \code{x = tits$coordx} and \code{y = tits$coordy}, but those defaults only worked if the data object \code{tits} was in the workspace. To run the code on page 263 of AHM1, you now need to specify the coordinates, eg, \code{map.Nmix.resi(fm5, x = tits$coordx, y = tits$coordy)}.
}

\examples{
# Not possible as would require a fitted model object. See section 6.9.3 of AHM1.
}

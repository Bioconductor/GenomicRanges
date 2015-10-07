\name{phicoef}

\alias{phicoef}

\title{Calculate the "phi coefficient" between two binary variables}

\description{
  The \code{phicoef} function calculates the "phi coefficient" between
  two binary variables.
}

\usage{
phicoef(x, y=NULL)
}

\arguments{
  \item{x, y}{
    Two logical vectors of the same length.
    If \code{y} is not supplied, \code{x} must be a 2x2 integer matrix
    (or an integer vector of length 4) representing the contingency table
    of two binary variables.
  }
}

\value{
  The "phi coefficient" between the two binary variables.
  This is a single numeric value ranging from -1 to +1.
}

\references{
  \url{http://en.wikipedia.org/wiki/Phi_coefficient}
}

\author{
  H. Pagès
}

\examples{
set.seed(33)
x <- sample(c(TRUE, FALSE), 100, replace=TRUE)
y <- sample(c(TRUE, FALSE), 100, replace=TRUE)
phicoef(x, y)
phicoef(rep(x, 10), c(rep(x, 9), y))

stopifnot(phicoef(table(x, y)) == phicoef(x, y))
stopifnot(phicoef(y, x) == phicoef(x, y))
stopifnot(phicoef(x, !y) == - phicoef(x, y))
stopifnot(phicoef(x, x) == 1)
}

\keyword{manip}

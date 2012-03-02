\name{strand-utils}

\alias{strand-utils}

\alias{strand,missing-method}
\alias{strand,character-method}
\alias{strand,factor-method}
\alias{strand,integer-method}
\alias{strand,logical-method}
\alias{strand,DataTable-method}

\title{Strand utilities}

\description{
  Five useful \code{strand} methods.
}

\usage{
strand(x)
}

\arguments{
  \item{x}{The object from which to obtain a strand factor, can be
    missing.}
}

\details{
  If \code{x} is missing, returns an empty factor with the
  standard levels that any strand factor should have: \code{+},
  \code{-}, and \code{*}.

  If \code{x} is a character vector or factor, it is coerced to a factor
  with the levels listed above.

  If \code{x} is an integer vector, it is coerced to a factor
  with the levels listed above. \code{1} and \code{-1}
  values in \code{x} are mapped to the \code{+} and \code{-}
  levels respectively. \code{NA}s in \code{x} produce \code{NA}s
  in the result.

  If \code{x} is a logical vector, it is coerced to a factor
  with the levels listed above. \code{FALSE} and \code{TRUE}
  values in \code{x} are mapped to the \code{+} and \code{-}
  levels respectively. \code{NA}s in \code{x} produce \code{NA}s
  in the result.

  If \code{x} inherits from \code{DataTable}, the \code{"strand"}
  column is returned as a factor with the levels listed above. If
  \code{x} has no \code{"strand"} column, this return value is
  populated with \code{NA}s.
}

\author{ Michael Lawrence }

\seealso{
  \code{\link[BiocGenerics]{strand}}
}

\examples{
strand()
strand(c("+", "-", NA, "*"))
strand(c(-1L, 1L, NA, -1L, NA))
strand(c(FALSE, FALSE, TRUE, NA, TRUE, FALSE))
}

\keyword{methods}

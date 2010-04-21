\name{strand}
\docType{methods}
\alias{strand-methods}

\alias{strand}
\alias{strand,missing-method}
\alias{strand,character-method}
\alias{strand,factor-method}
\alias{strand,integer-method}
\alias{strand,logical-method}
\alias{strand,DataTable-method}
\alias{strand<-}

\title{Accessing strand information}

\description{The \code{strand} generic is meant as an accessor for
  strand information. Four methods are defined by the GenomicRanges
  package, described below.
}

\usage{strand(x)}

\arguments{
  \item{x}{The object from which to obtain a strand factor, can be
    missing.}
}

\details{
  If \code{x} is missing, returns an empty factor with the
  standard levels that any strand factor should have: \code{+},
  \code{-}, and \code{*} (for either).

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

\examples{
strand()
strand(c("+", "-", NA, "*"))
strand(c(-1L, 1L, NA, -1L, NA))
strand(c(FALSE, FALSE, TRUE, NA, TRUE, FALSE))
}

\keyword{methods}

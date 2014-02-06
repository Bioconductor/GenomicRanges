\name{strand-utils}

\alias{strand-utils}

\alias{strand,missing-method}
\alias{strand,NULL-method}
\alias{strand,character-method}
\alias{strand,factor-method}
\alias{strand,integer-method}
\alias{strand,logical-method}
\alias{strand,Rle-method}
\alias{strand,DataTable-method}
\alias{strand<-,DataTable-method}

\title{Strand utilities}

\description{
  Some useful \code{strand} methods.
}

\usage{
\S4method{strand}{missing}(x)
\S4method{strand}{character}(x)
\S4method{strand}{factor}(x)
\S4method{strand}{integer}(x)
\S4method{strand}{logical}(x)
\S4method{strand}{Rle}(x)
\S4method{strand}{DataTable}(x)
\S4method{strand}{DataTable}(x) <- value
}

\arguments{
  \item{x}{The object from which to obtain a strand factor, can be
    missing.}
  \item{value}{Replacement value for the strand.}
}

\details{
  If \code{x} is missing, returns an empty factor with the
  "standard strand levels" i.e. \code{+}, \code{-}, and \code{*}.

  If \code{x} is a character vector or factor, it is coerced to a
  factor with the levels listed above. \code{NA} values in \code{x}
  are not accepted.

  If \code{x} is an integer vector, it is coerced to a factor
  with the levels listed above. \code{1}, \code{-1}, and \code{NA}
  values in \code{x} are mapped to the \code{+}, \code{-}, and
  \code{*} levels respectively.

  If \code{x} is a logical vector, it is coerced to a factor
  with the levels listed above. \code{FALSE}, \code{TRUE}, and \code{NA}
  values in \code{x} are mapped to the \code{+}, \code{-}, and
  \code{*} levels respectively.

  If \code{x} is a character-, factor-, integer-, or
  logical-\link[IRanges]{Rle}, it is transformed with
  \code{runValue(x) <- strand(runValue(x))} and returned.

  If \code{x} inherits from \code{DataTable}, the \code{"strand"}
  column is passed thru \code{strand()} and returned.
  If \code{x} has no \code{"strand"} column, this return value is
  populated with \code{*}s.
}

\value{
  A factor (or factor-\link[IRanges]{Rle}) with the "standard strand levels"
  (i.e. \code{+}, \code{-}, and \code{*}) and no NAs.
}

\author{Michael Lawrence and H. Pages}

\seealso{
  \code{\link[BiocGenerics]{strand}}
}

\examples{
strand()

x1 <- c("-", "*", "*", "+", "-", "*")
x2 <- factor(c("-", "-", "+", "-"))
x3 <- c(-1L, NA, NA, 1L, -1L, NA)
x4 <- c(TRUE, NA, NA, FALSE, TRUE, NA)

strand(x1)
strand(x2)
strand(x3)
strand(x4)

strand(Rle(x1))
strand(Rle(x2))
strand(Rle(x3))
strand(Rle(x4))

strand(DataFrame(score=2:-3))
strand(DataFrame(score=2:-3, strand=x3))
strand(DataFrame(score=2:-3, strand=Rle(x3)))

## Sanity checks:
target <- strand(x1)
stopifnot(identical(target, strand(x3)))
stopifnot(identical(target, strand(x4)))

stopifnot(identical(Rle(strand(x1)), strand(Rle(x1))))
stopifnot(identical(Rle(strand(x2)), strand(Rle(x2))))
stopifnot(identical(Rle(strand(x3)), strand(Rle(x3))))
stopifnot(identical(Rle(strand(x4)), strand(Rle(x4))))
}

\keyword{methods}

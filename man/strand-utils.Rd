\name{strand-utils}

\alias{strand-utils}

\alias{strand,missing-method}
\alias{strand,NULL-method}
\alias{strand,character-method}
\alias{strand,factor-method}
\alias{strand,integer-method}
\alias{strand,logical-method}
\alias{strand,Rle-method}
\alias{strand,RleList-method}
\alias{strand,DataFrame-method}
\alias{strand<-,DataFrame,ANY-method}

\alias{invertStrand,NULL-method}
\alias{invertStrand,character-method}
\alias{invertStrand,factor-method}
\alias{invertStrand,integer-method}
\alias{invertStrand,logical-method}
\alias{invertStrand,Rle-method}
\alias{invertStrand,RleList-method}

\title{Strand utilities}

\description{
  A bunch of useful \code{strand} and \code{invertStrand} methods.
}

\usage{
\S4method{strand}{missing}(x)
\S4method{strand}{character}(x)
\S4method{strand}{factor}(x)
\S4method{strand}{integer}(x)
\S4method{strand}{logical}(x)
\S4method{strand}{Rle}(x)
\S4method{strand}{RleList}(x)
\S4method{strand}{DataFrame}(x)
\S4method{strand}{DataFrame,ANY}(x) <- value

\S4method{invertStrand}{character}(x)
\S4method{invertStrand}{factor}(x)
\S4method{invertStrand}{integer}(x)
\S4method{invertStrand}{logical}(x)
\S4method{invertStrand}{Rle}(x)
\S4method{invertStrand}{RleList}(x)
}

\arguments{
  \item{x}{
    The object from which to obtain a \emph{strand factor}, \emph{strand
    factor \link[S4Vectors]{Rle}}, or \emph{strand factor \link{RleList}}
    object. Can be missing. See Details and Value sections below for more
    information.
  }
  \item{value}{
    Replacement value for the strand.
  }
}

\details{
  All the \code{strand} and \code{invertStrand} methods documented
  here return either a \emph{strand factor}, \emph{strand factor
  \link[S4Vectors]{Rle}}, or \emph{strand factor \link{RleList}} object.
  These are factor, factor-\link[S4Vectors]{Rle}, or factor-\link{RleList}
  objects containing the "standard strand levels" (i.e. \code{+}, \code{-},
  and \code{*}) and no NAs.
}

\value{
  All the \code{strand} and \code{invertStrand} methods documented here
  return an object that is \emph{parallel} to input object \code{x} when
  \code{x} is a character, factor, integer, logical, \link[S4Vectors]{Rle},
  or \link{RleList} object.

  For the \code{strand} methods:
  \itemize{
    \item If \code{x} is missing, returns an empty factor with the
          "standard strand levels" i.e. \code{+}, \code{-}, and \code{*}.

    \item If \code{x} is a character vector or factor, it is coerced to a
          factor with the levels listed above. \code{NA} values in \code{x}
          are not accepted.

    \item If \code{x} is an integer vector, it is coerced to a factor
          with the levels listed above. \code{1}, \code{-1}, and \code{NA}
          values in \code{x} are mapped to the \code{+}, \code{-}, and
          \code{*} levels respectively.

    \item If \code{x} is a logical vector, it is coerced to a factor
          with the levels listed above. \code{FALSE}, \code{TRUE}, and
          \code{NA} values in \code{x} are mapped to the \code{+}, \code{-},
          and \code{*} levels respectively.

    \item If \code{x} is a character-, factor-, integer-, or
          logical-\link[S4Vectors]{Rle}, it is transformed with
          \code{runValue(x) <- strand(runValue(x))} and returned.

    \item If \code{x} is an \link{RleList} object, each list element in
          \code{x} is transformed by calling \code{strand()} on it and
          the resulting \link{RleList} object is returned. More precisely
          the returned object is \code{endoapply(x, strand)}.
          Note that in addition to being \emph{parallel} to \code{x}, this
          object also has the same \emph{shape} as \code{x} (i.e. its list
          elements have the same lengths as in \code{x}).
    
    \item If \code{x} is a \code{DataFrame} object, the \code{"strand"}
          column is passed thru \code{strand()} and returned.
          If \code{x} has no \code{"strand"} column, this return value is
          populated with \code{*}s.
  }

  Each \code{invertStrand} method returns the same object as its corresponding
  \code{strand} method but with \code{"+"} and \code{"-"} switched.
}

\author{M. Lawrence and H. Pagès}

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
invertStrand(x1)
strand(x2)
invertStrand(x2)
strand(x3)
invertStrand(x3)
strand(x4)
invertStrand(x4)

strand(Rle(x1))
invertStrand(Rle(x1))
strand(Rle(x2))
invertStrand(Rle(x2))
strand(Rle(x3))
invertStrand(Rle(x3))
strand(Rle(x4))
invertStrand(Rle(x4))

x5 <- RleList(x1, character(0), as.character(x2))
strand(x5)
invertStrand(x5)

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

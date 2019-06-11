\name{GRangesFactor-class}
\docType{class}

\alias{class:GRangesFactor}
\alias{GRangesFactor-class}
\alias{GRangesFactor}

\alias{FactorToClass,GRanges-method}

\alias{seqnames,GRangesFactor-method}
\alias{start,GRangesFactor-method}
\alias{end,GRangesFactor-method}
\alias{width,GRangesFactor-method}
\alias{pos,GRangesFactor-method}
\alias{strand,GRangesFactor-method}
\alias{seqinfo,GRangesFactor-method}
\alias{granges,GRangesFactor-method}
\alias{ranges,GRangesFactor-method}

\alias{coerce,ANY,GRangesFactor-method}
\alias{coerce,Factor,GRanges-method}

\alias{show,GRangesFactor-method}

\title{GRangesFactor objects}

\description{
  A GRangesFactor object is a \link[S4Vectors]{Factor} derivative where the
  levels are a \link{GRanges} object.

  See \code{?\link[S4Vectors]{Factor}} and in the \pkg{S4Vectors} package
  for general information about \link[S4Vectors]{Factor} objects.
}

\usage{
GRangesFactor(x, levels=GRanges(), index=NULL, ...)  # constructor function
}

\arguments{
  \item{x, levels}{
    Like with the \code{\link[S4Vectors]{Factor}()} constructor function,
    at least one of \code{x} and \code{levels} must be specified.
    If \code{index} is \code{NULL}, both can be specified.

    When \code{levels} is specified, it must be a \link{GRanges} object (or
    derivative) with no duplicate ranges (i.e. \code{anyDuplicated(levels)}
    must return 0).

    When \code{x} and \code{levels} are both specified, they should
    both be \link{GRanges} objects or \link{GRanges} derivatives of
    the same class, and all the elements in \code{x} must be represented
    in \code{levels} (i.e. the integer vector returned by
    \code{match(x, levels)} should contain no \code{NA}s).
  }
  \item{index}{
    \code{NULL} or an integer (or numeric) vector of valid positive indices
    (no \code{NA}s) into \code{levels}.
  }
  \item{...}{
    Optional metadata columns.
  }
}

\details{
  Like with the \code{\link[S4Vectors]{Factor}()} constructor function,
  there are 4 different ways to use the \code{GRangesFactor()}
  constructor function. See Details section in the man page for
  \link[S4Vectors]{Factor} objects for more information.
}

\value{
  A GRangesFactor object.
}

\section{Accessors}{
  GRangesFactor objects support the accessors documented in the
  man page for \link[S4Vectors]{Factor} objects.

  In addition, the following getters are supported for convenience:
  \code{seqnames()}, \code{start()}, \code{end()}, \code{width()},
  \code{strand()}, \code{seqinfo()}, \code{granges()}, and \code{ranges()}.
  When called on GRangesFactor object \code{x}, they all behave as if they
  were called on \code{unfactor(x)}.
}

\section{Decoding a Factor}{
  Because a GRangesFactor object \code{x} is a \link[S4Vectors]{Factor}
  derivative, \code{unfactor(x)} can be used to \emph{decode} it.
  \code{unfactor(x)} returns an object of the same class as \code{levels(x)}
  (i.e. a \link{GRanges} object or derivative) and same length as \code{x}.

  See \code{?\link[S4Vectors]{unfactor}} for more information.
}

\section{Coercion}{
  GRangesFactor objects support the coercions documented in the
  man page for \link[S4Vectors]{Factor} objects.

  In addition, coercion back and forth between \link{GRanges} and
  GRangesFactor is supported via \code{as(x, "GRanges")} and
  \code{as(x, "GRangesFactor")}.
}

\section{Subsetting}{
  A GRangesFactor object can be subsetted with \code{[}, like a
  \link[S4Vectors]{Factor} object.
}

\section{Concatenation}{
  2 or more GRangesFactor objects can be concatenated with \code{c()}.
  The result of this concatenation is another GRangesFactor object.

  See Concatenation section in \code{?\link[S4Vectors]{Factor}}.
}

\section{Comparing & ordering}{
  See Comparing & Ordering section in \code{?\link[S4Vectors]{Factor}}.
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item \link{GRanges} objects.

    \item \link[S4Vectors]{Factor} objects in the \pkg{S4Vectors} package
          for the parent class of GRangesFactor.

    \item \code{\link[BiocGenerics]{anyDuplicated}} in the \pkg{BiocGenerics}
          package.
  }
}

\examples{
showClass("GRangesFactor")  # GRangesFactor extends Factor

## ---------------------------------------------------------------------
## CONSTRUCTOR & ACCESSORS
## ---------------------------------------------------------------------
set.seed(123)
ir0 <- IRanges(sample(5, 8, replace=TRUE), width=10, names=letters[1:8])
gr0 <- GRanges("chrA", ir0, ID=paste0("ID", 1:8))

## Use explicit levels:
gr1 <- GRanges("chrA", IRanges(1:6, width=10))
F1 <- Factor(gr0, levels=gr1)
F1
length(F1)
names(F1)
levels(F1)  # gr1
nlevels(F1)
as.integer(F1)  # encoding

## If we don't specify the levels, they'll be set to unique(gr0):
F2 <- Factor(gr0)
F2
length(F2)
names(F2)
levels(F2)  # unique(gr0)
nlevels(F2)
as.integer(F2)

## ---------------------------------------------------------------------
## DECODING
## ---------------------------------------------------------------------
unfactor(F1)

stopifnot(identical(gr0, unfactor(F1)))
stopifnot(identical(gr0, unfactor(F2)))

unfactor(F1, use.names=FALSE)
unfactor(F1, ignore.mcols=TRUE)

## ---------------------------------------------------------------------
## COERCION
## ---------------------------------------------------------------------
F2b <- as(gr0, "Factor")  # same as Factor(gr0)
stopifnot(identical(F2, F2b))

as.factor(F2)
as.factor(F1)

as.character(F1)  # same as unfactor(as.factor(F1)),
                  # and also same as as.character(unfactor(F1))

## On an ordinary factor 'f', 'as(f, "Factor")' and 'Factor(f)' are
## NOT the same:
f <- factor(sample(letters, 500, replace=TRUE), levels=letters)
as(f, "Factor")  # same levels as 'f'
Factor(f)        # levels **are** 'f'!

stopifnot(identical(f, as.factor(as(f, "Factor"))))

## ---------------------------------------------------------------------
## CONCATENATION
## ---------------------------------------------------------------------
gr3 <- GRanges("chrA", IRanges(c(5, 2, 8:6), width=10))
F3 <- Factor(levels=gr3, index=2:4)
F13 <- c(F1, F3)
F13
levels(F13)

stopifnot(identical(c(unfactor(F1), unfactor(F3)), unfactor(F13)))

## ---------------------------------------------------------------------
## COMPARING & ORDERING
## ---------------------------------------------------------------------
\dontrun{
  F1 == F2  # Error! (F1 and F2 don't have the same levels)
}

order(F1)  # same as order(as.integer(F1))
order(F2)  # same as order(as.integer(F2))

## The levels of the Factor influence the order of the table:
table(F1)
table(F2)
}
\keyword{methods}
\keyword{classes}

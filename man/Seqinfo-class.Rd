\name{Seqinfo-class}
\docType{class}

% Classes:
\alias{class:Seqinfo}
\alias{Seqinfo-class}
\alias{Seqinfo}

% Methods:
\alias{length,Seqinfo-method}
\alias{seqnames,Seqinfo-method}
\alias{seqnames<-,Seqinfo-method}
\alias{names,Seqinfo-method}
\alias{names<-,Seqinfo-method}
\alias{seqlengths,Seqinfo-method}
\alias{seqlengths<-,Seqinfo-method}
\alias{isCircular,Seqinfo-method}
\alias{isCircular<-,Seqinfo-method}
\alias{as.data.frame,Seqinfo-method}
\alias{show,Seqinfo-method}
\alias{merge,Seqinfo,missing-method}
\alias{merge,Seqinfo,NULL-method}
\alias{merge,NULL,Seqinfo-method}
\alias{merge,Seqinfo,Seqinfo-method}


\title{Seqinfo objects}

\description{
  A Seqinfo object is a table-like object that contains basic information
  about a set of genomic sequences. The table has 1 row per sequence and
  1 column per sequence attribute. Currently the only attributes are the
  length and circularity flag of the sequence but more attributes might
  be added in the future.
}

\details{
  Typically Seqinfo objects are not used directly but are part of
  higher level objects.
}

\section{Constructor}{
  \describe{
    \item{}{
        \code{Seqinfo(seqnames, seqlengths=NA, isCircular=NA)}:
        Creates a Seqinfo object.
    }
  }
}

\section{Accessor methods}{
  In the code snippets below, \code{x} is a Seqinfo object.

  \describe{
    \item{}{
      \code{length(x)}:
      Gets the number of sequences in \code{x}.
    }
    \item{}{
      \code{seqnames(x)}, \code{seqnames(x) <- value}:
      Gets/sets the names of the sequences in \code{x}.
      Those names must be non-NA, non-empty and unique.
      They are also called the \emph{keys} of the Seqinfo object.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
      Same as \code{seqnames(x)} and \code{seqnames(x) <- value}.
    }
    \item{}{
      \code{seqlengths(x)}, \code{seqlengths(x) <- value}:
      Gets/sets the lengths of the sequences in \code{x}.
    }
    \item{}{
      \code{isCircular(x)}, \code{isCircular(x) <- value}:
      Gets/sets the circularity flags of the sequences in \code{x}.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a Seqinfo object.

  \describe{
    \item{}{
      \code{as.data.frame(x)}:
      Turns \code{x} into a data frame.
    }
  }
}

\section{Combining Seqinfo objects}{
  There are no \code{c} or \code{rbind} method for Seqinfo objects.
  Both would be expected to just append the rows in \code{y} to the rows
  in \code{x} resulting in an object of length \code{length(x) + length(y)}.
  But that would tend to break the constraint that the seqnames of a Seqinfo
  object must be unique keys.

  So instead, a \code{merge} method is provided.

  In the code snippet below, \code{x} and \code{y} are Seqinfo objects.

  \describe{
    \item{}{
      \code{merge(x, y)}:
      Merge \code{x} and \code{y} into a single Seqinfo object where the
      keys (aka the seqnames) are \code{union(seqnames(x), seqnames(y))}.
      If a row in \code{y} has the same key as a row in \code{x}, and if
      the 2 rows contain compatible information (NA values are compatible
      with anything), then they are merged into a single row in the result.
      If they cannot be merged (because they contain different seqlengths
      and/or circularity flags), then an error is raised.
    }
  }
}

\author{H. Pages}

\seealso{
  \code{\link{seqinfo}}
}

\examples{
  x <- Seqinfo(seqnames=c("chr1", "chr2", "chr3", "chrM"),
               seqlengths=c(100, 200, NA, 15),
               isCircular=c(NA, FALSE, FALSE, TRUE))
  x
  y <- Seqinfo(seqnames=c("chr3", "chr4", "chrM"),
               seqlengths=c(300, NA, 15))
  y
  merge(x, y)  # rows for chr3 and chrM are merged

  ## Note that, strictly speaking, merging 2 Seqinfo objects is not
  ## a commutative operation, i.e., in general 'z1 <- merge(x, y)'
  ## is not identical to 'z2 <- merge(y, x)'. However 'z1' and 'z2'
  ## are guaranteed to contain the same information (i.e. the same
  ## rows, but typically not in the same order):
  merge(y, x)

  ## This contradicts what 'x' says about circularity of chr3 and chrM:
  isCircular(y)[c("chr3", "chrM")] <- c(TRUE, FALSE)
  y
  if (interactive()) {
    merge(x, y)  # raises an error
  }
}

\keyword{methods}
\keyword{classes}

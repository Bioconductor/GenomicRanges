\name{SeqInfo-class}
\docType{class}

% Classes:
\alias{class:SeqInfo}
\alias{SeqInfo-class}
\alias{SeqInfo}

% Methods:
\alias{length,SeqInfo-method}
\alias{seqnames,SeqInfo-method}
\alias{names,SeqInfo-method}
\alias{seqlengths,SeqInfo-method}
\alias{isCircular}
\alias{isCircular,SeqInfo-method}
\alias{as.data.frame,SeqInfo-method}
\alias{show,SeqInfo-method}


\title{SeqInfo objects}

\description{
  A SeqInfo object is a data.frame-like object that contains basic
  information about a set of genomic sequences. Currently only the
  length and circularity flag of each sequence is stored but more
  information might be added in the future.
}

\details{
  Typically SeqInfo objects are not used directly but are part of
  higher level objects.
}

\section{Constructor}{
  \describe{
    \item{}{
        \code{SeqInfo(seqnames, seqlengths=NA, isCircular=NA)}:
        Creates a SeqInfo object.
    }
  }
}

\section{Accessor methods}{
  In the code snippets below, \code{x} is a SeqInfo object.

  \describe{
    \item{}{
      \code{length(x)}:
      Gets the number of sequences in \code{x}.
    }
    \item{}{
      \code{seqnames(x)}, \code{names(x)}:
      Gets the names of the sequences in \code{x}.
    }
    \item{}{
      \code{seqlengths(x)}:
      Gets the lengths of the sequences in \code{x}.
    }
    \item{}{
      \code{isCircular(x)}:
      Gets the circularity flag of the sequences in \code{x}.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a SeqInfo object.

  \describe{
    \item{}{
      \code{as.data.frame(x)}:
      Turns \code{x} into a data frame.
    }
  }
}

\author{H. Pages}

\examples{
  ## Coming soon...
}

\keyword{methods}
\keyword{classes}

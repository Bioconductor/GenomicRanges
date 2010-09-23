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
\alias{isCircular}
\alias{isCircular<-}
\alias{isCircular,Seqinfo-method}
\alias{isCircular<-,Seqinfo-method}
\alias{as.data.frame,Seqinfo-method}
\alias{show,Seqinfo-method}


\title{Seqinfo objects}

\description{
  A Seqinfo object is a data.frame-like object that contains basic
  information about a set of genomic sequences. Currently only the
  length and circularity flag of each sequence is stored but more
  information might be added in the future.
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
      Gets/sets the circularity flag of the sequences in \code{x}.
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

\author{H. Pages}

\examples{
  ## Coming soon...
}

\keyword{methods}
\keyword{classes}

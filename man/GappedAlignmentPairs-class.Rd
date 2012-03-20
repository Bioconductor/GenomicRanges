\name{GappedAlignmentPairs-class}
\docType{class}

% Class:
\alias{class:GappedAlignmentPairs}
\alias{GappedAlignmentPairs-class}

% Constructors:
\alias{GappedAlignmentPairs}
\alias{readGappedAlignmentPairs}

% Accessors:
\alias{length,GappedAlignmentPairs-method}
\alias{names,GappedAlignmentPairs-method}
\alias{names<-,GappedAlignmentPairs-method}
\alias{seqnames,GappedAlignmentPairs-method}
\alias{strand,GappedAlignmentPairs-method}
\alias{elementMetadata,GappedAlignmentPairs-method}
\alias{elementMetadata<-,GappedAlignmentPairs-method}
\alias{seqinfo,GappedAlignmentPairs-method}
\alias{seqinfo<-,GappedAlignmentPairs-method}
\alias{first}
\alias{first,GappedAlignmentPairs-method}
\alias{last}
\alias{last,GappedAlignmentPairs-method}
\alias{left}
\alias{left,GappedAlignmentPairs-method}
\alias{right}
\alias{right,GappedAlignmentPairs-method}
\alias{isProperPair}
\alias{isProperPair,GappedAlignmentPairs-method}

% Vector methods:
\alias{[,GappedAlignmentPairs-method}

% List methods:
\alias{[[,GappedAlignmentPairs-method}
\alias{unlist,GappedAlignmentPairs-method}

% Coercion:
\alias{coerce,GappedAlignmentPairs,GRangesList-method}

% "show" method:
\alias{show,GappedAlignmentPairs-method}


\title{GappedAlignmentPairs objects}

\description{
  The GappedAlignmentPairs class is a container for "alignment pairs".
  All "alignment pairs" are made of a "first" alignment and a "last" alignment,
  and are typically representing the hits of a set of paired-end reads to
  a reference genome. In a given pair, the "first" alignment represents the
  hit of the first end of the read (aka "first segment in the template", using
  SAM Spec terminology), and the "last" alignment represents the hit of the
  second end of the read (aka "last segment in the template", using SAM Spec
  terminology).
}

\details{
  Typically, a GappedAlignmentPairs object will be created by loading
  records from a BAM (or SAM) file containing aligned paired-end reads.
  Each element in the resulting object will be obtained by pairing 2
  records. See the \code{GappedAlignmentPairs} function below for more
  information.
}

\section{Constructors}{
  \describe{
    \item{}{
      \code{readGappedAlignmentPairs(file, format="BAM", use.names=FALSE, ...)}:
      Read a file containing paired-end reads as a GappedAlignmentPairs
      object.
      By default (i.e. \code{use.names=FALSE}), the resulting object has no
      names. If \code{use.names} is \code{TRUE}, then the names are
      constructed from the query template names (QNAME field in a SAM/BAM
      file). Note that the 2 records in a pair of records have the same QNAME.

      Note that this function is just a front-end that delegates to the
      format-specific back-end function specified via the \code{format}
      argument. The \code{use.names} argument and any extra argument are
      passed to the back-end function.
      Only the BAM format is supported for now. Its back-end is the
      \code{\link[Rsamtools]{readBamGappedAlignmentPairs}} function
      defined in the Rsamtools package.
      See \code{?\link[Rsamtools]{readBamGappedAlignmentPairs}} for
      more information (you might need to install and load the Rsamtools
      package first).
    }
    \item{}{
      \code{GappedAlignmentPairs(first, last, isProperPair, names = NULL)}:
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{x} is a GappedAlignmentPairs object.

  \describe{
    \item{}{
      \code{length(x)}:
      Returns the number of alignment pairs in \code{x}.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
    }
    \item{}{
      \code{seqnames(x)}, \code{seqnames(x) <- value}:
    }
    \item{}{
      \code{strand(x)}, \code{strand(x) <- value}:
    }
    \item{}{
      \code{seqinfo(x)}, \code{seqinfo(x) <- value}:
      Gets or sets the information about the underlying sequences.
      \code{value} must be a \link{Seqinfo} object.
    }
    \item{}{
      \code{seqlevels(x)}, \code{seqlevels(x) <- value}:
      Gets or sets the sequence levels.
      \code{seqlevels(x)} is equivalent to \code{seqlevels(seqinfo(x))}
      or to \code{levels(rname(x))}, those 2 expressions being
      guaranteed to return identical character vectors on a
      GappedAlignmentPairs object. \code{value} must be a character vector
      with no NAs.
      See \code{?\link{seqlevels}} for more information.
    }
    \item{}{
      \code{seqlengths(x)}, \code{seqlengths(x) <- value}:
      Gets or sets the sequence lengths.
      \code{seqlengths(x)} is equivalent to \code{seqlengths(seqinfo(x))}.
      \code{value} can be a named non-negative integer or numeric vector
      eventually with NAs.
    }
    \item{}{
      \code{isCircular(x)}, \code{isCircular(x) <- value}:
      Gets or sets the circularity flags.
      \code{isCircular(x)} is equivalent to \code{isCircular(seqinfo(x))}.
      \code{value} must be a named logical vector eventually with NAs.
    }
    \item{}{
      \code{genome(x)}, \code{genome(x) <- value}:
      Gets or sets the genome identifier or assembly name for each sequence.
      \code{genome(x)} is equivalent to \code{genome(seqinfo(x))}.
      \code{value} must be a named character vector eventually with NAs.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a GappedAlignmentPairs object.

  \describe{
    \item{}{
      \code{as(x, "GRangesList")}:
    }
  }
}

\section{Subsetting and related operations}{
  In the code snippets below, \code{x} is a GappedAlignmentPairs object.

  \describe{
    \item{}{
      \code{x[i]}:
    }
  }
}

\author{
  H. Pages
}

\seealso{
}

\examples{
}

\keyword{methods}
\keyword{classes}

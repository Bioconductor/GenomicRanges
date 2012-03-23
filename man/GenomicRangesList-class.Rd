\name{GenomicRangesList-class}
\docType{class}

% Class:
\alias{class:GenomicRangesList}
\alias{class:SimpleGenomicRangesList}
\alias{GenomicRangesList-class}
\alias{SimpleGenomicRangesList-class}

% Constructors:
\alias{GenomicRangesList}

% Coercion
\alias{coerce,RangedDataList,GenomicRangesList-method}

\title{GenomicRangesList objects}

\description{
  A \code{GenomicRangesList} is a \code{\link[IRanges:List-class]{List}}
  of \code{\linkS4class{GenomicRanges}}. It is a virtual class;
  \code{SimpleGenomicRangesList} is the basic implementation.  The
  subclass \code{\linkS4class{GRangesList}} provides special behavior
  and is particularly efficient for storing a large number of elements.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GenomicRangesList(\dots)}: Constructs a
      \code{SimpleGenomicRangesList} with elements
      taken from the arguments in \dots. If the only argument is a list,
      the elements are taken from that list.
    }
  }
}

\section{Coercion}{
  \describe{
    \item{}{
      \code{as(from, "GenomicRangesList")}: Supported \code{from} types
      include:
      \describe{
        \item{RangedDataList}{Each element of \code{from} is coerced to
          a \code{GenomicRanges}.}
      }
    }
  }
}

\author{Michael Lawrence}

\seealso{
  \code{\linkS4class{GRangesList}}, which differs from
  \code{SimpleGenomicRangesList} in that the \code{GRangesList} treats
  its elements as single, compound ranges, particularly in overlap
  operations. \code{SimpleGenomicRangesList} is just a barebones list
  for now, without that compound semantic.
}

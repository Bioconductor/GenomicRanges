\name{GenomicRangesList-class}
\docType{class}

% Class:
\alias{class:GenomicRangesList}
\alias{GenomicRangesList-class}

% Constructors:
\alias{GenomicRangesList}

\title{GenomicRangesList objects}

\description{
  A \code{GenomicRangesList} is a
  \code{\link[IRanges:SimpleList]{SimpleList}} of
  \code{\linkS4class{GenomicRanges}} objects. Currently, this is more
  or less a mechanism for formally storing a list of such objects. There
  is no special behavior like that provided by \code{\linkS4class{GRangesList}}.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GenomicRangesList(\dots)}: Constructs a list with elements
      taken from the arguments in \dots. If the only argument is a list,
      the elements are taken from that list.
    }
  }
}

\author{Michael Lawrence}

\seealso{
  \code{\linkS4class{GRangesList}}, which differs from
  \code{GenomicRangesList} in that the \code{GRangesList} treats its
  elements as single, compound ranges, particularly in overlap
  operations. \code{GenomicRangesList} is just a barebones list for now,
  without that compound semantic.
}

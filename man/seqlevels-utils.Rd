\name{utils}
\docType{methods}
\alias{keepSeqlevels}
\alias{keepSeqlevels,GenomicRanges,GenomicRanges-method}
\alias{keepSeqlevels,GenomicRanges,GRangesList-method}
\alias{keepSeqlevels,GenomicRanges,GAlignments-method}
\alias{keepSeqlevels,GenomicRanges,character-method}
\alias{keepSeqlevels,GRangesList,GenomicRanges-method}
\alias{keepSeqlevels,GRangesList,GRangesList-method}
\alias{keepSeqlevels,GRangesList,GAlignments-method}
\alias{keepSeqlevels,GRangesList,character-method}
\alias{keepSeqlevels,GAlignments,GenomicRanges-method}
\alias{keepSeqlevels,GAlignments,GRangesList-method}
\alias{keepSeqlevels,GAlignments,GAlignments-method}
\alias{keepSeqlevels,GAlignments,character-method}

\alias{renameSeqlevels}
\alias{renameSeqlevels,GAlignments,character-method}
\alias{renameSeqlevels,GRangesList,character-method}
\alias{renameSeqlevels,GenomicRanges,character-method}

\title{
  (Deprecated) seqlevels utility functions
}

\description{
  Rename or subset the seqlevels in a \link{GenomicRanges}, \link{GRangesList} 
  or \link{GAlignments} object.
}

\usage{
  \S4method{keepSeqlevels}{GenomicRanges,character}(x, value, ...) 
  \S4method{renameSeqlevels}{GenomicRanges,character}(x, value, ...) 
}

\arguments{
  \item{x}{
    The \link{GenomicRanges}, \link{GRangesList} or \link{GAlignments} object 
    in which the seqlevels will be removed or renamed.
  }
  \item{value}{
    For \code{keepSeqlevels}, a \link{GRanges}, \link{GRangesList} or
    \code{GAlignments} or character vector. \code{x} is subset on the
    seqlevels in \code{value}.

    For \code{renameSeqlevels}, a named character vector where the names are the 
    \sQuote{old} and the values are the \sQuote{new} seqlevels. 
  }
  \item{...}{
    Arguments passed to other functions.
  }
}

\value{
  The \code{x} object with seqlevels removed or renamed. If \code{x} has
  no seqlevels (empty object) or no replacement values match the current
  seqlevels in \code{x} the unchanged \code{x} is returned.
}

\details{
  Many operations on \link{GRanges} objects require the seqlevels to match 
  before a comparison can be made (e.g., \code{findOverlaps(type="within")}).
  \code{keepSeqlevels} and \code{renameSeqlevels} are convenience functions 
  for subsetting and renaming the seqlevels of these objects. 

  \code{keepSeqlevels} subsets \code{x} based on the seqlevels provided in 
  \code{value}. If \code{value} does not match any seqlevels in \code{x} the 
  original \code{x} is returned. When \code{x} is a \link{GRangesList}, there 
  may be multiple chromosomes in a single list element. If not all chromosomes 
  are specified in \code{value}, a reduced list element is returned. All empty 
  list elements are dropped. See examples.

  \code{renameSeqlevels} renames the seqlevels in \code{x} to those provided 
  in \code{value}. \code{value} is a named character vector where the names are 
  the \sQuote{old} seqlevel names and the values are the \sQuote{new}. If 
  no names in \code{value} match the seqlevels in \code{x} the original \code{x} 
  is returned. 
}

\seealso{
  \code{\link{seqinfo}}
}

\examples{

## These functions are deprecated. Use seqlevels() instead. See the ?seqinfo 
## man page for examples of renaming and subsetting objects by seqlevels.

}

\author{Valerie Obenchain \email{vobencha@fhcrc.org}}

\keyword{methods}
\keyword{utilities}

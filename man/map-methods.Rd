\name{map-methods}

\alias{map,GenomicRanges,GRangesList-method}

\alias{coerce,RangesMapping,GenomicRanges-method}
\alias{granges,RangesMapping-method}

\title{Mapping ranges between sequences}
\description{
  A method for translating a set of input ranges through a
  \link{GRangesList} object. Returns a \link[IRanges]{RangesMapping}
  object.

  NOTE: The \code{\link[IRanges]{map}} generic function is defined and
  documented in the \pkg{IRanges} package.
  A method for translating a set of input ranges through a
  \link[GenomicAlignments]{GAlignments} object is defined and
  documented in the \pkg{GenomicAlignments} package.
}
\usage{
\S4method{map}{GenomicRanges,GRangesList}(from, to)
}
\arguments{
  \item{from}{The input ranges to map, usually
    a \code{\linkS4class{GenomicRanges}}}
  \item{to}{The alignment between the sequences in \code{from} and the
    sequences in the result.}
}
\value{
  An object of class \code{RangesMapping}. The \pkg{GenomicRanges} package
  provides some additional methods on this object:

  \item{}{\code{as(from, "GenomicRanges")}: Creates a
    \code{\linkS4class{GenomicRanges}} with \code{seqnames} and
    \code{ranges} from the \code{space} and \code{ranges} of
    \code{from}. The \code{hits} are coerced to a \code{DataFrame} and
    stored as the \code{values} of the result.
  }
  \item{}{\code{granges(x)}: Like the above, except returns just the range
    information as a \code{GRanges}, without the matching information.
  }
}
\details{
  Each element in \code{to} is taken to represent an alignment
  of a sequence on a genome. The typical case is a set of transcript
  models, as might be obtained via \code{GenomicFeatures::exonsBy}.
  The method translates the input ranges to be relative to the transcript
  start. This is useful, for example, when predicting coding consequences
  of changes to the genomic sequence.
}
\seealso{
  The \link[IRanges]{RangesMapping} class is the typical return value.
}

\author{M. Lawrence}


\name{map-methods}

\alias{map,GenomicRanges,GRangesList-method}
\alias{map,GenomicRanges,GappedAlignments-method}

\alias{coerce,RangesMapping,GenomicRanges-method}
\alias{granges,RangesMapping-method}

\title{Mapping ranges between sequences}
\description{
  The GenomicRanges package provides several methods for the
  \code{\link[IRanges]{map}} generic. They each translate a set of input
  ranges through a certain type of alignment and return a
  \code{\link[IRanges:RangesMapping-class]{RangesMapping}} object.
}
\usage{
\S4method{map}{GenomicRanges,GRangesList}(from, to)
\S4method{map}{GenomicRanges,GappedAlignments}(from, to)
}
\arguments{
  \item{from}{The input ranges to map, usually
    a \code{\linkS4class{GenomicRanges}}}
  \item{to}{The alignment between the sequences in \code{from} and the
    sequences in the result.}
}
\value{
  An object of class \code{RangesMapping}. The GenomicRanges package
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
  The methods currently depend on the type of \code{to}:
  \describe{
    \item{GRangesList}{Each element is taken to represent an alignment
      of a sequence on a genome. The typical case is a set of transcript
      models, as might be obtained via
      \code{GenomicFeatures::exonsBy}. The method translates the input
      ranges to be relative to the transcript start. This is useful, for
      example, when predicting coding consequences of changes to the
      genomic sequence.
    }
    \item{GappedAlignments}{Each element is taken to represent the
      alignment of a (read) sequence. The CIGAR string is used to
      translate the input ranges to be relative to the read start. This
      is useful, for example, when determining the cycle (read position)
      at which a particular genomic mismatch occurs.
    }
  }
}
\seealso{
  The \code{\link[IRanges:RangesMapping-class]{RangesMapping}} class is
  the typical return value.
}

\author{M. Lawrence}


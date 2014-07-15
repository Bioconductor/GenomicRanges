\name{map-methods}

\alias{map-methods}
\alias{map,GenomicRanges,GRangesList-method}

\title{Mapping ranges between sequences}

\description{
  \code{map} is defunct. Use \code{mapCoords} instead.

  A method for translating a set of input ranges through a
  \link{GRangesList} object. Returns a GRangesMapping object.

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
  An object of class GRangesMapping which contains both a
  \linkS4class{Hits} instance and a \linkS4class{GRanges} instance.

  \itemize{
    \item \code{hits}:
      Contains the matching (overlaps) from ranges in \code{from}
      to the ranges in \code{to}. Matching can be many-to-one or 
      one-to-many. See \code{?`{\link{Hits}}`} for details. 

    \item \code{GRanges}:
      Contains the mapped coordinates for ranges in \code{from} that 
      are `within` ranges in \code{to}. Length of the GRanges is the 
      same as the \code{Hits} object and contains one row for each overlap.
      Strand corresponds to the strand of the range in \code{to}.
  }
}

\details{
  Each element in \code{to} is taken to represent an alignment
  of a sequence on a genome. The typical case is a set of transcript
  models, as might be obtained via \code{GenomicFeatures::exonsBy}.
  The method translates the input ranges to be relative to the transcript
  start. This is useful, for example, when predicting coding consequences
  of changes to the genomic sequence. Note that output is produced only
  for ranges in \code{from} that fall completely `within` ranges in 
  \code{to}.
}

\seealso{
  The GRangesMapping class is the typical return value.
}

\examples{
## See ?mapCoords for examples
}

\author{M. Lawrence}

\name{map-methods}

\alias{map-methods}
\alias{map,GenomicRanges,GRangesList-method}

\title{Mapping ranges between sequences}

\description{
  A method for translating a set of input ranges through a
  \link{GRangesList} object. Returns a \linkS4class{GRangesMapping}
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
  An object of class \code{GRangesMapping} which contains both a
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
  The \linkS4class{GRangesMapping} class is the typical return value.
}

\examples{
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
gr <- GRanges("chr2L", IRanges(c(7500, 8400, 9000), width = 200,
              names = LETTERS[1:3]))

## Convert ranges in 'gr' to cds-centric coordinates.
cdsbytx <- cdsBy(txdb, "tx")
map1 <- map(gr, cdsbytx)
map1

## Mapped coordinates are returned for ranges completely 'within'
## coding regions. Range 'A' overlaps regions but is not 'within'.
findOverlaps(gr["A"], cdsbytx, type = "any")
findOverlaps(gr["A"], cdsbytx, type = "within")

## Convert to exon-centric coordinates.
exonsbytx <- exonsBy(txdb, "tx")
map2 <- map(gr, exonsbytx)

## Coercion of the GRangesMapping class to GenomicRanges produces 
## GRanges with the matching data (hits) as metadata columns.
map3 <- as(map2, "GenomicRanges")
map3

## Add the transcript names from 'exonsbytx'. 
mcols(map3)$tx_id <- names(exonsbytx)[map3$subjectHits]
map3
}

\author{M. Lawrence}

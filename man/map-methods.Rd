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
  An object of class \code{GRangesMapping}, which consists of a
  \linkS4class{Hits} instance and a \linkS4class{GRanges} instance.
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
gr <- GRanges("chr2L", IRanges(c(7500, 8400, 9000), width = 200))

## convert to cds-centric coordinates
cdsbytx <- cdsBy(txdb, "tx")
map1 <- map(gr, cdsbytx)
## Only ranges that fall completely 'within' a coding region are mapped.
hits(map1)
granges(map1)


## convert to exon-centric coordinates
exonsbytx <- exonsBy(txdb, "tx")
map2 <- map(gr, exonsbytx)
## A convienent form of the GRanges with matching data as 'mcols':
as(map2, "GenomicRanges")
}


\author{M. Lawrence}


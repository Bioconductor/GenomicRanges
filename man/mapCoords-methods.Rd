\name{mapCoords-methods}

\alias{mapCoords-methods}
\alias{mapCoords,GenomicRanges,GRangesList-method}

\title{Mapping ranges between sequences}

\description{
  A method for translating a set of input ranges through a
  \link{GRangesList} object. Returns a \linkS4class{GenomicRanges}
  object.

  NOTE: The \code{\link[IRanges]{mapCoords}} generic function is 
  defined and documented in the \pkg{IRanges} package. A method for 
  translating a set of input ranges through a
  \link[GenomicAlignments]{GAlignments} object is defined and
  documented in the \pkg{GenomicAlignments} package.
}

\usage{
\S4method{mapCoords}{GenomicRanges,GRangesList}(x, to, ..., 
          ignore.strand = FALSE, eltHits = FALSE)
}

\arguments{
  \item{x}{The input ranges to map, usually a \code{\linkS4class{GRanges}}.}
 
  \item{to}{The alignment between the sequences in \code{x} and the
    sequences in the result.}

  \item{ignore.strand}{Indicates if strand should be ignored in
    overlap operations.}

  \item{eltHits}{Indicates if element indices of the GRangesList
    hit in the overlap operation should be returned as a metadata 
    column. The `subjectHits` metadata column (always returned)
    is an index for the outer list elements; `eltHits` is and 
    index for the inner list elements.}

  \item{\dots}{Arguments passed to other methods.}
}

\value{
  A \code{GRanges} object of mapped coordinates with matching data
  as metadata columns (`queryHits` and `subjectHits`).

  Matching data are the result of calling \code{findOverlaps} with
  type `within` on  ranges in \code{x} (the query) and the 
  ranges in \code{to} (the subject). Matching can be many-to-one 
  or one-to-many; one row is reported for each match. 
}

\details{
  Each element in \code{to} is taken to represent an alignment
  of a sequence on a genome. The typical case is a set of transcript
  models, as might be obtained via \code{GenomicFeatures::exonsBy}.
  The method translates the input ranges to be relative to the transcript
  start. This is useful, for example, when predicting coding consequences
  of changes to the genomic sequence. Translated coordinates are only
  reported for ranges in \code{x} that fall completely `within` ranges 
  in \code{to}.
}

\seealso{
  The generic \link[IRanges]{mapCoords-methods} in the IRanges package.
}

\examples{
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
gr <- GRanges("chr2L", IRanges(c(7500, 8400, 9000), width = 200,
              names = LETTERS[1:3]))

## Convert ranges in 'gr' to cds-centric coordinates.
cdsbytx <- cdsBy(txdb, "tx")
map1 <- mapCoords(gr, cdsbytx)
map1

## Mapped coordinates are returned for ranges that fall completely 
## within ranges in 'to'. Range 'A' partially overlaps coding regions 
## from 3 different transcripts but is not 'within' any of them.
findOverlaps(gr["A"], cdsbytx, type = "any")
findOverlaps(gr["A"], cdsbytx, type = "within")

## Convert to exon-centric coordinates.
exonsbytx <- exonsBy(txdb, "tx")
map2 <- mapCoords(gr, exonsbytx)

## Use the mapping data to 
map3 <- as(map2, "GenomicRanges")
map3

## Add the transcript names from 'exonsbytx'. 
mcols(map3)$tx_id <- names(exonsbytx)[map3$subjectHits]
map3
}

\author{M. Lawrence}

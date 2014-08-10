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
          ignore.strand = FALSE, elt.loc = FALSE, elt.hits = FALSE)
}

\arguments{
  \item{x}{The input ranges to map, usually a \code{\linkS4class{GRanges}}.}
 
  \item{to}{The alignment between the sequences in \code{x} and the
    sequences in the result, usually a \code{\linkS4class{GRangesList}}.}

  \item{ignore.strand}{\code{logical}; When TRUE strand is ignored in 
    overlap operations.}

  \item{elt.loc}{\code{logical}; When TRUE the location in the
    individual list elements is returned as metadata column \code{eltLoc}. 
    For example, if \code{to} is a GRangesList of coding regions by transcript
    \code{eltLoc} is the position relative to each individual coding region.
    In contrast, the range in the output \code{GRanges} is the position 
    relative to the transcript, i.e., coding regions are concatenated.}

  \item{elt.hits}{\code{logical}; When TRUE the indices of the unlisted
    \code{to} hit by \code{x} in the overlap operation are returned
    as metadata column \code{eltHits}. Internal use only.}

  \item{\dots}{Arguments passed to other methods.}
}

\value{
  A \code{GRanges} object of mapped coordinates with matching data
  as metadata columns \code{queryHits} and \code{subjectHits}. 
  Matching data are the result of calling \code{findOverlaps} with
  type `within` on  ranges in \code{x} (the query) and the 
  ranges in \code{to} (the subject). Matching can be many-to-one 
  or one-to-many; one row is reported for each match.
 
  Optional metadata columns are \code{eltLoc} and \code{eltHits},
  see arguments section for details. The ranges in the output
  \code{GRanges} are position relative to the outer list element
  of \code{to}; all individual list elements are concatenated and
  counting starts at the 5' or 3' end depending on strand. In
  contrast, the \code{eltLoc} metadata column contains position
  relative to the individual list elements. 
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

\author{M. Lawrence and V. Obenchain \url{vobencha@fhcrc.org}}

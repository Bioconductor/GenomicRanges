\name{mapCoords-methods}

\alias{mapCoords}
\alias{mapCoords-methods}
\alias{mapCoords,GenomicRanges,GRangesList-method}
\alias{pmapCoords,GenomicRanges,GRangesList-method}


\title{Mapping ranges between sequences}

\description{
  A method for translating a set of input ranges through a
  \link{GRangesList} object. Returns a \linkS4class{GenomicRanges}
  object.

  The generics for \code{mapCoords} and \code{pmapCoords} are defined in the
  \pkg{IRanges} package. A method for translating a set of input ranges through
  a \link[GenomicAlignments]{GAlignments} object is defined and in the
  \pkg{GenomicAlignments} package.
}

\usage{
\S4method{mapCoords}{GenomicRanges,GRangesList}(x, to, ..., 
          ignore.strand = TRUE, elt.loc = FALSE, elt.hits = FALSE)

\S4method{pmapCoords}{GenomicRanges,GRangesList}(x, to, ..., 
          ignore.strand = TRUE, elt.loc = FALSE, elt.hits = FALSE)
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

  \item{\dots}{Arguments passed to other methods. Currently not used.}
}

\value{
  A \code{GRanges} object of mapped coordinates with matching data
  as metadata columns \code{queryHits} and \code{subjectHits}.

  Matching data are the result of calling \code{findOverlaps} with type `within`
  on  ranges in \code{x} (the query) and the ranges in \code{to} (the subject).
  In the case of \code{mapCoords} matching can be many-to-one or one-to-many;
  one row is reported for each match. For \code{pmapCoords} matching is
  one-to-one as the i-th element in \code{x} is only mapped to the i-th element
  in \code{to}.
 
  Optional metadata columns are \code{eltLoc} and \code{eltHits}, see arguments
  section for details. The ranges in the output \code{GRanges} are position
  relative to the outer list element of \code{to}; all individual list elements
  are concatenated and counting starts at the 5' or 3' end depending on strand.
  In contrast, the \code{eltLoc} metadata column contains position relative to
  the individual list elements. 
}

\details{
  Each element in \code{to} is taken to represent an alignment of a sequence on
  a genome. The typical case is a set of transcript models, as might be obtained
  via \code{GenomicFeatures::exonsBy}.  Each outer list element of the
  GRangesList represents a transcript while each each individual element is an
  exon in the transcript.

  \code{mapCoords} and \code{pmapCoords} translate the ranges in \code{x}
  relative to the transcript start (i.e., start of all ranges in \code{to}). The
  widths of the individual elements (exons in this example) are concatenated and
  counting starts at the 5' or 3' end depending on strand. Translated
  coordinates are only reported for ranges in \code{x} that fall completely
  `within` ranges in \code{to}.

  The transcript-centric coordinates are are useful, for example, when
  predicting coding consequences of changes to the genomic sequence.
  Additionally, the position of \code{x} relative to the individual element
  (i.e., exon) is returned as a metadata column when \code{elt.loc = TRUE}.

  \code{mapCoords} maps the i-th element in \code{x} to each element in
  \code{to} returning in a many-to-many mapping. In contrast, \code{pmapCoords}
  treats the two inputs as parallel vectors and maps the i-th element of
  \code{x} to the i-th element of \code{to} returning a maximum of one result
  per input element.
}

\seealso{
  \itemize{
    \item{The generic \link[IRanges]{mapCoords-methods} in the IRanges package.}
    \item{Additional methods in the GenomicAlignments package 
          \link[GenomicAlignments]{mapCoords-methods}.}
  }
}

\examples{
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
gr <- GRanges("chr2L", IRanges(c(7500, 8400, 9000), width = 200,
              names = LETTERS[1:3]))

## Coding region coordinates relative to transcript:
cdsbytx <- cdsBy(txdb, "tx")
map1 <- mapCoords(gr, cdsbytx)
map1

## Range A is not translated because it does not fall completely 
## 'within' any ranges.
findOverlaps(gr["A"], cdsbytx, type = "any")
findOverlaps(gr["A"], cdsbytx, type = "within")

## Exon region coordinates relative to trancript:
exonsbytx <- exonsBy(txdb, "tx")
map2 <- mapCoords(gr, exonsbytx)

## Add exon coordinates relative to the individual exons:
map3 <- mapCoords(gr, exonsbytx, elt.hits=TRUE, elt.loc=TRUE) 
map3

## Extract the individual exons hit:
unlist(exonsbytx, use.names=FALSE)[mcols(map3)$eltHits]
}

\author{M. Lawrence and V. Obenchain \url{vobencha@fhcrc.org}}

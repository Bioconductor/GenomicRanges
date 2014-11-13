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
\S4method{mapCoords}{GenomicRanges,GRangesList}(from, to, ..., 
          ignore.strand = TRUE, elt.hits = FALSE)

\S4method{pmapCoords}{GenomicRanges,GRangesList}(from, to, ..., 
          ignore.strand = TRUE, elt.hits = FALSE)
}

\arguments{
  \item{from}{The input ranges to map, usually a \code{\linkS4class{GRanges}}.}
 
  \item{to}{The alignment between the sequences in \code{from} and the
    sequences in the result, usually a \code{\linkS4class{GRangesList}}.}

  \item{ignore.strand}{\code{logical}; When TRUE strand is ignored in 
    overlap operations.}

  \item{elt.hits}{\code{logical}; When TRUE, the output includes a metadata 
    column, eltHits, with indices of the inner list elements of \code{to}
    hit by \code{from}. Useful for identifying elements of \code{to} hit
    by \code{from}. See examples.}

  \item{\dots}{Arguments passed to other methods. Currently not used.}
}

\value{
  A \code{GRanges} object of mapped coordinates with matching data as metadata
  columns \code{queryHits} and \code{subjectHits}.  The ranges in the output
  \code{GRanges} are position relative to the outer list element of \code{to};
  all individual list elements are concatenated and counting starts at the 5' or
  3' end depending on strand.

  Matching data are the result of calling \code{findOverlaps} with type `within`
  on  ranges in \code{from} (the query) and the ranges in \code{to} (the
  subject).  In the case of \code{mapCoords} matching can be many-to-one or
  one-to-many; one row is reported for each match. For \code{pmapCoords}
  matching is one-to-one as the i-th element in \code{from} is only mapped to
  the i-th element in \code{to}.
 
  When \code{elt.hits} is TRUE, the \code{eltHits} metadata column includes 
  the index of inner list elements in \code{to} hit by \code{from}. In
  some applications it may be useful to identify the exact list element that
  was overlapped. These elements can be extracted with the combination of
  \code{subjectHits} (outer list index) and \code{eltHits} (inner list index). 
}

\details{
  Each element in \code{to} is taken to represent an alignment of a sequence on
  a genome. The typical case is a set of transcript models, as might be obtained
  via \code{GenomicFeatures::exonsBy}.  Each outer list element of the
  GRangesList represents a transcript while each each individual element is an
  exon in the transcript.

  \code{mapCoords} and \code{pmapCoords} translate the ranges in \code{from}
  relative to the transcript start (i.e., start of all ranges in \code{to}). The
  widths of the individual elements (exons in this example) are concatenated and
  counting starts at the 5' or 3' end depending on strand. Translated
  coordinates are only reported for ranges in \code{from} that fall completely
  `within` ranges in \code{to}.

  The transcript-centric coordinates are are useful, for example, when
  predicting coding consequences of changes to the genomic sequence.

  \code{mapCoords} maps the i-th element in \code{from} to each element in
  \code{to} returning in a many-to-many mapping. In contrast, \code{pmapCoords}
  treats the two inputs as parallel vectors and maps the i-th element of
  \code{from} to the i-th element of \code{to} returning a maximum of one result
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
map1 <- mapCoords(gr, cdsbytx, ignore.strand = FALSE)
map1

## Range A is not translated because it does not fall completely 
## 'within' any ranges.
findOverlaps(gr["A"], cdsbytx, type = "any")
findOverlaps(gr["A"], cdsbytx, type = "within")

## Exon region coordinates relative to transcript:
exonsbytx <- exonsBy(txdb, "tx")
map2 <- mapCoords(gr, exonsbytx, ignore.strand = FALSE)

## Extracting elements from 'to':
##
## Some applications may require a closer look at the elements
## in 'to' used in the mapping. The subjectHits metadata column contains
## indices of the outer list elements hit. When elt.hits = TRUE the output 
## also includes an eltHits column which holds indices of the inner list
## elements. Together these can be used to extract the elements of 'to',
## which in this example are the exons.
map3 <- mapCoords(gr, exonsbytx, ignore.strand = FALSE, elt.hits=TRUE) 
map3

## The first 2 rows of 'map3' show the second range in 'from' hit list 
## elements 1 and 3.
map3[1:2, ]

## Extract list elements of 'to':
outerlist <- exonsbytx[mcols(map3[1:2, ])$subjectHits]

## Extract the exons from the list elements:
exons <- Map("[", outerlist, mcols(map3[1:2, ])$eltHits)
}

\author{M. Lawrence and V. Obenchain \url{vobencha@fhcrc.org}}

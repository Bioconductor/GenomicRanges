\name{GPos-class}
\docType{class}

% Classes

\alias{class:GPos}
\alias{GPos-class}
\alias{GPos}

\alias{class:UnstitchedGPos}
\alias{UnstitchedGPos-class}
\alias{UnstitchedGPos}

\alias{class:StitchedGPos}
\alias{StitchedGPos-class}
\alias{StitchedGPos}

% Accessors

\alias{pos,GPos-method}

% Coercion

\alias{coerce,ANY,UnstitchedGPos-method}
\alias{coerce,ANY,StitchedGPos-method}
\alias{coerce,ANY,GPos-method}
\alias{coerce,GRanges,UnstitchedGPos-method}
\alias{coerce,GRanges,StitchedGPos-method}
\alias{coerce,GRanges,GPos-method}
\alias{from_GPos_to_GRanges}
\alias{coerce,UnstitchedGPos,GRanges-method}
\alias{coerce,StitchedGPos,GRanges-method}
\alias{as.data.frame.GPos}
\alias{as.data.frame,GPos-method}

\alias{updateObject,GPos-method}

% Display

\alias{summary.GPos}
\alias{summary,GPos-method}
\alias{show,GPos-method}

\title{Memory-efficient representation of genomic positions}

\description{
  The GPos class is a container for storing a set of \emph{genomic positions}
  (a.k.a. \emph{genomic loci}).
  It exists in 2 flavors: UnstitchedGPos and StitchedGPos. Each flavor uses
  a particular internal representation:
  \itemize{
    \item In an UnstitchedGPos instance the positions are stored as an
          integer vector.
    \item In a StitchedGPos instance the positions are stored as an
          \link{IRanges} object where each range represents a run of
          \emph{consecutive positions} (i.e. a run of positions that
          are adjacent and in \emph{ascending order}). This storage is
          particularly memory-efficient when the vector of positions
          contains long runs of consecutive positions.
  }
  Because genomic positions can be seen as genomic ranges of width 1, the
  GPos class extends the \link{GenomicRanges} virtual class (via the
  \link{GRanges} class).
}

\usage{
## Constructor function
GPos(seqnames=NULL, pos=NULL, strand=NULL,
     ..., seqinfo=NULL, seqlengths=NULL, stitch=NA)
}

\arguments{
  \item{seqnames, strand, ..., seqinfo, seqlengths}{
    See documentation of the \code{\link{GRanges}()} constructor function
    for a description of these arguments.
  }
  \item{pos}{
    \code{NULL}, or an integer or numeric vector, or an \link[IRanges]{IRanges}
    or \link[IRanges]{IPos} object (or other \link[IRanges]{IntegerRanges}
    derivative).
    If not \code{NULL}, \code{GPos()} will try to turn it into an
    \link[IRanges]{IPos} derivative with \code{IPos(pos, stitch=stitch)}.

    When \code{pos} is an \link{IRanges} object (or other \link{IntegerRanges}
    derivative), each range in it is interpreted as a run of
    consecutive positions.
  }
  \item{stitch}{
    \code{TRUE}, \code{FALSE}, or \code{NA} (the default).

    Controls which internal representation should be used: StitchedGPos
    (when \code{stitch} is \code{TRUE}) or UnstitchedGPos (when \code{stitch}
    is \code{FALSE}).

    When \code{stitch} is \code{NA} (the default), which internal
    representation will be used depends on the flavour of the
    \link[IRanges]{IPos} derivative returned by \code{IPos(pos)}:
    UnstitchedGPos if \code{IPos(pos)} returns an
    \link[IRanges]{UnstitchedIPos} instance, and StitchedGPos if it
    returns a \link[IRanges]{StitchedIPos} instance.
  }
}

\details{
  Even though a \link{GRanges} object can be used for storing genomic
  positions, using a GPos object is more efficient. In particular the
  memory footprint of an UnstitchedGPos object is typically about half
  that of a \link{GRanges} object.

  OTOH the memory footprint of a StitchedGPos object can vary a lot but
  will never be worse than that of a \link{GRanges} object. However it
  will reduce dramatically if the vector of positions contains long runs
  of consecutive positions. In the worst case scenario (i.e. when the
  object contains no consecutive positions) its memory footprint will
  be the same as that of a \link{GRanges} object.

  Like for any \link[S4Vectors]{Vector} derivative, the length of a
  GPos object cannot exceed \code{.Machine$integer.max} (i.e. 2^31 on
  most platforms). \code{GPos()} will return an error if \code{pos}
  contains too many positions.
}

\value{
  An UnstitchedGPos or StitchedGPos object.
}

\section{Accessors}{

  \subsection{Getters}{
    GPos objects support the same set of getters as other \link{GenomicRanges}
    derivatives (i.e. \code{seqnames()}, \code{ranges()}, \code{start()},
    \code{end()}, \code{strand()}, \code{mcols()}, \code{seqinfo()},
    etc...), plus the \code{pos()} getter which is equivalent to
    \code{start()} or \code{end()}. See \code{?\link{GenomicRanges}} for the
    list of getters supported by \link{GenomicRanges} derivatives.

    Note that \code{ranges()} returns an \link[IRanges]{IPos} derivative
    instead of the \link[IRanges]{IRanges} object that one gets with other
    \link{GenomicRanges} derivatives. To get an \link[IRanges]{IRanges}
    object, you need to call \code{ranges()} again on this
    \link[IRanges]{IPos} derivative i.e. do \code{ranges(ranges(x))}
    on GPos object \code{x}.
  }

  \subsection{Setters}{
    Like \link{GRanges} objects, GPos derivatives support the following
    setters:
    \itemize{
      \item The \code{seqnames()} and \code{strand()} setters.
      \item The \code{names()}, \code{mcols()} and \code{metadata()} setters.
      \item The family of setters that operate on the \code{seqinfo}
            component of an object:
            \code{\link[GenomeInfoDb]{seqlevels}()},
            \code{\link[GenomeInfoDb]{seqlevelsStyle}()},
            \code{\link[GenomeInfoDb]{seqlengths}()},
            \code{\link[GenomeInfoDb]{isCircular}()},
            \code{\link[GenomeInfoDb]{genome}()},
            and \code{\link[GenomeInfoDb]{seqinfo}()}.
            These setters are defined and documented in the \pkg{GenomeInfoDb}
            package.
    }

    However, there is no \code{pos()} setter for GPos derivatives at the
    moment (although one might be added in the future).
  }
}

\section{Coercion}{
  From UnstitchedGPos to StitchedGPos and vice-versa: coercion back and
  forth between UnstitchedGPos and StitchedGPos is supported via
  \code{as(x, "StitchedGPos")} and \code{as(x, "UnstitchedGPos")}.
  This is the most efficient and recommended way to switch between the
  2 internal representations. Note that this switch can have dramatic
  consequences on memory usage so is for advanced users only.
  End users should almost never need to do this switch when following
  a typical workflow.

  From \link{GenomicRanges} to UnstitchedGPos, StitchedGPos, or GPos:
  A \link{GenomicRanges} derivative \code{x} in which all the ranges have
  a width of 1 can be coerced to an UnstitchedGPos or StitchedGPos object
  with \code{as(x, "UnstitchedGPos")} or \code{as(x, "StitchedGPos")},
  respectively.
  For convenience \code{as(x, "GPos")} is supported and is equivalent to
  \code{as(x, "UnstitchedGPos")}.

  From GPos to \link{GRanges}:
  A GPos derivative \code{x} can be coerced to a \link{GRanges} object
  with \code{as(x, "GRanges")}. However be aware that the resulting object
  can use thousands times (or more) memory than \code{x}!
  See "MEMORY USAGE" in the Examples section below.

  From GPos to ordinary R objects:
  Like with any other \link{GenomicRanges} derivative, \code{as.character()},
  \code{as.factor()}, and \code{as.data.frame()} work on a GPos derivative
  \code{x}. Note however that \code{as.data.frame(x)} returns a data frame
  with a \code{pos} column (containing \code{pos(x)}) instead of the
  \code{start}, \code{end}, and \code{width} columns that one gets with other
  \link{GenomicRanges} derivatives.
}

\section{Subsetting}{
  A GPos derivative can be subsetted exactly like a \link{GRanges} object.
}

\section{Concatenation}{
  GPos derivatives can be concatenated with \code{c()} or \code{append()}.
  See \code{?\link[S4Vectors]{c}} in the \pkg{S4Vectors} package for
  more information about concatenating Vector derivatives.
}

\section{Splitting and Relisting}{
  Like with any other \link{GRanges} object, \code{split()} and \code{relist()}
  work on a GPos derivative.
}

\note{
  Internal representation of GPos objects has changed in \pkg{GenomicRanges}
  1.29.10 (Bioc 3.6). Update any old object \code{x} with:
  \code{x <- updateObject(x, verbose=TRUE)}.
}

\author{
  Hervé Pagès; based on ideas borrowed from Georg Stricker
  \email{georg.stricker@in.tum.de} and Julien Gagneur
  \email{gagneur@in.tum.de}
}

\seealso{
  \itemize{
    \item The \link[IRanges]{IPos} class in the \pkg{IRanges}
          package for storing a set of \emph{integer positions} (i.e.
          integer ranges of width 1).

    \item The \link{GRanges} class for storing a set of \emph{genomic ranges}
          of arbitrary width.

    \item \link[GenomeInfoDb]{Seqinfo} objects and the
          \code{\link[GenomeInfoDb]{seqinfo}} accessor and family in the
          \pkg{GenomeInfoDb} package for accessing/modifying information
          about the underlying sequences of a \link{GenomicRanges} derivative.

    \item \link{GenomicRanges-comparison} for comparing and ordering genomic
          ranges and/or positions.

    \item \link[GenomicRanges]{findOverlaps-methods} for finding overlapping
          genomic ranges and/or positions.

    \item \link[GenomicRanges]{intra-range-methods} and
          \link[GenomicRanges]{inter-range-methods}
          for \emph{intra range} and \emph{inter range} transformations
          of \link{GenomicRanges} derivatives.

    \item \link[GenomicRanges]{coverage-methods} for computing the coverage
          of a set of genomic ranges and/or positions.

    \item \link[GenomicRanges]{nearest-methods} for finding the nearest
          genomic range/position neighbor.

    \item The \code{\link[BSgenome]{snpsBySeqname}},
          \code{\link[BSgenome]{snpsByOverlaps}}, and
          \code{\link[BSgenome]{snpsById}} methods for
          \link[BSgenome]{SNPlocs} objects defined in the \pkg{BSgenome}
          package for extractors that return a GPos derivative.

    \item \link[SummarizedExperiment]{SummarizedExperiment} objects and
          derivatives in the \pkg{SummarizedExperiment} package.
  }
}

\examples{
showClass("GPos")  # shows the known subclasses

## ---------------------------------------------------------------------
## BASIC EXAMPLES
## ---------------------------------------------------------------------

## Example 1:
gpos1a <- GPos(seqnames=Rle(c("chr1", "chr2", "chr1"), c(10, 6, 4)),
               pos=c(44:53, 5:10, 2:5))
gpos1a  # unstitched

length(gpos1a)
seqnames(gpos1a)
pos(gpos1a)  # same as 'start(gpos1a)' and 'end(gpos1a)'
strand(gpos1a)
as.character(gpos1a)
as.data.frame(gpos1a)
as(gpos1a, "GRanges")
as.data.frame(as(gpos1a, "GRanges"))
gpos1a[9:17]

gpos1b <- GPos(seqnames=Rle(c("chr1", "chr2", "chr1"), c(10, 6, 4)),
               pos=c(44:53, 5:10, 2:5), stitch=TRUE)
gpos1b  # stitched

## 'gpos1a' and 'gpos1b' are semantically equivalent, only their
## internal representations differ:
all(gpos1a == gpos1b)

gpos1c <- GPos(c("chr1:44-53", "chr2:5-10", "chr1:2-5"))
gpos1c  # stitched

identical(gpos1b, gpos1c)

## Example 2:
pos_runs <- GRanges("chrI", IRanges(c(1, 6, 12, 17), c(5, 10, 16, 20)),
                    strand=c("*", "-", "-", "+"))
gpos2 <- GPos(pos_runs)
gpos2  # stitched
strand(gpos2)

## Example 3:
gpos3A <- gpos3B <- GPos(c("chrI:1-1000", "chrI:1005-2000"))
npos <- length(gpos3A)

mcols(gpos3A)$sample <- Rle("sA")
sA_counts <- sample(10, npos, replace=TRUE)
mcols(gpos3A)$counts <- sA_counts

mcols(gpos3B)$sample <- Rle("sB")
sB_counts <- sample(10, npos, replace=TRUE)
mcols(gpos3B)$counts <- sB_counts

gpos3 <- c(gpos3A, gpos3B)
gpos3

## Example 4:
library(BSgenome.Scerevisiae.UCSC.sacCer2)
genome <- BSgenome.Scerevisiae.UCSC.sacCer2
gpos4 <- GPos(seqinfo(genome))
gpos4  # all the positions along the genome are represented
mcols(gpos4)$dna <- do.call("c", unname(as.list(genome)))
gpos4

## Note however that, like for any Vector derivative, the length of a
## GPos derivative cannot exceed '.Machine$integer.max' (i.e. 2^31 on
## most platforms) so the above only works with a "small" genome.
## For example it doesn't work with the Human genome:
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
\dontrun{
  GPos(seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene))  # error!
}

## You can use isSmallGenome() to check upfront whether the genome is
## "small" or not.
isSmallGenome(genome)  # TRUE
isSmallGenome(TxDb.Hsapiens.UCSC.hg19.knownGene)  # FALSE

## ---------------------------------------------------------------------
## MEMORY USAGE
## ---------------------------------------------------------------------

## Coercion to GRanges works...
gr4 <- as(gpos4, "GRanges")
gr4
## ... but is generally not a good idea:
object.size(gpos4)
object.size(gr4)     # 8 times bigger than the StitchedGPos object!

## Shuffling the order of the positions impacts memory usage:
gpos4r <- rev(gpos4)
object.size(gpos4r)  # significantly
gpos4s <- sample(gpos4)
object.size(gpos4s)  # even worse!

## If one anticipates a lot of shuffling of the genomic positions,
## then an UnstitchedGPos object should be used instead:
gpos4b <- as(gpos4, "UnstitchedGPos")
object.size(gpos4b)  # initial size is bigger than stitched version
object.size(rev(gpos4b))  # size didn't change
object.size(sample(gpos4b))  # size increased, but is still < stitched
                             # version

## AN IMPORTANT NOTE: In the worst situations, GPos still performs as
## good as a GRanges object.
object.size(as(gpos4r, "GRanges"))  # same size as 'gpos4r'
object.size(as(gpos4s, "GRanges"))  # same size as 'gpos4s'

## Best case scenario is when the object is strictly sorted (i.e.
## positions are in strict ascending order).
## This can be checked with:
is.unsorted(gpos4, strict=TRUE)  # 'gpos4' is strictly sorted

## ---------------------------------------------------------------------
## USING MEMORY-EFFICIENT METADATA COLUMNS
## ---------------------------------------------------------------------
## In order to keep memory usage as low as possible, it is recommended
## to use a memory-efficient representation of the metadata columns that
## we want to set on the object. Rle's are particularly well suited for
## this, especially if the metadata columns contain long runs of
## identical values. This is the case for example if we want to use a
## GPos object to represent the coverage of sequencing reads along a
## genome.

## Example 5:
library(pasillaBamSubset)
library(Rsamtools)  # for the BamFile() constructor function
bamfile1 <- BamFile(untreated1_chr4())
bamfile2 <- BamFile(untreated3_chr4())
gpos5 <- GPos(seqinfo(bamfile1))
library(GenomicAlignments)  # for "coverage" method for BamFile objects
cvg1 <- unlist(coverage(bamfile1), use.names=FALSE)
cvg2 <- unlist(coverage(bamfile2), use.names=FALSE)
mcols(gpos5) <- DataFrame(cvg1, cvg2)
gpos5

object.size(gpos5)  # lightweight

## Keep only the positions where coverage is at least 10 in one of the
## 2 samples:
gpos5[mcols(gpos5)$cvg1 >= 10 | mcols(gpos5)$cvg2 >= 10]

## ---------------------------------------------------------------------
## USING A GPos OBJECT IN A RangedSummarizedExperiment OBJECT
## ---------------------------------------------------------------------
## Because the GPos class extends the GenomicRanges virtual class, a
## GPos derivative can be used as the rowRanges component of a
## RangedSummarizedExperiment object.

## As a 1st example, we show how the counts for samples sA and sB in
## 'gpos3' can be stored in a SummarizedExperiment object where the rows
## correspond to unique genomic positions and the columns to samples:
library(SummarizedExperiment)
counts <- cbind(sA=sA_counts, sB=sB_counts)
mcols(gpos3A) <- NULL
rse3 <- SummarizedExperiment(list(counts=counts), rowRanges=gpos3A)
rse3
rowRanges(rse3)
head(assay(rse3))

## Finally we show how the coverage data from Example 5 can be easily
## stored in a lightweight SummarizedExperiment derivative:
cvg <- mcols(gpos5)
mcols(gpos5) <- NULL
rse5 <- SummarizedExperiment(list(cvg=cvg), rowRanges=gpos5)
rse5
rowRanges(rse5)
assay(rse5)

## Keep only the positions where coverage is at least 10 in one of the
## 2 samples:
rse5[assay(rse5)$cvg1 >= 10 | assay(rse5)$cvg2 >= 10]
}
\keyword{methods}
\keyword{classes}

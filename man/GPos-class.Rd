\name{GPos-class}
\docType{class}

\alias{class:GPos}
\alias{GPos-class}
\alias{GPos}

\alias{length,GPos-method}
\alias{names,GPos-method}
\alias{names<-,GPos-method}
\alias{seqnames,GPos-method}
\alias{pos}
\alias{pos,GPos-method}
\alias{start,GPos-method}
\alias{end,GPos-method}
\alias{width,GPos-method}
\alias{ranges,GPos-method}
\alias{strand,GPos-method}
\alias{seqinfo,GPos-method}
\alias{extractROWS,GPos-method}
\alias{show,GPos-method}
\alias{c,GPos-method}

\title{GPos objects}

\description{
  The GPos class is a container for storing a set of genomic \emph{positions},
  that is, genomic ranges of width 1. Even though a \link{GRanges} object
  could be used for that, using a GPos object is much more memory-efficient,
  especially when the object contains long runs of adjacent positions.
}

\usage{
## Constructor
GPos(pos_runs)
}

\arguments{
  \item{pos_runs}{
    A \link{GRanges} object (or any other \link{GenomicRanges} derivative)
    where each range is interpreted as a run of adjacent genomic positions.
    If \code{pos_runs} is not a \link{GenomicRanges} object,
    \code{GPos()} first tries to coerce it to one with
    \code{as(pos_runs, "GenomicRanges", strict=FALSE)}.
  }
}

\value{
  A GPos object.
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item \link{GRanges} objects.

    \item \link[SummarizedExperiment]{SummarizedExperiment} objects in the
          \pkg{SummarizedExperiment} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## BASIC EXAMPLES
## ---------------------------------------------------------------------

## Example 1:
gpos1 <- GPos(c("chr1:44-53", "chr1:5-10", "chr2:2-5"))
gpos1
pos(gpos1)  # same as 'start(gpos1)' and 'end(gpos1)'
gpos1[9:17]

## Example 2:
pos_runs <- GRanges("chrI", IRanges(c(1, 6, 12, 17), c(5, 10, 16, 20)),
                    strand=c("+", "-", "-", "+"))
gpos2 <- GPos(pos_runs)
gpos2

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
## GPos object cannot exceed '.Machine$integer.max' (i.e. 2^31 on most
## platforms) so the above only works with a "small" genome.
## For example it doesn't work with the Human genome:
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
\dontrun{
  GPos(seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene))  # error!
}

## You can use isSmallGenome() to check upfront whether the genome is
## "small" or not.
isSmallGenome(genome)
isSmallGenome(TxDb.Hsapiens.UCSC.hg19.knownGene)

## ---------------------------------------------------------------------
## MEMORY USAGE
## ---------------------------------------------------------------------

## Coercion to GRanges works...
gr4 <- as(gpos4, "GRanges")
gr4
## ... but is generally not a good idea:
object.size(gpos4)
object.size(gr4)  # 6951 times bigger than the GPos object!

## Shuffling the order of the positions impacts memory usage:
object.size(rev(gpos4))     # significantly
object.size(sample(gpos4))  # even worse!

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
cov1 <- unlist(coverage(bamfile1), use.names=FALSE)
cov2 <- unlist(coverage(bamfile2), use.names=FALSE)
mcols(gpos5) <- DataFrame(cov1, cov2)
gpos5

## Keep only the positions where coverage is at least 10 in one of the
## 2 samples:
gpos5[mcols(gpos5)$cov1 >= 10 | mcols(gpos5)$cov2 >= 10]

## ---------------------------------------------------------------------
## USING A GPos OBJECT IN A SummarizedExperiment OBJECT
## ---------------------------------------------------------------------
## Because the GPos class extends the GenomicRanges virtual class, a GPos
## object can be used as the rowRanges component of a SummarizedExperiment
## object.

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

## Finally we show how the coverage data from Example 5 can be
## efficiently stored in a SummarizedExperiment:
cov <- mcols(gpos5)
mcols(gpos5) <- NULL
rse5 <- SummarizedExperiment(list(cov=cov), rowRanges=gpos5)
rse5
rowRanges(rse5)
assay(rse5)

## Keep only the positions where coverage is at least 10 in one of the
## 2 samples:
rse5[assay(rse5)$cov1 >= 10 | assay(rse5)$cov2 >= 10]
}
\keyword{methods}
\keyword{classes}

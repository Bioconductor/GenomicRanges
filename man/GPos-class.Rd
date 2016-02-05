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
  that is, genomic ranges of width 1. It's more memory-efficient than
  the \link{GRanges} container when the object contains long runs of adjacent
  positions.
}

\usage{
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

\details{
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item \link{GRanges} objects.
  }
}

\examples{
## ---------------------------------------------------------------------
## BASIC EXAMPLES
## ---------------------------------------------------------------------

## Example 1
gpos1 <- GPos(c("chr1:44-53", "chr1:5-10", "chr2:2-5"))
gpos1
pos(gpos1)  # same as 'start(gpos1)' and 'end(gpos1)'
gpos1[9:17]

## Example 2
pos_runs <- GRanges("chrI", IRanges(c(1, 6, 12, 17), c(5, 10, 16, 20)),
                    strand=c("+", "-", "-", "+"))
gpos2 <- GPos(pos_runs)
gpos2

## Example 3
gpos3A <- gpos3B <- GPos(c("chrI:1-1000", "chrI:1005-2000"))
npos <- length(gpos3A)

mcols(gpos3A)$sample <- Rle("sA")
sA_counts <- sample(1:10, npos, replace=TRUE)
mcols(gpos3A)$counts <- sA_counts

mcols(gpos3B)$sample <- Rle("sB")
sB_counts <- sample(1:10, npos, replace=TRUE)
mcols(gpos3B)$counts <- sB_counts

gpos3 <- c(gpos3A, gpos3B)
gpos3

## Example 4
library(BSgenome.Scerevisiae.UCSC.sacCer2)
gpos4 <- GPos(seqinfo(BSgenome.Scerevisiae.UCSC.sacCer2))
gpos4  # all the positions along the genome are represented

## ---------------------------------------------------------------------
## MEMORY USAGE
## ---------------------------------------------------------------------

## Coercion to GRanges works...
gr4 <- as(gpos4, "GRanges")
gr4
## ... but is generally not a good idea:
object.size(gpos4)
object.size(gr4)  # 6951 times bigger than the GPos object!

## In order to keep memory usage as low as possible, it is highly
## recommended to use Rle objects only for the metadata columns:
npos <- length(gpos4)
counts <- rep(Rle(sample(150L, replace=TRUE, npos / 10)), each=10)
counts <- c(counts, Rle(0L, npos - length(counts)))
mcols(gpos4)$counts <- counts
gpos4

## ---------------------------------------------------------------------
## USING A GPos OBJECT IN A SummarizedExperiment OBJECT
## ---------------------------------------------------------------------
## A GPos object can be used as the rowRanges component of a
## SummarizedExperiment object.

## As a 1st example, we show how the counts for samples sA and sB in
## 'gpos3' can be stored in a SummarizedExperiment object where the rows
## correspond to unique genomic positions and the columns to samples:
library(SummarizedExperiment)
counts <- cbind(sA=sA_counts, sB=sB_counts)
mcols(gpos3A) <- NULL
rse <- SummarizedExperiment(list(counts=counts), rowRanges=gpos3A)
rse
rowRanges(rse)
head(assay(rse))

## If the number of genomic positions to store in the SummarizedExperiment
## object is big (e.g. millions), then it's necessary to use a
## memory-efficient representation of the assay data like sparse matrices,
## Rle matrices (not implemented yet), on-disk data, etc...
## Here we use the Matrix package to create 2 sparse matrices of
## dimensions "length of sacCer2 genome" x "nb of samples":
library(Matrix)
var1 <- Matrix(0, nrow=length(gpos4), ncol=3)  # 3 samples
colnames(var1) <- c("s1", "s2", "s3")
var2 <- var1
var1[sample(prod(dim(var1)), 5000000)] <-
                                sample(15L, replace=TRUE, 5000000)
var2[sample(prod(dim(var2)), 5000000)] <-
                                sample(15L, replace=TRUE, 5000000)

## Store them in a SummarizedExperiment object:
rse <- SummarizedExperiment(list(var1=var1, var2=var2), rowRanges=gpos4)
rse
rowRanges(rse)
head(assay(rse, 1))
head(assay(rse, 2))
}
\keyword{methods}
\keyword{classes}

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
  that is, genomic ranges of length 1. It's more memory-efficient than
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
## 1st example
pos_runs <- GRanges(c("chr1", "chr1", "chr2"),
                    IRanges(c(44, 5, 2), width=c(10, 6, 4)))
gpos1 <- GPos(pos_runs)
gpos1
pos(gpos1)  # same as 'start(gpos1)' and 'end(gpos1)'
gpos1[9:17]

## 2nd example
pos_runs <- GRanges("chrI", IRanges(c(1, 6, 12, 17), c(5, 10, 16, 20)),
                    strand=c("+", "-", "-", "+"))
gpos2 <- GPos(pos_runs)
gpos2

## 3rd example
pos_runs <- GRanges("chrI", IRanges(c(1, 1005), c(1000, 2000)))
gpos3A <- gpos3B <- GPos(pos_runs)
mcols(gpos3A)$sample <- Rle("sA")
mcols(gpos3B)$sample <- Rle("sB")
gpos3 <- c(gpos3A, gpos3B)

## 4th example
library(BSgenome.Scerevisiae.UCSC.sacCer2)
gpos4 <- GPos(seqinfo(BSgenome.Scerevisiae.UCSC.sacCer2))
gpos4  # all the positions along the genome are represented

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
}
\keyword{methods}
\keyword{classes}

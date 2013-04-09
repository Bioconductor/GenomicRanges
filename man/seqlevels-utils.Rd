\name{utils}
\docType{methods}
\alias{keepSeqlevels}
\alias{keepSeqlevels,GenomicRanges,GenomicRanges-method}
\alias{keepSeqlevels,GenomicRanges,GRangesList-method}
\alias{keepSeqlevels,GenomicRanges,GAlignments-method}
\alias{keepSeqlevels,GenomicRanges,character-method}
\alias{keepSeqlevels,GRangesList,GenomicRanges-method}
\alias{keepSeqlevels,GRangesList,GRangesList-method}
\alias{keepSeqlevels,GRangesList,GAlignments-method}
\alias{keepSeqlevels,GRangesList,character-method}
\alias{keepSeqlevels,GAlignments,GenomicRanges-method}
\alias{keepSeqlevels,GAlignments,GRangesList-method}
\alias{keepSeqlevels,GAlignments,GAlignments-method}
\alias{keepSeqlevels,GAlignments,character-method}

\alias{renameSeqlevels}
\alias{renameSeqlevels,GAlignments,character-method}
\alias{renameSeqlevels,GRangesList,character-method}
\alias{renameSeqlevels,GenomicRanges,character-method}

\title{
  seqlevels utility functions
}

\description{
  Rename or subset the seqlevels in a \link{GenomicRanges}, \link{GRangesList} 
  or \link{GAlignments} object.
}

\usage{
  \S4method{keepSeqlevels}{GenomicRanges,character}(x, value, ...) 
  \S4method{renameSeqlevels}{GenomicRanges,character}(x, value, ...) 
}

\arguments{
  \item{x}{
    The \link{GenomicRanges}, \link{GRangesList} or \link{GAlignments} object 
    in which the seqlevels will be removed or renamed.
  }
  \item{value}{
    For \code{keepSeqlevels}, a \link{GRanges}, \link{GRangesList} or
    \code{GAlignments} or character vector. \code{x} is subset on the
    seqlevels in \code{value}.

    For \code{renameSeqlevels}, a named character vector where the names are the 
    \sQuote{old} and the values are the \sQuote{new} seqlevels. 
  }
  \item{...}{
    Arguments passed to other functions.
  }
}

\value{
  The \code{x} object with seqlevels removed or renamed. If \code{x} has
  no seqlevels (empty object) or no replacement values match the current
  seqlevels in \code{x} the unchanged \code{x} is returned.
}

\details{
  Many operations on \link{GRanges} objects require the seqlevels to match 
  before a comparison can be made (e.g., \code{findOverlaps(type="within")}).
  \code{keepSeqlevels} and \code{renameSeqlevels} are convenience functions 
  for subsetting and renaming the seqlevels of these objects. 

  \code{keepSeqlevels} subsets \code{x} based on the seqlevels provided in 
  \code{value}. If \code{value} does not match any seqlevels in \code{x} the 
  original \code{x} is returned. When \code{x} is a \link{GRangesList}, there 
  may be multiple chromosomes in a single list element. If not all chromosomes 
  are specified in \code{value}, a reduced list element is returned. All empty 
  list elements are dropped. See examples.

  \code{renameSeqlevels} renames the seqlevels in \code{x} to those provided 
  in \code{value}. \code{value} is a named character vector where the names are 
  the \sQuote{old} seqlevel names and the values are the \sQuote{new}. If 
  no names in \code{value} match the seqlevels in \code{x} the original \code{x} 
  is returned. 
}

\seealso{
  \link[VariantAnnotation]{VCF}
}

\examples{
gr1 <- GRanges(seqnames = c("chr1", "chr2"),
               ranges = IRanges(c(7,13), width = 3),
               strand = c("+", "-"), score = 3:4, GC = c(0.3, 0.5))
gr2 <- GRanges(seqnames = c("chr1", "chr1", "chr2", "chr3", "chr3"),
               ranges = IRanges(c(1, 4, 8, 9, 16), width=5),
               strand = "-", score = c(3L, 2L, 5L, 6L, 2L), 
               GC = c(0.4, 0.1, 0.55, 0.20, 0.10))
gr3 <- GRanges(seqnames = c("CHROM4", "CHROM4"), 
               ranges = IRanges(c(20, 45), width=6), 
               strand = "+", score = c(2L, 5L), GC = c(0.30, 0.45))
 
## GRanges :
gr3_rename <- renameSeqlevels(gr3, c(CHROM4="chr4"))
gr3_rename
 
gr2_subset_chr <- keepSeqlevels(gr2, c("chr1", "chr2"))
gr2_subset_gr <- keepSeqlevels(gr2, gr1)
identical(gr2_subset_chr, gr2_subset_gr)
 
## GRangesList :
grl1 <- GRangesList("gr1" = gr1, "gr2" = gr2, "gr3" = gr3)
grl2 <- GRangesList("gr1" = gr1, "gr2" = gr2, "gr3" = gr3_rename)
grl1_rename <- renameSeqlevels(grl1, c(CHROM4="chr4"))
identical(grl1_rename, grl2)
 
grl1_subset <- keepSeqlevels(grl1, "chr3")
 
## GAlignments :
library(Rsamtools)
galn_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
galn <- readGAlignments(galn_file)
 
galn_rename <- renameSeqlevels(galn, c(seq2="chr2"))
galn_subset <- keepSeqlevels(galn_rename, gr1)
galn_subset

## See ?VCF for examples using renameSeqlevels and keepSeqlevels with 
## VCF class objects
}

\author{Valerie Obenchain \email{vobencha@fhcrc.org}}

\keyword{methods}
\keyword{utilities}


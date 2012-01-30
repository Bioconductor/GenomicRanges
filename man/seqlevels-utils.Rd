\name{utils}
\docType{methods}
\alias{keepSeqlevels}
\alias{keepSeqlevels,GenomicRanges,GenomicRanges-method}
\alias{keepSeqlevels,GenomicRanges,GRangesList-method}
\alias{keepSeqlevels,GenomicRanges,GappedAlignments-method}
\alias{keepSeqlevels,GenomicRanges,character-method}
\alias{keepSeqlevels,GRangesList,GenomicRanges-method}
\alias{keepSeqlevels,GRangesList,GRangesList-method}
\alias{keepSeqlevels,GRangesList,GappedAlignments-method}
\alias{keepSeqlevels,GRangesList,character-method}
\alias{keepSeqlevels,GappedAlignments,GenomicRanges-method}
\alias{keepSeqlevels,GappedAlignments,GRangesList-method}
\alias{keepSeqlevels,GappedAlignments,GappedAlignments-method}
\alias{keepSeqlevels,GappedAlignments,character-method}
\alias{keepSeqlevels,VCF,character-method}

\alias{renameSeqlevels}
\alias{renameSeqlevels,GappedAlignments,character-method}
\alias{renameSeqlevels,GRangesList,character-method}
\alias{renameSeqlevels,GenomicRanges,character-method}
\alias{renameSeqlevels,VCF,character-method}

\title{
  seqlevels utility functions
}

\description{
  Rename or subset the seqlevels in a \link{GenomicRanges}, \link{GRangesList} 
  or \link{GappedAlignments} object.
}

\usage{
  \S4method{keepSeqlevels}{GenomicRanges,GenomicRanges}(x, value, ...) 
  \S4method{keepSeqlevels}{GenomicRanges,GRangesList}(x, value, ...) 
  \S4method{keepSeqlevels}{GenomicRanges,GappedAlignments}(x, value, ...) 
  \S4method{keepSeqlevels}{GenomicRanges,character}(x, value, ...) 
  \S4method{keepSeqlevels}{GRangesList,GenomicRanges}(x, value, ...) 
  \S4method{keepSeqlevels}{GRangesList,GRangesList}(x, value, ...) 
  \S4method{keepSeqlevels}{GRangesList,GappedAlignments}(x, value, ...) 
  \S4method{keepSeqlevels}{GRangesList,character}(x, value, ...) 
  \S4method{keepSeqlevels}{GappedAlignments,GenomicRanges}(x, value, ...) 
  \S4method{keepSeqlevels}{GappedAlignments,GRangesList}(x, value, ...) 
  \S4method{keepSeqlevels}{GappedAlignments,GappedAlignments}(x, value, ...) 
  \S4method{keepSeqlevels}{GappedAlignments,character}(x, value, ...) 
  \S4method{keepSeqlevels}{VCF,character}(x, value, ...) 

  \S4method{renameSeqlevels}{GenomicRanges,character}(x, value, ...) 
  \S4method{renameSeqlevels}{GRangesList,character}(x, value, ...) 
  \S4method{renameSeqlevels}{GappedAlignments,character}(x, value, ...) 
  \S4method{renameSeqlevels}{VCF,character}(x, value, ...) 
}

\arguments{
  \item{x}{
    The \link{GenomicRanges}, \link{GRangesList} or \link{GappedAlignments} object 
    for which the seqlevels will be removed or renamed.
  }
  \item{value}{
    For \code{keepSeqlevels}, a \link{GRanges}, \link{GRangesList} or
    \code{GappedAlignments} or character vector. The seqlevels from 
    \code{value} object will be used to subset the seqlevels in the
    \code{x}. 

    For \code{renameSeqlevels}, a named character vector where the names are the 
    'old' seqlevels and the values are the 'new' seqlevels. seqlevels in 
    the \code{x} will be renamed from the 'old' to the 'new'. 
  }
  \item{...}{
    Arguments passed to other functions.
  }
}

\value{
  The \code{x} object with seqlevels removed or renamed.
}

\details{
  Many operations on \link{GRanges} objects require the seqlevels 
  to match before a comparison can be made 
  (e.g., \code{findOverlaps(type="within")}).
  \code{keepSeqlevels} and \code{renameSeqlevels} are convenience functions 
  for subsetting and renaming the seqlevels of these objects. 

  The \code{keepSeqlevels} function subsets the seqlevels of \code{x} 
  based on the seqlevels provided in \code{value}. seqlevels contained
  in \code{value} that are not in \code{x} are ignored. If
  none of the seqlevels in the \code{value} are present in
  \code{x} an error is thrown. 
  When \code{x} is a \link{GRangesList}, keepSeqlevels retains only
  the seqlevels specified in \code{value}. If a list element
  has multiple chromosomes but not all chromosomes are specified
  in \code{value}, a reduced list element is returned.
  All empty list elements are dropped. See examples.

  \code{renameSeqlevels} renames the seqlevels in \code{x} based on
  those provided in \code{value}. \code{value} is a named
  character vector where the names should match existing seqlevels in
  \code{x}. The matching seqlevels will be renamed to the corresponding
  list values. List names in \code{value} that are not in 
  \code{x} are ignored. An error is thrown if none of the list names in 
  \code{value} match the seqlevels in \code{x}

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

## GappedAlignments :
library(Rsamtools)
galn_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
galn <- readGappedAlignments(galn_file)

galn_rename <- renameSeqlevels(galn, c(seq2="chr2"))
galn_subset <- keepSeqlevels(galn_rename, gr1)
galn_subset
}


\author{Valerie Obenchain \email{vobencha@fhcrc.org}}

\keyword{methods}
\keyword{utilities}


\name{makeSeqnameIds}

\alias{makeSeqnameIds}
\alias{sortSeqlevels}

\title{Sort sequence levels in "natural" order}

\description{
  \code{sortSeqlevels} sorts sequence levels in "natural" order.

  \code{makeSeqnameIds} assigns a unique ID to each unique sequence name.
  The returned IDs span 1:N where N is the number of unique sequence names
  provided. \code{sortSeqlevels} is based on \code{makeSeqnameIds}.
}

\usage{
sortSeqlevels(seqlevels, X.is.sexchrom=NA)
makeSeqnameIds(seqnames, X.is.sexchrom=NA)
}

\arguments{
  \item{seqlevels}{
    A character vector containing sequence names. Typically obtained with
    \code{\link{seqlevels}()}.
  }
  \item{X.is.sexchrom}{
    A logical indicating whether X refers to the sexual chromosome 
    or to chromosome with Roman Numeral X.
  }
  \item{seqnames}{
    A character vector or factor containing sequence names.
  }
}

\value{
  For \code{sortSeqlevels}: A character level containing the sequence levels
  sorted in "natural" order.

  For \code{makeSeqnameIds}: An integer vector the same length as
  \code{seqnames}. The values span 1:N where N is the number of unique
  sequence names provided.
}

\author{
  H. Pages
}

\examples{
sortSeqlevels(c("11", "Y", "1", "10", "9", "M", "2"))

seqlevels <- c("chrXI", "chrY", "chrI", "chrX", "chrIX", "chrM", "chrII")
sortSeqlevels(seqlevels)
sortSeqlevels(seqlevels, X.is.sexchrom=TRUE)
sortSeqlevels(seqlevels, X.is.sexchrom=FALSE)

seqlevels <- c("chr2RHet", "chr4", "chrUextra", "chrYHet",
               "chrM", "chrXHet", "chr2LHet", "chrU",
               "chr3L", "chr3R", "chr2R", "chrX")
sortSeqlevels(seqlevels)

library(BSgenome.Scerevisiae.UCSC.sacCer2)
makeSeqnameIds(seqnames(Scerevisiae))
makeSeqnameIds(seqnames(Scerevisiae)[c(1:5,5:1)])
}

\keyword{manip}

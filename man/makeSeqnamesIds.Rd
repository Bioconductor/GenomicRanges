\name{makeSeqnameIds}

\alias{makeSeqnameIds}

\title{Assigns a unique ID to each unique sequence name}

\description{
    Assign a unique ID to each unique sequence name. The 
    returned IDs span 1:N where N is the number of unique
    sequence names provided.
}

\usage{
makeSeqnameIds(seqnames, X.is.sexchrom=NA)
}

\arguments{
  \item{seqnames}{
    A character vector or factor of sequence names.
  }
  \item{X.is.sexchrom}{
    A logical indicating whether X refers to the sexual chromosome 
    or to chromosome with Roman Numeral X.
  }
}

\value{
  An integer vector the same length as \code{seqnames}. The values
  span 1:N where N is the number of unique sequence names provided.
}

\author{
  H. Pages
}

\examples{
library(BSgenome.Scerevisiae.UCSC.sacCer2)
makeSeqnameIds(seqnames(Scerevisiae))
makeSeqnameIds(seqnames(Scerevisiae)[c(1:5,5:1)])
}

\keyword{manip}

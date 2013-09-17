\name{makeSeqnameIds}

\alias{makeSeqnameIds}

\title{Assign sequence IDs to sequence names}

\description{
  \code{makeSeqnameIds} assigns a unique ID to each unique sequence name
  in the input vector. The returned IDs span 1:N where N is the number of
  unique sequence names in the input vector.
}

\usage{
makeSeqnameIds(seqnames, X.is.sexchrom=NA)
}

\arguments{
  \item{seqnames}{
    A character vector or factor containing sequence names.
  }
  \item{X.is.sexchrom}{
    A logical indicating whether X refers to the sexual chromosome 
    or to chromosome with Roman Numeral X.
  }
}

\value{
  An integer vector of the same length as \code{seqnames}.
  The values in the vector span 1:N where N is the number of unique
  sequence names in the input vector.
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item \code{\link{sortSeqlevels}} for sorting the sequence levels
          of an object in "natural" order.
  }
}

\examples{
library(BSgenome.Scerevisiae.UCSC.sacCer2)
makeSeqnameIds(seqnames(Scerevisiae))
makeSeqnameIds(seqnames(Scerevisiae)[c(1:5,5:1)])
}

\keyword{manip}

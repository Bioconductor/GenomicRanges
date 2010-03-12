\name{seqnames}

\alias{seqnames}
\alias{seqnames<-}

\title{Accessing the sequence names stored in an object}

\description{The \code{seqnames} generic is meant as an accessor for
  the sequence names stored in an object.
}

\usage{seqnames(x)}

\arguments{
  \item{x}{The object from which to access the sequence names.}
}

\details{
  Various classes provide a \code{seqnames} method for getting or
  setting the sequence names of an object of the class.

  See for example \code{?`\link{seqnames,GRanges-method}`} for the
  \link{GRanges} class and \code{?`\link{seqnames,GRangesList-method}`}
  for the \link{GRangesList} class, both classes being defined in
  the GenomicRanges package.
  See also \code{?`\link[BSgenome]{seqnames,BSgenome-method}`} for the
  \link[BSgenome]{BSgenome} class defined in the BSgenome package (you
  might need to install and load the package first).

  The full list can be seen with \code{showMethods("seqnames")}
  (for the get methods) and \code{showMethods("seqnames<-")}
  (for the replacement methods).
}

\seealso{
  \link{seqnames,GRanges-method},
  \link{seqnames,GRangesList-method},
  \link[BSgenome]{seqnames,BSgenome-method}
}

\examples{
  showMethods("seqnames")
  showMethods("seqnames<-")

  ?`seqnames,GRanges-method`
}

\keyword{methods}

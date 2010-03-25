\name{seqnames}

\alias{seqnames}
\alias{seqnames<-}

\alias{seqlengths}
\alias{seqlengths<-}

\title{Accessing the sequence names and lengths stored in an object}

\description{The \code{seqnames} and \code{seqlengths} generics are meant
  as accessors for the sequence names and lengths stored in an object,
  respectively.
}

\usage{
  seqnames(x)
  seqlengths(x)
}

\arguments{
  \item{x}{The object from which to access the sequence information.}
}

\details{
  Various classes provide \code{seqnames} and \code{seqlengths} methods
  for getting or setting the sequence names and lengths of an object of
  the class, respectively.

  See for example \code{?`\link{seqnames,GRanges-method}`} for the
  \link{GRanges} class and \code{?`\link{seqnames,GRangesList-method}`}
  for the \link{GRangesList} class, both classes being defined in
  the GenomicRanges package.

  For more examples see
  \code{?`\link[BSgenome]{seqnames,BSgenome-method}`} for the
  \link[BSgenome]{BSgenome} class defined in the BSgenome package
  and \code{?`\link[GenomicFeatures]{seqnames,TranscriptDb-method}`}
  for the \link[GenomicFeatures]{TranscriptDb} class defined in the
  GenomicFeatures package (you might need to install and load these
  packages first).

  The full list can be seen with \code{showMethods("seqnames")} and
  \code{showMethods("seqlengths")} (for the get methods) and
  \code{showMethods("seqnames<-")} and \code{showMethods("seqlengths<-")}
  (for the replacement methods).
}

\seealso{
  \link{seqnames,GRanges-method},
  \link{seqnames,GRangesList-method},
  \link[BSgenome]{seqnames,BSgenome-method}
  \link[GenomicFeatures]{seqnames,TranscriptDb-method}
}

\examples{
  showMethods("seqnames")
  showMethods("seqnames<-")

  showMethods("seqlengths")
  showMethods("seqlengths<-")

  ?`seqnames,GRanges-method`
}

\keyword{methods}

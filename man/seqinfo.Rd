\name{seqinfo}

\alias{seqinfo}
\alias{seqinfo<-}
\alias{seqnames}
\alias{seqnames<-}
\alias{seqlengths}
\alias{seqlengths,ANY-method}
\alias{seqlengths<-}
\alias{isCircular}
\alias{isCircular,ANY-method}
\alias{isCircular<-}
\alias{isCircularWithKnownLength}
\alias{isCircularWithKnownLength,ANY-method}


\title{Accessing sequence information}

\description{
  A set of generic functions for getting/setting sequence information
  such as the sequence names, lengths and circularity flags.
}

\usage{
  seqinfo(x)
  seqinfo(x) <- value

  seqnames(x)
  seqnames(x) <- value

  seqlengths(x)
  seqlengths(x) <- value

  isCircular(x)
  isCircular(x) <- value

  isCircularWithKnownLength(x)
}

\arguments{
  \item{x}{The object from which to access the sequence information.}
  \item{value}{A vector containing the sequence information to store.}
}

\details{
  Various classes provide methods for those generic functions.

  See for example \code{?`\link{seqinfo,GRanges-method}`} for the
  \link{GRanges} class and \code{?`\link{seqinfo,GRangesList-method}`}
  for the \link{GRangesList} class, both classes being defined in
  the GenomicRanges package.

  For more examples see
  \code{?`\link[BSgenome]{seqnames,BSgenome-method}`} for the
  \link[BSgenome]{BSgenome} class defined in the BSgenome package
  and \code{?`\link[GenomicFeatures]{seqnames,TranscriptDb-method}`}
  for the \link[GenomicFeatures]{TranscriptDb} class defined in the
  GenomicFeatures package (you might need to install and load these
  packages first).

  The full list of methods defined for a given generic can be seen
  with e.g. \code{showMethods("seqinfo")} or
  \code{showMethods("seqnames")} (for the getters),
  and \code{showMethods("seqinfo<-")} or \code{showMethods("seqnames<-")}
  (for the setters aka \emph{replacement methods}). Note that this
  shows only methods defined in packages that are currently attached.

  \code{isCircularWithKnownLength(x)} is formally defined as
  \code{(isCircular(x) \%in\% TRUE) & !is.na(seqlengths(x))}.
}

\seealso{
  \link{Seqinfo-class},
  \link{seqinfo,GRanges-method},
  \link{seqinfo,GRangesList-method},
  \link[BSgenome]{seqnames,BSgenome-method}
  \link[GenomicFeatures]{seqnames,TranscriptDb-method}
}

\examples{
  showMethods("seqinfo")
  showMethods("seqinfo<-")

  showMethods("seqnames")
  showMethods("seqnames<-")

  if (interactive())
    ?`seqnames,GRanges-method`
}

\keyword{methods}

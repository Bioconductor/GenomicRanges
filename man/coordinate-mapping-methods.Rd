\name{coordinate-mapping-methods}
\alias{coordinate-mapping-methods}
\alias{coordinate-mapping}

\alias{mapToGenome}
\alias{pmapToGenome}

\alias{mapToTranscript}
\alias{pmapToTranscript}


\title{Map range coordinates between genome and transcript space}

\description{
  Coordinate mapping methods for \link{GenomicRanges} objects.

  Methods in this package are under construction. See the
  \code{?`coordinate-mapping-methods`} man page in \pkg{GenomicAlignments}
  for a description of existing methods.
}

\usage{
}

\arguments{
  \item{from}{
    A \code{\linkS4class{Ranges}} or \code{\linkS4class{GenomicRanges}}
    object of positions to be mapped.

    \itemize{
      \item{} For \code{mapToTranscript} \code{from} is in genome-based 
        coordinates 
      \item{} For \code{mapToGenome} \code{from} is in transcript-based 
        coordinates 
    }
  }
  \item{to}{
    A \code{\linkS4class{GenomicRanges}} object that represents the alignment 
    (aka mapping) between the sequences in \code{from} and the sequences in the
    result.

    The current implementation assumes that positions in \code{from} 
    occur within the \code{to} region.
  }
  \item{ignore.strand}{
    When TRUE, strand is ignored in overlap operations.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\details{
  See \code{?`coordinate-mapping-methods`} in 
  the \pkg{GenomicAlignments} package for existing methods.
}

\value{
  An object the same class as \code{from}.
}

\seealso{
  \itemize{
    \item \code{?`coordinate-mapping-methods`} in the \pkg{GenomicAlignments} 
      package for methods on \link{GAlignments} objects.
  }
}

\examples{
## See examples in GenomicAlignments package.
}

\keyword{methods}
\keyword{utilities}

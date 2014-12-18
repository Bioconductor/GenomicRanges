\name{coordinate-mapping-methods}
\alias{coordinate-mapping-methods}
\alias{coordinate-mapping}

%\alias{mapToGenome}
%\alias{mapToGenome,GenomicRanges,GenomicRanges}
%\alias{mapToGenome,GenomicRanges,GRangesList}
%\alias{pmapToGenome}
%\alias{pmapToGenome,GenomicRanges,GenomicRanges}
%\alias{pmapToGenome,GenomicRanges,GRangesList}

\alias{mapToTranscript}
\alias{mapToTranscript,GenomicRanges,GenomicRanges}
\alias{mapToTranscript,GenomicRanges,GRangesList}
\alias{pmapToTranscript}
\alias{pmapToTranscript,GenomicRanges,GenomicRanges}
\alias{pmapToTranscript,GenomicRanges,GRangesList}


\title{Map range coordinates between genome and transcript space}

\description{
  Coordinate mapping methods for \linkS4class{GenomicRanges} objects.

  Methods in this package are under construction. See the
  \code{?`coordinate-mapping-methods`} man page in \pkg{GenomicAlignments}
  for a description of existing methods.
}

\usage{
}

\arguments{
  \item{x}{
    A \code{\linkS4class{Ranges}} or \code{\linkS4class{GenomicRanges}}
    object of positions to be mapped.

    \itemize{
      \item{} For \code{mapToTranscript} \code{x} is in genome-based 
        coordinates 
      \item{} For \code{mapToGenome} \code{x} is in transcript-based 
        coordinates 
    }
  }
  \item{alignment}{
    A \code{\linkS4class{GenomicRanges}} object that represents the alignment 
    (aka mapping) between the sequences in \code{x} and the sequences in the
    result.

    The current implementation assumes that positions in \code{x} 
    occur within the \code{alignment} region.
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

  Methods where \code{alignment} is GRangesList assumes all inner list elements
  in the same outer list element are of the same strand. 

}

\value{
  An object the same class as \code{x}.
}

\seealso{
  \itemize{
    \item \code{?`coordinate-mapping-methods`} in the \pkg{GenomicAlignments} 
      package for methods on \linkS4class{GAlignments} objects.
  }
}

\examples{
}

\keyword{methods}
\keyword{utilities}

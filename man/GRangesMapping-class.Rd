\name{GRangesMapping-class}
\docType{class}
\alias{GRangesMapping-class}
\alias{GRangesMapping}

% accessors
\alias{hits}
\alias{granges,GRangesMapping-method}
\alias{dim,GRangesMapping-method}
\alias{length,GRangesMapping-method}
\alias{queryHits,GRangesMapping-method}
\alias{subjectHits,GRangesMapping-method}

% coercion
\alias{coerce,GRangesMapping,RangedData-method}
\alias{coerce,GRangesMapping,GenomicRanges-method}


\title{GRangesMapping objects}

\description{
  The GRangesMapping class is a container for both the genomic 
  locations and the matching results output from the \code{map} 
  function.
}

\details{
  The \code{map} generic converts a set of ranges to the equivalent
  ranges on another sequence, through some sort of alignment between
  sequences, and outputs a \code{GRangesMapping} object. There are two
  components of that object: the transformed GRanges, and the
  \code{hits}, a \code{\linkS4class{Hits}} object of the same length
  that matches each input range to a destination sequence (useful when
  the alignment is one/many to many).
}

\section{\code{GRangesMapping} Accessors}{
  \describe{
    \item{}{\code{granges(x)}: Gets the mapped GRanges.}
    \item{}{\code{hits(x)}: Gets the matching between the input ranges
      and the destination sequences (of which there may be more than
      one).}
    \item{}{\code{dim(x)}: Same as \code{dim(hits(x))}.}
    \item{}{\code{length(x)}: Same as \code{length(hits(x))}.}
    \item{}{\code{subjectHits(x)}: Same as \code{subjectHits(hits(x))}.}
    \item{}{\code{queryHits(x)}: Same as \code{queryHits(hits(x))}.}
  }
}

\section{\code{GRangesMapping} Coercion}{
  \describe{
    \item{}{\code{as(from, "RangedData")}: Converts
      a \code{GRangesMapping} into
      a \code{\linkS4class{RangedData}}. The \code{ranges}/\code{space}
      in the \code{RangedData} are the \code{ranges}/\code{space}
      of \code{from}, and the \code{values} result from the coercion of
      the \code{hits} to a \code{DataFrame}.
    }
    \item{}{\code{as(from, "GenomicRanges")}: Returns the 
      \code{\linkS4class{GRanges}} with the \code{hits}
      (coerced to a \code{DataFrame}) as the \code{mcols}.
    }
  }
}

\author{ Michael Lawrence }

\seealso{
  Methods on the generic \code{map}, which generates an instance of
  this class, are defined in other packages, like GenomicRanges.
}

\examples{
# See ?map-methods in GenomicRanges or GenomicAlignments for examples.
}

\keyword{methods}
\keyword{classes}

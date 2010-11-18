\name{seqinfo}

\alias{seqinfo}
\alias{seqinfo<-}
\alias{seqnames}
\alias{seqnames<-}
\alias{seqlevels}
\alias{seqlevels,ANY-method}
\alias{seqlevels<-}
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
  from/on an object.
}

\usage{
  seqinfo(x)
  seqinfo(x) <- value

  seqnames(x)
  seqnames(x) <- value

  seqlevels(x)
  seqlevels(x) <- value

  seqlengths(x)
  seqlengths(x) <- value

  isCircular(x)
  isCircular(x) <- value

  isCircularWithKnownLength(x)
}

\arguments{
  \item{x}{The object from/on which to get/set the sequence information.}
  \item{value}{Typically a \link{Seqinfo} object for the \code{seqinfo}
  setter. A vector containing the sequence information to store for the
  other setters.}
}

\details{
  Various classes implement methods for those generic functions.

  The \link{Seqinfo} class plays a central role for those generics because:
  \itemize{
  \item It implements methods for all the above generics except for
        \code{seqinfo}. That is, the \code{seqnames}, \code{seqlevels},
        \code{seqlengths}, \code{isCircular} and
        \code{isCircularWithKnownLength} getters are defined
        for \link{Seqinfo} objects.
  \item Classes that implement the \code{seqinfo} getter are typically
        expected to return a \link{Seqinfo} object.
  \item Default \code{seqlevels}, \code{seqlengths}, \code{isCircular}
        and \code{isCircularWithKnownLength} getters are provided.
        By default, \code{seqlevels(x)} does \code{seqlevels(seqinfo(x))},
        \code{seqlengths(x)} does \code{seqlengths(seqinfo(x))},
        \code{isCircular(x)} does \code{isCircular(seqinfo(x))},
        and \code{isCircularWithKnownLength(x)} does
        \code{isCircularWithKnownLength(seqinfo(x))}.
        So any class with a \code{seqinfo} getter that returns a
        \link{Seqinfo} object also has all those getters working
        out-of-the-box.
        See the \link{GRanges}, \link{GRangesList} and \link{GappedAlignments}
        classes as examples of such classes (those 3 classes are defined
        in the GenomicRanges package).
  }
}

\note{
  The full list of methods defined for a given generic can
  be seen with e.g. \code{showMethods("seqinfo")} or
  \code{showMethods("seqnames")} (for the getters),
  and \code{showMethods("seqinfo<-")} or \code{showMethods("seqnames<-")}
  (for the setters aka \emph{replacement methods}).
  Please be aware that this shows only methods defined in packages
  that are currently attached.
}

\seealso{
  \link{Seqinfo-class},
  \link{GRanges-class},
  \link{GRangesList-class},
  \link{GappedAlignments-class}
}

\examples{
  showMethods("seqinfo")
  showMethods("seqinfo<-")

  showMethods("seqnames")
  showMethods("seqnames<-")

  if (interactive())
    ?`GRanges-class`
}

\keyword{methods}

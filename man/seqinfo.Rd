\name{seqinfo}

\alias{seqinfo}
\alias{seqinfo<-}
\alias{seqnames}
\alias{seqnames<-}
\alias{seqlevels}
\alias{seqlevels,ANY-method}
\alias{seqlevels<-}
\alias{seqlevels<-,ANY-method}
\alias{seqlengths}
\alias{seqlengths,ANY-method}
\alias{seqlengths<-}
\alias{seqlengths<-,ANY-method}
\alias{isCircular}
\alias{isCircular,ANY-method}
\alias{isCircular<-}
\alias{isCircular<-,ANY-method}
\alias{genome}
\alias{genome,ANY-method}
\alias{genome<-}
\alias{genome<-,ANY-method}
\alias{seqnameStyle}
\alias{seqnameStyle,ANY-method}
\alias{seqnameStyle<-}
\alias{seqnameStyle<-,ANY-method}

\alias{seqinfo,List-method}
\alias{seqinfo<-,List-method}
\alias{seqinfo,RangedData-method}
\alias{seqinfo<-,RangedData-method}
\alias{seqinfo,RangesList-method}
\alias{seqnames,RangesList-method}
\alias{seqnames,RangedData-method}

\title{Accessing sequence information}

\description{
  A set of generic functions for getting/setting sequence information
  from/on an object.
}

\usage{
seqinfo(x)
seqinfo(x, new2old=NULL, force=FALSE) <- value

seqnames(x)
seqnames(x) <- value

seqlevels(x)
seqlevels(x, force=FALSE) <- value

seqlengths(x)
seqlengths(x) <- value

isCircular(x)
isCircular(x) <- value

genome(x)
genome(x) <- value

seqnameStyle(x)
seqnameStyle(x) <- value
}

\arguments{
  \item{x}{
    The object from/on which to get/set the sequence information.
  }
  \item{new2old}{
    The \code{new2old} argument allows the user to rename, drop, add and/or
    reorder the "sequence levels" in \code{x}.

    \code{new2old} can be \code{NULL} or an integer vector with one element
    per row in \link{Seqinfo} object \code{value} (i.e. \code{new2old} and
    \code{value} must have the same length) describing how the "new" sequence
    levels should be mapped to the "old" sequence levels, that is, how the
    rows in \code{value} should be mapped to the rows in \code{seqinfo(x)}.
    The values in \code{new2old} must be >= 1 and <= \code{length(seqinfo(x))}.
    \code{NA}s are allowed and indicate sequence levels that are being added.
    Old sequence levels that are not represented in \code{new2old} will be
    dropped, but this will fail if those levels are in use (e.g. if \code{x}
    is a \link{GRanges} object with ranges defined on those sequence levels)
    unless \code{force=TRUE} is used (see below).

    If \code{new2old=NULL}, then sequence levels can only be added to the
    existing ones, that is, \code{value} must have at least as many rows
    as \code{seqinfo(x)} (i.e. \code{length(values) >= length(seqinfo(x))})
    and also \code{seqlevels(values)[seq_len(length(seqlevels(x)))]} must be
    identical to \code{seqlevels(x)}.
  }
  \item{force}{
    Force dropping sequence levels currently in use. This is achieved by
    dropping the elements in \code{x} where those levels are used (hence
    typically reducing the length of \code{x}).
  }
  \item{value}{
    Typically a \link{Seqinfo} object for the \code{seqinfo} setter.

    Either a named or unnamed character vector for the \code{seqlevels}
    setter.

    A vector containing the sequence information to store for the other
    setters.
  }
}

\details{
  Various classes implement methods for those generic functions.

  The \link{Seqinfo} class plays a central role for those generics because:

  \itemize{
  \item It has methods for all those generics (except \code{seqinfo}).
        That is, the \code{seqnames}, \code{seqlevels}, \code{seqlengths},
        \code{isCircular}, \code{genome} and \code{seqnameStyle} getters
        and setters are defined for \link{Seqinfo} objects.

  \item For classes that implement it, the \code{seqinfo} getter should
        return a \link{Seqinfo} object.

  \item Default \code{seqlevels}, \code{seqlengths}, \code{isCircular},
        \code{genome}, and \code{seqnameStyle} getters and setters are
        provided.
        By default, \code{seqlevels(x)} does \code{seqlevels(seqinfo(x))},
        \code{seqlengths(x)} does \code{seqlengths(seqinfo(x))},
        \code{isCircular(x)} does \code{isCircular(seqinfo(x))},
        \code{genome(x)} does \code{genome(seqinfo(x))},
        and \code{seqnameStyle(x)} does \code{seqnameStyle(seqinfo(x))}.
        So any class with a \code{seqinfo} getter will have all the above
        getters work out-of-the-box. If, in addition, the class defines
        a \code{seqinfo} setter, then all the corresponding setters will
        also work out-of-the-box.

        See the \link{GRanges}, \link{GRangesList}, \link{GappedAlignments},
        and \link{GappedAlignmentPairs} classes for examples of classes that
        define the \code{seqinfo} getter and setter (those 4 classes are
        defined in the GenomicRanges package).

        See the \link[GenomicFeatures]{TranscriptDb} class (defined in the
        GenomicFeatures package) for an example of a class that defines only
        the \code{seqinfo} getter (no setter).
      }

      The GenomicRanges package defines \code{seqinfo} and
      \code{seqinfo<-} methods for these low-level IRanges data
      structures: \code{List}, \code{RangesList} and
      \code{RangedData}. Those objects do not have the means to formally
      store sequence information. Thus, the wrappers simply store the
      \code{Seqinfo} object within \code{metadata(x)}. Initially, the
      metadata is empty, so there is some effort to generate a
      reasonable default \code{Seqinfo}. The names of any \code{List}
      are taken as the \code{seqnames}, and the \code{universe} of
      \code{RangesList} or \code{RangedData} is taken as the
      \code{genome}.
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
  \itemize{
    \item \link{Seqinfo-class}.
    \item \link{GRanges-class}.
    \item \link{GRangesList-class}.
    \item \link{GappedAlignments-class}.
    \item \link{GappedAlignmentPairs-class}.
    \item \link[GenomicFeatures]{TranscriptDb-class}.
  }
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

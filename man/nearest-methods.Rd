\name{nearest-methods}

\alias{nearest-methods}

\alias{precede,GenomicRanges,GenomicRanges-method}
\alias{precede,GenomicRanges,missing-method}
\alias{follow,GenomicRanges,GenomicRanges-method}
\alias{follow,GenomicRanges,missing-method}
\alias{nearest,GenomicRanges,GenomicRanges-method}
\alias{nearest,GenomicRanges,missing-method}
\alias{distance,GenomicRanges,GenomicRanges-method}
\alias{distanceToNearest,GenomicRanges,GenomicRanges-method}
\alias{distanceToNearest,GenomicRanges,missing-method}

\title{Finding the nearest genomic range neighbor}

\description{
  The \code{nearest}, \code{precede} and \code{follow} methods find
  nearest neighbors between \link{GenomicRanges} objects (i.e. objects
  that belong to the \link{GenomicRanges} class or one of its subclasses,
  this includes for example \link{GRanges} objects).
}

\usage{
\S4method{precede}{GenomicRanges,GenomicRanges}(x, subject,
        select=c("arbitrary", "all"), ignore.strand=FALSE, ...)

\S4method{precede}{GenomicRanges,missing}(x, subject,
        select=c("arbitrary", "all"), ignore.strand=FALSE, ...)

\S4method{follow}{GenomicRanges,GenomicRanges}(x, subject,
        select=c("arbitrary", "all"), ignore.strand=FALSE, ...)

\S4method{follow}{GenomicRanges,missing}(x, subject,
        select=c("arbitrary", "all"), ignore.strand=FALSE, ...)

\S4method{nearest}{GenomicRanges,GenomicRanges}(x, subject,
        select=c("arbitrary", "all"), ignore.strand=FALSE, ...)

\S4method{nearest}{GenomicRanges,missing}(x, subject,
        select=c("arbitrary", "all"), ignore.strand=FALSE, ...)

\S4method{distance}{GenomicRanges,GenomicRanges}(x, y, ignore.strand=FALSE, ...)

\S4method{distanceToNearest}{GenomicRanges,GenomicRanges}(x, subject,
        ignore.strand=FALSE, ...)
\S4method{distanceToNearest}{GenomicRanges,missing}(x, subject,
        ignore.strand=FALSE, ...)
}

\arguments{
  \item{x, subject}{\link{GenomicRanges} objects.}
  \item{select}{See \code{?`\link[IRanges]{nearest-methods}`}.}
  \item{ignore.strand}{
    \code{TRUE} or \code{FALSE}. Whether the strand of the input ranges
    should be ignored or not. See details below.
  }
  \item{...}{Additional arguments for methods.}
}

\details{
\code{precede} identifies which subject(s) the query precedes. Returns the
index of the range in \code{subject} that is directly preceded by the range
in \code{x}. When \code{ignore.strand=TRUE}, the strand for both \code{x}
and \code{subject} are set to \code{"+"}. When \code{select="arbitrary"}
an integer vector is returned with a single match per \code{x}. If no match
is found an \code{NA} is returned. When \code{select="all"} a
\link[IRanges]{Hits} object is returned with all matches for \code{x}.
If \code{x} does not have a match in \code{subject} the \code{x} is not
included in the \code{Hits} object. Overlapping ranges are excluded.
Matching by strand :
\itemize{
  \item{}{x on + strand can match to ranges on both + and 
    * strands. In the case of a tie the first range by order 
    is chosen.
  } 
  \item{}{x on - strand can match to ranges on both - and 
    * strands. In the case of a tie the first range by order 
    is chosen.
  } 
  \item{}{x on * strand can match to ranges on any of +, 
    - or * strands. In the case of a tie the first range by
    order is chosen.
  }
}

\code{follow} identifies which subject(s) the query follows. Returns the
index of the interval in \code{subject} that is directly followed by the
range in \code{x}. When \code{ignore.strand=TRUE}, both \code{x} and
\code{subject} strand are set to \code{"+"}.  When \code{select="arbitrary"}
an integer vector is returned with a single match per \code{x}. If no match
is found an \code{NA} is returned. When \code{select="all"} a
\link[IRanges]{Hits} object is returned with all matches for \code{x}.
If \code{x} does not have a match in \code{subject} the \code{x} is not
included in the \code{Hits} object. Overlapping ranges are excluded. See
the \code{precede} above for details of strand matching.

\code{nearest} does conventional nearest neighbor finding. It returns an
integer vector containing the index of the nearest neightbor range in
\code{subject} for each range in \code{x}. If there is no nearest neighbor
\code{NA} is returned. If \code{ignore.strand=TRUE}, both \code{x} and 
\code{subject} strand are set to \code{"+"}. \code{nearest} calls 
\code{precede} and \code{follow} and thus the strand matching for \code{*} 
follows the conventions documented under those functions. 

\code{distance} calculates the number of positions separating two features.
The value is zero if the features overlap and \code{NA} if the features
are on different sequences, or different strands (if \code{ignore.strand} is
\code{FALSE}).

\code{distanceToNearest} returns the distance for each range in \code{x}
to its nearest neighbor in the \code{subject}. If \code{ignore.strand=TRUE},
both \code{x} and \code{subject} strand are set to \code{"+"}. 
}

\author{P. Aboyoun and V. Obenchain <vobencha@fhcrc.org>}

\seealso{
  \itemize{
    \item The \link{GenomicRanges} and \link{GRanges} classes.
    \item The \link[IRanges]{Ranges} and \link[IRanges]{Hits} classes in the
          IRanges package.
    \item The \link[IRanges]{nearest-methods} man page in the IRanges package.
    \item \link[GenomicRanges]{findOverlaps-methods} for finding just the
          overlapping ranges.
  }
}

\examples{
## ------------------------------------------
## precede() and follow()
## ------------------------------------------
query <- GRanges("A", IRanges(c(5, 20), width=1), strand="+")
subject <- GRanges("A", IRanges(rep(c(10, 15), 2), width=1),
                        strand=c("+", "+", "-", "-"))
precede(query, subject)
follow(query, subject)

strand(query) <- "-"
precede(query, subject)
follow(query, subject)

## ties choose first in order
query <- GRanges("A", IRanges(10, width=1), c("+", "-", "*"))
subject <- GRanges("A", IRanges(c(5, 5, 5, 15, 15, 15), width=1),
                        rep(c("+", "-", "*"), 2))
precede(query, subject)
precede(query, rev(subject))

## ------------------------------------------
## nearest()
## ------------------------------------------
## When multiple ranges overlap an "arbitrary" range is chosen
query <- GRanges("A", IRanges(5, 15))
subject <- GRanges("A", IRanges(c(1, 15), c(5, 19)))
nearest(query, subject)

## select="all" returns all hits
nearest(query, subject, select="all")

## Ranges in 'x' will self-select when 'subject' is present
query <- GRanges("A", IRanges(c(1, 10), width=5))
nearest(query, query)

## Ranges in 'x' will not self-select when 'subject' is missing
nearest(query)
}
 
\keyword{utilities}

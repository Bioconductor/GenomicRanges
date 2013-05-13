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
  The \code{nearest}, \code{precede}, \code{follow}, \code{distance}
  and \code{distanceToNearest} methods for \code{\linkS4class{GenomicRanges}}
  objects and subclasses.
}

\usage{
\S4method{precede}{GenomicRanges,GenomicRanges}(x, subject, select=c("arbitrary", "all"), ignore.strand=FALSE, ...)
\S4method{precede}{GenomicRanges,missing}(x, subject, select=c("arbitrary", "all"), ignore.strand=FALSE, ...)

\S4method{follow}{GenomicRanges,GenomicRanges}(x, subject, select=c("arbitrary", "all"), ignore.strand=FALSE, ...)
\S4method{follow}{GenomicRanges,missing}(x, subject, select=c("arbitrary", "all"), ignore.strand=FALSE, ...)

\S4method{nearest}{GenomicRanges,GenomicRanges}(x, subject, select=c("arbitrary", "all"), ignore.strand=FALSE, ...)
\S4method{nearest}{GenomicRanges,missing}(x, subject, select=c("arbitrary", "all"), ignore.strand=FALSE, ...)

\S4method{distanceToNearest}{GenomicRanges,GenomicRanges}(x, subject, ignore.strand=FALSE, ...)
\S4method{distanceToNearest}{GenomicRanges,missing}(x, subject, ignore.strand=FALSE, ...)

\S4method{distance}{GenomicRanges,GenomicRanges}(x, y, ignore.strand=FALSE, ...)
}

\arguments{
  \item{x}{The query \link{GenomicRanges} instance.
  }
  \item{subject}{The subject \link{GenomicRanges} instance
    within which the nearest neighbors are found. Can be missing,
    in which case \code{x} is also the subject. 
  }
  \item{y}{For the \code{distance} method, a \code{GRanges} instance.
    Cannot be missing. If \code{x} and \code{y} are not the same length, 
    the shortest will be recycled to match the length of the longest.
  }
  \item{select}{Logic for handling ties. By default, all methods
    select a single interval (arbitrary for \code{nearest},
    the first by order in \code{subject} for \code{precede}, and the
    last for \code{follow}). 

    When \code{select="all"} a \link[IRanges]{Hits} object is returned with 
    all matches for \code{x}. If \code{x} does not have a match in 
    \code{subject} the \code{x} is not included in the \code{Hits} object. 
  }
  \item{ignore.strand}{
    A \code{logical} indicating if the strand of the input ranges
    should be ignored. When \code{TRUE}, strand is set to \code{'+'}.
  }
  \item{...}{Additional arguments for methods.}
}

\details{
  \itemize{
    \item{nearest: }{
      Performs conventional nearest neighbor finding.
      Returns an integer vector containing the index of the nearest neighbor 
      range in \code{subject} for each range in \code{x}. If there is no 
      nearest neighbor \code{NA} is returned. For details of the algorithm
      see the man page in IRanges, ?\code{nearest}. 
    }
    \item{precede: }{
      For each range in \code{x}, \code{precede} returns
      the index of the range in \code{subject} that is directly 
      preceded by the range in \code{x}. Overlapping ranges are excluded.
      \code{NA} is returned when there are no qualifying ranges in
      \code{subject}.
    }
    \item{follow: }{
      The opposite of \code{precede}, \code{follow} returns
      the index of the range in \code{subject} that is directly followed by the
      range in \code{x}. Overlapping ranges are excluded. \code{NA} is returned 
      when there are no qualifying ranges in \code{subject}.
    }
    \item{Orientation and Strand: }{
      The relevant orientation for \code{precede} and \code{follow} 
      is 5' to 3', consistent with the direction of translation. 
      Because positional numbering along a chromosome is from left to 
      right and transcription takes place from 5' to 3', \code{precede} and 
      \code{follow} can appear to have `opposite' behavior on the \code{+}
      and \code{-} strand. Using positions 5 and 6 as an example, 5 precedes 
      6 on the \code{+} strand but follows 6 on the \code{-} strand.

      A range with strand \code{*} can be compared to ranges on either the
      \code{+} or \code{-} strand. Below we outline the priority when ranges 
      on multiple strands are compared. When \code{ignore.strand=TRUE} all 
      ranges are treated as if on the \code{+} strand. 
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
    }
    \item{distanceToNearest: }{Returns the distance for each range in \code{x}
      to its nearest neighbor in the \code{subject}. 
    }
    \item{distance: }{
      Returns the distance for each range in \code{x} to the range in \code{y}.
      The behavior of \code{distance} has changed in Bioconductor 2.12.
      See the man page \code{?distance} in IRanges for details.  
    }
  }
}

\value{
  For \code{nearest}, \code{precede} and \code{follow}, an integer
  vector of indices in \code{subject}, or a \code{\linkS4class{Hits}} if 
  \code{select="all"}.

  For \code{distanceToNearest}, a \code{Hits} object with a column for
  the \code{query} index (queryHits), \code{subject} index (subjectHits)
  and the \code{distance} between the pair.

  For \code{distance}, an integer vector of distances between the ranges
  in \code{x} and \code{y}.
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
  ## -----------------------------------------------------------
  ## precede() and follow()
  ## -----------------------------------------------------------
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
 
  ## ignore.strand=TRUE treats all ranges as '+'
  precede(query[1], subject[4:6], select="all", ignore.strand=FALSE)
  precede(query[1], subject[4:6], select="all", ignore.strand=TRUE)


  ## -----------------------------------------------------------
  ## nearest()
  ## -----------------------------------------------------------
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

  ## -----------------------------------------------------------
  ## distance(), distanceToNearest()
  ## -----------------------------------------------------------
  ## Adjacent, overlap, separated by 1
  query <- GRanges("A", IRanges(c(1, 2, 10), c(5, 8, 11)))
  subject <- GRanges("A", IRanges(c(6, 5, 13), c(10, 10, 15)))
  distance(query, subject)

  ## recycling
  distance(query[1], subject)

  ## zero-width ranges
  zw <- GRanges("A", IRanges(4,3))
  stopifnot(distance(zw, GRanges("A", IRanges(3,4))) == 0L)
  sapply(-3:3, function(i) 
      distance(shift(zw, i), GRanges("A", IRanges(4,3))))

  query <- GRanges(c("A", "B"), IRanges(c(1, 5), width=1))
  distanceToNearest(query, subject)

  ## distance() with GRanges and TranscriptDb see ?TranscriptDb
  ## man page in GenomicFeatures package.
}
 
\keyword{utilities}

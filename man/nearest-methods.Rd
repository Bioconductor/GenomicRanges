\name{nearest-methods}

\alias{nearest-methods}

\alias{precede}
\alias{precede,GenomicRanges,GenomicRanges-method}
\alias{precede,GenomicRanges,missing-method}
\alias{follow}
\alias{follow,GenomicRanges,GenomicRanges-method}
\alias{follow,GenomicRanges,missing-method}
\alias{nearest}
\alias{nearest,GenomicRanges,GenomicRanges-method}
\alias{nearest,GenomicRanges,missing-method}
\alias{nearestKNeighbors}
\alias{nearestKNeighbors,GenomicRanges,GenomicRanges-method}
\alias{nearestKNeighbors,GenomicRanges,missing-method}
\alias{distance}
\alias{distance,GenomicRanges,GenomicRanges-method}
\alias{distanceToNearest}
\alias{distanceToNearest,GenomicRanges,GenomicRanges-method}
\alias{distanceToNearest,GenomicRanges,missing-method}

\title{Finding the nearest genomic range/position neighbor}

\description{
  The \code{nearest}, \code{precede}, \code{follow}, \code{distance},
  \code{nearestKNeighbors}, and \code{distanceToNearest} methods for
  \code{\linkS4class{GenomicRanges}} objects and subclasses.
}

\usage{
\S4method{precede}{GenomicRanges,GenomicRanges}(x, subject, 
    select=c("first", "all"), ignore.strand=FALSE)
\S4method{precede}{GenomicRanges,missing}(x, subject, 
    select=c("first", "all"), ignore.strand=FALSE)

\S4method{follow}{GenomicRanges,GenomicRanges}(x, subject, 
    select=c("last", "all"), ignore.strand=FALSE)
\S4method{follow}{GenomicRanges,missing}(x, subject, 
    select=c("last", "all"), ignore.strand=FALSE)

\S4method{nearest}{GenomicRanges,GenomicRanges}(x, subject, 
    select=c("arbitrary", "all"), ignore.strand=FALSE)
\S4method{nearest}{GenomicRanges,missing}(x, subject, 
    select=c("arbitrary", "all"), ignore.strand=FALSE)

\S4method{nearestKNeighbors}{GenomicRanges,GenomicRanges}(x, subject, k=1L,
    select=c("arbitrary", "all"), ignore.strand=FALSE)
\S4method{nearestKNeighbors}{GenomicRanges,missing}(x, subject, k=1L,
    select=c("arbitrary", "all"), ignore.strand=FALSE)

\S4method{distanceToNearest}{GenomicRanges,GenomicRanges}(x, subject, 
    ignore.strand=FALSE, ...)
\S4method{distanceToNearest}{GenomicRanges,missing}(x, subject, 
    ignore.strand=FALSE, ...)

\S4method{distance}{GenomicRanges,GenomicRanges}(x, y, 
    ignore.strand=FALSE, ...)
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
  \item{k}{For the \code{nearestKNeighbors} method, an \code{integer} declaring how 
    many nearest neighbors to find.
  }
  \item{select}{Logic for handling ties. By default, all methods
    select a single interval (arbitrary for \code{nearest},
    the first by order in \code{subject} for \code{precede}, and the
    last for \code{follow}). 

    When \code{select="all"} a \link[S4Vectors]{Hits} object is returned with 
    all matches for \code{x}.
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
      see the man page in the \pkg{IRanges} package (\code{?nearest}). 
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
    \item{nearestKNeighbors: }{
      Performs conventional k-nearest neighbor finding.
      Returns an \link[IRanges]{IntegerList} containing the index of the
      k-nearest neighbors in \code{subject} for each range in \code{x}. If there
      is no nearest neighbor \code{NA} is returned. If \code{select="all"} is
      specified, ties will be included in the resulting 
      \link[IRanges]{IntegerList}.
    }
    \item{Orientation and strand for \code{precede} and \code{follow}: }{
      Orientation is 5' to 3', consistent with the direction of translation. 
      Because positional numbering along a chromosome is from left to 
      right and transcription takes place from 5' to 3', \code{precede} and 
      \code{follow} can appear to have `opposite' behavior on the \code{+}
      and \code{-} strand. Using positions 5 and 6 as an example, 5 precedes 
      6 on the \code{+} strand but follows 6 on the \code{-} strand.

      The table below outlines the orientation when ranges on different 
      strands are compared. In general, a feature on \code{*} is considered
      to belong to both strands. The single exception is when both \code{x} 
      and \code{subject} are \code{*} in which case both are treated as \code{+}.
\preformatted{
       x  |  subject  |  orientation 
     -----+-----------+----------------
a)     +  |  +        |  ---> 
b)     +  |  -        |  NA
c)     +  |  *        |  --->
d)     -  |  +        |  NA
e)     -  |  -        |  <---
f)     -  |  *        |  <---
g)     *  |  +        |  --->
h)     *  |  -        |  <---
i)     *  |  *        |  --->  (the only situation where * arbitrarily means +)
}

    }
    \item{distanceToNearest: }{Returns the distance for each range in \code{x}
      to its nearest neighbor in the \code{subject}. 
    }
    \item{distance: }{
      Returns the distance for each range in \code{x} to the range in \code{y}.
      The behavior of \code{distance} has changed in Bioconductor 2.12.
      See the man page \code{?distance} in the \pkg{IRanges} package for
      details.  
    }
  }
}

\value{
  For \code{nearest}, \code{precede} and \code{follow}, an integer
  vector of indices in \code{subject}, or a \link[S4Vectors]{Hits} if 
  \code{select="all"}.

  For \code{nearestKNeighbors}, an \link[IRanges]{IntegerList} of vertices in
  \code{subject}.

  For \code{distanceToNearest}, a \link[S4Vectors]{Hits} object with a
  column for the \code{query} index (queryHits), \code{subject} index
  (subjectHits) and the \code{distance} between the pair.

  For \code{distance}, an integer vector of distances between the ranges
  in \code{x} and \code{y}.
}

\author{P. Aboyoun and V. Obenchain}

\seealso{
  \itemize{
    \item The \link{GenomicRanges} and \link{GRanges} classes.
    \item The \link[IRanges]{IntegerRanges} class in the \pkg{IRanges} package.
    \item The \link[S4Vectors]{Hits} class in the \pkg{S4Vectors} package.
    \item The \link[IRanges]{nearest-methods} man page in the \pkg{IRanges}
          package.
    \item \link[GenomicRanges]{findOverlaps-methods} for finding just the
          overlapping ranges.
    \item The \link[GenomicFeatures]{nearest-methods} man page in the 
          \pkg{GenomicFeatures} package.
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
## nearestKNeighbors()
## -----------------------------------------------------------
## Without an argument, k defaults to 1
query <- GRanges("A", IRanges(c(2, 5), c(8, 15)))
subject <- GRanges("A", IRanges(c(1, 4, 10, 15), c(5, 7, 12, 19)))
nearestKNeighbors(query, subject)

## Return multiple neighbors with k > 1
nearestKNeighbors(query, subject, k=3)

## select="all" returns all hits
nearestKNeighbors(query, subject, select="all")

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

## distance() with GRanges and TxDb see the 
## ?'distance,GenomicRanges,TxDb-method' man 
## page in the GenomicFeatures package.
}
 
\keyword{utilities}

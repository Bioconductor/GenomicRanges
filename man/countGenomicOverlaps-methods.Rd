\name{countGenomicOverlaps}

\alias{countGenomicOverlaps}
\alias{countGenomicOverlaps,GRangesList,GRangesList-method}
\alias{countGenomicOverlaps,GenomicRanges,GenomicRanges-method}
\alias{countGenomicOverlaps,GenomicRanges,GRangesList-method}
\alias{countGenomicOverlaps,GRangesList,GenomicRanges-method}
\alias{countGenomicOverlaps,GappedAlignments,GenomicRanges-method}
\alias{countGenomicOverlaps,GappedAlignments,GRangesList-method}


\title{Count Read Hits in Genomic Features} 

\description{
  Count read hits per exon or transcript and resolve multi-hit reads. 
}

\usage{
  \S4method{countGenomicOverlaps}{GRangesList,GRangesList}(
    query, subject, 
    type = c("any", "start", "end", "within", "equal"),
    resolution = c("none", "divide", "uniqueDisjoint"),
    ignore.strand = FALSE, splitreads = TRUE, ...) 
}

\arguments{
  \item{query}{
    A \link{GRangesList}, \link{GRanges}, or \link{GappedAlignments}
    object. The \code{query} is intended to be a
    \code{\link{GRangesList}} where each list element is a read. When
    the length of a list element is greater than 1, it is assumed to
    represent a spliced read.  Both \link{GRanges} and
    \link{GappedAlignments} are coerced to a \link{GRangesList} with
    each row as a single list element. If the cigar in the
    \link{GappedAlignments} has gaps, the read will be represented in
    the corresponding \link{GRangesList} object as a split read (i.e.,
    list element of length 2 or greater).
  }
  \item{subject}{
    A \code{\link{GRangesList}}, or a 
    \link{GRanges} object. 
    The subject is expected to be a list of genomic features, 
    specifically genes, with each row representing a feature
    where the feature could be an exon or transcript. If a 
    \link{GRanges}
    object is provided, it will be coerced to a 
    \link{GRangesList} 
    with each original range representing as a single list element
    (i.e., as a single gene). 
  }
  \item{type}{
    See \code{\link[IRanges]{findOverlaps}} in the IRanges package for
    a description of this argument.
  }
  \item{resolution}{
    A \code{character(1)} string of "none", "divide", or "uniqueDisjoint".
    These rule sets are used to distribute read hits when 
    multiple subjects are hit by the same query.

    \itemize{
      \item "none" : No conflict resolution is performed. All queries that
            hit more than 1 subject are dropped.
      \item "divide" : The hit from a single query is divided equally among 
            all subjects that were hit. If a query hit 4 subjects 
            each subject is assigned 1/4 of a hit.
      \item "uniqueDisjoint" : Subjects hit by a common query are 
            partitioned into disjoint intervals. Any regions that are shared 
            between the subjects are discarded. If the read overlaps one of 
            these remaining unique disjoint regions the hit is assigned to 
            that feature. If the read overlaps both or none of the regions, 
            no hit is assigned. Therefore, unlike the \code{divide} option, 
            \code{uniqueDisjoint} does not resolve multi-hit conflict in 
            all situations.
   }
  }
  \item{ignore.strand}{
    A logical value indicating if strand should be considered when matching.
  }
  \item{splitreads}{
    A logical value indicating if split reads should be included.
  }
  \item{...}{Additional arguments, perhaps used by methods defined on
  this generic.
  }
}

\details{
  The \code{countGenomicOverlaps} methods use the \code{findOverlaps} 
  function in conjunction with a resolution method to identify overlaps
  and resolve queries that match multiple subjects.
  The usual \code{type} argument of \code{findOverlaps}
  is used to specify the type of overlap. The \code{resolution}
  argument is used to select a method to resolve the conflict
  when a query hits more than 1 subject. Here the
  term `hit' means an overlap identified by \code{findOverlaps}.
 
  The primary difference in the handling of split reads vs
  simple reads (i.e., no gap in the CIGAR) is the portion of
  the read hit each split read fragment has to contribute.
  All reads, whether simple or split, have an overall value
  of 1 to contribute to a subject they hit. In the case of the
  split reads, this value is further divided by the number of
  fragments in the read. For example, if a split read has 3
  fragments (i.e., two gaps in the CIGAR) each
  fragment has a value of 1/3 to contribute to the subject
  they hit. As with the simple reads, depending upon the
  \code{resolution} chosen the value may be divided, fully 
  assigned or discarded.

  More detailed examples can be found in the \code{countGenomicOverlaps}
  vignette. 

}

\value{
  A \link{GRangesList} object with an
  additional metadata column specifying the number of hits.
}

\author{Valerie Obenchain and Martin Morgan}

\examples{
rng1 <- function(s, w)
GRanges(seq="chr1", IRanges(s, width=w), strand="+")

rng2 <- function(s, w)
GRanges(seq="chr2", IRanges(s, width=w), strand="+")

subj <- GRangesList(A=rng1(1000, 500),
                    B=rng2(2000, 900),
                    C=rng1(c(3000, 3600), c(500, 300)),
                    D=rng2(c(7000, 7500), c(600, 300)),
                    E1=rng1(4000, 500), E2=rng1(c(4300, 4500), c(400, 400)),
                    F=rng2(3000, 500),
                    G=rng1(c(5000, 5600), c(500, 300)),
                    H1=rng1(6000, 500), H2=rng1(6600, 400))

query <- GRangesList(a=rng1(1400, 500),
                     b=rng2(2700, 100),
                     c=rng1(3400, 300),
                     d=rng2(7100, 600),
                     e=rng1(4200, 500),
                     f=rng2(c(3100, 3300), 50),
                     g=rng1(c(5400, 5600), 50),
                     h=rng1(c(6400, 6600), 50))

## Overlap type = "any"
none <- countGenomicOverlaps(query, subj, type="any", resolution="none")
divide <- countGenomicOverlaps(query, subj, type="any", resolution="divide")
uniqueDisjoint <- countGenomicOverlaps(query, subj, type="any", 
    resolution="uniqueDisjoint")
data.frame(none = values(unlist(none))[["hits"]], 
           divide = values(unlist(divide))[["hits"]],
           uniqDisj = values(unlist(uniqueDisjoint))[["hits"]])

## Split read with 4 fragments where
## 2 fragments hit the same subject, 
## 1 fragment hits two different subjects,
## 1 fragment hits a single subject 
splitreads <- GRangesList(c(rng1(c(3000, 3200, 4000), 100), rng1(5400, 300)))
split_none <- countGenomicOverlaps(splitreads, subj, type="any")
split_divide <- countGenomicOverlaps(splitreads, subj, type="any", 
    resolution="divide")

data.frame(none = values(unlist(split_none))[["hits"]],
           divide = values(unlist(split_divide))[["hits"]])

}

\keyword{methods}
\keyword{utilities}


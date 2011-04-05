\name{countGenomicOverlaps}

\alias{countGenomicOverlaps,GRangesList,GRangesList-method}
\alias{countGenomicOverlaps,GenomicRanges,GenomicRanges-method}
\alias{countGenomicOverlaps,GenomicRanges,GRangesList-method}
\alias{countGenomicOverlaps,GRangesList,GenomicRanges-method}
\alias{countGenomicOverlaps,GappedAlignments,GenomicRanges-method}
\alias{countGenomicOverlaps,GappedAlignments,GRangesList-method}


\title{Count Read Hits in Genomic Features} 

\description{
  Count read hits per exon and resolve multi-hit reads. 
}

\usage{
  \S4method{countGenomicOverlaps}{GRangesList,GRangesList}(
    query, subject, 
    type = c("any", "start", "end", "within", "equal"),
    resolution = c("none", "divide", "uniqueDisjoint"),
    junctionwt = 0.5, ignore.strand = FALSE, ...) 
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
    list element of length 2).
  }
  \item{subject}{
    A \code{\link{GRangesList}}, or a 
    \link{GRanges} object. 
    The subject is expected to be a list of genomic features, 
    specifically genes, with each row representing an exon. If a 
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
    A \code{character(1)} string of "none", "divide", "uniqueDisjoint", or "coverage".
    These rule sets are used to distribute read hits when 
    multiple subjects are hit by the same query.

    \itemize{
      \item "none" : No conflict resolution is peformed. All queries that
            hit more than 1 subject are dropped. No split reads are counted.
      \item "divide" : The hit from a single query is divided equally amoung 
            all subjects that were hit. If a query hit 4 subjects 
            each subject is assigned 1/4. 
      \item "uniqueDisjoint" : Subjects hit by a common query are 
            partitioned into disjoint intervals. Any regions that are shared 
            between the subjects are discarded. If the read overlaps one of 
            these remaining unique disjoint regions the hit is assigned to 
            that exon. If the read overlaps both or none of the regions, 
            no hit is assigned. Therefore, unlike the \code{'divide'} option, 
            \code{'uniqueDisjoint'} does not resolve multi-hit conflict in 
            all situations.
   }
  }
  \item{junctionwt}{
    A numeric value representing the hit weight for junction reads. A junction
    read is defined as a split read query that hits 2 exons either in the
    same gene or in different genes. The default is 1. 
  }
  \item{ignore.strand}{
    A logical value indicating if strand should be considered when matching.
  }
  \item{...}{Additional arguments, perhaps used by methods defined on
  this generic.
  }
}

\details{
  The \code{countGenomicOverlaps} methods use the \code{findOverlaps} 
  function in conjunction with defined rule sets to identify overlaps
  and resolve queries that match multiple subjects.
  The usual \code{type} argument of \code{findOverlaps}
  is used to specify the type of overlap. the \code{resolve}
  argument is used to select a method to resolve the conflict
  when a query hits more than 1 subject. Here the
  term `hit` means an overlap identified by \code{findOverlaps}.
 
  Split reads are handled separately from the single and multiple 
  hit reads. If both portions of a split read hit the same subject, the
  subject is assigned the hit. When 2 different subjects are hit
  each is assigned the default \code{junctionwt} value of 0.5. 
  The user can modify the value assigned to such junction reads with the
  \code{junctionwt} argument. 

}

\value{
  A \link{GRangesList} object with an
  additional metadata column specifying the number of hits.
}

\author{Valerie Obenchain and Martin Morgan}

\examples{
rng <- function(s, w)
GRanges(seq="chr1", IRanges(s, width=w), strand="+")

rng2 <- function(s, w)
GRanges(seq="chr2", IRanges(s, width=w), strand="+")

subj <- GRangesList(A=rng(1000, 500),
                    B=rng2(2000, 900),
                    C=rng(c(3000, 3600), c(500, 300)),
                    D=rng2(c(7000, 7500), c(600, 300)),
                    E1=rng(4000, 500), E2=rng(c(4300, 4500), c(400, 400)),
                    F=rng2(3000, 500),
                    G=rng(c(5000, 5600), c(500, 300)),
                    H1=rng(6000, 500), H2=rng(6600, 400))

query <- GRangesList(a=rng(1400, 500),
                     b=rng2(2700, 100),
                     c=rng(3400, 300),
                     d=rng2(7100, 600),
                     e=rng(4200, 500),
                     f=rng2(c(3100, 3300), 50),
                     g=rng(c(5400, 5600), 50),
                     h=rng(c(6400, 6600), 50))

## overlap type="any"
a1 <- GenomicRanges:::countGenomicOverlaps(query, subj, type="any", resolution="none")
a2 <- GenomicRanges:::countGenomicOverlaps(query, subj, type="any", resolution="none",
                           junctionwt=5)
a3 <- GenomicRanges:::countGenomicOverlaps(query, subj, type="any", resolution="divide")
a4 <- GenomicRanges:::countGenomicOverlaps(query, subj, type="any", resolution="uniqueDisjoint")

res_any <- data.frame(none = values(unlist(a1))[["hits"]],
                      noneWithjctwt = values(unlist(a2))[["hits"]],
                      divide = values(unlist(a3))[["hits"]],
                      uniqDisj = values(unlist(a4))[["hits"]])

}

\keyword{methods}
\keyword{utilities}


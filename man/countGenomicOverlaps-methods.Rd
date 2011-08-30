\name{countGenomicOverlaps}

\alias{countGenomicOverlaps}
\alias{countGenomicOverlaps,GRangesList,GRangesList-method}
\alias{countGenomicOverlaps,GenomicRanges,GenomicRanges-method}
\alias{countGenomicOverlaps,GenomicRanges,GRangesList-method}
\alias{countGenomicOverlaps,GRangesList,GenomicRanges-method}
\alias{countGenomicOverlaps,GenomicRanges,GappedAlignments-method}
\alias{countGenomicOverlaps,GRangesList,GappedAlignments-method}


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
    A \link{GRangesList}, or a \link{GRanges} of genomic features.
    These are the annotations that define the genomic regions and will often 
    be the result of calling "exonsBy" or "transcriptsBy" on a 
    \link[GenomicFeatures]{TranscriptDb} object. If a \link{GRangesList} is 
    provided, each top level of the list represents a "super" such as a gene 
    and each row is a "sub" such as an exon or transcript.  When \code{query}
    is a \link{GRanges} all rows are considered to be of the same level
    (e.g., all genes, all exons or all transcripts).
  }
  \item{subject}{
    A \link{GRangesList}, \link{GRanges}, or \link{GappedAlignments}
    representing the data (e.g., reads). 

    List structures as the \code{subject} are used to represent reads with multiple parts
    (i.e., gaps in the CIGAR). When a \link{GappedAlignments} is provided 
    it is coerced to a \link{GRangesList} object. If any of the reads in the
    \link{GappedAlignments} have gaps, the corresponding 
    \link{GRangesList} will have multiple elements for that top level list. When
    \code{subject} is a \link{GRanges}, it is assumed that all reads are
    simple and do not have multiple parts. 
  }
  \item{type}{
    See \code{\link[IRanges]{findOverlaps}} in the IRanges package for
    a description of this argument.
  }
  \item{resolution}{
    A \code{character(1)} string of "none", "divide", or "uniqueDisjoint".
    These rule sets are used to distribute read hits when 
    multiple queries are hit by the same subject.

    \itemize{
      \item "none" : No conflict resolution is performed. All subjects that
            hit more than 1 query are dropped.
      \item "divide" : The hit from a single subject is divided equally among 
            all queries that were hit. If a subject hit 4 queries 
            each query is assigned 1/4 of a hit.
      \item "uniqueDisjoint" : Queries hit by a common subject are 
            partitioned into disjoint intervals. Any regions that are shared 
            between the queries are discarded. If the read overlaps one of 
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
  and resolve subjects (reads) that match multiple queries (annotation regions).
  The usual \code{type} argument of \code{findOverlaps}
  is used to specify the type of overlap. The \code{resolution}
  argument is used to select a method to resolve the conflict
  when a subject hits more than 1 query. Here the
  term `hit' means an overlap identified by \code{findOverlaps}.
 
  The primary difference in the handling of split reads vs
  simple reads (i.e., no gap in the CIGAR) is the portion of
  the read hit each split read fragment has to contribute.
  All reads, whether simple or split, have an overall value
  of 1 to contribute to a query they hit. In the case of the
  split reads, this value is further divided by the number of
  fragments in the read. For example, if a split read has 3
  fragments (i.e., two gaps in the CIGAR) each
  fragment has a value of 1/3 to contribute to the query 
  they hit. As with the simple reads, depending upon the
  \code{resolution} chosen the value may be divided, fully 
  assigned or discarded.

  More detailed examples can be found in the \code{countGenomicOverlaps}
  vignette. 

}

\value{
  A \link{SummarizedExperiment} object. The counts are held in the \code{assays}
  slot, the \code{subject} is in the \code{rowData} slot and the \code{colData}
  slot will either be \code{NULL} or hold any metadata that was present in the
  \code{query}.
}

\author{Valerie Obenchain and Martin Morgan}

\examples{
rng1 <- function(s, w)
GRanges(seq="chr1", IRanges(s, width=w), strand="+")

rng2 <- function(s, w)
GRanges(seq="chr2", IRanges(s, width=w), strand="+")

query <- GRangesList(A=rng1(1000, 500),
                    B=rng2(2000, 900),
                    C=rng1(c(3000, 3600), c(500, 300)),
                    D=rng2(c(7000, 7500), c(600, 300)),
                    E1=rng1(4000, 500), E2=rng1(c(4300, 4500), c(400, 400)),
                    F=rng2(3000, 500),
                    G=rng1(c(5000, 5600), c(500, 300)),
                    H1=rng1(6000, 500), H2=rng1(6600, 400))

subj <- GRangesList(a=rng1(1400, 500),
                     b=rng2(2700, 100),
                     c=rng1(3400, 300),
                     d=rng2(7100, 600),
                     e=rng1(4200, 500),
                     f=rng2(c(3100, 3300), 50),
                     g=rng1(c(5400, 5600), 50),
                     h=rng1(c(6400, 6600), 50))

## Overlap type = "any"
none <- countGenomicOverlaps(query, subj, 
                             type="any", resolution="none")
divide <- countGenomicOverlaps(query, subj, 
                               type="any", resolution="divide")
uniqueDisjoint <- countGenomicOverlaps(query, subj, type="any", 
                                       resolution="uniqueDisjoint")
data.frame(none = assays(none)$counts, 
           divide = assays(divide)$counts,
           uniqDisj = assays(uniqueDisjoint)$counts)


## Split read with 4 fragments :
splitreads <- GRangesList(c(rng1(c(3000, 3200, 4000), 100), rng1(5400, 300)))
## Unlist both the splitreads and the query to see 
## - read fragments 1 and 2 both hit query 3
## - read fragment 3 hits query 7
## - read fragment 4 hits query 11 and 12 
findOverlaps(unlist(query), unlist(splitreads))

## Use countGenomicOverlaps to avoid double counting.
## Because this read has 4 parts each part contributes a count of 1/4.
## When resolution="none" only reads that hit a single region are counted. 
split_none <- countGenomicOverlaps(query, splitreads, type="any",
                                   resolution="none")
## When resolution="divide" all reads are counted by dividing their count 
## evenly between the regions they hit. Region 3 of the query was hit
## by two reads each contributing a count of 1/4. Region 7 was hit
## by one read contributing a count of 1/4. Regions 11 and 12 were both
## hit by the same read resulting in having to share (i.e., "divide") the 
## single 1/4 hit read 4 had to contribute.
split_divide <- countGenomicOverlaps(query, splitreads, 
                                     type="any", resolution="divide")

data.frame(none = assays(split_none)$counts,
           divide = assays(split_divide)$counts)

}

\keyword{methods}
\keyword{utilities}


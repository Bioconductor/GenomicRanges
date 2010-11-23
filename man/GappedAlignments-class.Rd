\name{GappedAlignments-class}
\docType{class}

% Class:
\alias{class:GappedAlignments}
\alias{GappedAlignments-class}
\alias{GappedAlignments}

% Constructors:
\alias{readGappedAlignments}

% Accessors:
\alias{length,GappedAlignments-method}
\alias{rname}
\alias{rname,GappedAlignments-method}
\alias{seqnames,GappedAlignments-method}
\alias{rname<-}
\alias{rname<-,GappedAlignments-method}
\alias{seqnames<-,GappedAlignments-method}
\alias{strand,GappedAlignments-method}
\alias{cigar}
\alias{cigar,GappedAlignments-method}
\alias{qwidth}
\alias{qwidth,GappedAlignments-method}
\alias{grglist}
\alias{grglist,GappedAlignments-method}
\alias{granges}
\alias{granges,GappedAlignments-method}
\alias{grg}
\alias{grg,GappedAlignments-method}
\alias{rglist}
\alias{rglist,GappedAlignments-method}
\alias{ranges,GappedAlignments-method}
\alias{start,GappedAlignments-method}
\alias{end,GappedAlignments-method}
\alias{width,GappedAlignments-method}
\alias{ngap,GappedAlignments-method}
\alias{seqinfo,GappedAlignments-method}
\alias{seqinfo<-,GappedAlignments-method}
\alias{seqlevels<-,GappedAlignments-method}

% Coercion:
\alias{as.data.frame,GappedAlignments-method}
\alias{coerce,GappedAlignments,GRangesList-method}
\alias{coerce,GappedAlignments,GRanges-method}
\alias{coerce,GappedAlignments,RangesList-method}
\alias{coerce,GappedAlignments,Ranges-method}

% Subsetting:
\alias{[,GappedAlignments-method}

% "show" method:
\alias{show,GappedAlignments-method}

% Other methods:
\alias{qnarrow}
\alias{qnarrow,GappedAlignments-method}
\alias{narrow,GappedAlignments-method}
\alias{pintersect,GappedAlignments,GRanges-method}
\alias{pintersect,GRanges,GappedAlignments-method}
\alias{coverage,GappedAlignments-method}
\alias{findOverlaps,GappedAlignments,ANY-method}
\alias{findOverlaps,ANY,GappedAlignments-method}
\alias{findOverlaps,GappedAlignments,GappedAlignments-method}
\alias{countOverlaps,GappedAlignments,ANY-method}
\alias{countOverlaps,ANY,GappedAlignments-method}
\alias{countOverlaps,GappedAlignments,GappedAlignments-method}
\alias{subsetByOverlaps,GappedAlignments,ANY-method}
\alias{subsetByOverlaps,ANY,GappedAlignments-method}
\alias{subsetByOverlaps,GappedAlignments,GappedAlignments-method}
\alias{match,GappedAlignments,ANY-method}
\alias{match,ANY,GappedAlignments-method}
\alias{match,GappedAlignments,GappedAlignments-method}
\alias{\%in\%,GappedAlignments,ANY-method}
\alias{\%in\%,ANY,GappedAlignments-method}
\alias{\%in\%,GappedAlignments,GappedAlignments-method}


\title{GappedAlignments objects}

\description{
  The GappedAlignments class is a simple container which
  purpose is to store a set of alignments that will hold just enough
  information for supporting the operations described below.

  WARNING! This is work-in-progress. Expect frequent changes in
  functionalities.
}

\details{
  A GappedAlignments object is a vector-like object where each element
  describes an alignment i.e. how a given sequence (called "query"
  or "read", typically short) aligns to a reference sequence (typically
  long).

  Most of the time, a GappedAlignments object will be created by loading
  records from a BAM (or SAM) file and each element in the resulting
  object will correspond to a record. BAM/SAM records generally contain
  a lot of information but only part of that information is loaded
  in the GappedAlignments object. In particular, we discard the query
  sequences (SEQ field), the query ids (QNAME field), the query qualities
  (QUAL), the mapping qualities (MAPQ) and any other information that
  is not needed in order to support the operations or methods described
  below.

  This means that multi-reads (i.e. reads with multiple hits in the
  reference) won't receive any special treatment i.e. the various SAM/BAM
  records corresponding to a multi-read will show up in the GappedAlignments
  object as if they were coming from different/unrelated queries.
  Also paired-end reads will be treated as single-end reads and the
  pairing information will be lost.

  Each element of a GappedAlignments object consists of:
  \itemize{
    \item The name of the reference sequence. (This is the RNAME field
          in a SAM/BAM record.)
    \item The strand in the reference sequence to which the query is
          aligned. (This information is stored in the FLAG field in a
          SAM/BAM record.)
    \item The CIGAR string in the "Extended CIGAR format" (see the SAM
          Format Specifications for the details).
    \item The 1-based leftmost position/coordinate of the clipped query
          relative to the reference sequence. We will refer to it as
          the "start" of the query. (This is the POS field in a SAM/BAM
          record.)
    \item The 1-based rightmost position/coordinate of the clipped query
          relative to the reference sequence. We will refer to it as
          the "end" of the query. (This is NOT explicitly stored in a
          SAM/BAM record but can be inferred from the POS and CIGAR fields.)
          Note that all positions/coordinates are always relative to
          the first base at the 5' end of the plus strand of the reference
          sequence, even when the query is aligned to the minus strand.
    \item The genomic intervals between the "start" and "end" of the query
          that are "covered" by the alignment. Saying that the full
          [start,end] interval is covered is the same as saying that the
          alignment has no gap (no N in the CIGAR). It is then considered
          a simple alignment. Note that a simple alignment can have
          mismatches or deletions (in the reference). In other words, a
          deletion, encoded with a D, is NOT considered a gap.
  }

  Note that the last 2 items are not expicitly stored in the GappedAlignments
  object: they are inferred on-the-fly from the CIGAR and the "start".

  The rest of this man page will focus on describing how to:
  \itemize{
    \item Access the information stored in a GappedAlignments object
          in a way that is independent from how the data are actually
          stored internally.
    \item How to create and manipulate a GappedAlignments object.
  }
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{readGappedAlignments(file, format="BAM", ...)}:
      Read a file as a GappedAlignments object.
      The function is just a front-end that delegates to a format-specific
      back-end function (any extra argument is passed to the back-end
      function).
      Only the BAM format is supported for now. Its back-end is the
      \code{\link[Rsamtools]{readBamGappedAlignments}} function defined
      in the Rsamtools package.
      See \code{?\link[Rsamtools]{readBamGappedAlignments}} for
      more information (you might need to install and load the package
      first).
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{x} is a GappedAlignments object.

  \describe{
    \item{}{
      \code{length(x)}:
      Returns the number of alignments in \code{x}.
    }
    \item{}{
      \code{rname(x)}, \code{rname(x) <- value}:
      Gets or sets the name of the reference sequences for each alignment
      in \code{x}.
      \code{value} can be an \link[IRanges]{Rle} object, a character vector,
      or a factor.
    }
    \item{}{
      \code{seqnames(x)}, \code{seqnames(x) <- value}:
      Same as \code{rname(x)} and \code{rname(x) <- value}.
    }
    \item{}{
      \code{strand(x)}:
      Returns a character factor of length \code{length(x)}
      (with levels +, - and *) containing the strand in the reference
      sequence to which the query is aligned.
    }
    \item{}{
      \code{cigar(x)}:
      Returns a character vector of length \code{length(x)}
      containing the CIGAR string for each alignment.
    }
    \item{}{
      \code{qwidth(x)}:
      Returns an integer vector of length \code{length(x)}
      containing the length of the query *after* hard clipping
      (i.e. the length of the query sequence that is stored in
      the corresponding SAM/BAM record).
    }
    \item{}{
      \code{grglist(x)}, \code{granges(x)}, \code{rglist(x)}, \code{ranges(x)}:
      Returns either a \link{GRangesList} object, or a \link{GRanges} object,
      or a \link[IRanges]{RangesList} object, or a \link[IRanges]{Ranges}
      object of length \code{length(x)} where each element represents the
      regions in the reference to which a query is aligned.
      See Details section above for more information.
      More precisely, the \link[IRanges]{RangesList} object returned by
      \code{rglist(x)} is a \link[IRanges]{CompressedNormalIRangesList} object,
      and the \link[IRanges]{Ranges} object returned by \code{ranges(x)} is
      an \link[IRanges]{IRanges} object.
    }
    \item{}{
      \code{start(x)}, \code{end(x)}:
      Returns an integer vector of length \code{length(x)}
      containing the "start" and "end" (respectively) of the query
      for each alignment. See Details section above for the exact
      definitions of the "start" and "end" of a query.
      Note that \code{start(x)} and \code{end(x)} are equivalent
      to \code{start(granges(x))} and \code{end(granges(x))},
      respectively (or, alternatively, to \code{min(rglist(x))} and
      \code{max(rglist(x))}, respectively).
    }
    \item{}{
      \code{width(x)}:
      Equivalent to \code{width(granges(x))} (or, alternatively, to
      \code{end(x) - start(x) + 1L}).
      Note that this is generally different from \code{qwidth(x)}
      except for alignments with a trivial CIGAR string (i.e. a
      string of the form \code{"<n>M"} where <n> is a number).
    }
    \item{}{
      \code{ngap(x)}:
      Returns an integer vector of length \code{length(x)}
      containing the number of gaps for each alignment.
      Equivalent to \code{elementLengths(rglist(x)) - 1L}.
    }
    \item{}{
      \code{seqinfo(x)}, \code{seqinfo(x) <- value}:
      Gets or sets the information about the underlying sequences.
      \code{value} must be a \link{Seqinfo} object.
    }
    \item{}{
      \code{seqlevels(x)}, \code{seqlevels(x) <- value}:
      Gets or sets the sequence levels.
      \code{seqlevels(x)} is equivalent to \code{seqlevels(seqinfo(x))}
      or to \code{levels(rname(x))}, those 2 expressions being
      guaranteed to return indentical character vectors on a GappedAlignments
      object. \code{value} must be a character vector with no NAs.
    }
    \item{}{
      \code{seqlengths(x)}, \code{seqlengths(x) <- value}:
      Gets or sets the sequence lengths.
      \code{seqlengths(x)} is equivalent to \code{seqlengths(seqinfo(x))}.
      \code{value} can be a named non-negative integer or numeric vector
      eventually with NAs.
    }
    \item{}{
      \code{isCircular(x)}, \code{isCircular(x) <- value}:
      Gets or sets the circularity flags.
      \code{isCircular(x)} is equivalent to \code{isCircular(seqinfo(x))}.
      \code{value} must be a named logical vector eventually with NAs.
    }
  }
}

\section{Subsetting and related operations}{
  In the code snippets below, \code{x} is a GappedAlignments object.

  \describe{
    \item{}{
      \code{x[i]}:
      Returns a new GappedAlignments object made of the selected
      alignments. \code{i} can be a numeric or logical vector.
    }
  }
}

\section{Other methods}{

  \describe{
    \item{}{
      \code{qnarrow(x, start=NA, end=NA, width=NA)}:
      \code{x} is a GappedAlignments object.
      Returns a new GappedAlignments object of the same length as \code{x}
      describing how the narrowed query sequences align to the reference.
      The \code{start}/\code{end}/\code{width} arguments describe how
      to narrow the query sequences. They must be vectors of integers.
      NAs and negative values are accepted and "solved" according to the
      rules of the SEW (Start/End/Width) interface (see
      \code{?\link[IRanges]{solveUserSEW}} for the details).
    }
    \item{}{
      \code{narrow(x, start=NA, end=NA, width=NA)}:
      \code{x} is a GappedAlignments object.
      Returns a new GappedAlignments object of the same length as \code{x}
      describing the narrowed alignments. Unlike with \code{qnarrow}
      now the \code{start}/\code{end}/\code{width} arguments describe
      the narrowing on the reference side, not the query side.
      Like with \code{qnarrow}, they must be vectors of integers.
      NAs and negative values are accepted and "solved" according to the
      rules of the SEW (Start/End/Width) interface (see
      \code{?\link[IRanges]{solveUserSEW}} for the details).
    }
    \item{}{
      \code{pintersect(x, y)}:
      Either \code{x} is a GappedAlignments object and \code{y} is a
      GRanges object or \code{x} is a GRanges object and \code{y} is a
      GappedAlignments object.
      Returns a new GappedAlignments object of the same length as the 
      GappedAlignments input arguments. Like with \code{narrow}, the
      resulting "parallel" intersection is with respect to the reference.
    }
    \item{}{
      \code{coverage(x, shift=0L, width=NULL, weight=1L)}:
      \code{x} is a GappedAlignments object.
      Returns a named \link[IRanges]{RleList} object with one element
      (integer-Rle) per unique reference sequence. Each element represents
      \code{x}'s coverage of the corresponding reference sequence, that is,
      how many times each nucleotide position in the sequence is covered
      by the alignments in \code{x}.
      Note that the semantic of the \code{coverage} method for
      GappedAlignments objects is different from the semantic of the method
      for \link[IRanges]{Ranges} objects (the latter returns a single
      integer-Rle object representing the coverage of all ranges
      relatively to a unique imaginary reference sequence).
    }
    \item{}{
      \code{findOverlaps(query, subject, ...)},
      \code{countOverlaps(query, subject, ...)},
      \code{subsetByOverlaps(query, subject, ...)},
      \code{match(x, table, nomatch=NA_integer_, incomparables=NULL)},
      \code{x \%in\% table}:
      \code{query} or \code{subject} or both are GappedAlignments objects.
      \code{findOverlaps(query, subject, ...)} is equivalent to
      \code{findOverlaps(grglist(query), subject, ...)} when
      \code{query} is a GappedAlignments object, or to
      \code{findOverlaps(query, grglist(subject), ...)} when
      \code{subject} is a GappedAlignments object, or to
      \code{findOverlaps(grglist(query), grglist(subject), ...)} when
      both are GappedAlignments objects.
      The same apply to \code{countOverlaps(query, subject, ...)} and
      \code{subsetByOverlaps(query, subject, ...)}.
      See \code{?`\link{findOverlaps,GRangesList,GenomicRanges-method}`},
      \code{?`\link{countOverlaps,GRangesList,GenomicRanges-method}`} and
      \code{?`\link{subsetByOverlaps,GRangesList,GenomicRanges-method}`}
      for more information (in particular for descriptions of the
      extra arguments and the returned object).
    }
  }
}

\references{
  \url{http://samtools.sourceforge.net/}
}

\author{
  H. Pages and P. Aboyoun
}

\seealso{
  \code{\link[Rsamtools]{readBamGappedAlignments}},
  \link{GRangesList-class},
  \link{GRanges-class},
  \link[IRanges]{CompressedNormalIRangesList-class},
  \link[IRanges]{IRanges-class},
  \code{\link{seqinfo}},
  \link{Seqinfo-class},
  \code{\link[IRanges]{coverage}},
  \link[IRanges]{RleList-class},
  \link{pintersect,GRanges,GRanges-method},
  \link{findOverlaps,GRangesList,GenomicRanges-method},
  \link{countOverlaps,GRangesList,GenomicRanges-method},
  \link{subsetByOverlaps,GRangesList,GenomicRanges-method}
}

\examples{
library(Rsamtools)  # the toy file below is there
aln1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
aln1 <- readGappedAlignments(aln1_file)
aln1

## ---------------------------------------------------------------------
## A. BASIC MANIPULATION
## ---------------------------------------------------------------------
length(aln1)
head(aln1)
head(rname(aln1))
seqlevels(aln1)

## Rename the reference sequences:
rname(aln1) <- sub("seq", "chr", rname(aln1))
seqlevels(aln1)

head(strand(aln1))
head(cigar(aln1))
head(qwidth(aln1))
table(qwidth(aln1))

grglist(aln1)  # a GRangesList object
granges(aln1)  # a GRanges object
rglist(aln1)   # a CompressedNormalIRangesList object
ranges(aln1)   # an IRanges object
stopifnot(identical(elementLengths(grglist(aln1)), elementLengths(rglist(aln1))))

head(start(aln1))
head(end(aln1))
head(width(aln1))
head(ngap(aln1))

## ---------------------------------------------------------------------
## B. SUBSETTING
## ---------------------------------------------------------------------
aln1[strand(aln1) == "-"]
aln1[grep("I", cigar(aln1), fixed=TRUE)]
aln1[grep("N", cigar(aln1), fixed=TRUE)]  # no gaps

## A confirmation that all the queries map to the reference with no
## gaps:
stopifnot(all(ngap(aln1) == 0))

## Different ways to subset:
aln1[6]             # a GappedAlignments object of length 1
grglist(aln1)[[6]]  # a GRanges object of length 1
rglist(aln1)[[6]]   # a NormalIRanges object of length 1

## Ds are NOT gaps:
ii <- grep("D", cigar(aln1), fixed=TRUE)
aln1[ii]
ngap(aln1[ii])
grglist(aln1[ii])

## qwidth() vs width():
aln1[qwidth(aln1) != width(aln1)]

## This MUST return an empty object:
aln1[cigar(aln1) == "35M" & qwidth(aln1) != 35]
## but this doesn't have too:
aln1[cigar(aln1) != "35M" & qwidth(aln1) == 35]

## ---------------------------------------------------------------------
## C. qnarrow()/narrow()
## ---------------------------------------------------------------------
## Note that there is no difference between qnarrow() and narrow() when
## all the alignments are simple and with no indels.

## This trims 3 nucleotides on the left and 5 nucleotides on the right
## of each alignment:
qnarrow(aln1, start=4, end=-6)
## Note that the 'start' and 'end' arguments specify what part of each
## query sequence should be kept (negative values being relative to the
## right end of the query sequence), not what part should be trimmed.

## Trimming on the left doesn't change the "end" of the queries.
qnarrow(aln1, start=21)
stopifnot(identical(end(qnarrow(aln1, start=21)), end(aln1)))

## ---------------------------------------------------------------------
## D. coverage()
## ---------------------------------------------------------------------
coverage(aln1)

## ---------------------------------------------------------------------
## E. findOverlaps()/countOverlaps()
## ---------------------------------------------------------------------
findOverlaps(aln1, granges(aln1)[1])
sum(countOverlaps(aln1, granges(aln1)[1]))
subsetByOverlaps(aln1, granges(aln1)[1])
table(match(aln1, granges(aln1)[1]), useNA = "ifany")
table(aln1 \%in\% granges(aln1)[1])
}

\keyword{methods}
\keyword{classes}

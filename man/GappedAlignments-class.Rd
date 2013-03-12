\name{GappedAlignments-class}
\docType{class}

% Class:
\alias{class:GappedAlignments}
\alias{GappedAlignments-class}

% Constructors:
\alias{GappedAlignments}
\alias{updateObject,GappedAlignments-method}
\alias{readGappedAlignments}

% Accessors:
\alias{length,GappedAlignments-method}
\alias{names,GappedAlignments-method}
\alias{seqnames,GappedAlignments-method}
\alias{rname}
\alias{rname,GappedAlignments-method}
\alias{strand,GappedAlignments-method}
\alias{names<-,GappedAlignments-method}
\alias{seqnames<-,GappedAlignments-method}
\alias{rname<-}
\alias{rname<-,GappedAlignments-method}
\alias{strand<-,GappedAlignments-method}
\alias{cigar}
\alias{cigar,GappedAlignments-method}
\alias{qwidth}
\alias{qwidth,GappedAlignments-method}
\alias{start,GappedAlignments-method}
\alias{end,GappedAlignments-method}
\alias{width,GappedAlignments-method}
\alias{ngap,GappedAlignments-method}
\alias{elementMetadata<-,GappedAlignments-method}
\alias{seqinfo,GappedAlignments-method}
\alias{seqinfo<-,GappedAlignments-method}

% Coercion:
\alias{grglist}
\alias{granges}
\alias{introns}
\alias{rglist}
\alias{grglist,GappedAlignments-method}
\alias{granges,GappedAlignments-method}
\alias{introns,GappedAlignments-method}
\alias{rglist,GappedAlignments-method}
\alias{ranges,GappedAlignments-method}
\alias{coerce,GappedAlignments,GRangesList-method}
\alias{coerce,GappedAlignments,GRanges-method}
\alias{coerce,GappedAlignments,RangesList-method}
\alias{coerce,GappedAlignments,Ranges-method}
\alias{as.data.frame,GappedAlignments-method}

% Subsetting:
\alias{[,GappedAlignments-method}

% Combining:
\alias{c,GappedAlignments-method}

% "show" method:
\alias{show,GappedAlignments-method}

% Other methods:
\alias{qnarrow}
\alias{qnarrow,GappedAlignments-method}
\alias{narrow,GappedAlignments-method}
\alias{splitAsListReturnedClass,GappedAlignments-method}


\title{GappedAlignments objects}

\description{
  The GappedAlignments class is a simple container which
  purpose is to store a set of alignments that will hold just enough
  information for supporting the operations described below.
}

\details{
  A GappedAlignments object is a vector-like object where each element
  describes an alignment i.e. how a given sequence (called "query"
  or "read", typically short) aligns to a reference sequence (typically
  long).

  Typically, a GappedAlignments object will be created by loading
  records from a BAM (or SAM) file and each element in the resulting
  object will correspond to a record. BAM/SAM records generally contain
  a lot of information but only part of that information is loaded
  in the GappedAlignments object. In particular, we discard the query
  sequences (SEQ field), the query qualities (QUAL), the mapping qualities
  (MAPQ) and any other information that is not needed in order to support
  the operations or methods described below.

  This means that multi-reads (i.e. reads with multiple hits in the
  reference) won't receive any special treatment i.e. the various SAM/BAM
  records corresponding to a multi-read will show up in the GappedAlignments
  object as if they were coming from different/unrelated queries.
  Also paired-end reads will be treated as single-end reads and the
  pairing information will be lost (see \code{?\link{GappedAlignmentPairs}}
  for how to handle aligned paired-end reads).

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

  Optionally, a GappedAlignments object can have names (accessed thru the
  \code{\link[base]{names}} generic function) which will be coming from
  the QNAME field of the SAM/BAM records.

  The rest of this man page will focus on describing how to:
  \itemize{
    \item Access the information stored in a GappedAlignments object
          in a way that is independent from how the data are actually
          stored internally.
    \item How to create and manipulate a GappedAlignments object.
  }
}

\section{Constructors}{
  \describe{
    \item{}{
      \code{readGappedAlignments(file, format="BAM", use.names=FALSE, ...)}:
      Read a file containing aligned reads as a GappedAlignments object.
      By default (i.e. \code{use.names=FALSE}), the resulting object has no
      names. If \code{use.names} is \code{TRUE}, then the names are
      constructed from the query template names (QNAME field in a SAM/BAM
      file).

      Note that this function is just a front-end that delegates to the
      format-specific back-end function specified via the \code{format}
      argument. The \code{use.names} argument and any extra argument are
      passed to the back-end function.
      Only the BAM format is supported for now. Its back-end is the
      \code{\link[Rsamtools]{readBamGappedAlignments}} function defined
      in the Rsamtools package.
      See \code{?\link[Rsamtools]{readBamGappedAlignments}} for
      more information (you might need to install and load the Rsamtools
      package first).
    }
    \item{}{
      \code{GappedAlignments(seqnames=Rle(factor()), pos=integer(0),
                             cigar=character(0),
                             strand=NULL, names=NULL, seqlengths=NULL, ...)}:
      Low-level GappedAlignments constructor. Generally not used directly.
      Named arguments in \code{...} are used as metadata columns.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{x} is a GappedAlignments object.

  \describe{
    \item{}{
      \code{length(x)}:
      Return the number of alignments in \code{x}.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
      Get or set the names of \code{x}.
      See \code{readGappedAlignments} above for how to automatically extract
      and set the names from the file to read.
    }
    \item{}{
      \code{seqnames(x)}, \code{seqnames(x) <- value}:
      Get or set the name of the reference sequence for each alignment
      in \code{x} (see Details section above for more information about
      the RNAME field of a SAM/BAM file).
      \code{value} can be a factor, or a 'factor' \link[IRanges]{Rle},
      or a character vector.
    }
    \item{}{
      \code{rname(x)}, \code{rname(x) <- value}:
      Same as \code{seqnames(x)} and \code{seqnames(x) <- value}.
    }
    \item{}{
      \code{strand(x)}, \code{strand(x) <- value}:
      Get or set the strand for each alignment in \code{x} (see Details
      section above for more information about the strand of an alignment).
      \code{value} can be a factor (with levels +, - and *), or a 'factor'
      \link[IRanges]{Rle}, or a character vector.
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
      Returns an integer vector of the same length as \code{x} containing
      the number of gaps (i.e. N operations in the CIGAR) for each alignment.
      Equivalent to \code{unname(elementLengths(rglist(x))) - 1L}.
    }
    \item{}{
      \code{seqinfo(x)}, \code{seqinfo(x) <- value}:
      Get or set the information about the underlying sequences.
      \code{value} must be a \link{Seqinfo} object.
    }
    \item{}{
      \code{seqlevels(x)}, \code{seqlevels(x) <- value}:
      Get or set the sequence levels.
      \code{seqlevels(x)} is equivalent to \code{seqlevels(seqinfo(x))}
      or to \code{levels(seqnames(x))}, those 2 expressions being
      guaranteed to return identical character vectors on a GappedAlignments
      object. \code{value} must be a character vector with no NAs.
      See \code{?\link{seqlevels}} for more information.
    }
    \item{}{
      \code{seqlengths(x)}, \code{seqlengths(x) <- value}:
      Get or set the sequence lengths.
      \code{seqlengths(x)} is equivalent to \code{seqlengths(seqinfo(x))}.
      \code{value} can be a named non-negative integer or numeric vector
      eventually with NAs.
    }
    \item{}{
      \code{isCircular(x)}, \code{isCircular(x) <- value}:
      Get or set the circularity flags.
      \code{isCircular(x)} is equivalent to \code{isCircular(seqinfo(x))}.
      \code{value} must be a named logical vector eventually with NAs.
    }
    \item{}{
      \code{genome(x)}, \code{genome(x) <- value}:
      Get or set the genome identifier or assembly name for each sequence.
      \code{genome(x)} is equivalent to \code{genome(seqinfo(x))}.
      \code{value} must be a named character vector eventually with NAs.
    }
    \item{}{
      \code{seqnameStyle(x)}:
      Get or set the seqname style for \code{x}.
      Note that this information is not stored in \code{x} but inferred
      by looking up \code{seqnames(x)} against a seqname style database
      stored in the seqnames.db metadata package (required).
      \code{seqnameStyle(x)} is equivalent to \code{seqnameStyle(seqinfo(x))}
      and can return more than 1 seqname style (with a warning)
      in case the style cannot be determined unambiguously.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a GappedAlignments object.

  \describe{
    \item{}{
      \code{grglist(x, order.as.in.query=FALSE,
                       drop.D.ranges=FALSE)},
      \code{rglist(x, order.as.in.query=FALSE,
                      drop.D.ranges=FALSE)}:

      Return either a \link{GRangesList} or a \link[IRanges]{RangesList}
      object of length \code{length(x)} where the i-th element represents
      the ranges (with respect to the reference) of the i-th alignment in
      \code{x}.

      More precisely, the \link[IRanges]{RangesList} object returned by
      \code{rglist(x)} is a \link[IRanges]{CompressedIRangesList} object.

      The \code{order.as.in.query} toggle affects the order of the ranges
      \emph{within} each top-level element of the returned object.

      If \code{FALSE} (the default), then the ranges are ordered from 5' to 3'
      in elements associated with the plus strand (i.e. corresponding to
      alignments located on the plus strand), and from 3' to 5' in elements
      associated with the minus strand. So, whatever the strand is, the ranges
      are in ascending order (i.e. left-to-right).

      If \code{TRUE}, then the order of the ranges in elements associated
      with the \emph{minus} strand is reversed. So they end up being
      ordered from 5' to 3' too, which means that they are now in decending
      order (i.e. right-to-left). It also means that, when
      \code{order.as.in.query=TRUE} is used, the ranges are
      \emph{always} ordered consistently with the original "query template",
      that is, in the order defined by walking the "query template" from the
      beginning to the end.
      
      If \code{drop.D.ranges} is \code{TRUE}, then deletions (D operations
      in the CIGAR) are treated like gaps (N operations in the CIGAR),
      that is, the ranges corresponding to deletions are dropped.

      See Details section above for more information.
    }
    \item{}{
      \code{granges(x)}, \code{ranges(x)}:
      Return either a \link{GRanges} or a \link[IRanges]{Ranges}
      object of length \code{length(x)} where each element represents the
      regions in the reference to which a query is aligned.

      More precisely, the \link[IRanges]{Ranges} object returned by
      \code{ranges(x)} is an \link[IRanges]{IRanges} object.
    }
    \item{}{
      \code{introns(x)}: Extract the gaps (i.e. N operations in the CIGAR)
      as a \link{GRangesList} object of the same length as \code{x}.
      Equivalent to:
      \preformatted{
    psetdiff(granges(x), grglist(x, order.as.in.query=TRUE))
      }
    }
    \item{}{
      \code{as(x, "GRangesList")}, \code{as(x, "GRanges")},
      \code{as(x, "RangesList")}, \code{as(x, "Ranges")}:
      Alternate ways of doing \code{grglist(x)}, \code{granges(x)},
      \code{rglist(x)}, \code{ranges(x)}, respectively.
    }
  }
}

\section{Subsetting and related operations}{
  In the code snippets below, \code{x} is a GappedAlignments object.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new GappedAlignments object made of the selected
      alignments. \code{i} can be a numeric or logical vector.
    }
  }
}

\section{Combining}{ 
  \describe{
    \item{}{
      \code{c(...)}:
      Concatenates the GappedAlignment objects in \code{...}.
    }
  }
}

\section{Other methods}{

  \describe{
    \item{}{
      \code{qnarrow(x, start=NA, end=NA, width=NA)}:
      \code{x} is a GappedAlignments object.
      Return a new GappedAlignments object of the same length as \code{x}
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
      Return a new GappedAlignments object of the same length as \code{x}
      describing the narrowed alignments. Unlike with \code{qnarrow}
      now the \code{start}/\code{end}/\code{width} arguments describe
      the narrowing on the reference side, not the query side.
      Like with \code{qnarrow}, they must be vectors of integers.
      NAs and negative values are accepted and "solved" according to the
      rules of the SEW (Start/End/Width) interface (see
      \code{?\link[IRanges]{solveUserSEW}} for the details).
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
  \itemize{
    \item \link{GappedAlignmentPairs-class}.
    \item \code{\link[Rsamtools]{readBamGappedAlignments}}.
    \item \link{GRangesList-class}.
    \item \link{GRanges-class}.
    \item \link{findOverlaps-methods}.
    \item \link{coverage-methods}.
    \item \code{\link{seqinfo}}.
    \item \link[IRanges]{CompressedIRangesList-class}.
    \item \link{setops-methods}.
  }
}

\examples{
library(Rsamtools)  # for ScanBamParam() and the ex1.bam file
ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
gal <- readGappedAlignments(ex1_file, param=ScanBamParam(what="flag"))
gal

## ---------------------------------------------------------------------
## A. BASIC MANIPULATION
## ---------------------------------------------------------------------
length(gal)
head(gal)
names(gal)  # no names by default
seqnames(gal)
strand(gal)
head(cigar(gal))
head(qwidth(gal))
table(qwidth(gal))
head(start(gal))
head(end(gal))
head(width(gal))
head(ngap(gal))
seqlevels(gal)

## Rename the reference sequences:
seqlevels(gal) <- sub("seq", "chr", seqlevels(gal))
seqlevels(gal)

grglist(gal)  # a GRangesList object
stopifnot(identical(unname(elementLengths(grglist(gal))), ngap(gal) + 1L))
granges(gal)  # a GRanges object
rglist(gal)   # a CompressedIRangesList object
stopifnot(identical(unname(elementLengths(rglist(gal))), ngap(gal) + 1L))
ranges(gal)   # an IRanges object
introns(gal)  # a GRangesList object
stopifnot(identical(unname(elementLengths(introns(gal))), ngap(gal)))

## ---------------------------------------------------------------------
## B. SUBSETTING
## ---------------------------------------------------------------------
gal[strand(gal) == "-"]
gal[grep("I", cigar(gal), fixed=TRUE)]
gal[grep("N", cigar(gal), fixed=TRUE)]  # no gaps

## A confirmation that all the queries map to the reference with no
## gaps:
stopifnot(all(ngap(gal) == 0))

## Different ways to subset:
gal[6]             # a GappedAlignments object of length 1
grglist(gal)[[6]]  # a GRanges object of length 1
rglist(gal)[[6]]   # a NormalIRanges object of length 1

## D operations are NOT gaps:
ii <- grep("D", cigar(gal), fixed=TRUE)
gal[ii]
ngap(gal[ii])
grglist(gal[ii])

## qwidth() vs width():
gal[qwidth(gal) != width(gal)]

## This MUST return an empty object:
gal[cigar(gal) == "35M" & qwidth(gal) != 35]
## but this doesn't have too:
gal[cigar(gal) != "35M" & qwidth(gal) == 35]

## ---------------------------------------------------------------------
## C. qnarrow()/narrow()
## ---------------------------------------------------------------------
## Note that there is no difference between qnarrow() and narrow() when
## all the alignments are simple and with no indels.

## This trims 3 nucleotides on the left and 5 nucleotides on the right
## of each alignment:
qnarrow(gal, start=4, end=-6)
## Note that the 'start' and 'end' arguments specify what part of each
## query sequence should be kept (negative values being relative to the
## right end of the query sequence), not what part should be trimmed.

## Trimming on the left doesn't change the "end" of the queries.
qnarrow(gal, start=21)
stopifnot(identical(end(qnarrow(gal, start=21)), end(gal)))
}

\keyword{methods}
\keyword{classes}

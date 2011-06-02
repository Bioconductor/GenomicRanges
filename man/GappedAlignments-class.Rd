\name{GappedAlignments-class}
\docType{class}

% Class:
\alias{class:GappedAlignments}
\alias{GappedAlignments-class}
\alias{GappedAlignments}

% Constructors:
\alias{updateObject,GappedAlignments-method}
\alias{readGappedAlignments}

% Accessors:
\alias{length,GappedAlignments-method}
\alias{names,GappedAlignments-method}
\alias{rname}
\alias{rname,GappedAlignments-method}
\alias{seqnames,GappedAlignments-method}
\alias{strand,GappedAlignments-method}
\alias{rname<-}
\alias{rname<-,GappedAlignments-method}
\alias{seqnames<-,GappedAlignments-method}
\alias{strand<-,GappedAlignments-method}
\alias{cigar}
\alias{cigar,GappedAlignments-method}
\alias{qwidth}
\alias{qwidth,GappedAlignments-method}
\alias{grglist}
\alias{grglist,GappedAlignments-method}
\alias{granges}
\alias{granges,GappedAlignments-method}
\alias{grg}
\alias{rglist}
\alias{rglist,GappedAlignments-method}
\alias{ranges,GappedAlignments-method}
\alias{start,GappedAlignments-method}
\alias{end,GappedAlignments-method}
\alias{width,GappedAlignments-method}
\alias{ngap,GappedAlignments-method}
\alias{elementMetadata,GappedAlignments-method}
\alias{seqinfo,GappedAlignments-method}
\alias{seqinfo<-,GappedAlignments-method}

% Coercion:
\alias{as.data.frame,GappedAlignments-method}
\alias{coerce,GappedAlignments,GRangesList-method}
\alias{coerce,GappedAlignments,GRanges-method}
\alias{coerce,GappedAlignments,RangesList-method}
\alias{coerce,GappedAlignments,Ranges-method}

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

  Most of the time, a GappedAlignments object will be created by loading
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

\section{Constructor}{
  \describe{
    \item{}{
      \code{readGappedAlignments(file, format="BAM", use.names=FALSE, ...)}:
      Read a file as a GappedAlignments object.
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
      \code{names(x)}:
      Gets the names of \code{x}. See \code{readGappedAlignments} above
      for how the names are set.
    }
    \item{}{
      \code{rname(x)}, \code{rname(x) <- value}:
      Gets or sets the name of the reference sequence for each alignment
      in \code{x} (see Details section above for more information about
      the RNAME field of a SAM/BAM file).
      \code{value} can be a factor, or a 'factor' \link[IRanges]{Rle},
      or a character vector.
    }
    \item{}{
      \code{seqnames(x)}, \code{seqnames(x) <- value}:
      Same as \code{rname(x)} and \code{rname(x) <- value}.
    }
    \item{}{
      \code{strand(x)}, \code{strand(x) <- value}:
      Gets or sets the strand for each alignment in \code{x} (see Details
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
      \code{grglist(x, drop.D.ranges = FALSE)}, \code{granges(x)},
      \code{rglist(x, drop.D.ranges = FALSE)}, \code{ranges(x)}:
      Returns either a \link{GRangesList} object, or a \link{GRanges} object,
      or a \link[IRanges]{RangesList} object, or a \link[IRanges]{Ranges}
      object of length \code{length(x)} where each element represents the
      regions in the reference to which a query is
      aligned. If \code{drop.D.ranges} is \code{TRUE} for
      either \code{grglist} or \code{rglist}, the ranges corresponding
      to deletions in the CIGAR string are dropped, i.e., they are not
      considered part of the alignment but are treated like the N (intron)
      CIGAR element.
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
  \code{\link{seqinfo}},
  \link[IRanges]{CompressedNormalIRangesList-class},
  \link[IRanges]{IRanges-class},
  \link{coverage-methods},
  \link{setops-methods},
  \link{findOverlaps-methods}
}

\examples{
library(Rsamtools)  # because file ex1.bam is in this package
galn_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
galn <- readGappedAlignments(galn_file)
galn

## ---------------------------------------------------------------------
## A. BASIC MANIPULATION
## ---------------------------------------------------------------------
length(galn)
head(galn)
names(galn)  # no names by default
head(rname(galn))
seqlevels(galn)

## Rename the reference sequences:
seqlevels(galn) <- sub("seq", "chr", seqlevels(galn))
seqlevels(galn)

head(strand(galn))
head(cigar(galn))
head(qwidth(galn))
table(qwidth(galn))

grglist(galn)  # a GRangesList object
granges(galn)  # a GRanges object
rglist(galn)   # a CompressedNormalIRangesList object
ranges(galn)   # an IRanges object
stopifnot(identical(elementLengths(grglist(galn)), elementLengths(rglist(galn))))

head(start(galn))
head(end(galn))
head(width(galn))
head(ngap(galn))

## ---------------------------------------------------------------------
## B. SUBSETTING
## ---------------------------------------------------------------------
galn[strand(galn) == "-"]
galn[grep("I", cigar(galn), fixed=TRUE)]
galn[grep("N", cigar(galn), fixed=TRUE)]  # no gaps

## A confirmation that all the queries map to the reference with no
## gaps:
stopifnot(all(ngap(galn) == 0))

## Different ways to subset:
galn[6]             # a GappedAlignments object of length 1
grglist(galn)[[6]]  # a GRanges object of length 1
rglist(galn)[[6]]   # a NormalIRanges object of length 1

## Ds are NOT gaps:
ii <- grep("D", cigar(galn), fixed=TRUE)
galn[ii]
ngap(galn[ii])
grglist(galn[ii])

## qwidth() vs width():
galn[qwidth(galn) != width(galn)]

## This MUST return an empty object:
galn[cigar(galn) == "35M" & qwidth(galn) != 35]
## but this doesn't have too:
galn[cigar(galn) != "35M" & qwidth(galn) == 35]

## ---------------------------------------------------------------------
## C. qnarrow()/narrow()
## ---------------------------------------------------------------------
## Note that there is no difference between qnarrow() and narrow() when
## all the alignments are simple and with no indels.

## This trims 3 nucleotides on the left and 5 nucleotides on the right
## of each alignment:
qnarrow(galn, start=4, end=-6)
## Note that the 'start' and 'end' arguments specify what part of each
## query sequence should be kept (negative values being relative to the
## right end of the query sequence), not what part should be trimmed.

## Trimming on the left doesn't change the "end" of the queries.
qnarrow(galn, start=21)
stopifnot(identical(end(qnarrow(galn, start=21)), end(galn)))
}

\keyword{methods}
\keyword{classes}

\name{GappedAlignmentPairs-class}
\docType{class}

% Class:
\alias{class:GappedAlignmentPairs}
\alias{GappedAlignmentPairs-class}

% Constructors:
\alias{GappedAlignmentPairs}
\alias{readGappedAlignmentPairs}

% Accessors:
\alias{length,GappedAlignmentPairs-method}
\alias{names,GappedAlignmentPairs-method}
\alias{names<-,GappedAlignmentPairs-method}
\alias{first}
\alias{first,GappedAlignmentPairs-method}
\alias{last}
\alias{last,GappedAlignmentPairs-method}
\alias{left}
\alias{left,GappedAlignmentPairs-method}
\alias{right}
\alias{right,GappedAlignmentPairs-method}
\alias{seqnames,GappedAlignmentPairs-method}
\alias{strand,GappedAlignmentPairs-method}
\alias{ngap,GappedAlignmentPairs-method}
\alias{isProperPair}
\alias{isProperPair,GappedAlignmentPairs-method}
\alias{elementMetadata<-,GappedAlignmentPairs-method}
\alias{seqinfo,GappedAlignmentPairs-method}
\alias{seqinfo<-,GappedAlignmentPairs-method}

% Vector methods:
\alias{[,GappedAlignmentPairs-method}

% List methods:
\alias{[[,GappedAlignmentPairs,ANY,ANY-method}
\alias{unlist,GappedAlignmentPairs-method}

% Coercion:
\alias{grglist,GappedAlignmentPairs-method}
\alias{introns,GappedAlignmentPairs-method}
\alias{coerce,GappedAlignmentPairs,GRangesList-method}

% "show" method:
\alias{show,GappedAlignmentPairs-method}


\title{GappedAlignmentPairs objects}

\description{
  The GappedAlignmentPairs class is a container for "alignment pairs".
}

\details{
  A GappedAlignmentPairs object is a list-like object where each element
  describes an "alignment pair".

  An "alignment pair" is made of a "first" and a "last" alignment,
  and is formally represented by a \link{GappedAlignments} object of
  length 2. It is typically representing a hit of a paired-end read to
  the reference genome that was used by the aligner. More precisely,
  in a given pair, the "first" alignment represents the hit of the first
  end of the read (aka "first segment in the template", using SAM Spec
  terminology), and the "last" alignment represents the hit of the second
  end of the read (aka "last segment in the template", using SAM Spec
  terminology).

  In general, a GappedAlignmentPairs object will be created by loading
  records from a BAM (or SAM) file containing aligned paired-end reads,
  using the \code{readGappedAlignmentPairs} function (see below).
  Each element in the returned object will be obtained by pairing 2
  records.
}

\section{Constructors}{
  \describe{
    \item{}{
      \code{readGappedAlignmentPairs(file, format="BAM", use.names=FALSE, ...)}:
      Read a file containing paired-end reads as a GappedAlignmentPairs
      object.
      By default (i.e. \code{use.names=FALSE}), the resulting object has no
      names. If \code{use.names} is \code{TRUE}, then the names are
      constructed from the query template names (QNAME field in a SAM/BAM
      file). Note that the 2 records in a pair of records have the same QNAME.

      Note that this function is just a front-end that delegates to the
      format-specific back-end function specified via the \code{format}
      argument. The \code{use.names} argument and any extra argument are
      passed to the back-end function.
      Only the BAM format is supported for now. Its back-end is the
      \code{\link[Rsamtools]{readBamGappedAlignmentPairs}} function
      defined in the Rsamtools package.
      See \code{?\link[Rsamtools]{readBamGappedAlignmentPairs}} for
      more information (you might need to install and load the Rsamtools
      package first).
    }
    \item{}{
      \code{GappedAlignmentPairs(first, last, isProperPair, names=NULL)}:
      Low-level GappedAlignmentPairs constructor. Generally not used directly.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{x} is a GappedAlignmentPairs object.

  \describe{
    \item{}{
      \code{length(x)}:
      Return the number of alignment pairs in \code{x}.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
      Get or sets the names of \code{x}.
      See \code{readGappedAlignmentPairs} above for how to automatically
      extract and set the names from the file to read.
    }
    \item{}{
      \code{first(x, invert.strand=FALSE)},
      \code{last(x, invert.strand=FALSE)}:
      Get the "first" or "last" alignment for each alignment pair in
      \code{x}.
      The result is a \link{GappedAlignments} object of the same length
      as \code{x}.
      If \code{invert.strand=TRUE}, then the strand is inverted on-the-fly,
      i.e. "+" becomes "-", "-" becomes "+", and "*" remains unchanged.
    }
    \item{}{
      \code{left(x)}:
      Get the "left" alignment for each alignment pair in \code{x}.
      By definition, the "left" alignment in a pair is the alignment that
      is on the + strand. If this is the "first" alignment, then it's returned
      as-is by \code{left(x)}, but if this is the "last" alignment, then it's
      returned by \code{left(x)} with the strand inverted.
    }
    \item{}{
      \code{right(x)}:
      Get the "right" alignment for each alignment pair in \code{x}.
      By definition, the "right" alignment in a pair is the alignment that
      is on the - strand. If this is the "first" alignment, then it's returned
      as-is by \code{right(x)}, but if this is the "last" alignment, then it's
      returned by \code{right(x)} with the strand inverted.
    }
    \item{}{
      \code{seqnames(x)}:
      Get the name of the reference sequence for each alignment pair
      in \code{x}. This comes from the RNAME field of the BAM file and
      has the same value for the 2 records in a pair
      (\code{\link[Rsamtools]{makeGappedAlignmentPairs}}, the function
      used by \code{\link[Rsamtools]{readBamGappedAlignmentPairs}} for
      doing the pairing, rejects pairs with incompatible RNAME values).
    }
    \item{}{
      \code{strand(x)}:
      Get the strand for each alignment pair in \code{x}. By definition
      the strand of an alignment pair is the strand of the \emph{"first"}
      alignment in the pair. In a GappedAlignmentPairs object, the strand
      of the "last" alignment in a pair is \emph{always} the opposite of the
      strand of the "first" alignment
      (\code{\link[Rsamtools]{makeGappedAlignmentPairs}}, the function
      used by \code{\link[Rsamtools]{readBamGappedAlignmentPairs}} for
      doing the pairing, rejects pairs where the "first" and "last"
      alignments are on the same strand).
    }
    \item{}{
      \code{ngap(x)}:
      Equivalent to \code{ngap(first(x)) + ngap(last(x))}.
    }
    \item{}{
      \code{isProperPair(x)}:
      Get the "isProperPair" flag bit (bit 0x2 in SAM Spec) set by
      the aligner for each alignment pair in \code{x}.
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
      guaranteed to return identical character vectors on a
      GappedAlignmentPairs object. \code{value} must be a character vector
      with no NAs.
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
      List the matching seqname styles for \code{x}.
      \code{seqnameStyle(x)} is equivalent to \code{seqnameStyle(seqinfo(x))}.
      Note that this information is not stored in \code{x} but inferred
      by looking up \code{seqlevels(x)} against a seqname style database
      stored in the seqnames.db metadata package (required).
    }
  }
}

\section{Vector methods}{
  In the code snippets below, \code{x} is a GappedAlignmentPairs object.

  \describe{
    \item{}{
      \code{x[i]}:
      Return a new GappedAlignmentPairs object made of the selected
      alignment pairs. 
    }
  }
}

\section{List methods}{
  In the code snippets below, \code{x} is a GappedAlignmentPairs object.

  \describe{
    \item{}{
      \code{x[[i]]}:
      Extract the i-th alignment pair as a \link{GappedAlignments} object
      of length 2. As expected \code{x[[i]][1]} and \code{x[[i]][2]} are
      respectively the "first" and "last" alignments in the pair.
    }
    \item{}{
      \code{unlist(x, use.names=TRUE)}:
      Return the \link{GappedAlignments} object conceptually defined
      by \code{c(x[[1]], x[[2]], ..., x[[length(x)]])}.
      \code{use.names} determines whether \code{x} names should be
      propagated to the result or not.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a GappedAlignmentPairs object.

  \describe{
    \item{}{
      \code{grglist(x, order.as.in.query=FALSE, drop.D.ranges=FALSE)}:

      Return a \link{GRangesList} object of length \code{length(x)}
      where the i-th element represents the ranges (with respect to the
      reference) of the i-th alignment pair in \code{x}.

      IMPORTANT: The strand of the ranges coming from the "last" alignment
      in the pair is \emph{always} inverted.

      The \code{order.as.in.query} toggle affects the order of the ranges
      \emph{within} each top-level element of the returned object.

      If \code{FALSE} (the default), then the "left" ranges are placed before
      the "right" ranges, and, within each left or right group, are ordered
      from 5' to 3' in elements associated with the plus strand and from 3'
      to 5' in elements associated with the minus strand.
      More formally, the i-th element in the returned \link{GRangesList}
      object can be defined as \code{c(grl1[[i]], grl2[[i]])}, where
      \code{grl1} is \code{grglist(left(x))} and \code{grl2} is
      \code{grglist(right(x))}.

      If \code{TRUE}, then the "first" ranges are placed before the "last"
      ranges, and, within each first or last group, are \emph{always}
      ordered from 5' to 3', whatever the strand is.
      More formally, the i-th element in the returned \link{GRangesList}
      object can be defined as \code{c(grl1[[i]], grl2[[i]])}, where
      \code{grl1} is \code{grglist(first(x),
                                   order.as.in.query=TRUE)}
      and
      \code{grl2} is \code{grglist(last(x, invert.strand=TRUE),
                                   order.as.in.query=TRUE)}.

      Note that the relationship between the 2 \link{GRangesList} objects
      obtained with \code{order.as.in.query} being respectively
      \code{FALSE} or \code{TRUE} is simpler than it sounds: the only
      difference is that the order of the ranges in elements associated
      with the \emph{minus} strand is reversed.

      Finally note that, in the latter, the ranges are \emph{always} ordered
      consistently with the original "query template", that is, in the order
      defined by walking the "query template" from the beginning to the end.

      If \code{drop.D.ranges} is \code{TRUE}, then deletions (Ds in the
      CIGAR) are treated like gaps (Ns in the CIGAR), that is, the ranges
      corresponding to deletions are dropped.
    }
    \item{}{
      \code{introns(x)}: Extract the gaps (i.e. N operations in the CIGAR)
      of the "first" and "last" alignments of each pair as a
      \link{GRangesList} object of the same length as \code{x}.
      Equivalent to (but faster than):
      \preformatted{
    introns1 <- introns(first(x))
    introns2 <- introns(last(x, invert.strand=TRUE))
    mendoapply(c, introns1, introns2)
      }
    }
    \item{}{
      \code{as(x, "GRangesList")}:
      An alternate way of doing \code{grglist(x)}.
    }
  }
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item \link{GappedAlignments-class}.
    \item \code{\link[Rsamtools]{readBamGappedAlignmentPairs}}.
    \item \code{\link[Rsamtools]{makeGappedAlignmentPairs}}.
    \item \link{GRangesList-class}.
    \item \link{GRanges-class}.
    \item \link{findOverlaps-methods}.
    \item \link{coverage-methods}.
    \item \code{\link{seqinfo}}.
  }
}

\examples{
ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
galp <- readGappedAlignmentPairs(ex1_file, use.names=TRUE)
galp

length(galp)
head(galp)
head(names(galp))
first(galp)
last(galp)
last(galp, invert.strand=TRUE)
left(galp)
right(galp)
seqnames(galp)
strand(galp)
head(ngap(galp))
table(isProperPair(galp))
seqlevels(galp)

## Rename the reference sequences:
seqlevels(galp) <- sub("seq", "chr", seqlevels(galp))
seqlevels(galp)

galp[[1]]
unlist(galp)

grglist(galp)  # a GRangesList object
grglist(galp, order.as.in.query=TRUE)
stopifnot(identical(unname(elementLengths(grglist(galp))), ngap(galp) + 2L))
introns(galp)  # a GRangesList object
stopifnot(identical(unname(elementLengths(introns(galp))), ngap(galp)))
}

\keyword{methods}
\keyword{classes}

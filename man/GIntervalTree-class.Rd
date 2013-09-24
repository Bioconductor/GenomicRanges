\name{GIntervalTree-class}
\docType{class}

% Class:
\alias{class:GIntervalTree}
\alias{GIntervalTree-class}
\alias{GIntervalTree}

% Constructors:
\alias{GIntervalTree}
%\alias{updateObject,GRanges-method}

% Coercion:
\alias{coerce,GRanges,GIntervalTree-method}
\alias{coerce,GIntervalTree,GRanges-method}

% Accessors:
\alias{seqnames,GIntervalTree-method}
\alias{ranges,GIntervalTree-method}
\alias{strand,GIntervalTree-method}
\alias{elementMetadata,GIntervalTree-method}
\alias{names,GIntervalTree-method}
\alias{seqinfo,GIntervalTree-method}
\alias{score,GIntervalTree-method}

% Ranges methods:
\alias{start,GIntervalTree-method}
\alias{end,GIntervalTree-method}
\alias{width,GIntervalTree-method}

% Vector methods:
\alias{[,GIntervalTree-method}
\alias{show,GIntervalTree-method}

\title{GIntervalTree objects}

\description{
  The GRanges class is a container for the genomic locations and
  their associated annotations.
}

\details{
  A common type of query that arises when working with intervals is
  finding which intervals in one set overlap those in another. An
  efficient family of algorithms for answering such queries is known as
  the Interval Tree. The \code{GIntervalTree} class implements persistent
  Interval Trees for efficient querying of genomic intervals. It uses
  the \link[IRanges]{IntervalForest} class to store a set of trees,
  one for each \code{seqlevel} in a \code{\link{GRanges}} object.

  The simplest approach for finding overlaps is to call the
  \code{\link{findOverlaps}} function on a \link{Ranges} or other object
  with range information. See the man page of \code{\link{findOverlaps}}
  for how to use this and other related functions. A \code{GIntervalTree} 
  object is a derivative of \link{GenomicRanges} and
  stores its genomic ranges as a set of trees (a forest, with one tree per 
  \code{seqlevel}) that is optimized for overlap queries.
  Thus, for repeated queries against the same subject, it is more
  efficient to create a \code{GIntervalTree} once for the subject using
  the constructor described below and then perform the queries against
  the \code{GIntervalTree} instance.

  Like its \link{GenomicRanges} parent class, the GIntervalTree class stores 
  the sequences of genomic locations and associated
  annotations. Each element in the sequence is comprised of a sequence name,
  an interval, a \link{strand}, and optional metadata columns (e.g. score, GC
  content, etc.). This information is stored in four components as in 
  \link{GenomicRanges}, but two of these components are treated in a specific way:
  \describe{
    \item{\code{ranges}}{an \link[IRanges]{IntervalForest} object containing
      the ranges stored as a set of interval trees.}
    \item{\code{seqnames}}{these are not stored directly in this class, but are obtained from
    the \code{partitioning} component of the \link[IRanges]{IntervalForest} object stored
    in \code{ranges}.}
  }

  Note that \code{GIntervalTree} objects are not supported for \link{GRanges} objects
  with circular genomes.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GIntervalTree(x)}: Creates a \code{GIntervalTree} object from a \link{GRanges}
            object.
      \describe{
        \item{\code{x}}{a \link{GRanges} object containing the genomic ranges.}
      }
    }
  }
}

\section{Coercion}{
  \describe{
    \item{}{
      \code{as(from, "GIntervalTree")}: Creates a \code{GIntervalTree} object from a
      \link{GRanges} object.
      \code{as(from, "GRanges")}: Creates a \link{GRanges} object from an
      \code{GIntervalTree} object
    }
  }
}

\section{Accessors}{
  In the following code snippets, \code{x} is a GIntervalTree object.

  \describe{
    \item{}{
      \code{length(x)}:
      Get the number of elements.
    }
    \item{}{
      \code{seqnames(x)}:
      Get the sequence names.
    }
    \item{}{
      \code{ranges(x)}:
      Get the ranges as an \link[IRanges]{IRanges} object. This is for consistency
      with the \code{ranges} accessor for \link{GRanges} objects. To access the underlying
      \link[IRanges]{IntervalForest} object use the \code{obj@ranges} form.
    }
    \item{}{
      \code{strand(x)}:
      Get the strand. 
    }
    \item{}{
      \code{mcols(x, use.names=FALSE)}, \code{mcols(x) <- value}:
      Get or set the metadata columns.
      If \code{use.names=TRUE} and the metadata columns are not \code{NULL},
      then the names of \code{x} are propagated as the row names of the
      returned \link{DataFrame} object.
      When setting the metadata columns, the supplied value must be \code{NULL}
      or a data.frame-like object (i.e. \link{DataTable} or data.frame)
      object holding element-wise metadata.
    }
    \item{}{
       \code{elementMetadata(x)}, \code{elementMetadata(x) <- value},
       \code{values(x)}, \code{values(x) <- value}:
       Alternatives to \code{mcols} functions. Their use is discouraged.
    }
    \item{}{
      \code{seqinfo(x)}:
      Get or set the information about the underlying sequences.
      \code{value} must be a \link{Seqinfo} object.
    }
    \item{}{
      \code{seqlevels(x)}:
      Get the sequence levels. These are stored in the \code{partition}
      slot of the underlying \link[IRanges]{IntervalForest} object.
    }
    \item{}{
      \code{seqlengths(x)}:
      Get the sequence lengths.
    }
    \item{}{
      \code{isCircular(x)}:
      Get the circularity flags. Note that GIntervalTree objects are not supported
      for circular genomes.
    }
    \item{}{
      \code{genome(x)}:
      Get or the genome identifier or assembly name for each sequence.
    }
    \item{}{
      \code{seqnameStyle(x)}:
      Get the seqname style for \code{x}.
    }
    \item{}{
      \code{score(x)}: Get the \dQuote{score} column from the element
      metadata, if any.
    }
  }
}

\section{Ranges methods}{
  In the following code snippets, \code{x} is a GIntervalTree object.

  \describe{
    \item{}{
      \code{start(x)}:
      Get \code{start(ranges(x))}.
    }
    \item{}{
      \code{end(x)}:
      Get \code{end(ranges(x))}.
    }
    \item{}{
      \code{width(x)}:
      Get \code{width(ranges(x))}.
    }
  }
}

\section{Subsetting}{
  In the code snippets below, \code{x} is a GIntervalTree object.

  \describe{
    \item{}{
      \code{x[i, j]}:
      Get elements \code{i} with optional metadata columns
      \code{mcols(x)[,j]}, where \code{i} can be missing; an NA-free
      logical, numeric, or character vector; or a 'logical' Rle object.
    }
    }
}

\section{Other methods}{

  \describe{
    \item{}{
      \code{show(x)}:
      By default the \code{show} method displays 5 head and 5 tail
      elements. This can be changed by setting the global options
      \code{showHeadLines} and \code{showTailLines}. If the object
      length is less than (or equal to) the sum of these 2 options
      plus 1, then the full object is displayed.
    }
  }
}

\author{Hector Corrada Bravo, P. Aboyoun}

\seealso{
  \code{\link{seqinfo}},
  \link[IRanges]{IntervalForest},
  \link[IRanges]{IntervalTree},
  \link[GenomicRanges]{findOverlaps-methods},
}

\examples{
seqinfo <- Seqinfo(paste0("chr", 1:3), c(1000, 2000, 1500), NA, "mock1")
gr <-
  GRanges(seqnames =
          Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
          ranges = IRanges(
            1:10, width = 10:1, names = head(letters,10)),
          strand = Rle(
            strand(c("-", "+", "*", "+", "-")),
            c(1, 2, 2, 3, 2)),
          score = 1:10,
          GC = seq(1, 0, length=10),
          seqinfo=seqinfo)
tree <- GIntervalTree(gr)
tree

## Summarizing elements
table(seqnames(tree))
sum(width(tree))
summary(mcols(tree)[,"score"])

## find Overlaps
subject <-
    GRanges(seqnames =
              Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
            ranges =
              IRanges(1:10, width = 10:1, names = head(letters,10)),
            strand =
              Rle(strand(c("-", "+", "*", "+", "-")),
                  c(1, 2, 2, 3, 2)),
            score = 1:10,
            GC = seq(1, 0, length=10))
  query <-
    GRanges(seqnames = "chr2", ranges = IRanges(4:3, 6),
            strand = "+", score = 5:4, GC = 0.45)
  
  stree <- GIntervalTree(subject)
  findOverlaps(query, stree)
  countOverlaps(query, stree)
}

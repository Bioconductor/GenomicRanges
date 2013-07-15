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
\alias{length,GIntervalTree-method}
\alias{[,GIntervalTree-method}

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
    \item{\code{seqnames}}{in this class, these are not stored directly, but are stored
    as the \code{partition} component of the \link[IRanges]{IntervalForest} object stored
    in \code{ranges}.
      containing the sequence names.}
  }

  Note that \code{GIntervalTree} objects are not supported for \link{GRanges} objects
  with circular genomes.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GIntervalTree(x)}: Creates a GIntervalTree object from a \link{GRanges}
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
      \code{as(from, "GIntervalTree")}: Creates a GIntervalTree object from a
      \link{GRanges} object.
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
    \item{}{
      \code{head(x, n = 6L)}:
      If \code{n} is non-negative, returns the first n elements of the
      GIntervalTree object.
      If \code{n} is negative, returns all but the last \code{abs(n)} elements
      of the GIntervalTree object.
    }
    % \item{}{
    %   \code{rep(x, times, length.out, each)}:
    %   Repeats the values in \code{x} through one of the following conventions:
    %   \describe{
    %     \item{\code{times}}{Vector giving the number of times to repeat each
    %       element if of length \code{length(x)}, or to repeat the whole vector
    %       if of length 1.}
    %     \item{\code{length.out}}{Non-negative integer. The desired length of
    %       the output vector.}
    %     \item{\code{each}}{Non-negative integer.  Each element of \code{x} is
    %       repeated \code{each} times.}
    %   }
    % }
    % \item{}{
    %   \code{seqselect(x, start=NULL, end=NULL, width=NULL)}:
    %   Similar to \code{window}, except that multiple consecutive subsequences
    %   can be requested for concatenation. As such two of the three \code{start},
    %   \code{end}, and \code{width} arguments can be used to specify the
    %   consecutive subsequences. Alternatively, \code{start} can take a Ranges
    %   object or something that can be converted to a Ranges object like an
    %   integer vector, logical vector or logical Rle. If the concatenation of
    %   the consecutive subsequences is undesirable, consider using
    %   \code{\link{Views}}.
    % }
    % \item{}{
    %   \code{seqselect(x, start=NULL, end=NULL, width=NULL) <- value}:
    %   Similar to \code{window<-}, except that multiple consecutive subsequences
    %   can be replaced with a \code{value} whose length is a divisor of the
    %   number of elements it is replacing. As such two of the three \code{start},
    %   \code{end}, and \code{width} arguments can be used to specify the
    %   consecutive subsequences. Alternatively, \code{start} can take a Ranges
    %   object or something that can be converted to a Ranges object like an
    %   integer vector, logical vector or logical Rle.
    % }
    % \item{}{
    %   \code{subset(x, subset)}:
    %   Returns a new object of the same class as \code{x} made of the subset
    %   using logical vector \code{subset}, where missing values are taken as
    %   \code{FALSE}.
    % }
    \item{}{
      \code{tail(x, n = 6L)}:
      If \code{n} is non-negative, returns the last n elements of the
      GRanges object.
      If \code{n} is negative, returns all but the first \code{abs(n)} elements
      of the GRanges object.
    }
  %   \item{}{
  %     \code{window(x, start = NA, end = NA, width = NA, frequency = NULL, delta = NULL, ...)}:
  %     Extracts the subsequence window from the GRanges object using:
  %     \describe{
  %       \item{\code{start}, \code{end}, \code{width}}{The start, end, or width
  %         of the window. Two of the three are required.}
  %       \item{\code{frequency}, \code{delta}}{Optional arguments that specify
  %         the sampling frequency and increment within the window.}
  %     }
  %     In general, this is more efficient than using \code{"["} operator.
  %   }
  %   \item{}{
  %     \code{window(x, start = NA, end = NA, width = NA, keepLength = TRUE) <- value}:
  %     Replaces the subsequence window specified on the left (i.e. the subsequence
  %     in \code{x} specified by \code{start}, \code{end} and \code{width})
  %     by \code{value}.
  %     \code{value} must either be of class \code{class(x)}, belong to a subclass
  %     of \code{class(x)}, be coercible to \code{class(x)}, or be \code{NULL}.
  %     If \code{keepLength} is \code{TRUE}, the elements of \code{value} are
  %     repeated to create a GRanges object with the same number of elements
  %     as the width of the subsequence window it is replacing.
  %     If \code{keepLength} is \code{FALSE}, this replacement method can modify
  %     the length of \code{x}, depending on how the length of the left
  %     subsequence window compares to the length of \code{value}.
  %   }
  %   \item{}{
  %     \code{x$name}, \code{x$name <- value}:
  %     Shortcuts for \code{mcols(x)$name} and \code{mcols(x)$name <- value},
  %     respectively. Provided as a convenience, for GRanges objects *only*,
  %     and as the result of strong popular demand.
  %     Note that those methods are not consistent with the other \code{$}
  %     and \code{$<-} methods in the IRanges/GenomicRanges infrastructure,
  %     and might confuse some users by making them believe that a GRanges
  %     object can be manipulated as a data.frame-like object.
  %     Therefore we recommend using them only interactively, and we discourage
  %     their use in scripts or packages. For the latter, use
  %     \code{mcols(x)$name} and \code{mcols(x)$name <- value}, instead
  %     of \code{x$name} and \code{x$name <- value}, respectively.
  %   }
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

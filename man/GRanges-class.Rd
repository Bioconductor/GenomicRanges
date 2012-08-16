\name{GRanges-class}
\docType{class}

% Class:
\alias{class:GenomicRanges}
\alias{GenomicRanges-class}
\alias{GenomicRanges}
\alias{class:GRanges}
\alias{GRanges-class}
\alias{GRanges}
\alias{GenomicRangesORGRangesList-class}
\alias{GenomicRangesORmissing-class}

% Constructors:
\alias{GRanges}
\alias{updateObject,GRanges-method}

% Coercion:
\alias{coerce,RangedData,GRanges-method}
\alias{coerce,GenomicRanges,RangedData-method}
\alias{coerce,RangesList,GRanges-method}
\alias{coerce,GenomicRanges,RangesList-method}
\alias{coerce,GenomicRanges,GappedAlignments-method}
\alias{coerce,RleList,GRanges-method}
\alias{coerce,RleViewsList,GRanges-method}
\alias{as.data.frame,GenomicRanges-method}

% Accessors:
\alias{seqnames,GRanges-method}
\alias{seqnames<-,GenomicRanges-method}
\alias{seqnameStyle,GRanges-method}
\alias{ranges,GRanges-method}
\alias{ranges<-,GenomicRanges-method}
\alias{strand,GRanges-method}
\alias{strand<-,GenomicRanges-method}
\alias{elementMetadata,GenomicRanges-method}
\alias{elementMetadata<-,GenomicRanges-method}
\alias{names,GenomicRanges-method}
\alias{names<-,GenomicRanges-method}
\alias{seqinfo,GRanges-method}
\alias{seqinfo<-,GenomicRanges-method}
\alias{score,GenomicRanges-method}

% Ranges methods:
\alias{start,GenomicRanges-method}
\alias{start<-,GenomicRanges-method}
\alias{end,GenomicRanges-method}
\alias{end<-,GenomicRanges-method}
\alias{width,GenomicRanges-method}
\alias{width<-,GenomicRanges-method}

% Vector methods:
\alias{length,GenomicRanges-method}
\alias{[,GenomicRanges-method}
\alias{[<-,GenomicRanges,ANY,ANY,ANY-method}
\alias{c,GenomicRanges-method}
\alias{seqselect,GenomicRanges-method}
\alias{seqselect<-,GenomicRanges-method}
\alias{window,GenomicRanges-method}

% $ and $<- methods:
\alias{$,GenomicRanges-method}
\alias{$<-,GenomicRanges-method}

% "show" method:
\alias{show,GenomicRanges-method}

\title{GRanges objects}

\description{
  The GRanges class is a container for the genomic locations and
  their associated annotations.
}

\details{
  The GRanges class stores the sequences of genomic locations and associated
  annotations. Each element in the sequence is comprised of a sequence name,
  an interval, a \link{strand}, and optional metadata columns (e.g. score, GC
  content, etc.). This information is stored in four components:
  \describe{
    \item{\code{seqnames}}{a 'factor' \link[IRanges]{Rle} object
      containing the sequence names.}
    \item{\code{ranges}}{an \link[IRanges]{IRanges} object containing
      the ranges.}
    \item{\code{strand}}{a 'factor' \link[IRanges]{Rle} object containing
      the \link{strand} information.}
    \item{\code{mcols}}{a \link[IRanges]{DataFrame} object
      containing the metadata columns. Columns cannot be named
      \code{"seqnames"}, \code{"ranges"}, \code{"strand"},
      \code{"seqlevels"}, \code{"seqlengths"}, \code{"isCircular"},
      \code{"genome"}, \code{"start"}, \code{"end"}, \code{"width"},
      or \code{"element"}.}
    \item{\code{seqinfo}}{a \link{Seqinfo} object containing information
      about the set of genomic sequences present in the GRanges object.}
  }
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GRanges(seqnames = Rle(), ranges = IRanges(),
                    strand = Rle("*", length(seqnames)),
                    ..., seqinfo)}:
      Creates a GRanges object.
      \describe{
        \item{\code{seqnames}}{Rle object, character vector, or factor
          containing the sequence names.}
        \item{\code{ranges}}{IRanges object containing the ranges.}
        \item{\code{strand}}{Rle object, character vector, or factor
          containing the strand information.}
        \item{\code{seqinfo}}{a Seqinfo object containing allowed
          sequence names and lengths (or NA) for each
          \code{level(seqnames)}.}
        \item{\code{\ldots}}{Optional metadata columns.
          These columns cannot be named
          \code{"start"}, \code{"end"}, \code{"width"}, or
          \code{"element"}. A named integer vector \code{"seqlength"}
	  can be used instead of \code{seqinfo}.}
      }
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a GRanges object.

  \describe{
    \item{}{
      \code{as(from, "GRanges")}: Creates a GRanges object from a
      RangedData, RangesList, RleList or RleViewsList object.
    }
    \item{}{
      \code{as(from, "RangedData")}:
      Creates a RangedData object from a GRanges
      object. The \code{strand} and metadata columns become columns
      in the result. The \code{seqlengths(from)}, \code{isCircular(from)},
      and \code{genome(from)} vectors are stored in the metadata columns
      of \code{ranges(rd)}.
    }
    \item{}{
      \code{as(from, "RangesList")}:
      Creates a RangesList object from a GRanges
      object. The \code{strand} and metadata columns become \emph{inner}
      metadata columns (i.e. metadata columns on the ranges).
      The \code{seqlengths(from)}, \code{isCircular(from)}, and
      \code{genome(from)} vectors become the metadata columns.
    }
    \item{}{
      \code{as(from, "GappedAlignments")}:
      Creates a GappedAlignments object from a GRanges object. The metadata
      columns are propagated. cigar values are created from the sequence
      width unless a "cigar" metadata column already exists in \code{from}.
    }
    \item{}{
      \code{as.data.frame(x, row.names = NULL, optional = FALSE, ...)}:
      Creates a data.frame with columns \code{seqnames} (factor),
      \code{start} (integer), \code{end} (integer), \code{width} (integer),
      \code{strand} (factor), as well as the additional metadata columns
      stored in \code{mcols(x)}. Pass an explicit
      \code{stringsAsFactors=TRUE/FALSE} argument via \code{\ldots} to
      override the default conversions for the metadata columns in
      \code{mcols(x)}.
    }
  }
}

\section{Accessors}{
  In the following code snippets, \code{x} is a GRanges object.

  \describe{
    \item{}{
      \code{length(x)}:
      Gets the number of elements.
    }
    \item{}{
      \code{seqnames(x)}, \code{seqnames(x) <- value}:
      Gets or sets the sequence names.
      \code{value} can be an \link[IRanges]{Rle} object, a character vector,
      or a factor.
    }
    \item{}{
      \code{seqnameStyle(x)}:
      Lists the matching seqname styles for \code{x}.
    }
    \item{}{
      \code{ranges(x)}, \code{ranges(x) <- value}:
      Gets or sets the ranges. \code{value} can be a Ranges object.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
      Gets or sets the names of the elements.
    }
    \item{}{
      \code{strand(x)}, \code{strand(x) <- value}:
      Gets or sets the strand. \code{value} can be an Rle object, character
      vector, or factor.
    }
    \item{}{
      \code{mcols(x)}, \code{mcols(x) <- value}:
      Gets or sets the metadata columns.
      \code{value} can be \code{NULL}, or a data.frame-like object (i.e.
      \link{DataFrame} or data.frame) holding element-wise metadata.
    }
    \item{}{
       \code{elementMetadata(x)}, \code{elementMetadata(x) <- value},
       \code{values(x)}, \code{values(x) <- value}:
       Alternatives to \code{mcols} functions. Their use is discouraged.
    }
    \item{}{
      \code{seqinfo(x)}, \code{seqinfo(x) <- value}:
      Gets or sets the information about the underlying sequences.
      \code{value} must be a \link{Seqinfo} object.
    }
    \item{}{
      \code{seqlevels(x)}, \code{seqlevels(x, force=FALSE) <- value}:
      Gets or sets the sequence levels.
      \code{seqlevels(x)} is equivalent to \code{seqlevels(seqinfo(x))}
      or to \code{levels(seqnames(x))}, those 2 expressions being
      guaranteed to return identical character vectors on a GRanges object.
      \code{value} must be a character vector with no NAs.
      See \code{?\link{seqlevels}} for more information.
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
    \item{}{
      \code{genome(x)}, \code{genome(x) <- value}:
      Gets or sets the genome identifier or assembly name for each sequence.
      \code{genome(x)} is equivalent to \code{genome(seqinfo(x))}.
      \code{value} must be a named character vector eventually with NAs.
    }
    \item{}{
      \code{score(x)}: Gets the \dQuote{score} column from the element
      metadata, if any.
    }
  }
}

\section{Ranges methods}{
  In the following code snippets, \code{x} is a GRanges object.

  \describe{
    \item{}{
      \code{start(x)}, \code{start(x) <- value}:
      Gets or sets \code{start(ranges(x))}.
    }
    \item{}{
      \code{end(x)}, \code{end(x) <- value}:
      Gets or sets \code{end(ranges(x))}.
    }
    \item{}{
      \code{width(x)}, \code{width(x) <- value}:
      Gets or sets \code{width(ranges(x))}.
    }
  }
}

\section{Splitting and Combining}{
  In the code snippets below, \code{x} is a GRanges object.

  \describe{
    \item{}{
      \code{append(x, values, after = length(x))}:
      Inserts the \code{values} into \code{x} at the position given by
      \code{after}, where \code{x} and \code{values} are of the same
      class.
    }
    \item{}{
      \code{c(x, ...)}:
      Combines \code{x} and the GRanges objects in \code{...} together.
      Any object in \code{...} must belong to the same class as \code{x},
      or to one of its subclasses, or must be \code{NULL}.
      The result is an object of the same class as \code{x}.
    }
    \item{}{
      \code{c(x, ..., ignore.mcols=FALSE)}
      If the \code{GRanges} objects have metadata columns (represented as one
      \link{DataFrame} per object), each such \link{DataFrame} must have the
      same columns in order to combine successfully. In order to circumvent
      this restraint, you can pass in an \code{ignore.mcols=TRUE} argument
      which will combine all the objects into one and drop all of their
      metadata columns.
    }
    \item{}{
      \code{split(x, f, drop=FALSE)}:
      Splits \code{x} according to \code{f} to create a
      \link{GRangesList} object.
      If \code{f} is a list-like object then \code{drop} is ignored
      and \code{f} is treated as if it was
      \code{rep(seq_len(length(f)), sapply(f, length))},
      so the returned object has the same shape as \code{f} (it also
      receives the names of \code{f}).
      Otherwise, if \code{f} is not a list-like object, empty list
      elements are removed from the returned object if \code{drop} is
      \code{TRUE}.
    }
  }
}

\section{Subsetting}{
  In the code snippets below, \code{x} is a GRanges object.

  \describe{
    \item{}{
      \code{x[i, j]}, \code{x[i, j] <- value}:
      Gets or sets elements \code{i} with optional metadata columns
      \code{mcols(x)[,j]}, where \code{i} can be missing; an NA-free
      logical, numeric, or character vector; or a 'logical' Rle object.
    }
    \item{}{
      \code{x[i,j] <- value}:
      Replaces elements \code{i} and optional metadata columns \code{j}
      with \code{value}.
    }
    \item{}{
      \code{head(x, n = 6L)}:
      If \code{n} is non-negative, returns the first n elements of the
      GRanges object.
      If \code{n} is negative, returns all but the last \code{abs(n)} elements
      of the GRanges object.
    }
    \item{}{
      \code{rep(x, times, length.out, each)}:
      Repeats the values in \code{x} through one of the following conventions:
      \describe{
        \item{\code{times}}{Vector giving the number of times to repeat each
          element if of length \code{length(x)}, or to repeat the whole vector
          if of length 1.}
        \item{\code{length.out}}{Non-negative integer. The desired length of
          the output vector.}
        \item{\code{each}}{Non-negative integer.  Each element of \code{x} is
          repeated \code{each} times.}
      }
    }
    \item{}{
      \code{seqselect(x, start=NULL, end=NULL, width=NULL)}:
      Similar to \code{window}, except that multiple consecutive subsequences
      can be requested for concatenation. As such two of the three \code{start},
      \code{end}, and \code{width} arguments can be used to specify the
      consecutive subsequences. Alternatively, \code{start} can take a Ranges
      object or something that can be converted to a Ranges object like an
      integer vector, logical vector or logical Rle. If the concatenation of
      the consecutive subsequences is undesirable, consider using
      \code{\link{Views}}.
    }
    \item{}{
      \code{seqselect(x, start=NULL, end=NULL, width=NULL) <- value}:
      Similar to \code{window<-}, except that multiple consecutive subsequences
      can be replaced with a \code{value} whose length is a divisor of the
      number of elements it is replacing. As such two of the three \code{start},
      \code{end}, and \code{width} arguments can be used to specify the
      consecutive subsequences. Alternatively, \code{start} can take a Ranges
      object or something that can be converted to a Ranges object like an
      integer vector, logical vector or logical Rle.
    }
    \item{}{
      \code{subset(x, subset)}:
      Returns a new object of the same class as \code{x} made of the subset
      using logical vector \code{subset}, where missing values are taken as
      \code{FALSE}.
    }
    \item{}{
      \code{tail(x, n = 6L)}:
      If \code{n} is non-negative, returns the last n elements of the
      GRanges object.
      If \code{n} is negative, returns all but the first \code{abs(n)} elements
      of the GRanges object.
    }
    \item{}{
      \code{window(x, start = NA, end = NA, width = NA, frequency = NULL, delta = NULL, ...)}:
      Extracts the subsequence window from the GRanges object using:
      \describe{
        \item{\code{start}, \code{end}, \code{width}}{The start, end, or width
          of the window. Two of the three are required.}
        \item{\code{frequency}, \code{delta}}{Optional arguments that specify
          the sampling frequency and increment within the window.}
      }
      In general, this is more efficient than using \code{"["} operator.
    }
    \item{}{
      \code{window(x, start = NA, end = NA, width = NA, keepLength = TRUE) <- value}:
      Replaces the subsequence window specified on the left (i.e. the subsequence
      in \code{x} specified by \code{start}, \code{end} and \code{width})
      by \code{value}.
      \code{value} must either be of class \code{class(x)}, belong to a subclass
      of \code{class(x)}, be coercible to \code{class(x)}, or be \code{NULL}.
      If \code{keepLength} is \code{TRUE}, the elements of \code{value} are
      repeated to create a GRanges object with the same number of elements
      as the width of the subsequence window it is replacing.
      If \code{keepLength} is \code{FALSE}, this replacement method can modify
      the length of \code{x}, depending on how the length of the left
      subsequence window compares to the length of \code{value}.
    }
    \item{}{
      \code{x$name}, \code{x$name <- value}:
      Shortcuts for \code{mcols(x)$name} and \code{mcols(x)$name <- value},
      respectively. Provided as a convenience, for GRanges objects *only*,
      and as the result of strong popular demand.
      Note that those methods are not consistent with the other \code{$}
      and \code{$<-} methods in the IRanges/GenomicRanges infrastructure,
      and might confuse some users by making them believe that a GRanges
      object can be manipulated as a data.frame-like object.
      Therefore we recommend using them only interactively, and we discourage
      their use in scripts or packages. For the latter, use
      \code{mcols(x)$name} and \code{mcols(x)$name <- value}, instead
      of \code{x$name} and \code{x$name <- value}, respectively.
    }
  }
}

\author{P. Aboyoun}

\seealso{
  \link{GRangesList-class},
  \code{\link{seqinfo}},
  \link[IRanges]{Vector-class},
  \link[IRanges]{Ranges-class},
  \link[IRanges]{Rle-class},
  \link[IRanges]{DataFrame-class},
  \link{GenomicRanges-utils},
  \link{coverage-methods},
  \link{setops-methods},
  \link{findOverlaps-methods}
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
  gr

  ## Summarizing elements
  table(seqnames(gr))
  sum(width(gr))
  summary(mcols(gr)[,"score"])

  ## Renaming the underlying sequences
  seqlevels(gr)
  seqlevels(gr) <- sub("chr", "Chrom", seqlevels(gr))
  gr
  seqlevels(gr) <- sub("Chrom", "chr", seqlevels(gr)) # revert

  ## Combining objects
  gr2 <- GRanges(seqnames=Rle(c('chr1', 'chr2', 'chr3'), c(3, 3, 4)),
                 IRanges(1:10, width=5), strand='-',
                 score=101:110, GC = runif(10),
                 seqinfo=seqinfo)
  gr3 <- GRanges(seqnames=Rle(c('chr1', 'chr2', 'chr3'), c(3, 4, 3)),
                 IRanges(101:110, width=10), strand='-',
                 score=21:30,
                 seqinfo=seqinfo)
  some.gr <- c(gr, gr2)
 
  ## all.gr <- c(gr, gr2, gr3) ## (This would fail)
  all.gr <- c(gr, gr2, gr3, ignore.mcols=TRUE)
}

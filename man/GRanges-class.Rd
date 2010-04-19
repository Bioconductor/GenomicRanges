\name{GRanges-class}
\docType{class}

% Class
\alias{class:GRanges}
\alias{GRanges-class}
\alias{GRanges}

% Constructor-like function:
\alias{GRanges}

% Coercion methods:
\alias{coerce,RangedData,GRanges-method}
\alias{as.data.frame,GRanges-method}

% Accessor methods:
\alias{seqnames,GRanges-method}
\alias{seqnames<-,GRanges-method}
\alias{ranges,GRanges-method}
\alias{ranges<-,GRanges-method}
\alias{strand,GRanges-method}
\alias{strand<-,GRanges-method}
\alias{seqlengths,GRanges-method}
\alias{seqlengths<-,GRanges-method}
\alias{elementMetadata,GRanges-method}
\alias{elementMetadata<-,GRanges-method}
\alias{names,GRanges-method}
\alias{names<-,GRanges-method}

% Ranges methods:
\alias{start,GRanges-method}
\alias{start<-,GRanges-method}
\alias{end,GRanges-method}
\alias{end<-,GRanges-method}
\alias{width,GRanges-method}
\alias{width<-,GRanges-method}
\alias{flank,GRanges-method}
\alias{resize,GRanges-method}
\alias{shift,GRanges-method}
\alias{disjoin,GRanges-method}
\alias{gaps,GRanges-method}
\alias{range,GRanges-method}
\alias{reduce,GRanges-method}
\alias{coverage,GRanges-method}

% Sequence methods:
\alias{[,GRanges-method}
\alias{[<-,GRanges-method}
\alias{c,GRanges-method}
\alias{length,GRanges-method}
\alias{rev,GRanges-method}
\alias{seqselect,GRanges-method}
\alias{seqselect<-,GRanges-method}
\alias{split,GRanges-method}
\alias{window,GRanges-method}

% show method:
\alias{show,GRanges-method}


\title{GRanges objects}

\description{
  The GRanges class is a container for the genomic locations and
  their associated annotations.
}

\details{
  The GRanges class stores the sequences of genomic locations and associated
  annotations. Each element in the sequence is comprised of a sequence name,
  an interval, a \link{strand}, and optional element metadata (e.g. score, GC
  content, etc.). This information is stored in four slots:
  \describe{
    \item{\code{seqnames}}{a 'factor' \link[IRanges]{Rle} object
      containing the sequence names.}
    \item{\code{ranges}}{an \link[IRanges]{IRanges} object containing
      the ranges.}
    \item{\code{strand}}{a 'factor' \link[IRanges]{Rle} object containing
      the \link{strand} information.}
    \item{\code{elementMetadata}}{a \link[IRanges]{DataFrame} object
      containing the annotation columns. Columns cannot be named
      \code{"seqnames"}, \code{"ranges"}, \code{"strand"}, \code{"seqlengths"},
      \code{"start"}, \code{"end"}, \code{"width"}, or \code{"element"}.}
  }
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GRanges(seqnames = Rle(), ranges = IRanges(),
                    strand = Rle("*", length(seqnames)),
                    ...,
                    seqlengths =
                    structure(rep(NA_integer_, length(levels(seqnames))),
                              names = levels(seqnames)))}:
      Creates a GRanges object.
      \describe{
        \item{\code{seqnames}}{Rle object, character vector, or factor
          containing the sequence names.}
        \item{\code{ranges}}{IRanges object containing the ranges.}
        \item{\code{strand}}{Rle object, character vector, or factor
          containing the strand information.}
        \item{\code{seqlengths}}{a named integer vector containing the
          sequence lengths for each \code{level(seqnames)}.}
        \item{\code{\ldots}}{Optional annotation columns for the
          \code{elementMetadata} slot. These columns cannot be named
          \code{"start"}, \code{"end"}, \code{"width"}, or \code{"element"}.}
      }
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a GRanges object.

  \describe{
    \item{}{
      \code{as(from, "GRanges")}:
      Creates a GRanges object from a RangedData object.
    }
    \item{}{
      \code{as.data.frame(x, row.names = NULL, optional = FALSE)}:
      Creates a data.frame with columns \code{seqnames} (factor),
      \code{start} (integer), \code{end} (integer), \code{width} (integer),
      \code{strand} (factor), as well as the additional columns stored in
      \code{elementMetadata(x)}.
    }
  }
}

\section{Accessors}{
  In the following code snippets, \code{x} is a GRanges object.

  \describe{
    \item{}{
      \code{seqnames(x)}, \code{seqnames(x) <- value}:
      Gets or sets the sequence names. \code{value} can be an Rle object,
      character vector, or factor.
    }
    \item{}{
      \code{ranges(x)}, \code{ranges(x) <- value}:
      Gets or sets the ranges. \code{value} can be a Ranges object.
    }
    \item{}{
      \code{strand(x)}, \code{strand(x) <- value}:
      Gets or sets the strand. \code{value} can be an Rle object, character
      vector, or factor.
    }
    \item{}{
      \code{seqlengths(x)}, \code{seqlengths(x) <- value}:
      Gets or sets the seqlengths. \code{value} can be a named non-negative
      integer or numeric vector.
    }
    \item{}{
      \code{elementMetadata(x)}, \code{elementMetadata(x) <- value}:
      Gets or sets the optional data columns. \code{value} can be a DataFrame,
      data.frame object, or NULL.
    }
    \item{}{
      \code{values(x), values(x) <- value}:
      Alternative to \code{elementMetadata} functions.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
      Gets or sets the names of the elements.
    }
    \item{}{
      \code{length(x)}:
      Gets the number of elements.
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
    \item{}{
      \code{flank(x, width, start = TRUE, both = FALSE, use.names = TRUE)}:
      Returns a new GRanges object containing intervals of width \code{width}
      that flank the intervals in \code{x}.
      The \code{start} argument takes a logical indicating whether \code{x}
      should be flanked at the "start" (\code{TRUE}) or the "end"
      (\code{FALSE}), which for \code{strand(x) != "-"} is \code{start(x)} and
      \code{end(x)} respectively and for \code{strand(x) == "-"} is 
      code{end(x)} and \code{start(x)} respectively.
      The \code{both} argument takes a single logical value indicating whether
      the flanking region \code{width} positions extends \emph{into} the range.
      If \code{both = TRUE}, the resulting range thus straddles the end point,
      with \code{width} positions on either side.
    }
    \item{}{
      \code{resize(x, width, use.names = TRUE)}:
      Returns a new GRanges object containing intervals that have been
      resized to width \code{width} based on the \code{strand(x)} values.
      Elements where \code{strand(x) == "+"} are anchored at \code{start(x)},
      elements where \code{strand(x) == "-"} are anchored at the \code{end(x)},
      and elements where \code{strand = "*"} are anchored at
      \code{(end(x) - start(x))\%/\%2}. The \code{use.names} argument
      determines whether or not to keep the names on the ranges.
    }
    \item{}{
      \code{shift(x, shift, use.names = TRUE)}:
      Returns a new GRanges object containing intervals with start and
      end values that have been shifted by integer vector \code{shift}.
      The \code{use.names} argument determines whether or not to
      keep the names on the ranges.
    }
    \item{}{
      \code{disjoin(x)}:
      Returns a new GRanges object containing disjoint ranges for each
      distinct (seqname, strand) pairing. The names (\code{names(x)}) and
      the columns in \code{x} are dropped.
    }
    \item{}{
      \code{gaps(x, start = 1L, end = seqlengths(x))}:
      Returns a new GRanges object containing complemented ranges for each
      distinct (seqname, strand) pairing. The names (\code{names(x)}) and
      the columns in \code{x} are dropped. See
      \code{?\link[IRanges:RangesList-utils]{gaps}} for more information
      about range complements and for a description of the optional arguments.
    }
    \item{}{
      \code{range(x, ...)}:
      Returns a new GRanges object containing range bounds for each
      distinct (seqname, strand) pairing. The names (\code{names(x)}) and
      the columns in \code{x} are dropped.
    }
    \item{}{
      \code{reduce(x, drop.empty.ranges = FALSE, min.gapwidth = 1L)}:
      Returns a new GRanges object containing reduced ranges for each
      distinct (seqname, strand) pairing. The names (\code{names(x)}) and
      the columns in \code{x} are dropped. See
      \code{?\link[IRanges:RangesList-utils]{reduce}} for more information
      about range reduction and for a description of the optional arguments.
    }
    \item{}{
      \code{coverage(x, shift = list(0L), width = as.list(seqlengths(x)), weight = list(1L))}:
      Returns a named \link[IRanges]{RleList} object with one element
      ('integer' Rle) per unique sequence name representing how many
      times each position in the sequence is covered by the intervals in
      \code{x}. The \code{shift}, \code{width}, and \code{weight} arguments
      take list arguments, possibly named by the unique sequence names in
      \code{x}, whose elements are passed into the \code{coverage} method for
      \link{IRanges} object. See \code{?\link[IRanges]{coverage}} for more
      information on these optional arguments.
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
      \code{split(x, f = seq_len(length(x)), drop = FALSE)}:
      Splits \code{x} into a \link{GRangesList}, according to \code{f},
      dropping elements corresponding to unrepresented levels if \code{drop}
      is \code{TRUE}. Split factor \code{f} defaults to splitting each
      element of \code{x} into a separate element in the resulting
      \link{GRangesList} object.
    }
  }
}

\section{Subsetting}{
  In the code snippets below, \code{x} is a GRanges object.

  \describe{
    \item{}{
      \code{x[i, j]}, \code{x[i, j] <- value}:
      Gets or sets elements \code{i} with optional elementMetadata columns
      \code{elementMetadata(x)[,j]}, where \code{i} can be missing; an NA-free
      logical, numeric, or character vector; or a 'logical' Rle object.
    }
    \item{}{
      \code{x[i,j] <- value}:
      Replaces elements \code{i} and optional elementMetadata columns \code{j}
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
      \code{rev(x)}:
      Returns a new object of the same class as \code{x} made of the original
      elements in the reverse order.
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
      can be replaced by a \code{value} whose length is a divisor of the number
      of elements it is replacing. As such two of the three \code{start},
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
  }
}

\author{P. Aboyoun}

\seealso{
  \link{GRangesList},
  \link[IRanges]{Sequence},
  \link[IRanges]{Ranges},
  \link[IRanges]{Rle},
  \link[IRanges]{DataFrame}
}

\examples{
  gr <-
    GRanges(seqnames =
            Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
            ranges =
            IRanges(1:10, width = 10:1, names = head(letters,10)),
            strand =
            Rle(strand(c("-", "+", "*", "+", "-")),
                c(1, 2, 2, 3, 2)),
            score = 1:10,
            GC = seq(1, 0, length=10))
  gr

  # Summarizing elements
  table(seqnames(gr))
  sum(width(gr))
  summary(elementMetadata(gr)[,"score"]) # or values(gr)
  coverage(gr)

  # Changing sequence name
  unique(seqnames(gr))
  seqnames(gr) <- sub("chr", "Chrom", seqnames(gr))
  gr

  # Intra-interval operations
  flank(gr, 10)
  resize(gr, 10)
  shift(gr, 1)

  # Inter-interval operations
  disjoin(gr)
  gaps(gr, start = 1, end = 10)
  range(gr)
  reduce(gr)
}

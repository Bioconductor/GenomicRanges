\name{GRangesList-class}
\docType{class}

% Class
\alias{class:GRangesList}
\alias{GRangesList-class}
\alias{GRangesList}

% Constructor-like function:
\alias{GRangesList}

% Coercion methods:
\alias{coerce,RangedDataList,GRangesList-method}
\alias{as.data.frame,GRangesList-method}
\alias{unlist,GRangesList-method}

% Accessor methods:
\alias{seqnames,GRangesList-method}
\alias{seqnames<-,GRangesList-method}
\alias{ranges,GRangesList-method}
\alias{ranges<-,GRangesList-method}
\alias{strand,GRangesList-method}
\alias{strand<-,GRangesList-method}
\alias{seqlengths,GRangesList-method}
\alias{seqlengths<-,GRangesList-method}
\alias{elementMetadata,GRangesList-method}
\alias{elementMetadata<-,GRangesList-method}

% RangesList methods:
\alias{start,GRangesList-method}
\alias{start<-,GRangesList-method}
\alias{end,GRangesList-method}
\alias{end<-,GRangesList-method}
\alias{width,GRangesList-method}
\alias{width<-,GRangesList-method}
\alias{shift,GRangesList-method}
\alias{coverage,GRangesList-method}

% Sequence methods:
\alias{[,GRangesList-method}
\alias{[<-,GRangesList-method}

% show method:
\alias{show,GRangesList-method}


\title{GRangesList objects}

\description{
  The GRangesList class is a container for storing a collection
  of GRanges objects.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GRangesList(...)}:
      Creates a GRangesList object using GRanges objects
      supplied in \code{\ldots}.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a GRangesList object.

  \describe{
    \item{}{
      \code{as(from, "GRangesList")}:
      Creates a GRangesList object from a RangedDataList object.
    }
    \item{}{
      \code{as.data.frame(x, row.names = NULL, optional = FALSE)}:
      Creates a data.frame with columns \code{element} (character),
      \code{seqnames} (factor), \code{start} (integer),
      \code{end} (integer), \code{width} (integer), \code{strand} (factor),
      as well as the additional columns stored in
      \code{elementMetadata(unlist(x))}.
    }
    \item{}{\code{as.list(x, use.names = TRUE)}:
      Creates a list containing the elements of \code{x}.
    }
  }
}

\section{Accessors}{
  In the following code snippets, \code{x} is a GRanges object.
  
  \describe{
    \item{}{
      \code{seqnames(x)}, \code{seqnames(x) <- value}:
      Gets or sets the sequence names in the form of an RleList. \code{value}
      can be an RleList or CharacterList.
    }
    \item{}{
      \code{ranges(x)}, \code{ranges(x) <- value}:
      Gets or sets the ranges in the form of a CompressedIRangesList.
      \code{value} can be a RangesList object.
    }
    \item{}{
      \code{strand(x)}, \code{strand(x) <- value}:
      Gets or sets the strand in the form of an RleList. \code{value} can be
      an RleList or CharacterList object.
    }
    \item{}{
      \code{seqlengths(x)}, \code{seqlengths(x) <- value}:
      Gets or sets the seqlengths. \code{value} can be a named non-negative
      integer or numeric vector.
    }
  }
}

\section{List methods}{
  In the following code snippets, \code{x} is a GRangesList object.

  \describe{
    \item{}{
      \code{length(x)}:
      Gets the number of elements.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
      Gets or sets the names of the elements.
    }
    \item{}{
      \code{elementLengths(x)}:
      Gets the \code{length} of each of the elements.
    }
    \item{}{
      \code{isEmpty(x)}:
      Returns a logical indicating either if the GRangesList has no
      elements or if all its elements are empty.
    }
  }
}

\section{RangesList methods}{
  In the following code snippets, \code{x} is a GRangesList object.

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
      \code{shift(x, shift, use.names=TRUE)}:
      Returns a new GRangesList object containing intervals with start and
      end values that have been shifted by integer vector \code{shift}.
      The \code{use.names} argument determines whether or not to
      keep the names on the ranges.
    }
    \item{}{
      \code{coverage(x, shift = list(0L), width = list(NULL), weight = list(1L))}:
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

\section{Combining}{
  In the code snippets below, \code{x} is a GRangesList object.

  \describe{
    \item{}{\code{append(x, values, after = length(x))}:
      Inserts the \code{values} into \code{x} at the position given by
      \code{after}, where \code{x} and \code{values} are of the same
      class.
    }
    \item{}{
      \code{c(x, ...)}:
      Combines \code{x} and the GRangesList objects in \code{...}
      together. Any object in \code{...} must belong to the same class
      as \code{x}, or to one of its subclasses, or must be \code{NULL}.
      The result is an object of the same class as \code{x}.
    }
    \item{}{\code{unlist(x, recursive = TRUE, use.names = TRUE)}:
      Concatenates the elements of \code{x} into a single GRanges
      object.
    }
  }
}

\section{Subsetting}{
  In the following code snippets, \code{x} is a GRangesList object.

  \describe{
    \item{}{
      \code{x[i, j]}, \code{x[i, j] <- value}:
      Gets or sets elements \code{i} with optional values columns
      \code{values(x)[,j]}, where \code{i} can be missing; an NA-free
      logical, numeric, or character vector; a 'logical' Rle object, or
      an AtomicList object.
    }
    \item{}{
      \code{x[[i]]}, \code{x[[i]] <- value}:
      Gets or sets element \code{i}, where \code{i} is a numeric or character
      vector of length 1.
    }
    \item{}{
      \code{x$name}, \code{x$name <- value}:
      Gets or sets element \code{name}, where \code{name} is a name or character
      vector of length 1.
    }
    \item{}{
      \code{head(x, n = 6L)}:
      If \code{n} is non-negative, returns the first n elements of the
      GRangesList object.
      If \code{n} is negative, returns all but the last \code{abs(n)} elements
      of the GRangesList object.
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

\section{Looping}{
  In the code snippets below, \code{x} is a GRangesList object.

  \describe{
    \item{}{
      \code{endoapply(X, FUN, ...)}:
      Similar to \code{\link[base]{lapply}}, but performs an endomorphism,
      i.e. returns an object of \code{class(X)}.
    }
    \item{}{
      \code{lapply(X, FUN, ...)}:
      Like the standard \code{\link[base]{lapply}} function defined in the
      base package, the \code{lapply} method for GRangesList objects
      returns a list of the same length as \code{X}, with each element being
      the result of applying \code{FUN} to the corresponding element of
      \code{X}.
    }
    \item{}{
      \code{Map(f, ...)}:
      Applies a function to the corresponding elements of given
      GRangesList objects.
    }
    \item{}{
      \code{mapply(FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE)}:
      Like the standard \code{\link[base]{mapply}} function defined in the
      base package, the \code{mapply} method for GRangesList objects is a
      multivariate version of \code{sapply}.
    }
    \item{}{
      \code{mendoapply(FUN, ..., MoreArgs = NULL)}:
      Similar to \code{\link[base]{mapply}}, but performs an endomorphism
      across multiple objects, i.e. returns an object of
      \code{class(list(...)[[1]])}.
    }
    \item{}{
      \code{Reduce(f, x, init, right = FALSE, accumulate = FALSE)}:
      Uses a binary function to successively combine the elements of \code{x}
      and a possibly given initial value.
      \describe{
        \item{f}{A binary argument function.}
        \item{init}{An R object of the same kind as the elements of \code{x}.}
        \item{right}{A logical indicating whether to proceed from left to right
          (default) or from right to left.}
        \item{nomatch}{The value to be returned in the case when "no match" (no
          element satisfying the predicate) is found.}
      }
    }
    \item{}{
      \code{sapply(X, FUN, ..., simplify=TRUE, USE.NAMES=TRUE)}:
      Like the standard \code{\link[base:lapply]{sapply}} function defined in
      the base package, the \code{sapply} method for GRangesList objects
      is a user-friendly version of \code{lapply} by default returning a vector
      or matrix if appropriate.
    }
  }
}

\author{P. Aboyoun}

\seealso{
  \link{GRanges},
  \link[IRanges]{Sequence},
  \link[IRanges]{RangesList},
  \link[IRanges]{RleList},
  \link[IRanges]{DataFrameList}
}

\examples{
  ## Construction using GRangesList
  gr1 <-
    GRanges(seqnames = "chr2", ranges = IRanges(3, 6),
            strand = "+", score = 5L, GC = 0.45)
  gr2 <-
    GRanges(seqnames = c("chr1", "chr1"),
            ranges = IRanges(c(7,13), width = 3),
            strand = c("+", "-"), score = 3:4, GC = c(0.3, 0.5))
  gr3 <-
    GRanges(seqnames = c("chr1", "chr2"),
            ranges = IRanges(c(1, 4), c(3, 9)),
            strand = c("-", "-"), score = c(6L, 2L), GC = c(0.4, 0.1))
  grl <- GRangesList("gr1" = gr1, "gr2" = gr2, "gr3" = gr3)
  grl

  # Summarizing elements
  elementLengths(grl)
  table(seqnames(grl))
  coverage(grl)

  # Extracting subsets
  grl[seqnames(grl) == "chr1", ]
  grl[seqnames(grl) == "chr1" & strand(grl) == "+", ]

  # Changing sequence name
  unique(seqnames(grl))
  seqnames(grl) <- sub("chr", "Chrom", seqnames(grl))
  grl
}

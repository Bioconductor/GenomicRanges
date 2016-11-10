\name{GRangesList-class}
\docType{class}

% Class
\alias{class:GRangesList}
\alias{GRangesList-class}
\alias{GRangesList}

% Constructors:
\alias{GRangesList}
\alias{makeGRangesListFromFeatureFragments}
\alias{updateObject,GRangesList-method}

% Accessor methods:
\alias{seqnames,GRangesList-method}
\alias{seqnames<-,GRangesList-method}
\alias{ranges,GRangesList-method}
\alias{ranges<-,GRangesList-method}
\alias{start,GRangesList-method}
\alias{start<-,GRangesList-method}
\alias{end,GRangesList-method}
\alias{end<-,GRangesList-method}
\alias{width,GRangesList-method}
\alias{width<-,GRangesList-method}
\alias{strand,GRangesList-method}
\alias{strand<-,GRangesList,ANY-method}
\alias{strand<-,GRangesList,character-method}
\alias{elementMetadata,GRangesList-method}
\alias{elementMetadata<-,GRangesList-method}
\alias{seqinfo,GRangesList-method}
\alias{seqinfo<-,GRangesList-method}
\alias{score,GRangesList-method}
\alias{score<-,GRangesList-method}

% Coercion:
\alias{coerce,GRangesList,CompressedIRangesList-method}
\alias{coerce,GRangesList,IRangesList-method}
\alias{coerce,GRangesList,RangesList-method}
\alias{coerce,GRanges,GRangesList-method}

% Subsetting:
\alias{[,GRangesList,ANY-method}
\alias{[<-,GRangesList-method}
\alias{[<-,GRangesList,ANY-method}
\alias{[<-,GRangesList,ANY,ANY,ANY-method}
\alias{[[<-,GRangesList-method}
\alias{[[<-,GRangesList,ANY-method}
\alias{[[<-,GRangesList,ANY,ANY-method}

% extractList() and family:
\alias{relistToClass,GRanges-method}

% show method:
\alias{show,GRangesList-method}


\title{GRangesList objects}

\description{
  The GRangesList class is a container for storing a collection of
  GRanges objects. It is derived from GenomicRangesList.
}

\section{Constructors}{
  \describe{
    \item{}{
      \code{GRangesList(...)}:
      Creates a GRangesList object using GRanges objects
      supplied in \code{\ldots}.
    }
    \item{}{
      \code{makeGRangesListFromFeatureFragments(seqnames=Rle(factor()),
            fragmentStarts=list(), fragmentEnds=list(), fragmentWidths=list(),
            strand=character(0), sep=",")}:
      Constructs a GRangesList object from a list of fragmented features.
      See the Examples section below.
    }
  }
}

\section{Accessors}{
  In the following code snippets, \code{x} is a GRanges object.
  
  \describe{
    \item{}{
      \code{length(x)}:
      Get the number of list elements.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
      Get or set the names on \code{x}.
    }
    \item{}{
      \code{elementNROWS(x)}:
      Get a vector of the \code{length} of each of the list element.
    }
    \item{}{
      \code{isEmpty(x)}:
      Returns a logical indicating either if the GRangesList has no
      elements or if all its elements are empty.
    }
    \item{}{
      \code{seqnames(x)}, \code{seqnames(x) <- value}:
      Get or set the sequence names in the form of an RleList. \code{value}
      can be an RleList or CharacterList object.
    }
    \item{}{
      \code{ranges(x, use.mcols=FALSE)}, \code{ranges(x) <- value}:
      Get or set the ranges in the form of a CompressedIRangesList.
      \code{value} can be a RangesList object.
    }
    \item{}{
      \code{start(x)}, \code{start(x) <- value}:
      Get or set \code{start(ranges(x))}.
    }
    \item{}{
      \code{end(x)}, \code{end(x) <- value}:
      Get or set \code{end(ranges(x))}.
    }
    \item{}{
      \code{width(x)}, \code{width(x) <- value}:
      Get or set \code{width(ranges(x))}.
    }
    \item{}{
      \code{strand(x)}, \code{strand(x) <- value}:
      Get or set the strand in the form of an RleList. \code{value} can be
      an RleList, CharacterList or single character. \code{value} as a
      single character converts all ranges in \code{x} to the same
      \code{value}; for selective strand conversion (i.e., mixed \dQuote{+}
      and \dQuote{-}) use RleList or CharacterList. 
    }
    \item{}{
      \code{mcols(x, use.names=FALSE)}, \code{mcols(x) <- value}:
      Get or set the metadata columns.
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
      Get or set the information about the underlying sequences.
      \code{value} must be a \link{Seqinfo} object.
    }
    \item{}{
      \code{seqlevels(x)},
      \code{seqlevels(x, pruning.mode=c("error", "coarse", "fine", "tidy")) <- value}:
      Get or set the sequence levels.
      \code{seqlevels(x)} is equivalent to \code{seqlevels(seqinfo(x))}
      or to \code{levels(seqnames(x))}, those 2 expressions being
      guaranteed to return identical character vectors on a GRangesList
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
      \code{seqlevelsStyle(x)}, \code{seqlevelsStyle(x) <- value}:
      Get or set the seqname style for \code{x}.
      See the \link[GenomeInfoDb]{seqlevelsStyle} generic getter and setter
      in the \pkg{GenomeInfoDb} package for more information.
    }
    \item{}{
      \code{score(x), score(x) <- value}: Get or set the \dQuote{score}
      metadata column.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a GRangesList object.

  \describe{
    \item{}{
      \code{as.data.frame(x, row.names = NULL, optional = FALSE,
         ..., value.name = "value", use.outer.mcols = FALSE, 
         group_name.as.factor = FALSE)}: 
      Coerces \code{x} to a \code{data.frame}. See as.data.frame on the 
      \code{List} man page for details (?\code{List}). 
    }
    \item{}{\code{as.list(x, use.names = TRUE)}:
      Creates a list containing the elements of \code{x}.
    }
    \item{}{\code{as(x, "IRangesList")}:
      Turns \code{x} into an \link[IRanges]{IRangesList} object.
    }
  }
}

\section{Subsetting}{
  In the following code snippets, \code{x} is a GRangesList object.

  \describe{
    \item{}{
      \code{x[i, j]}, \code{x[i, j] <- value}:
      Get or set elements \code{i} with optional metadata columns
      \code{mcols(x)[,j]}, where \code{i} can be missing; an NA-free
      logical, numeric, or character vector; a 'logical' Rle object, or
      an AtomicList object.
    }
    \item{}{
      \code{x[[i]]}, \code{x[[i]] <- value}:
      Get or set element \code{i}, where \code{i} is a numeric or character
      vector of length 1.
    }
    \item{}{
      \code{x$name}, \code{x$name <- value}:
      Get or set element \code{name}, where \code{name} is a name or character
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
  }
}

\section{Combining}{
  In the code snippets below, \code{x} is a GRangesList object.

  \describe{
    \item{}{
      \code{c(x, ...)}:
      Combines \code{x} and the GRangesList objects in \code{...}
      together. Any object in \code{...} must belong to the same class
      as \code{x}, or to one of its subclasses, or must be \code{NULL}.
      The result is an object of the same class as \code{x}.
    }
    \item{}{
      \code{append(x, values, after = length(x))}:
      Inserts the \code{values} into \code{x} at the position given by
      \code{after}, where \code{x} and \code{values} are of the same
      class.
    }
    \item{}{\code{unlist(x, recursive = TRUE, use.names = TRUE)}:
      Concatenates the elements of \code{x} into a single GRanges
      object.
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

\author{P. Aboyoun & H. Pagès}

\seealso{
  \link{GRanges-class},
  \code{\link[GenomeInfoDb]{seqinfo}},
  \link[S4Vectors]{Vector-class},
  \link[IRanges]{RangesList-class},
  \link[IRanges]{RleList-class},
  \link[IRanges]{DataFrameList-class},
  \link[GenomicRanges]{intra-range-methods},
  \link[GenomicRanges]{inter-range-methods},
  \link{coverage-methods},
  \link{setops-methods},
  \link{findOverlaps-methods}
}

\examples{
## Construction with GRangesList():
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

## Summarizing elements:
elementNROWS(grl)
table(seqnames(grl))

## Extracting subsets:
grl[seqnames(grl) == "chr1", ]
grl[seqnames(grl) == "chr1" & strand(grl) == "+", ]

## Renaming the underlying sequences:
seqlevels(grl)
seqlevels(grl) <- sub("chr", "Chrom", seqlevels(grl))
grl

## Coerce to IRangesList (seqnames and strand information is lost):
as(grl, "IRangesList")

## isDisjoint():
isDisjoint(grl)

## disjoin():
disjoin(grl)  # metadata columns and order NOT preserved

## Construction with makeGRangesListFromFeatureFragments():
filepath <- system.file("extdata", "feature_frags.txt",
                        package="GenomicRanges")
featfrags <- read.table(filepath, header=TRUE, stringsAsFactors=FALSE) 
grl2 <- with(featfrags,
             makeGRangesListFromFeatureFragments(seqnames=targetName,
                                                 fragmentStarts=targetStart,
                                                 fragmentWidths=blockSizes,
                                                 strand=strand))
names(grl2) <- featfrags$RefSeqID
grl2
}

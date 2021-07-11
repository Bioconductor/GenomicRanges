\name{GRanges-class}
\docType{class}

\alias{class:IRanges_OR_IPos}
\alias{IRanges_OR_IPos-class}
\alias{IRanges_OR_IPos}

\alias{class:GenomicRanges}
\alias{GenomicRanges-class}
\alias{GenomicRanges}

\alias{class:GenomicPos}
\alias{GenomicPos-class}
\alias{GenomicPos}

\alias{GenomicRanges_OR_missing-class}

\alias{class:GRanges}
\alias{GRanges-class}
\alias{GRanges}

\alias{parallel_slot_names,GRanges-method}

\alias{updateObject,GRanges-method}

\alias{update,GRanges-method}

% Constructors:
\alias{GRanges}

% Accessors:
\alias{length,GenomicRanges-method}
\alias{seqnames,GRanges-method}
\alias{seqnames<-,GenomicRanges-method}
\alias{ranges,GRanges-method}
\alias{ranges<-,GenomicRanges-method}
\alias{start,GenomicRanges-method}
\alias{start<-,GenomicRanges-method}
\alias{end<-,GenomicRanges-method}
\alias{width,GenomicRanges-method}
\alias{width<-,GenomicRanges-method}
\alias{strand,GRanges-method}
\alias{strand<-,GenomicRanges,ANY-method}
\alias{names,GenomicRanges-method}
\alias{names<-,GenomicRanges-method}
\alias{$,GenomicRanges-method}
\alias{$<-,GenomicRanges-method}
\alias{seqinfo,GRanges-method}
\alias{seqinfo,List-method}
\alias{seqinfo<-,GenomicRanges-method}
\alias{seqinfo<-,List-method}
\alias{score,GenomicRanges-method}
\alias{score<-,GenomicRanges-method}
\alias{granges,GenomicRanges-method}

% Coercion:
\alias{coerce,GenomicRanges,GRanges-method}
\alias{as.character,GenomicRanges-method}
\alias{as.factor,GenomicRanges-method}
\alias{as.data.frame,GenomicRanges-method}
\alias{coerce,GenomicRanges,CompressedIRangesList-method}
\alias{coerce,GenomicRanges,IRangesList-method}
\alias{coerce,GenomicRanges,IntegerRangesList-method}
\alias{coerce,character,GRanges-method}
\alias{coerce,factor,GRanges-method}
\alias{coerce,IntegerRangesList,GRanges-method}
\alias{coerce,Seqinfo,GRanges-method}
\alias{coerce,Seqinfo,IntegerRangesList-method}
\alias{coerce,GenomicRanges,Grouping-method}
\alias{coerce,ANY,GenomicRanges-method}

% Subsetting:
\alias{[,list_OR_List,GenomicRanges-method}
\alias{getListElement,GenomicRanges-method}

% Display:
\alias{summary.GenomicRanges}
\alias{summary,GenomicRanges-method}
\alias{show,GenomicRanges-method}

% Concatenation:
\alias{bindROWS,GenomicRanges-method}

\title{GRanges objects}

\description{
  The GRanges class is a container for the genomic locations and
  their associated annotations.
}

\details{
  GRanges is a vector of genomic locations and associated
  annotations. Each element in the vector is comprised of a sequence name,
  an interval, a \link{strand}, and optional metadata columns (e.g. score, GC
  content, etc.). This information is stored in four components:
  \describe{
    \item{\code{seqnames}}{a 'factor' \link[S4Vectors]{Rle} object
      containing the sequence names.}
    \item{\code{ranges}}{an \link[IRanges]{IRanges} object containing
      the ranges.}
    \item{\code{strand}}{a 'factor' \link[S4Vectors]{Rle} object containing
      the \link{strand} information.}
    \item{\code{mcols}}{a \link[S4Vectors]{DataFrame} object
      containing the metadata columns. Columns cannot be named
      \code{"seqnames"}, \code{"ranges"}, \code{"strand"},
      \code{"seqlevels"}, \code{"seqlengths"}, \code{"isCircular"},
      \code{"start"}, \code{"end"}, \code{"width"}, or \code{"element"}.}
    \item{\code{seqinfo}}{a \link{Seqinfo} object containing information
      about the set of genomic sequences present in the GRanges object.}
  }
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GRanges(seqnames=NULL, ranges=NULL, strand=NULL,
                    ..., seqinfo=NULL, seqlengths=NULL)}:
      Creates a GRanges object.
      \describe{
        \item{\code{seqnames}}{
          \code{NULL}, or an \link[S4Vectors]{Rle} object, character vector,
          or factor containing the sequence names.
        }
        \item{\code{ranges}}{
          \code{NULL}, or an \link[IRanges]{IRanges} object containing the
          ranges.
        }
        \item{\code{strand}}{
          \code{NULL}, or an \link[S4Vectors]{Rle} object, character vector,
          or factor containing the strand information.
        }
        \item{\code{...}}{
          Metadata columns to set on the GRanges object. All the metadata
          columns must be vector-like objects of the same length as the object
          to construct. They cannot be named \code{"start"}, \code{"end"},
          \code{"width"}, or \code{"element"}.
        }
        \item{\code{seqinfo}}{
          Either \code{NULL}, or a \link[GenomeInfoDb]{Seqinfo} object,
          or a character vector of unique sequence names (a.k.a.
          \emph{seqlevels}), or a named numeric vector of sequence lengths.
          When not \code{NULL}, \code{seqinfo} must be compatible with the
          sequence names in \code{seqnames}, that is, it must have one entry
          for each unique sequence name in \code{seqnames}. Note that it can
          have additional entries i.e. entries for seqlevels not present
          in \code{seqnames}.
        }
        \item{\code{seqlengths}}{
          \code{NULL}, or an integer vector named with \code{levels(seqnames)}
          and containing the lengths (or NA) for each level in
          \code{levels(seqnames)}.
        }
      }
      If \code{ranges} is not supplied and/or NULL then the constructor
      proceeds in 2 steps:
      \enumerate{
        \item An initial GRanges object is created with
              \code{as(seqnames, "GRanges")}.
        \item Then this GRanges object is updated according to whatever
              non-NULL remaining arguments were passed to the call to
              \code{GRanges()}.
      }
      As a consequence of this behavior, \code{GRanges(x)} is equivalent to
      \code{as(x, "GRanges")}.
    }
  }
}

\section{Accessors}{
  In the following code snippets, \code{x} is a GRanges object.

  \describe{
    \item{}{
      \code{length(x)}:
      Get the number of elements.
    }
    \item{}{
      \code{seqnames(x)}, \code{seqnames(x) <- value}:
      Get or set the sequence names.
      \code{value} can be an \link[S4Vectors]{Rle} object, a character vector,
      or a factor.
    }
    \item{}{
      \code{ranges(x)}, \code{ranges(x) <- value}:
      Get or set the ranges. \code{value} can be an
      \link[IRanges]{IntegerRanges} object.
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
      Get or set the strand. \code{value} can be an Rle object, character
      vector, or factor.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
      Get or set the names of the elements.
    }
    \item{}{
      \code{mcols(x, use.names=FALSE)}, \code{mcols(x) <- value}:
      Get or set the metadata columns.
      If \code{use.names=TRUE} and the metadata columns are not \code{NULL},
      then the names of \code{x} are propagated as the row names of the
      returned \link{DataFrame} object.
      When setting the metadata columns, the supplied value must be \code{NULL}
      or a data-frame-like object (i.e. \link{DataFrame} or data.frame)
      holding element-wise metadata.
    }
    \item{}{
       \code{elementMetadata(x)}, \code{elementMetadata(x) <- value},
       \code{values(x)}, \code{values(x) <- value}:
       Alternatives to \code{mcols} functions. Their use is discouraged.
    }
    \item{}{
      \code{seqinfo(x)}, \code{seqinfo(x) <- value}:
      Get or set the information about the underlying sequences.
      \code{value} must be a \link[GenomeInfoDb]{Seqinfo} object.
    }
    \item{}{
      \code{seqlevels(x)},
      \code{seqlevels(x, pruning.mode=c("error", "coarse", "fine", "tidy")) <- value}:
      Get or set the sequence levels.
      \code{seqlevels(x)} is equivalent to \code{seqlevels(seqinfo(x))}
      or to \code{levels(seqnames(x))}, those 2 expressions being
      guaranteed to return identical character vectors on a GRanges object.
      \code{value} must be a character vector with no NAs.
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
      column from the element metadata.
    }
    \item{}{
      \code{granges(x, use.names=FALSE, use.mcols=FALSE)}: Squeeze the genomic
      ranges out of \link{GenomicRanges} object \code{x} and return them in a
      GRanges object \emph{parallel} to \code{x} (i.e. same length as \code{x}).
      If \code{use.mcols} is \code{TRUE}, the metadata columns are propagated.
      If \code{x} is a \link{GenomicRanges} derivative with \emph{extra column
      slots}, these will be propagated as metadata columns on the returned
      GRanges object.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a GRanges object.

  \describe{
    \item{}{
      \code{as(from, "GRanges")}: Creates a GRanges object from a character
      vector, a factor, or \link[IRanges]{IntegerRangesList} object.

      When \code{from} is a character vector (or a factor), each element
      in it must represent a genomic range in format \code{chr1:2501-2800}
      (unstranded range) or \code{chr1:2501-2800:+} (stranded range).
      \code{..} is also supported as a separator between the start and end
      positions. Strand can be \code{+}, \code{-}, \code{*}, or missing.
      The names on \code{from} are propagated to the returned GRanges object.
      See \code{as.character()} and \code{as.factor()} below for the
      reverse transformations.

      Coercing a data.frame or DataFrame into a GRanges object is also
      supported. See \code{\link{makeGRangesFromDataFrame}} for the details.
    }
    \item{}{
      \code{as(from, "IntegerRangesList")}:
      Creates a \link[IRanges]{IntegerRangesList} object from a GRanges
      object. The \code{strand} and metadata columns become \emph{inner}
      metadata columns (i.e. metadata columns on the ranges).
      The \code{seqlengths(from)}, \code{isCircular(from)}, and
      \code{genome(from)} vectors become the metadata columns.
    }
    \item{}{
      \code{as.character(x, ignore.strand=FALSE)}:
      Turn GRanges object \code{x} into a character vector where each
      range in \code{x} is represented by a string in format
      \code{chr1:2501-2800:+}. If \code{ignore.strand} is TRUE or if
      \emph{all} the ranges in \code{x} are unstranded (i.e. their strand
      is set to \code{*}), then all the strings in the output are in
      format \code{chr1:2501-2800}.

      The names on \code{x} are propagated to the returned character vector.
      Its metadata (\code{metadata(x)}) and metadata columns (\code{mcols(x)})
      are ignored.

      See \code{as(from, "GRanges")} above for the reverse transformation.
    }
    \item{}{
      \code{as.factor(x)}:
      Equivalent to
\preformatted{  factor(as.character(x), levels=as.character(sort(unique(x))))
}
      See \code{as(from, "GRanges")} above for the reverse transformation.

      Note that \code{table(x)} is supported on a GRanges object. It is
      equivalent to, but much faster than, \code{table(as.factor(x))}.
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
    \item{}{
      \code{as(from, "Grouping")}: Creates a
      \code{\link[IRanges]{ManyToOneGrouping}} object that groups
      \code{from} by seqname, strand, start and end (same as the default
      sort order). This makes it convenient, for example, to aggregate a
      GenomicRanges object by range.
    }
  }

  In the code snippets below, \code{x} is a \link[GenomeInfoDb]{Seqinfo}
  object.

  \describe{
    \item{}{
      \code{as(x, "GRanges")}, \code{as(x, "GenomicRanges")},
      \code{as(x, "IntegerRangesList")}: Turns \link[GenomeInfoDb]{Seqinfo}
      object \code{x} (with no \code{NA} lengths) into a GRanges or
      \link[IRanges]{IntegerRangesList}.
    }
  }
}

\section{Subsetting}{
  In the code snippets below, \code{x} is a GRanges object.
  \describe{
    \item{}{
      \code{x[i]}:
      Return a new GRanges object made of the elements selected by \code{i}.
    }
    \item{}{
      \code{x[i, j]}:
      Like the above, but allow the user to conveniently subset the metadata
      columns thru \code{j}.
    }
    \item{}{
      \code{x[i] <- value}:
      Replacement version of \code{x[i]}.
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

  See \code{?`\link[S4Vectors]{[}`} in the \pkg{S4Vectors} package for more
  information about subsetting Vector derivatives and for an important note
  about the \code{x[i, j]} form.

  Note that a GRanges object can be used as a subscript to subset a
  list-like object that has names on it. In that case, the names on the
  list-like object are interpreted as sequence names.
  In the code snippets below, \code{x} is a list or \link{List} object with
  names on it, and the subscript \code{gr} is a GRanges object with all its
  seqnames being valid \code{x} names.

  \describe{
    \item{}{
      \code{x[gr]}:
      Return an object of the same class as \code{x} and \emph{parallel}
      to \code{gr}. More precisely, it's conceptually doing:
\preformatted{  lapply(gr, function(gr1) x[[seqnames(gr1)]][ranges(gr1)])
}
    }
  }
}

\section{Concatenation}{
  \describe{
    \item{}{
      \code{c(x, ..., ignore.mcols=FALSE)}:
      Concatenate GRanges object \code{x} and the GRanges objects in
      \code{...} together.
      See \code{?\link[S4Vectors]{c}} in the \pkg{S4Vectors} package
      for more information about concatenating Vector derivatives.
    }
  }
}

\section{Splitting}{
  \describe{
    \item{}{
      \code{split(x, f, drop=FALSE)}:
      Splits GRanges object \code{x} according to \code{f} to create a
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

\section{Displaying}{
  In the code snippets below, \code{x} is a GRanges object.

  \describe{
    \item{}{
      \code{show(x)}:
      By default the \code{show} method displays 5 head and 5 tail
      elements. This can be changed by setting the global options
      \code{showHeadLines} and \code{showTailLines}. If the object
      length is less than (or equal to) the sum of these 2 options
      plus 1, then the full object is displayed.
      Note that these options also affect the display of
      \link[GenomicAlignments]{GAlignments} and
      \link[GenomicAlignments]{GAlignmentPairs} objects (defined in
      the \pkg{GenomicAlignments} package), as well as other objects
      defined in the \pkg{IRanges} and \pkg{Biostrings} packages (e.g.
      \link[IRanges]{IRanges} and \link[Biostrings]{DNAStringSet} objects).
    }
  }
}

\author{P. Aboyoun and H. Pagès}

\seealso{
  \itemize{
    \item The \link[IRanges]{IRanges} class in the \pkg{IRanges}
          package for storing a set of \emph{integer ranges}.

    \item The \link{GPos} class for representing a set of \emph{genomic
          positions} (i.e. \emph{genomic ranges} of width 1, a.k.a.
          \emph{genomic loci}).

    \item \code{\link{makeGRangesFromDataFrame}} for making a GRanges object
          from a data.frame or \link[S4Vectors]{DataFrame} object.

    \item \link[GenomeInfoDb]{Seqinfo} objects and the
          \code{\link[GenomeInfoDb]{seqinfo}} accessor and family in the
          \pkg{GenomeInfoDb} package for accessing/modifying information
          about the underlying sequences of a \link{GenomicRanges} derivative.

    \item \link{GenomicRanges-comparison} for comparing and ordering genomic
          ranges and/or positions.

    \item \link[GenomicRanges]{findOverlaps-methods} for finding overlapping
          genomic ranges and/or positions.

    \item \link[GenomicRanges]{intra-range-methods} and
          \link[GenomicRanges]{inter-range-methods}
          for \emph{intra range} and \emph{inter range} transformations
          of \link{GenomicRanges} derivatives.

    \item \link[GenomicRanges]{coverage-methods} for computing the coverage
          of a set of genomic ranges and/or positions.

    \item \link[GenomicRanges]{setops-methods} for set operations on GRanges
          objects.

    \item \link[GenomicRanges]{nearest-methods} for finding the nearest
          genomic range/position neighbor.

    \item \code{\link{absoluteRanges}} for transforming genomic ranges into
          \emph{absolute} ranges (i.e. into ranges on the sequence obtained
          by virtually concatenating all the sequences in a genome).

    \item \code{\link{tileGenome}} for putting tiles on a genome.

    \item \link{genomicvars} for manipulating genomic variables.

    \item \link{GRangesList} objects.

    \item \link[S4Vectors]{Vector}, \link[S4Vectors]{Rle}, and
          \link[S4Vectors]{DataFrame} objects in the \pkg{S4Vectors} package.
  }
}

\examples{
showClass("GRanges")  # shows the known subclasses

## ---------------------------------------------------------------------
## CONSTRUCTION
## ---------------------------------------------------------------------
## Specifying the bare minimum i.e. seqnames and ranges only. The
## GRanges object will have no names, no strand information, and no
## metadata columns:
gr0 <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
               IRanges(1:10, width=10:1))
gr0

## Specifying names, strand, metadata columns. They can be set on an
## existing object:
names(gr0) <- head(letters, 10)
strand(gr0) <- Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2))
mcols(gr0)$score <- 1:10
mcols(gr0)$GC <- seq(1, 0, length=10)
gr0

## ... or specified at construction time:
gr <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
              IRanges(1:10, width=10:1, names=head(letters, 10)),
              Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
              score=1:10, GC=seq(1, 0, length=10))
stopifnot(identical(gr0, gr))

## Specifying the seqinfo. It can be set on an existing object:
seqinfo <- Seqinfo(paste0("chr", 1:3), c(1000, 2000, 1500), NA, "mock1")
seqinfo(gr0) <- merge(seqinfo(gr0), seqinfo)
seqlevels(gr0) <- seqlevels(seqinfo)

## ... or specified at construction time:
gr <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
              IRanges(1:10, width=10:1, names=head(letters, 10)),
              Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
              score=1:10, GC=seq(1, 0, length=10),
              seqinfo=seqinfo)
stopifnot(identical(gr0, gr))

## ---------------------------------------------------------------------
## COERCION
## ---------------------------------------------------------------------
## From GRanges:
as.character(gr)
as.factor(gr)
as.data.frame(gr)

## From character to GRanges:
x1 <- "chr2:56-125"
as(x1, "GRanges")
as(rep(x1, 4), "GRanges")
x2 <- c(A=x1, B="chr1:25-30:-")
as(x2, "GRanges")

## From data.frame to GRanges:
df <- data.frame(chrom="chr2", start=11:15, end=20:24)
gr3 <- as(df, "GRanges")

## Alternatively, coercion to GRanges can be done by just calling the
## GRanges() constructor on the object to coerce:
gr1 <- GRanges(x1)  # same as as(x1, "GRanges")
gr2 <- GRanges(x2)  # same as as(x2, "GRanges")
gr3 <- GRanges(df)  # same as as(df, "GRanges")

## Sanity checks:
stopifnot(identical(as(x1, "GRanges"), gr1))
stopifnot(identical(as(x2, "GRanges"), gr2))
stopifnot(identical(as(df, "GRanges"), gr3))

## ---------------------------------------------------------------------
## SUMMARIZING ELEMENTS
## ---------------------------------------------------------------------
table(seqnames(gr))
table(strand(gr))
sum(width(gr))
table(gr)
summary(mcols(gr)[,"score"])

## The number of lines displayed in the 'show' method are controlled
## with two global options:
longGR  <- sample(gr, 25, replace=TRUE)
longGR
options(showHeadLines=7)
options(showTailLines=2)
longGR

## Revert to default values
options(showHeadLines=NULL)
options(showTailLines=NULL)

## ---------------------------------------------------------------------
## INVERTING THE STRAND
## ---------------------------------------------------------------------
invertStrand(gr)

## ---------------------------------------------------------------------
## RENAMING THE UNDERLYING SEQUENCES
## ---------------------------------------------------------------------
seqlevels(gr)
seqlevels(gr) <- sub("chr", "Chrom", seqlevels(gr))
gr
seqlevels(gr) <- sub("Chrom", "chr", seqlevels(gr)) # revert

## ---------------------------------------------------------------------
## COMBINING OBJECTS
## ---------------------------------------------------------------------
gr2 <- GRanges(seqnames=Rle(c('chr1', 'chr2', 'chr3'), c(3, 3, 4)),
               IRanges(1:10, width=5),
               strand='-',
               score=101:110, GC=runif(10),
               seqinfo=seqinfo)
gr3 <- GRanges(seqnames=Rle(c('chr1', 'chr2', 'chr3'), c(3, 4, 3)),
               IRanges(101:110, width=10),
               strand='-',
               score=21:30,
               seqinfo=seqinfo)
some.gr <- c(gr, gr2)

c(gr, gr2, gr3)
c(gr, gr2, gr3, ignore.mcols=TRUE)

## ---------------------------------------------------------------------
## USING A GRANGES OBJECT AS A SUBSCRIPT TO SUBSET ANOTHER OBJECT
## ---------------------------------------------------------------------
## Subsetting *by* a GRanges subscript is supported only if the object
## to subset is a named list-like object:
x <- RleList(chr1=101:120, chr2=2:-8, chr3=31:40)
x[gr]
}
\keyword{methods}
\keyword{classes}

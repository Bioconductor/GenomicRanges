\name{SummarizedExperiment-class}
\docType{class}

% Class
\alias{SummarizedExperiment-class}
\alias{SummarizedExperiment}
\alias{Assays-class}

% Constructor-like function:
\alias{SummarizedExperiment}
\alias{SummarizedExperiment,missing-method}
\alias{SummarizedExperiment,matrix-method}
\alias{SummarizedExperiment,list-method}
\alias{SummarizedExperiment,SimpleList-method}

% Coercion methods:
\alias{updateObject,SummarizedExperiment-method}

% Accessor methods:
\alias{assays}
\alias{assays,SummarizedExperiment-method}
\alias{assays<-}
\alias{assays<-,SummarizedExperiment,list-method}
\alias{assays<-,SummarizedExperiment,SimpleList-method}
\alias{assay}
\alias{assay,SummarizedExperiment,ANY-method}
\alias{assay,SummarizedExperiment,character-method}
\alias{assay,SummarizedExperiment,numeric-method}
\alias{assay,SummarizedExperiment,missing-method}
\alias{assay<-}
\alias{assay<-,SummarizedExperiment,character,matrix-method}
\alias{assay<-,SummarizedExperiment,numeric,matrix-method}
\alias{assay<-,SummarizedExperiment,missing,matrix-method}
\alias{rowData}
\alias{rowData,SummarizedExperiment-method}
\alias{rowData<-}
\alias{rowData<-,SummarizedExperiment,GenomicRanges-method}
\alias{rowData<-,SummarizedExperiment,GRangesList-method}
\alias{colData}
\alias{colData,SummarizedExperiment-method}
\alias{colData<-}
\alias{colData<-,SummarizedExperiment,DataFrame-method}
\alias{exptData}
\alias{exptData,SummarizedExperiment-method}
\alias{exptData<-}
\alias{exptData<-,SummarizedExperiment,SimpleList-method}
\alias{exptData<-,SummarizedExperiment,list-method}
\alias{dim,SummarizedExperiment-method}
\alias{dimnames,SummarizedExperiment-method}
\alias{dimnames<-,SummarizedExperiment,NULL-method}
\alias{dimnames<-,SummarizedExperiment,list-method}

% Subset methods
\alias{[,SummarizedExperiment-method}
\alias{[<-,SummarizedExperiment,ANY,ANY,SummarizedExperiment-method}
\alias{$,SummarizedExperiment-method}
\alias{$<-,SummarizedExperiment,ANY-method}
\alias{$<-,SummarizedExperiment-method}
\alias{[[,SummarizedExperiment,ANY,missing-method}
\alias{[[<-,SummarizedExperiment,ANY,missing,ANY-method}
\alias{[[<-,SummarizedExperiment,ANY,missing-method}

% Combine methods 
\alias{cbind,SummarizedExperiment-method}
\alias{rbind,SummarizedExperiment-method}

% GenomicRanges compatibility methods
\alias{compare,ANY,SummarizedExperiment-method}
\alias{compare,SummarizedExperiment,ANY-method}
\alias{compare,SummarizedExperiment,SummarizedExperiment-method}
\alias{countOverlaps,Vector,SummarizedExperiment-method}
\alias{countOverlaps,SummarizedExperiment,Vector-method}
\alias{countOverlaps,SummarizedExperiment,SummarizedExperiment-method}
\alias{coverage,SummarizedExperiment-method}
\alias{disjointBins,SummarizedExperiment-method}
\alias{distance,ANY,SummarizedExperiment-method}
\alias{distance,SummarizedExperiment,ANY-method}
\alias{distance,SummarizedExperiment,SummarizedExperiment-method}
\alias{distanceToNearest,ANY,SummarizedExperiment-method}
\alias{distanceToNearest,SummarizedExperiment,ANY-method}
\alias{distanceToNearest,SummarizedExperiment,SummarizedExperiment-method}
\alias{duplicated,SummarizedExperiment-method}
\alias{elementMetadata,SummarizedExperiment-method}
\alias{elementMetadata<-,SummarizedExperiment-method}
\alias{end,SummarizedExperiment-method}
\alias{end<-,SummarizedExperiment-method}
\alias{findOverlaps,Vector,SummarizedExperiment-method}
\alias{findOverlaps,SummarizedExperiment,Vector-method}
\alias{findOverlaps,SummarizedExperiment,SummarizedExperiment-method}
\alias{flank,SummarizedExperiment-method}
\alias{follow,ANY,SummarizedExperiment-method}
\alias{follow,SummarizedExperiment,ANY-method}
\alias{follow,SummarizedExperiment,SummarizedExperiment-method}
\alias{granges,SummarizedExperiment-method}
\alias{isDisjoint,SummarizedExperiment-method}
\alias{match,Vector,SummarizedExperiment-method}
\alias{match,SummarizedExperiment,Vector-method}
\alias{match,SummarizedExperiment,SummarizedExperiment-method}
\alias{\%in\%,Vector,SummarizedExperiment-method}
\alias{\%in\%,SummarizedExperiment,Vector-method}
\alias{\%in\%,SummarizedExperiment,SummarizedExperiment-method}
\alias{mcols,SummarizedExperiment-method}
\alias{mcols<-,SummarizedExperiment-method}
\alias{narrow,SummarizedExperiment-method}
\alias{nearest,ANY,SummarizedExperiment-method}
\alias{nearest,SummarizedExperiment,ANY-method}
\alias{nearest,SummarizedExperiment,SummarizedExperiment-method}
\alias{nearest,SummarizedExperiment,missing-method}
\alias{order,SummarizedExperiment-method}
\alias{overlapsAny,SummarizedExperiment,SummarizedExperiment-method}
\alias{overlapsAny,SummarizedExperiment,Vector-method}
\alias{overlapsAny,Vector,SummarizedExperiment-method}
\alias{precede,ANY,SummarizedExperiment-method}
\alias{precede,SummarizedExperiment,ANY-method}
\alias{precede,SummarizedExperiment,SummarizedExperiment-method}
\alias{ranges,SummarizedExperiment-method}
\alias{ranges<-,SummarizedExperiment-method}
\alias{rank,SummarizedExperiment-method}
\alias{resize,SummarizedExperiment-method}
\alias{restrict,SummarizedExperiment-method}
\alias{seqinfo,SummarizedExperiment-method}
\alias{seqinfo<-,SummarizedExperiment-method}
\alias{seqnames,SummarizedExperiment-method}
\alias{shift,SummarizedExperiment-method}
\alias{sort,SummarizedExperiment-method}
\alias{splitAsListReturnedClass,GRanges-method}
\alias{start,SummarizedExperiment-method}
\alias{start<-,SummarizedExperiment-method}
\alias{strand,SummarizedExperiment-method}
\alias{strand<-,SummarizedExperiment-method}
\alias{subsetByOverlaps,SummarizedExperiment,SummarizedExperiment-method}
\alias{subsetByOverlaps,SummarizedExperiment,Vector-method}
\alias{subsetByOverlaps,Vector,SummarizedExperiment-method}
\alias{values,SummarizedExperiment-method}
\alias{values<-,SummarizedExperiment-method}
\alias{width,SummarizedExperiment-method}
\alias{width<-,SummarizedExperiment-method}

% show method:
\alias{show,SummarizedExperiment-method}

\title{SummarizedExperiment instances}

\description{

  The SummarizedExperiment class is an eSet-like container where rows
  represent ranges of interest (as a \code{\linkS4class{GRanges}-class})
  and columns represent samples (with sample data summarized as a
  \code{\linkS4class{DataFrame}-class}). A \code{SummarizedExperiment}
  contains one or more assays, each represented by a matrix of numeric
  or other mode.

}

\usage{

## Constructors

SummarizedExperiment(assays, ...)
\S4method{SummarizedExperiment}{SimpleList}(assays, rowData = GRangesList(),
    colData = DataFrame(), exptData = SimpleList(), ...,
    verbose = FALSE)
\S4method{SummarizedExperiment}{missing}(assays, ...)
\S4method{SummarizedExperiment}{list}(assays, ...)
\S4method{SummarizedExperiment}{matrix}(assays, ...)

## Accessors

assays(x, ..., withDimnames=TRUE)
assays(x, ...) <- value
assay(x, i, ...)
assay(x, i, ...) <- value
rowData(x, ...)
rowData(x, ...) <- value
colData(x, ...)
colData(x, ...) <- value
exptData(x, ...)
exptData(x, ...) <- value
\S4method{dim}{SummarizedExperiment}(x)
\S4method{dimnames}{SummarizedExperiment}(x)
\S4method{dimnames}{SummarizedExperiment,NULL}(x) <- value
\S4method{dimnames}{SummarizedExperiment,list}(x) <- value

## colData access

\S4method{$}{SummarizedExperiment}(x, name)
\S4method{$}{SummarizedExperiment,ANY}(x, name) <- value
\S4method{[[}{SummarizedExperiment,ANY,missing}(x, i, j, ...)
\S4method{[[}{SummarizedExperiment,ANY,missing,ANY}(x, i, j, ...) <- value

## rowData access
## see 'GRanges compatibility', below

## Subsetting

\S4method{[}{SummarizedExperiment}(x, i, j, ..., drop=TRUE)
\S4method{[}{SummarizedExperiment,ANY,ANY,SummarizedExperiment}(x, i, j) <- value

## Combining 

\S4method{cbind}{SummarizedExperiment}(..., deparse.level=1)
\S4method{rbind}{SummarizedExperiment}(..., deparse.level=1)

## Coercion

\S4method{updateObject}{SummarizedExperiment}(object, ..., verbose=FALSE)

}

\arguments{

  \item{assays}{A \code{list} or \code{SimpleList} of matrix elements,
    or a \code{matrix}. Each element of the list must have the same
    dimensions, and dimension names (if present) must be consistent
    across elements and with the row names of \code{rowData},
    \code{colData}.}

  \item{rowData}{A \code{GRanges} or \code{GRangesList} instance describing 
    the ranges of interest. Row names, if present, become the row names of the
    \code{SummarizedExperiment}. The length of the \code{GRanges} or the
    \code{GRangesList} must equal the number of rows of the 
    matrices in \code{assays}.}

  \item{colData}{An optional \code{DataFrame} describing the
    samples. Row names, if present, become the column names of the
    \code{SummarizedExperiment}.}

  \item{exptData}{An optional \code{SimpleList} of arbitrary content
    describing the overall experiment.}

  \item{...}{For \code{SummarizedExperiment}, S4 methods \code{list}
    and \code{matrix}, arguments identical to those of the
    \code{SimpleList} method.

    For \code{assay}, \code{...} may contain \code{withNames}, which is
    forwarded to \code{assays}.

    For \code{cbind, rbind}, \code{...} contains \code{SummarizedExperiment}
    objects to be combined.

    For other accessors, ignored.

  }

  \item{verbose}{A \code{logical(1)} indicating whether messages
    about data coercion during construction should be printed.}

  \item{x, object}{An instance of \code{SummarizedExperiment}-class.}

  \item{i, j}{For \code{assay}, \code{assay<-}, \code{i} is a integer or
    numeric scalar; see \sQuote{Details} for additional constraints.

    For \code{[,SummarizedExperiment}, \code{[,SummarizedExperiment<-},
    \code{i}, \code{j} are instances that can act to subset the
    underlying \code{rowData}, \code{colData}, and \code{matrix}
    elements of \code{assays}.

    For \code{[[,SummarizedExperiment},
    \code{[[<-,SummarizedExperiment}, \code{i} is a scalar index (e.g.,
    \code{character(1)} or \code{integer(1)}) into a column of
    \code{colData}.

  }

  \item{name}{A symbol representing the name of a column of
  \code{colData}.}

  \item{withDimnames}{A \code{logical(1)}, indicating whether dimnames
    should be applied to extracted assay elements.}

  \item{drop}{A \code{logical(1)}, ignored by these methods.}

  \item{value}{An instance of a class specified in the S4 method
    signature or as outlined in \sQuote{Details}.}

  \item{deparse.level}{See \code{?base::\link[base]{cbind}} for a description of
    this argument.}

}

\details{

  The \code{SummarizedExperiment} class is meant for numeric and other
  data types derived from a sequencing experiment. The structure is
  rectangular, like an \code{eSet} in \pkg{Biobase}.

  The rows of a \code{SummarizedExperiment} instance represent ranges
  (in genomic coordinates) of interest. The ranges of interest are
  described by a \code{\linkS4class{GRanges}-class} or a 
  \code{\linkS4class{GRangesList}-class} instance, accessible
  using the \code{rowData} function, described below. The \code{GRanges}
  and \code{GRangesList} classes contains sequence (e.g., chromosome) name, 
  genomic coordinates, and strand information. Each range can be annotated 
  with additional data; this data might be used to describe the range 
  (analogous to annotations associated with genes in a \code{eSet}) or to 
  summarize results (e.g., statistics of differential abundance) relevant to the
  range. Rows may or may not have row names; they often will not.

  Each column of a \code{SummarizedExperiment} instance represents a
  sample.  Information about the samples are stored in a
  \code{\linkS4class{DataFrame}-class}, accessible using the function
  \code{colData}, described below. The \code{DataFrame} must have as
  many rows as there are columns in the \code{SummarizedExperiment},
  with each row of the \code{DataFrame} providing information on the
  sample in the corresponding column of the
  \code{SummarizedExperiment}. Columns of the \code{DataFrame} represent
  different sample attributes, e.g., tissue of origin, etc. Columns of
  the \code{DataFrame} can themselves be annotated (via the
  \code{\link{mcols}} function) in a fashion similar to the
  \code{varMetadata} facilities of the \code{eSet} class. Column names
  typically provide a short identifier unique to each sample.

  A \code{SummarizedExperiment} can also contain information about the
  overall experiment, for instance the lab in which it was conducted,
  the publications with which it is associated, etc. This information is
  stored as a \code{\linkS4class{SimpleList}-class}, accessible using
  the \code{exptData} function. The form of the data associated with the
  experiment is left to the discretion of the user.

  The \code{SummarizedExperiment} is appropriate for matrix-like
  data. The data are accessed using the \code{assays} function,
  described below. This returns a \code{SimpleList}-class instance. Each
  element of the list must itself be a matrix (of any mode) and must
  have dimensions that are the same as the dimensions of the
  \code{SummarizedExperiment} in which they are stored. Row and column
  names of each matrix must either be \code{NULL} or match those of the
  \code{SummarizedExperiment} during construction. It is convenient for
  the elements of \code{SimpleList} of assays to be named.

  The \code{SummarizedExperiment} class has the following slots; this
  detail of class structure is not relevant to the user.
  \describe{
    \item{\code{exptData}}{A \link{SimpleList}-class instance containing
      information about the overall experiment.}
    \item{\code{rowData}}{A \link{GRanges}-class instance defining the
      ranges of interest and associated metadata.}
    \item{\code{colData}}{A \link{DataFrame}-class instance describing
      the samples and associated metadata.}
    \item{\code{assays}}{A \link{SimpleList}-class instance, each
      element of which is a matrix summarizing data associated with the
      corresponding range and sample.}
  }
}

\section{Constructor}{

  Instances are constructed using the \code{SummarizedExperiment}
  function with arguments outlined above.

}

\section{Coercion}{

  Package version 1.9.59 introduced a new way of representing
  \sQuote{assays}. If you have a serialized instance \code{x} of a
  \code{SummarizedExperiment} (e.g., from using the \code{save} function
  with a version of GenomicRanges prior to 1.9.59), it should be updated
  by invoking \code{x <- updateObject(x)}.

}

\section{Accessors}{

  In the following code snippets, \code{x} is a
  \code{SummarizedExperiment} instance.

  \describe{

    \item{\code{assays(x)}, \code{assays(x) <- value}:}{Get or set the
      assays. \code{value} is a \code{list} or \code{SimpleList}, each
      element of which is a matrix with the same dimensions as
      \code{x}.}

    \item{\code{assay(x, i)}, \code{assay(x, i) <- value}:}{A convenient
      alternative (to \code{assays(x)[[i]]}, \code{assays(x)[[i]] <-
      value}) to get or set the \code{i}th (default first) assay
      element. \code{value} must be a matrix of the same dimension as
      \code{x}, and with dimension names \code{NULL} or consistent with
      those of \code{x}.}

    \item{\code{rowData(x)}, \code{rowData(x) <- value}:}{Get or set the
      row data. \code{value} is a \code{GenomicRanges} instance. Row
      names of \code{value} must be NULL or consistent with the existing
      row names of \code{x}.}
   
    \item{\code{colData(x)}, \code{colData(x) <- value}:}{Get or set the
      column data. \code{value} is a \code{DataFrame} instance. Row
      names of \code{value} must be NULL or consistent with the existing
      column names of \code{x}.}

    \item{\code{exptData(x)}, \code{exptData(x) <- value}:}{Get or set
      the experiment data. \code{value} is a \code{list} or
      \code{SimpleList} instance, with arbitrary content.}

    \item{\code{dim(x)}:}{Get the dimensions (ranges x samples) of the
      \code{SummarizedExperiment}.}

    \item{\code{dimnames(x)}, \code{dimnames(x) <- value}:}{Get or set
      the dimension names. \code{value} is usually a list of length 2,
      containing elements that are either \code{NULL} or vectors of
      appropriate length for the corresponding dimension. \code{value}
      can be \code{NULL}, which removes dimension names. This method
      implies that \code{rownames}, \code{rownames<-}, \code{colnames},
      and \code{colnames<-} are all available.}

  }
}


\section{GRanges compatibility (rowData access)}{

  Many \code{\link{GRanges}-class} operations are supported on
  \sQuote{SummarizedExperiment} and derived instances, using
  \code{rowData}.

  Supported operations include: \code{\link{compare}},
  \code{\link{countOverlaps}}, \code{\link{coverage}},
  \code{\link{disjointBins}}, \code{\link{distance}},
  \code{\link{distanceToNearest}}, \code{\link{duplicated}},
  \code{\link{end}}, \code{\link{end<-}}, \code{\link{findOverlaps}},
  \code{\link{flank}}, \code{\link{follow}}, \code{\link{granges}},
  \code{\link{isDisjoint}}, \code{\link{match}}, \code{\link{mcols}},
  \code{\link{mcols<-}}, \code{\link{narrow}}, \code{\link{nearest}},
  \code{\link{order}}, \code{\link{overlapsAny}}, \code{\link{precede}},
  \code{\link{ranges}},
  \code{\link{ranges<-}}, \code{\link{rank}}, \code{\link{resize}},
  \code{\link{restrict}}, \code{\link{seqinfo}},
  \code{\link{seqinfo<-}}, \code{\link{seqnames}}, \code{\link{shift}},
  \code{\link{sort}}, \code{splitAsListReturnedClass},
  \code{\link{start}}, \code{\link{start<-}},
  \code{\link{strand}}, \code{\link{strand<-}},
  \code{\link{subsetByOverlaps}}, \code{\link{width}},
  \code{\link{width<-}}.

  Not all \code{\link{GRanges}-class} operations are supported, because
  they do not make sense for \sQuote{SummarizedExperiment} objects
  (e.g., length, name, as.data.frame, c, splitAsList), involve
  non-trivial combination or splitting of rows (e.g., disjoin, gaps,
  reduce, unique), or have not yet been implemented (Ops,
  map, seqselect, seqselect<-, window, window<-).

}

\section{Subsetting}{

  In the code snippets below, \code{x} is a \code{SummarizedExperiment}
  instance.

  \describe{

    \item{\code{x[i,j]}, \code{x[i,j] <- value}:}{Create or replace a
      subset of \code{x}. \code{i}, \code{j} can be \code{numeric},
      \code{logical}, \code{character}, or \code{missing}. \code{value}
      must be a \code{SummarizedExperiment} instance with dimensions,
      dimension names, and assay elements consistent with the subset
      \code{x[i,j]} being replaced.}

  }

  Additional subsetting accessors provide convenient access to
  \code{colData} columns

  \describe{

    \item{\code{x$name}, \code{x$name <- value}}{Access or replace
    column \code{name} in \code{x}.}

    \item{\code{x[[i, ...]]}, \code{x[[i, ...]] <- value}}{Access or
    replace column \code{i} in \code{x}.}

  }

}

\section{Combining}{

  In the code snippets below, \code{...} are \code{SummarizedExperiment}
  instances to be combined.

  \describe{

    \item{\code{cbind(...)}, \code{rbind(...)}:}{
      \code{cbind} combines objects with identical ranges (\code{rowData}) 
      but different samples (columns in \code{assays}). Objects with a 
      different number of rows or different rownames will not be
      combined.

      \code{rbind} combines objects with different ranges (\code{rowData}) 
      and the same subjects (columns in \code{assays}). Objects with a
      different number of columns or different colnames will not be
      combined.

      Currently \code{exptData} from all objects are combined into a
      single \code{SimpleList} without checking for data duplication.

    }

  }

}

\section{Implementation and Extension}{

  This section contains advanced material meant for package developers.

  \code{SummarizedExperiment} is implemented as an S4 class, and can be
  extended in the usual way, using
  \code{contains="SummarizedExperiment"} in the new class definition.

  In addition, the representation of the \code{assays} slot of
  \code{SummarizedExperiment} is as a virtual class \code{Assays}. This
  allows derived classes (\code{contains="Assays"}) to easily implement
  alternative requirements for the assays, e.g., backed by file-based
  storage like NetCDF or the \code{ff} package, while re-using the
  existing \code{SummarizedExperiment} class without modification. The
  requirements on \code{Assays} are list-like semantics (e.g.,
  \code{sapply}, \code{[[} subsetting, \code{names}) with elements
  having matrix- or array-like semantics (e.g., \code{dim},
  \code{dimnames}). These requirements can be made more precise if
  developers express interest.

  The current \code{assays} slot is implemented as a reference class
  that has copy-on-change semantics. This means that modifying non-assay
  slots does not copy the (large) assay data, and at the same time the
  user is not surprised by reference-based semantics. Updates to
  non-assay slots are very fast; updating the assays slot itself can be
  5x or more faster than with an S4 instance in the slot.

  In a little more detail, a small reference class hierarchy (not
  exported from the GenomicRanges name space) defines a reference class
  \code{ShallowData} with a single field \code{data} of type \code{ANY},
  and a derived class \code{ShallowSimpleListAssays} that specializes
  the type of \code{data} as \code{SimpleList}, and
  \code{contains=c("ShallowData", "Assays")}. The assays slot contains
  an instance of \code{ShallowSimpleListAssays}. Invoking
  \code{assays()} on a \code{SummarizedExperiment} re-dispatches from
  the \code{assays} slot to retrieve the \code{SimpleList} from the
  field of the reference class. This was achieved by implementing a
  generic (not exported) \code{value(x, name, ...)}, with a method
  implemented on \code{SummarizedExperiment} that retrieves a slot when
  \code{name} is a slot containing an S4 object in \code{x}, and a field
  when \code{name} is a slot containing a \code{ShallowData} instance in
  \code{x}. Copy-on-change semantics is maintained by implementing the
  \code{clone} method (\code{clone} methods are supposed to do a deep
  copy, \code{update} methods a shallow copy; the \code{clone} generic
  is introduced, and not exported, in the GenomicRanges package). The
  \sQuote{getter} and \sQuote{setter} code for methods implemented on
  \code{SummarizedExperiment} use \code{value} for slot access, and
  \code{clone} for replacement. This makes it easy to implement
  \code{ShallowData} instances for other slots if the need arises.

}

\author{Martin Morgan, \url{mtmorgan@fhcrc.org}}

\seealso{
  \link{GRanges},
  \link[IRanges]{DataFrame},
  \link[IRanges]{SimpleList},
}

\examples{
  nrows <- 200; ncols <- 6
  counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
  rowData <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE))
  colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                       row.names=LETTERS[1:6])
  sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
                 rowData=rowData, colData=colData)
  sset
  assays(sset) <- endoapply(assays(sset), asinh)
  head(assay(sset))

  sset[, sset$Treatment == "ChIP"]

  ## cbind is appropriate for combining objects with the same ranges and 
  ## different samples.
  se1 <- sset
  se2 <- se1[,1:3]
  colnames(se2) <- letters[1:ncol(se2)] 
  cmb1 <- cbind(se1, se2)

  ## rbind is appropriate for combining objects with the same samples and
  ## different ranges.
  se1 <- sset
  se2 <- se1[1:50,]
  rownames(se2) <- letters[1:nrow(se2)] 
  cmb2 <- rbind(se1, se2)

}

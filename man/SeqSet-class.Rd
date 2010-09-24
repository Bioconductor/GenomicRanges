\name{SeqSet-class}
\docType{class}

% Class
\alias{SeqSet-class}
\alias{SeqSet}

% Constructor-like function:
\alias{SeqSet}

% Coercion methods:

% Accessor methods:
\alias{assays,SeqSet-method}
\alias{assays<-,SeqSet,list-method}
\alias{assays<-,SeqSet,SimpleList-method}
\alias{assay,SeqSet,ANY-method}
\alias{assay,SeqSet,character-method}
\alias{assay,SeqSet,numeric-method}
\alias{assay,SeqSet,missing-method}
\alias{assay<-,SeqSet-method}
\alias{assay<-,SeqSet,character-method}
\alias{assay<-,SeqSet,numeric-method}
\alias{assay<-,SeqSet,missing-method}
\alias{rowData,SeqSet-method}
\alias{rowData<-,SeqSet,GenomicRanges-method}
\alias{colData,SeqSet-method}
\alias{colData<-,SeqSet,DataFrame-method}
\alias{exptData,SeqSet-method}
\alias{exptData<-,SeqSet,SimpleList-method}
\alias{exptData<-,SeqSet,list-method}
\alias{dim,SeqSet-method}
\alias{dimnames,SeqSet-method}
\alias{dimnames<-,SeqSet,NULL-method}
\alias{dimnames<-,SeqSet,list-method}

% Subset methods
\alias{[,SeqSet-method}
\alias{[<-,SeqSet,ANY,ANY,SeqSet-method}

% show method:
\alias{show,SeqSet-method}

\title{SeqSet objects}

\description{

  The SeqSet class is an eSet-like container where rows represent ranges
  of interest (as a \code{\linkS4class{GRanges}-class}) and columns
  represent samples (with sample data summarized as a
  \code{\linkS4class{DataFrame}-class}). A \code{SeqSet} contains one or
  more assays, each represented by a matrix of numeric or other mode.

}

\details{

  The \code{SeqSet} class is meant for numeric and other data types
  derived from a sequencing experiment. The structure is rectangular,
  like an \code{eSet} in \pkg{Biobase}.

  The rows of a \code{SeqSet} instance represent ranges (in genomic
  coordinates) of interest. The ranges of interest are described by a
  \code{\linkS4class{GRanges}-class} instance, accessible using the
  \code{rowData} function, described below. The \code{GRanges} class
  contains sequence (e.g., chromosome) name, genomic coordinates, and
  strand information. Each range can be annotated with additional data;
  this data might be used to describe the range (analogous to
  annotations associated with genes in a \code{eSet}) or to summarize
  results (e.g., statistics of differential abundance) relevant to the
  range. Rows may or may not have row names; they often will not.

  Each column of a \code{SeqSet} instance represents a sample.
  Information about the samples are stored in a
  \code{\linkS4class{DataFrame}-class}, accessible using the function
  \code{colData}, described below. The \code{DataFrame} must have as
  many rows as there are columns in the \code{SeqSet}, with each row of
  the \code{DataFrame} providing information on the sample in the
  corresponding column of the \code{SeqSet}. Columns of the
  \code{DataFrame} represent different sample attributes, e.g., tissue
  of origin, etc. Columns of the \code{DataFrame} can themselves be
  annotated (via the \code{\link{values}} function) in a fashion similar
  to the \code{varMetadata} facilities of the \code{eSet} class. Column
  names typically provide a short identifier unique to each sample.

  A \code{SeqSet} can also contain information about the overall
  experiment, for instance the lab in which it was conducted, the
  publications with which it is associated, etc. This information is
  stored as a \code{\linkS4class{SimpleList}-class}, accessible using
  the \code{exptData} function. The form of the data associated with the
  experiment is left to the discretion of the user.

  The \code{SeqSet} is appropriate for matrix-like data. The data are
  accessed using the \code{assays} function, described below. This
  returns a \code{SimpleList}-class instance. Each element of the list
  must itself be a matrix (of any mode) and must have dimensions that
  are the same as the dimensions of the \code{SeqSet} in which they are
  stored. Row and column names of each matrix must either be \code{NULL}
  or match those of the \code{SeqSet} during construction. It is
  convenient for the elements of \code{SimpleList} of assays to be
  named.

  The \code{SeqSet} class has the following slots; this detail of class
  structure is not relevant to the user.
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
  \describe{
    \item{\code{SeqSet(assays=SimpleList(), rowData=GRanges(),
	colData=DataFrame(), exptData=SimpleList(), ...)}:}{

      Creates a \code{SeqSet} instance, with arguments as follows:
      \describe{

	\item{\code{assays}:}{A \code{SimpleList} of matrix elements. Each
	  element must have the same dimensions, and dimension names (if
	  present) must be consistent across elements and with the row
	  names of \code{rowData}, \code{colData}.}

	\item{\code{rowData}:}{A \code{GRanges} describing the ranges of
	  interest. Row names, if present, become the row names of the
	  \code{SeqSet}.}

	\item{\code{colData}:}{A \code{DataFrame} describing the
	  samples. Row names, if present, become the column names of the
	  \code{SeqSet}.}

	\item{\code{exptData}:}{A \code{SimpleList} of arbitrary content
	  describing the overall experiment.}

      }
    }
  }
}

% \section{Coercion}{}

\section{Accessors}{
  In the following code snippets, \code{x} is a \code{SeqSet} object.

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
      \code{SeqSet}.}

    \item{\code{dimnames(x)}, \code{dimnames(x) <- value}:}{Get or set
      the dimension names. \code{value} is usually a list of length 2,
      containing elements that are either \code{NULL} or vectors of
      appropriate length for the corresponding dimension. \code{value}
      can be \code{NULL}, which removes dimension names. This method
      implies that \code{rownames}, \code{rownames<-}, \code{colnames},
      and \code{colnames<-} are all available.}

  }
}

\section{Subsetting}{
  In the code snippets below, \code{x} is a \code{SeqSet} object.

  \describe{

    \item{\code{x[i,j]}, \code{x[i,j] <- value}:}{Create or replace a
      subset of \code{x}. \code{i}, \code{j} can be \code{numeric},
      \code{logical}, \code{character}, or \code{missing}. \code{value}
      must be a \code{SeqSet} instance with dimensions, dimension names,
      and assay elements consistent with the subset \code{x[i,j]} being
      replaced.}
  }
}

\author{Martin Morgan}

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
  sset <- SeqSet(assays=SimpleList(counts=counts),
                 rowData=rowData, colData=colData)
  sset
  assays(sset) <- endoapply(assays(sset), asinh)
  head(assay(sset))
}

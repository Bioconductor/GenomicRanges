\name{makeGRangesListFromDataFrame}

\alias{makeGRangesListFromDataFrame}

\title{Make a GRangesList object from a data.frame or DataFrame}

\description{
    \code{makeGRangesListFromDataFrame} extends the
    \link[GenomicRanges]{makeGRangesFromDataFrame} functionality from
    \code{GenomicRanges}. It can take a data-frame-like object as input
    and tries to automatically find the columns that describe the genomic
    ranges. It returns a \link{GRangesList} object. This is different from the
    \code{makeGRangesFromDataFrame} function by requiring a
    \code{split.field}. The \code{split.field} acts like the
    "f" argument in the \code{\link[S4Vectors]{split}} function. This factor
    must be of the same length as the number of rows in the \code{DataFrame}
    argument. The \code{split.field} may also be a character vector.
}

\usage{
makeGRangesListFromDataFrame(df,
                             split.field = NULL,
                             names.field = NULL,
                             ...)
}
\arguments{
  \item{df}{
      A \code{DataFrame} or \code{data.frame} class object
  }
  \item{split.field}{
      A character string of a recognized column name in \code{df} that contains
      the grouping. This column defines how the rows of \code{df} are split and
      is typically a \code{factor} or \code{character} vector.  When
      \code{split.field} is not provided the \code{df} will be split by the
      number of rows.
  }
  \item{names.field}{
      An optional single \code{character} string indicating the name of the
      column in \code{df} that designates the names for the ranges in the
      elements of the \link{GRangesList}.
  }
  \item{...}{
      Additional arguments passed on to \link{makeGRangesFromDataFrame}
  }
}
\value{
  A \linkS4class{GRangesList} of the same length as the number of levels or
  unique character strings in the \code{df} column indicated by
  \code{split.field}. When \code{split.field} is not provided the \code{df}
  is split by row and the resulting \linkS4class{GRangesList} has the 
  same length as nrow(df).

  Names on the individual ranges are taken from the \code{names.field}
  argument. Names on the outer list elements of the \linkS4class{GRangesList} 
  are propagated from \code{split.field}.
}

\author{
  M. Ramos
}

\seealso{
  \itemize{
      \item \link{makeGRangesFromDataFrame}
  }
}

\examples{
## ---------------------------------------------------------------------
## BASIC EXAMPLES
## ---------------------------------------------------------------------

df <- data.frame(chr="chr1", start=11:15, end=12:16,
                 strand=c("+","-","+","*","."), score=1:5,
                 specimen = c("a", "a", "b", "b", "c"),
                 gene_symbols = paste0("GENE", letters[1:5]))
df

grl <- makeGRangesListFromDataFrame(df, split.field = "specimen",
                                    names.field = "gene_symbols")
grl
names(grl)

## Keep metadata columns
makeGRangesListFromDataFrame(df, split.field = "specimen",
                             keep.extra.columns = TRUE)
}

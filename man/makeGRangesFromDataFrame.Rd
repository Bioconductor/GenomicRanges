\name{makeGRangesFromDataFrame}

\alias{makeGRangesFromDataFrame}

\alias{coerce,data.frame,GRanges-method}
\alias{coerce,DataFrame,GRanges-method}

\title{Make a GRanges object from a data.frame or DataFrame}

\description{
  \code{makeGRangesFromDataFrame} finds the fields in the input that describe
  genomic ranges and returns them as a \link{GRanges} object.

  For convenience, coercing a data.frame or \link[S4Vectors]{DataFrame}
  \code{df} into a \link{GRanges} object is supported and does
  \code{makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)}
}

\usage{
makeGRangesFromDataFrame(df,
    keep.extra.columns=FALSE,
    ignore.strand=FALSE,
    seqinfo=NULL,
    seqnames.field=c("seqnames", "chr", "chrom"),
    start.field=c("start", "chromStart"),
    end.field=c("end", "chromEnd", "stop", "chromStop"),
    strand.field="strand",
    starts.in.df.are.0based=FALSE)
}

\arguments{
  \item{df}{
    A data.frame or \link[S4Vectors]{DataFrame} object.
  }
  \item{keep.extra.columns}{
    \code{TRUE} or \code{FALSE} (the default).
    If \code{TRUE}, then the columns in \code{df} that are not used to form
    the genomic ranges returned in the \link{GRanges} object will be stored
    in it as metadata columns. Otherwise, they will be ignored.
  }
  \item{ignore.strand}{
    \code{TRUE} or \code{FALSE} (the default).
    If \code{TRUE}, then the strand of the returned \link{GRanges} object
    will be set to \code{"*"}.
  }
  \item{seqinfo}{
    Either \code{NULL}, or a \link{Seqinfo} object, or a character vector
    of seqlevels, or a named numeric vector of sequence lengths.
    When not \code{NULL}, it must be compatible with the genomic ranges
    in \code{df} i.e. it must include at least the sequence levels
    represented in \code{df}.
  }
  \item{seqnames.field}{
    A character vector of recognized names for the column in \code{df}
    that contains the chromosome name (a.k.a. sequence name) associated
    with each genomic range.
    Only the first name in \code{seqnames.field} that is found
    in \code{colnames(df)} will be used.
    If no one is found, then an error is raised.
  }
  \item{start.field}{
    A character vector of recognized names for the column in \code{df}
    that contains the start positions of the genomic ranges.
    Only the first name in \code{start.field} that is found
    in \code{colnames(df)} will be used.
    If no one is found, then an error is raised.
  }
  \item{end.field}{
    A character vector of recognized names for the column in \code{df}
    that contains the end positions of the genomic ranges.
    Only the first name in \code{start.field} that is found
    in \code{colnames(df)} will be used.
    If no one is found, then an error is raised.
  }
  \item{strand.field}{
    A character vector of recognized names for the column in \code{df}
    that contains the strand associated with each genomic range.
    Only the first name in \code{strand.field} that is found
    in \code{colnames(df)} will be used.
    If no one is found or if \code{ignore.strand} is \code{TRUE},
    then the strand of the returned \link{GRanges} object will be
    set to \code{"*"}.
  }
  \item{starts.in.df.are.0based}{
    \code{TRUE} or \code{FALSE} (the default).
    If \code{TRUE}, then the start positions of the genomic ranges in
    \code{df} are considered to be 0-based and are converted to 1-based
    in the returned \link{GRanges} object.
    This feature is intended to make it more convenient to handle input
    that contains data obtained from resources using the "0-based
    start" convention. A notorious example of such resource is the UCSC
    Table Browser (\url{http://genome.ucsc.edu/cgi-bin/hgTables}).
  }
}

\value{
  A \link{GRanges} object with one element per row in the input.

  If the \code{seqinfo} argument was supplied, the returned object will
  have exactly the seqlevels specified in \code{seqinfo} and in the same
  order.

  If \code{df} has non-automatic row names (i.e. \code{rownames(df)} is
  not \code{NULL} or \code{seq_len(nrow(df))}), then they will be used to
  set the names of the returned \link{GRanges} object.
}

\author{
  H. Pages, based on a proposal by Kasper Daniel Hansen
}

\seealso{
  \itemize{
    \item \link{GRanges} objects.

    \item \link{Seqinfo} objects.

    \item The \code{\link{makeGRangesListFromFeatureFragments}} function
          for making a \link{GRangesList} object from a list of fragmented
          features.

    \item The \code{\link[rtracklayer]{getTable}} function in the
          \pkg{rtracklayer} package for an R interface to the UCSC
          Table Browser.

    \item \link[S4Vectors]{DataFrame} objects in the \pkg{IRanges} package.
  }
}

\examples{
df <- data.frame(chr="chr1", start=11:13, end=12:14,
                 strand=c("+","-","+"), score=1:3)

makeGRangesFromDataFrame(df)
gr <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
gr2 <- as(df, "GRanges")  # equivalent to the above
stopifnot(identical(gr, gr2))

makeGRangesFromDataFrame(df, ignore.strand=TRUE)
makeGRangesFromDataFrame(df, keep.extra.columns=TRUE,
                             ignore.strand=TRUE)

makeGRangesFromDataFrame(df, seqinfo=paste0("chr", 4:1))
makeGRangesFromDataFrame(df, seqinfo=c(chrM=NA, chr1=500, chrX=100))
makeGRangesFromDataFrame(df, seqinfo=Seqinfo(paste0("chr", 4:1)))

if (require(rtracklayer)) {
  session <- browserSession()
  genome(session) <- "sacCer2"
  query <- ucscTableQuery(session, "Most Conserved")
  df <- getTable(query)

  ## A common pitfall is to forget that the UCSC Table Browser uses the
  ## "0-based start" convention:
  gr0 <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
  head(gr0)
  min(start(gr0))

  ## The start positions need to be converted into 1-based positions,
  ## to adhere to the convention used in Bioconductor:
  gr1 <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE,
                                  starts.in.df.are.0based=TRUE)
  head(gr1)
}
}

\keyword{manip}

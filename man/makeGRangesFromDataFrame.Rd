\name{makeGRangesFromDataFrame}

\alias{makeGRangesFromDataFrame}
\alias{coerce,data.frame,GRanges-method}
\alias{coerce,DataFrame,GRanges-method}

\alias{makeGPosFromDataFrame}
\alias{coerce,data.frame,GPos-method}
\alias{coerce,DataFrame,GPos-method}

\title{Make a GRanges object from a data.frame or DataFrame}

\description{
  \code{makeGRangesFromDataFrame} and \code{makeGPosFromDataFrame} both
  take a data-frame-like object as input and try to automatically find the
  columns that describe genomic ranges (for \code{makeGRangesFromDataFrame})
  or genomic positions (for \code{makeGPosFromDataFrame}).
  If successful, they return them in a \link{GRanges} or \link{GPos} object.

  The two functions are also the workhorses behind the coercion methods from
  data.frame (or \link[S4Vectors]{DataFrame}) to \link{GRanges} or \link{GPos}.
}

\usage{
makeGRangesFromDataFrame(df, keep.extra.columns=FALSE,
        ignore.strand=FALSE, seqinfo=NULL,
        seqnames.field=c("seqnames", "seqname",
                         "chromosome", "chrom", "chr", "chromosome_name",
                         "seqid"),
        start.field="start",
        end.field=c("end", "stop"),
        strand.field="strand",
        starts.in.df.are.0based=FALSE,
        na.rm=FALSE)

makeGPosFromDataFrame(df, keep.extra.columns=FALSE,
        ignore.strand=FALSE, seqinfo=NULL,
        seqnames.field=c("seqnames", "seqname",
                         "chromosome", "chrom", "chr", "chromosome_name",
                         "seqid"),
        pos.field=c("pos", "position", "positions"),
        strand.field="strand",
        na.rm=FALSE)
}

\arguments{
  \item{df}{
    A data.frame or \link[S4Vectors]{DataFrame} object. If not, then
    the function first tries to turn \code{df} into a data frame with
    \code{as.data.frame(df)}.
  }
  \item{keep.extra.columns}{
    \code{TRUE} or \code{FALSE} (the default).
    If \code{TRUE}, then the columns in \code{df} that are not used to form
    the genomic ranges or genomic positions of the returned \link{GRanges}
    or \link{GPos} object are returned as metadata columns on the object.
    Otherwise, they are ignored.

    Note that if \code{df} has a \code{width} column, then
    \code{makeGRangesFromDataFrame} will always ignore it.
  }
  \item{ignore.strand}{
    \code{TRUE} or \code{FALSE} (the default).
    If \code{TRUE}, then the strand of the returned \link{GRanges}
    or \link{GPos} object is set to \code{"*"}.
  }
  \item{seqinfo}{
    Either \code{NULL}, or a \link[Seqinfo]{Seqinfo} object,
    or a character vector of unique sequence names (a.k.a. \emph{seqlevels}),
    or a named numeric vector of sequence lengths.
    When not \code{NULL}, \code{seqinfo} must be compatible with the genomic
    ranges in \code{df}, that is, it must have one entry for each unique
    sequence name represented in \code{df}. Note that it can have additional
    entries i.e. entries for seqlevels not represented in \code{df}.
  }
  \item{seqnames.field}{
    A character vector of recognized names for the column in \code{df}
    that contains the chromosome name (a.k.a. sequence name) associated
    with each genomic range.
    Only the first name in \code{seqnames.field} that is found
    in \code{colnames(df)} is used.
    If no one is found, then an error is raised.
  }
  \item{start.field}{
    A character vector of recognized names for the column in \code{df}
    that contains the start positions of the genomic ranges.
    Only the first name in \code{start.field} that is found
    in \code{colnames(df)} is used.
    If no one is found, then an error is raised.
  }
  \item{end.field}{
    A character vector of recognized names for the column in \code{df}
    that contains the end positions of the genomic ranges.
    Only the first name in \code{start.field} that is found
    in \code{colnames(df)} is used.
    If no one is found, then an error is raised.
  }
  \item{pos.field}{
    A character vector of recognized names for the column in \code{df}
    that contains the genomic positions to use to make the \link{GPos} object.
    Only the first name in \code{pos.field} that is found
    in \code{colnames(df)} is used.
    If no one is found, then an error is raised.
  }
  \item{strand.field}{
    A character vector of recognized names for the column in \code{df}
    that contains the strand associated with each genomic range.
    Only the first name in \code{strand.field} that is found
    in \code{colnames(df)} is used.
    If no one is found or if \code{ignore.strand} is \code{TRUE},
    then the strand of the returned \link{GRanges} object is
    set to \code{"*"}.
  }
  \item{starts.in.df.are.0based}{
    \code{TRUE} or \code{FALSE} (the default).
    If \code{TRUE}, then the start positions of the genomic ranges in
    \code{df} are considered to be \emph{0-based} and are converted to
    \emph{1-based} in the returned \link{GRanges} object.
    This feature is intended to make it more convenient to handle input
    that contains data obtained from resources using the "0-based
    start" convention. A notorious example of such resource is the UCSC
    Table Browser (\url{http://genome.ucsc.edu/cgi-bin/hgTables}).
  }
  \item{na.rm}{
    \code{TRUE} or \code{FALSE} (the default).

    If \code{TRUE}, then rows in the \code{df} with missing start, end,
    or position values (i.e. the value is \code{NA}) are ignored.
    Otherwise, they trigger an error.
  }
}

\value{
  A \link{GRanges} or \link{GPos} object.

  If \code{na.rm} is set to \code{FALSE} (the default), then the
  returned object is guaranteed to have one element per row in the input.
  However, if \code{na.rm} is set to \code{TRUE}, then the length of the
  returned object can be less than \code{nrow(df)}.

  If the \code{seqinfo} argument was supplied, the returned object will
  have exactly the seqlevels specified in \code{seqinfo} and in the same
  order. Otherwise, the seqlevels are ordered according to the output of
  the \code{\link[Seqinfo]{rankSeqlevels}} function (except if
  \code{df} contains the seqnames in the form of a factor or factor-Rle,
  in which case the \code{levels()} of the seqnames become the seqlevels
  of the returned object as-is, that is, without any re-ordering).

  If \code{df} has non-automatic row names (i.e. \code{rownames(df)} is
  not \code{NULL} and is not \code{seq_len(nrow(df))}), then they will be
  used to set names on the returned \link{GRanges} or \link{GPos} object.
}

\note{
  Coercing a data.frame or \link[S4Vectors]{DataFrame} \code{df}
  to \link{GRanges} (with \code{as(df, "GRanges")}), or calling
  \code{GRanges(df)}, are both equivalent to calling
  \code{makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)}.

  Coercing a data.frame or \link[S4Vectors]{DataFrame} \code{df}
  to \link{GPos} (with \code{as(df, "GPos")}) is equivalent to calling
  \code{makeGPosFromDataFrame(df, keep.extra.columns=TRUE)}.
}

\author{
  H. Pagès, based on a proposal by Kasper Daniel Hansen
}

\seealso{
  \itemize{
    \item \link{GRanges} and \link{GPos} objects.

    \item \link[Seqinfo]{Seqinfo} objects and the
          \code{\link[Seqinfo]{rankSeqlevels}} function in the
          \pkg{Seqinfo} package.

    \item The \code{\link{makeGRangesListFromFeatureFragments}} function
          for making a \link{GRangesList} object from a list of fragmented
          features.

    \item \code{\link[IRanges]{makeIRangesFromDataFrame}} in the
          \pkg{IRanges} package to make an \link[IRanges]{IRanges}
          object from a data.frame or \link[S4Vectors]{DataFrame}.

    \item The \code{\link[rtracklayer]{getTable}} function in the
          \pkg{rtracklayer} package for an R interface to the UCSC
          Table Browser.

    \item \link[S4Vectors]{DataFrame} objects in the \pkg{S4Vectors} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## BASIC makeGRangesFromDataFrame() EXAMPLES
## ---------------------------------------------------------------------

df <- data.frame(chr="chr1", start=11:15, end=12:16,
                 strand=c("+","-","+","*","."), score=1:5)
df
makeGRangesFromDataFrame(df)  # strand value "." is replaced with "*"

## With NAs in range-defining columns:
df$start[5] <- df$end[2] <- NA
df
#makeGRangesFromDataFrame(df)  # error!
makeGRangesFromDataFrame(df, na.rm=TRUE)  # rows with NAs got dropped

## The strand column is optional:
df <- data.frame(chr="chr1", start=11:15, end=12:16, score=1:5)
makeGRangesFromDataFrame(df)

gr <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
gr2 <- as(df, "GRanges")  # equivalent to the above
stopifnot(identical(gr, gr2))
gr2 <- GRanges(df)        # equivalent to the above
stopifnot(identical(gr, gr2))

makeGRangesFromDataFrame(df, ignore.strand=TRUE)
makeGRangesFromDataFrame(df, keep.extra.columns=TRUE,
                             ignore.strand=TRUE)

makeGRangesFromDataFrame(df, seqinfo=paste0("chr", 4:1))
makeGRangesFromDataFrame(df, seqinfo=c(chrM=NA, chr1=500, chrX=100))
makeGRangesFromDataFrame(df, seqinfo=Seqinfo(paste0("chr", 4:1)))

## ---------------------------------------------------------------------
## BASIC makeGPosFromDataFrame() EXAMPLES
## ---------------------------------------------------------------------

df <- data.frame(
    rsids=    c("rs7927381", "rs79273813", "rs79273815",
                "rs79273819", "rs79273817", "rs79273811", "rs79273814"),
    seqnames= c("11", "1", "11", "1", "5", "8", "8"),
    positions=c("67579271", "179599437", "70746832",
                "5066529", "107006951", "95077418", "43582071"),
    alleles=  c("T/C", "A/G", "C/G, C/T", "G/A, G/C", "T/C, T/G",
                "G/A, G/T", "G/A, G/T"),
    genes=    c(NA, "TDRD5", "SHANK2", NA, "LOC102467213", "NDUFAF6", NA),
    row.names=letters[1:7]
)

makeGPosFromDataFrame(df)
makeGPosFromDataFrame(df, keep.extra.columns=TRUE)
as(df, "GPos")  # same as above

## With missing positions:
df$positions[2L] <- NA
df
#makeGPosFromDataFrame(df)  # error!
makeGPosFromDataFrame(df, na.rm=TRUE)  # rows with NAs got dropped

## ---------------------------------------------------------------------
## ABOUT AUTOMATIC DETECTION OF THE seqnames/start/end/strand COLUMNS
## ---------------------------------------------------------------------

## Automatic detection of the seqnames/start/end/strand columns is
## case insensitive:
df <- data.frame(ChRoM="chr1", StarT=11:15, stoP=12:16,
                 STRAND=c("+","-","+","*","."), score=1:5)
makeGRangesFromDataFrame(df)

## It also ignores a common prefix between the start and end columns:
df <- data.frame(seqnames="chr1", tx_start=11:15, tx_end=12:16,
                 strand=c("+","-","+","*","."), score=1:5)
makeGRangesFromDataFrame(df)

## The common prefix between the start and end columns is used to
## disambiguate between more than one seqnames column:
df <- data.frame(chrom="chr1", tx_start=11:15, tx_end=12:16,
                 tx_chr="chr2", score=1:5)
makeGRangesFromDataFrame(df)

## ---------------------------------------------------------------------
## 0-BASED VS 1-BASED START POSITIONS
## ---------------------------------------------------------------------

if (require(rtracklayer)) {
  session <- browserSession()
  genome(session) <- "sacCer2"
  query <- ucscTableQuery(session, "Assembly")
  df <- getTable(query)
  head(df)

  ## A common pitfall is to forget that the UCSC Table Browser uses the
  ## "0-based start" convention:
  gr0 <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE,
                                      start.field="chromStart",
                                      end.field="chromEnd")
  head(gr0)

  ## The start positions need to be converted into 1-based positions,
  ## to adhere to the convention used in Bioconductor:
  gr1 <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE,
                                      start.field="chromStart",
                                      end.field="chromEnd",
                                      starts.in.df.are.0based=TRUE)
  head(gr1)
}
}

\keyword{manip}

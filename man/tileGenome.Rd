\name{tileGenome}

\alias{tileGenome}

\title{Put (virtual) tiles on a given genome}

\description{
  \code{tileGenome} returns a set of genomic regions that form a
  partitioning of the specified genome. Each region is called a "tile".
}

\usage{
tileGenome(seqlengths, ntile, tilewidth, cut.last.tile.in.chrom=FALSE)
}

\arguments{
  \item{seqlengths}{
    Either a named numeric vector of chromosome lengths or a \link{Seqinfo}
    object. More precisely, if a named numeric vector, it must have a length
    >= 1, cannot contain NAs or negative values, and cannot have duplicated
    names. If a \link{Seqinfo} object, then it's first replaced with the
    vector of sequence lengths stored in the object (extracted from the object
    with the \code{\link{seqlengths}} getter), then the restrictions described
    previously apply to this vector.
  }
  \item{ntile}{
    The number of tiles to generate.
  }
  \item{tilewidth}{
    The desired tile width. The effective tile width might be slightly
    different but is guaranteed to never be more than the desired width.
  }
  \item{cut.last.tile.in.chrom}{
    Whether or not to cut the last tile in each chromosome.
    This is set to \code{FALSE} by default.
    Can be set to \code{TRUE} only when \code{tilewidth} is specified.
    In that case, a tile will never overlap with more than 1 chromosome
    and a \link{GRanges} object is returned with one element (i.e. one
    genomic range) per tile. 
  }
}

\value{
  If \code{cut.last.tile.in.chrom} is \code{FALSE} (the default),
  a \link{GRangesList} object with one list element per tile, each of
  them containing a number of genomic ranges equal to the number of
  chromosomes it overlaps with. Note that when the tiles are small (i.e.
  much smaller than the chromosomes), most of them only overlap with a
  single chromosome.

  A \link{GRanges} object if \code{cut.last.tile.in.chrom} is \code{TRUE},
  with one element (i.e. one genomic range) per tile.
}

\author{
  H. Pages, based on a proposal by Martin Morgan
}

\seealso{
  \itemize{
    \item \link{GRangesList} and \link{GRanges} objects.

    \item \link{Seqinfo} objects and the \code{\link{seqlengths}} getter.

    \item \link[IRanges]{IntegerList} objects.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. WITH A TOY GENOME
## ---------------------------------------------------------------------

seqlengths <- c(chr1=60, chr2=20, chr3=25)

## Create 5 tiles:
tiles <- tileGenome(seqlengths, ntile=5)
tiles
elementLengths(tiles)  # tiles 3 and 4 contain 2 ranges

width(tiles)
## Use sum() on this IntegerList object to get the effective tile
## widths:
sum(width(tiles))  # each tile covers exactly 21 genomic positions

## Create 9 tiles:
tiles <- tileGenome(seqlengths, ntile=9)
elementLengths(tiles)  # tiles 6 and 7 contain 2 ranges

table(sum(width(tiles)))  # some tiles cover 12 genomic positions,
                          # others 11

## Specify the tile width:
tiles <- tileGenome(seqlengths, tilewidth=20)
length(tiles)  # 6 tiles
table(sum(width(tiles)))  # effective tile width is <= specified

## Specify the tile width and cut the last tile in each chromosome:
tiles <- tileGenome(seqlengths, tilewidth=24,
                    cut.last.tile.in.chrom=TRUE)
tiles
width(tiles)  # each tile covers exactly 24 genomic positions, except
              # the last tile in each chromosome

## Partition a genome by chromosome ("natural partitioning"):
tiles <- tileGenome(seqlengths, tilewidth=max(seqlengths),
                    cut.last.tile.in.chrom=TRUE)
tiles  # one tile per chromosome

## sanity check
stopifnot(all.equal(setNames(end(tiles), seqnames(tiles)), seqlengths))

## ---------------------------------------------------------------------
## B. WITH A REAL GENOME
## ---------------------------------------------------------------------

library(BSgenome.Scerevisiae.UCSC.sacCer2)
tiles <- tileGenome(seqinfo(Scerevisiae), ntile=20)
tiles

tiles <- tileGenome(seqinfo(Scerevisiae), tilewidth=50000,
                    cut.last.tile.in.chrom=TRUE)
tiles

## ---------------------------------------------------------------------
## C. AN APPLICATION: COMPUTE THE AVERAGE PER BIN OF NUMERICAL VARIABLES
##    DEFINED ALONG A GENOME
## ---------------------------------------------------------------------

## The metadata columns of a GRanges object can be used to store
## variables defined along a genome. Let's create such a GRanges:

library(BSgenome.Scerevisiae.UCSC.sacCer2)

x1 <- GRanges("chrI",
              IRanges(c(1, 221, 251), c(150, 295, 360)),
              score=c(0.4, 8, -10),
              seqinfo=seqinfo(Scerevisiae))
x1

## The score metadata column can be seen as a variable defined along the
## Scerevisiae genome. For some applications, there is sometimes the
## need to compute the average score for each genomic region in a set of
## predefined fixed-width regions (sometimes called "bins").
## Let's use tileGenome() to create such a set of bins:

bins1 <- tileGenome(seqinfo(x1), tilewidth=100,
                    cut.last.tile.in.chrom=TRUE)

## The function we're going to use to compute the average score per bin
## can be defined as follow:
##   'x': a GRanges object with metadata columns representing numerical
##        variables defined along the genome 'x' is based on (i.e. the
##        genome described by 'seqinfo(x)').
##   'mcolnames': the column names of the numerical variables in
##        'mcols(x)' for which to compute the average per bin.
##   'bins': a GRanges object representing a partitioning of the genome
##        'x' is based on. Typically obtained by calling tileGenome()
##        with 'cut.last.tile.in.chrom=TRUE' on 'seqinfo(x)'.
averagePerBin <- function(x, mcolnames, bins)
{
    stopifnot(identical(seqinfo(x), seqinfo(bins)))
    if (is.null(mcolnames))
        return(bins)
    bins_per_chrom <- as(bins, "RangesList")
    averageMCol <- function(numvar)
    {
        cvg <- coverage(x, weight=numvar)
        views_list <- RleViewsList(
            lapply(names(cvg),
                   function(seqname)
                       Views(cvg[[seqname]], bins_per_chrom[[seqname]])))
        unlist(viewMeans(views_list), use.names=FALSE)
    }
    mcols(bins) <- DataFrame(lapply(mcols(x), averageMCol))
    bins
}

## Compute the average score per bin:

average_per_bin1 <- averagePerBin(x1, "score", bins1)
average_per_bin1

## Here is another example with more than one numerical variables:

makeRandomGRanges <- function(seqinfo, density, width)
{
    ans_seqlevels <- seqnames(seqinfo)
    ans_start <- lapply(ans_seqlevels,
                        function(name) {
                            seqlength <- seqlengths(seqinfo)[name]
                            sample(seqlength - width + 1L,
                                   seqlength * density,
                                   replace=TRUE)
                        })
    ans_seqnames <- Rle(factor(ans_seqlevels, levels=ans_seqlevels),
                        elementLengths(ans_start))
    ans_ranges <- IRanges(start=unlist(ans_start), width=width)
    GRanges(seqnames=ans_seqnames, ranges=ans_ranges, seqinfo=seqinfo)
}

set.seed(22)
x2 <- makeRandomGRanges(seqinfo(Scerevisiae), 0.01, 20)
mcols(x2) <- DataFrame(score=runif(length(x2)),
                       quality=sample(40, length(x2), replace=TRUE))
x2

bins2 <- tileGenome(seqinfo(x2), tilewidth=50000,
                    cut.last.tile.in.chrom=TRUE)

average_per_bin2 <- averagePerBin(x2, c("score", "quality"), bins2)
average_per_bin2
}

\keyword{manip}

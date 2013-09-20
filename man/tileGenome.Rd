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
## C. AN APPLICATION: COMPUTE THE BINNED AVERAGE OF A NUMERICAL VARIABLE
##    DEFINED ALONG A GENOME
## ---------------------------------------------------------------------

## 1. When the variable is stored in a named RleList object
## --------------------------------------------------------

## In Bioconductor, a variable defined along a genome is typically
## represented as a named RleList object with one list element per
## chromosome. Let's create such a variable:

library(BSgenome.Scerevisiae.UCSC.sacCer2)
set.seed(22)
my_var1 <- RleList(
    lapply(seqlengths(Scerevisiae),
           function(len) Rle(sample(-10:10, len, replace=TRUE))),
    compress=FALSE)
my_var1

## In some applications, there is sometimes the need to compute the
## average of 'my_var1' for each genomic region in a set of predefined
## fixed-width regions (sometimes called "bins"). Let's use
## tileGenome() to create such a set of bins:

bins1 <- tileGenome(seqinfo(Scerevisiae), tilewidth=100,
                    cut.last.tile.in.chrom=TRUE)

## We define the following function to compute the binned average of a
## numerical variable defined along a genome.
## Arguments:
##   'bins': a GRanges object representing the genomic bins.
##        Typically obtained by calling tileGenome() with
##        'cut.last.tile.in.chrom=TRUE'.
##   'numvar': a named RleList object representing a numerical
##        variable defined along the genome covered by 'bins', which
##        is the genome described by 'seqinfo(bins)'.
##   'mcolname': the name to give to the metadata column that will
##        contain the binned average in the returned object.
## Returns 'bins' with an additional metadata column named 'mcolname'
## containing the binned average.
binnedAverage <- function(bins, numvar, mcolname)
{
    stopifnot(is(bins, "GRanges"))
    stopifnot(is(numvar, "RleList"))
    stopifnot(identical(seqlevels(bins), names(numvar)))
    bins_per_chrom <- split(bins, seqnames(bins))
    means_list <- lapply(names(numvar),
        function(seqname) {
            views <- Views(numvar[[seqname]],
                           ranges(bins_per_chrom[[seqname]]))
            viewMeans(views)
        })
    new_mcol <- unsplit(means_list, as.factor(seqnames(bins)))
    mcols(bins)[[mcolname]] <- new_mcol
    bins
}

## Compute the binned average for 'my_var1':

bins1 <- binnedAverage(bins1, my_var1, "binned_var1")
bins1

## 2. When the variable is stored in a metadata column of a disjoint
##    GRanges object
## -----------------------------------------------------------------

## A GRanges object is said to be disjoint if it contains ranges
## that do not overlap with each other. This can be tested with the
## isDisjoint() function. For example, the GRanges object returned
## by tileGenome() is always guaranteed to be disjoint:

isDisjoint(bins1)

## In addition to named RleList objects, the metadata columns of a
## disjoint GRanges object can also be seen as variables defined
## along a genome. An obvious example is the "binned_var1" metadata
## column in 'bins1'. Another example is the "score" metadata column
## in the following GRanges object:

x2 <- GRanges("chrI",
              IRanges(c(1, 211, 291), c(150, 285, 377)),
              score=c(0.4, 8, -10),
              seqinfo=seqinfo(Scerevisiae))
x2

## Turning the score variable into a named RleList representation
## can be done by computing the weighted coverage of 'x2':

score <- coverage(x2, weight="score")
score

## Now we can pass 'score' to binnedAverage() to compute the average
## score per bin:

bins1 <- binnedAverage(bins1, score, "binned_score")
bins1

## With bigger bins:

bins2 <- tileGenome(seqinfo(x2), tilewidth=50000,
                    cut.last.tile.in.chrom=TRUE)
bins2 <- binnedAverage(bins2, score, "binned_score")
bins2

## Note that the binned variables in 'bins1' and 'bins2' can be
## turned back into named RleList objects:

binned_var1 <- coverage(bins1, weight="binned_var1")
stopifnot(all.equal(mean(binned_var1), mean(my_var1)))

binned_score <- coverage(bins2, weight="binned_score")
stopifnot(all.equal(mean(binned_score), mean(score)))

## Of course the "binned" variables are much more compact in memory
## than the original variables (they contain much less runs):

object.size(binned_var1)
object.size(my_var1)
}

\keyword{manip}

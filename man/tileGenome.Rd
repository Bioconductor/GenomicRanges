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

  If \code{cut.last.tile.in.chrom} is \code{TRUE}, a \link{GRanges}
  object with one element (i.e. one genomic range) per tile.
}

\author{
  H. Pagès, based on a proposal by M. Morgan
}

\seealso{
  \itemize{
    \item \link{genomicvars} for an example of how to compute the binned
          average of a numerical variable defined along a genome.

    \item \link{GRangesList} and \link{GRanges} objects.

    \item \link{Seqinfo} objects and the \code{\link{seqlengths}} getter.

    \item \link[IRanges]{IntegerList} objects.

    \item \link[IRanges]{Views} objects.
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
elementNROWS(tiles)  # tiles 3 and 4 contain 2 ranges

width(tiles)
## Use sum() on this IntegerList object to get the effective tile
## widths:
sum(width(tiles))  # each tile covers exactly 21 genomic positions

## Create 9 tiles:
tiles <- tileGenome(seqlengths, ntile=9)
elementNROWS(tiles)  # tiles 6 and 7 contain 2 ranges

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

## See '?genomicvars' for an example of how to compute the binned
## average of a numerical variable defined along a genome.
}

\keyword{manip}

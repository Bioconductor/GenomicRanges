\name{tile-methods}
\alias{tile-methods}

\alias{tile,GenomicRanges-method}

\title{Tile a GenomicRanges object}

\description{
  \code{\link[IRanges]{tile}} method for
  \link{GenomicRanges}. Partitions each range into a set of
  tiles. Tiles are defined in terms of their number or width.
}

\usage{
\S4method{tile}{GenomicRanges}(x, n, width)
}

\arguments{
  \item{x}{
    A \link{GenomicRanges} object, like a \code{GRanges}.
  }
  \item{n}{
    The number of tiles to generate.
    See \code{?\link[IRanges]{tile}} in the \pkg{IRanges} package
    for more information about this argument.
  }
  \item{width}{
    The (maximum) width of each tile.
    See \code{?\link[IRanges]{tile}} in the \pkg{IRanges} package
    for more information about this argument.
  }
}

\details{
  Splits \code{x} into a \code{GRangesList}, each element of which corresponds
  to a tile, or partition, of \code{x}. Specify the tile geometry with either
  \code{n} or \code{width} (not both). Passing \code{n} creates \code{n} tiles
  of approximately equal width, truncated by sequence end, while passing
  \code{width} tiles the region with ranges of the given width, again truncated
  by sequence end.
}

\value{
  A \code{GRangesList} object, each element of which corresponds to a tile.
}

\author{M. Lawrence}

\seealso{
  \code{\link[IRanges]{tile}} in the \pkg{IRanges} package.
}

\examples{
gr <- GRanges(
        seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges=IRanges(1:10, end=11),
        strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        seqlengths=c(chr1=11, chr2=12, chr3=13))

# split every range in half
tiles <- tile(gr, n = 2L)
stopifnot(all(elementNROWS(tiles) == 2L))

# split ranges into subranges of width 2
# odd width ranges must contain one subrange of width 1
tiles <- tile(gr, width = 2L)
stopifnot(all(all(width(tiles) \%in\% c(1L, 2L))))
}

\keyword{methods}
\keyword{utilities}

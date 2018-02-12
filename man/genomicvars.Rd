\name{genomicvars}

\alias{genomicvars}
\alias{genomicvariables}

\alias{coerce,RleList,GRanges-method}
\alias{coerce,RleViewsList,GRanges-method}
\alias{bindAsGRanges}
\alias{mcolAsRleList}
\alias{binnedAverage}


\title{Manipulating genomic variables}

\description{
  A \emph{genomic variable} is a variable defined along a genome.
  Here are 2 ways a genomic variable is generally represented in
  Bioconductor:
  \enumerate{
    \item as a named \link[IRanges]{RleList} object with one list element
          per chromosome;
    \item as a metadata column on a \emph{disjoint} GRanges object.
  }
  This man page documents tools for switching from one form to the other.
}

\usage{
bindAsGRanges(...)
mcolAsRleList(x, varname)
binnedAverage(bins, numvar, varname, na.rm=FALSE)
}

\arguments{
  \item{...}{
    One or more genomic variables in the form of named RleList objects.
  }
  \item{x}{
    A \emph{disjoint} \link{GRanges} object with metadata columns on it.
    A \link{GRanges} object is said to be \emph{disjoint} if it contains
    ranges that do not overlap with each other. This can be tested with
    \code{isDisjoint}. See \code{?`\link{isDisjoint,GenomicRanges-method}`}
    for more information about the \code{isDisjoint} method for
    \link{GRanges} objects.
  }
  \item{varname}{
    The name of the genomic variable.

    For \code{mcolAsRleList} this must be the name of the metadata column
    on \code{x} to be turned into an \link[IRanges]{RleList} object.

    For \code{binnedAverage} this will be the name of the metadata column
    that contains the binned average in the returned object.
  }
  \item{bins}{
    A \link{GRanges} object representing the genomic bins. Typically
    obtained by calling \code{\link{tileGenome}} with
    \code{cut.last.tile.in.chrom=TRUE}.
  }
  \item{numvar}{
    A named \link[IRanges]{RleList} object representing a numerical variable
    defined along the genome covered by \code{bins} (which is the genome
    described by \code{seqinfo(bins)}).
  }
  \item{na.rm}{
    A logical value indicating whether \code{NA} values should be
    stripped before the average is computed.
  }
}

\details{
  \code{bindAsGRanges} allows to switch the representation of one or
  more genomic variables from the \emph{named RleList} form to the
  \emph{metadata column on a disjoint GRanges object} form by binding
  the supplied named RleList objects together and putting them
  on the same \link{GRanges} object. This transformation is lossless.

  \code{mcolAsRleList} performs the opposite transformation and is also
  lossless (however the circularity flags and genome information in
  \code{seqinfo(x)} won't propagate). It works for any metadata column on
  \code{x} that can be put in \link[S4Vectors]{Rle} form i.e. that is an
  atomic vector or a factor.

  \code{binnedAverage} computes the binned average of a numerical variable
  defined along a genome.
}

\value{
  For \code{bindAsGRanges}: a \link{GRanges} object with 1 metadata column
  per supplied genomic variable.

  For \code{mcolAsRleList}: a named \link[IRanges]{RleList} object with
  1 list element per seqlevel in \code{x}.

  For \code{binnedAverage}: input \link{GRanges} object \code{bins} with
  an additional metadata column named \code{varname} containing the binned
  average.
}

\author{
  H. Pagès
}

\seealso{
  \itemize{
    \item \link[IRanges]{RleList} objects in the \pkg{IRanges} package.

    \item \link{coverage,GenomicRanges-method} for computing the coverage
          of a \link{GRanges} object.

    \item The \code{\link{tileGenome}} function for putting tiles on a
          genome.

    \item \link{GRanges} objects and \link{isDisjoint,GenomicRanges-method}
          for the \code{isDisjoint} method for \link{GenomicRanges} objects.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. TWO WAYS TO REPRESENT A GENOMIC VARIABLE
## -----------------------------------------------------------------

## 1) As a named RleList object
## ----------------------------
## Let's create a genomic variable in the "named RleList" form:
library(BSgenome.Scerevisiae.UCSC.sacCer2)
set.seed(55)
my_var <- RleList(
    lapply(seqlengths(Scerevisiae),
        function(seqlen) {
            tmp <- sample(50L, seqlen, replace=TRUE) %/% 50L
            Rle(cumsum(tmp - rev(tmp)))
        }
    ),
    compress=FALSE)
my_var

## 2) As a metadata column on a disjoint GRanges object
## ----------------------------------------------------
gr1 <- bindAsGRanges(my_var=my_var)
gr1

gr2 <- GRanges(c("chrI:1-150",
                 "chrI:211-285",
                 "chrI:291-377",
                 "chrV:51-60"),
               score=c(0.4, 8, -10, 2.2),
               id=letters[1:4],
               seqinfo=seqinfo(Scerevisiae))
gr2

## Going back to the "named RleList" form:
mcolAsRleList(gr1, "my_var")
score <- mcolAsRleList(gr2, "score")
score
id <- mcolAsRleList(gr2, "id")
id
bindAsGRanges(score=score, id=id)

## Bind 'my_var', 'score', and 'id' together:
gr3 <- bindAsGRanges(my_var=my_var, score=score, id=id)

## Sanity checks:
stopifnot(identical(my_var, mcolAsRleList(gr3, "my_var")))
stopifnot(identical(score, mcolAsRleList(gr3, "score")))
stopifnot(identical(id, mcolAsRleList(gr3, "id")))
gr2b <- bindAsGRanges(score=score, id=id)
seqinfo(gr2b) <- seqinfo(gr2)
stopifnot(identical(gr2, gr2b))

## ---------------------------------------------------------------------
## B. BIND TOGETHER THE COVERAGES OF SEVERAL BAM FILES
## ---------------------------------------------------------------------

library(pasillaBamSubset)
library(GenomicAlignments)
untreated1_cvg <- coverage(BamFile(untreated1_chr4()))
untreated3_cvg <- coverage(BamFile(untreated3_chr4()))
all_cvg <- bindAsGRanges(untreated1=untreated1_cvg,
                         untreated3=untreated3_cvg)

## Keep regions with coverage:
all_cvg[with(mcols(all_cvg), untreated1 + untreated3 >= 1)]

## Plot the coverage profiles with the Gviz package:
library(Gviz)
plotNumvars <- function(numvars, region, name="numvars", ...)
{
    stopifnot(is(numvars, "GRanges"))
    stopifnot(is(region, "GRanges"), length(region) == 1L)
    gtrack <- GenomeAxisTrack()
    dtrack <- DataTrack(numvars,
                        chromosome=as.character(seqnames(region)),
                        name=name,
                        groups=colnames(mcols(numvars)), type="l", ...)
    plotTracks(list(gtrack, dtrack), from=start(region), to=end(region))
}
plotNumvars(all_cvg, GRanges("chr4:1-25000"),
            "coverage", col=c("red", "blue"))
plotNumvars(all_cvg, GRanges("chr4:1.03e6-1.08e6"),
            "coverage", col=c("red", "blue"))

## Sanity checks:
stopifnot(identical(untreated1_cvg, mcolAsRleList(all_cvg, "untreated1")))
stopifnot(identical(untreated3_cvg, mcolAsRleList(all_cvg, "untreated3")))

## ---------------------------------------------------------------------
## C. COMPUTE THE BINNED AVERAGE OF A NUMERICAL VARIABLE DEFINED ALONG A
##    GENOME
## ---------------------------------------------------------------------

## In some applications (e.g. visualization), there is the need to compute
## the average of a genomic variable for a set of predefined fixed-width
## regions (sometimes called "bins").
## Let's use tileGenome() to create such a set of bins:
bins1 <- tileGenome(seqinfo(Scerevisiae), tilewidth=100,
                    cut.last.tile.in.chrom=TRUE)

## Compute the binned average for 'my_var' and 'score':
bins1 <- binnedAverage(bins1, my_var, "binned_var")
bins1
bins1 <- binnedAverage(bins1, score, "binned_score")
bins1

## Binned average in "named RleList" form:
binned_var1 <- mcolAsRleList(bins1, "binned_var")
binned_var1
stopifnot(all.equal(mean(my_var), mean(binned_var1)))  # sanity check

mcolAsRleList(bins1, "binned_score")

## With bigger bins:
bins2 <- tileGenome(seqinfo(Scerevisiae), tilewidth=50000,
                    cut.last.tile.in.chrom=TRUE)
bins2 <- binnedAverage(bins2, my_var, "binned_var")
bins2 <- binnedAverage(bins2, score, "binned_score")
bins2

binned_var2 <- mcolAsRleList(bins2, "binned_var")
binned_var2
stopifnot(all.equal(mean(my_var), mean(binned_var2)))  # sanity check

mcolAsRleList(bins2, "binned_score")

## Not surprisingly, the "binned" variables are much more compact in
## memory than the original variables (they contain much less runs):
object.size(my_var)
object.size(binned_var1)
object.size(binned_var2)

## ---------------------------------------------------------------------
## D. SANITY CHECKS
## ---------------------------------------------------------------------

bins3 <- tileGenome(c(chr1=10, chr2=8), tilewidth=5,
                    cut.last.tile.in.chrom=TRUE)
my_var3 <- RleList(chr1=Rle(c(1:3, NA, 5:7)), chr2=Rle(c(-3, NA, -3, NaN)))
bins3 <- binnedAverage(bins3, my_var3, "binned_var3", na.rm=TRUE)
binned_var3 <- mcols(bins3)$binned_var3
stopifnot(
  identical(mean(my_var3$chr1[1:5], na.rm=TRUE),
            binned_var3[1]),
  identical(mean(c(my_var3$chr1, 0, 0, 0)[6:10], na.rm=TRUE),
            binned_var3[2]),
  #identical(mean(c(my_var3$chr2, 0), na.rm=TRUE),
  #          binned_var3[3]),
  identical(0, binned_var3[4])
)
}

\keyword{manip}

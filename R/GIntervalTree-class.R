#' GIntervalTree class
#' 
#' Defines persistent interval trees for GRanges objects.
#' 
#' 
#' @name GIntervalTree-class
#' @family GIntervalTree
#' 
#' @exportClass GIntervalTree
#' @import GenomicRanges
#' @import BiocGenerics
setClass("GIntervalTree",
         contains="GenomicRanges",
         representation(
           ranges="IntervalForest",
           strand="Rle",
           elementMetadata="DataFrame",
           seqinfo="Seqinfo"),
         prototype(
           strand=Rle(strand()))
)

.valid.GIntervalTree.length <- function(x) {
  n <- length(x@ranges)
  if ((length(strand(x)) != n)
      || (nrow(mcols(x)) != n))
    return("slot lengths are not all equal")
  NULL
}

.valid.GIntervalTree.ranges <- function(x) {
  if (class(x@ranges) != "IntervalForest")
    return("'ranges(x)' must be a IntervalForest instance")
  NULL
}

.valid.GIntervalTree <- function(x) {
  c(.valid.GIntervalTree.length(x),
    .valid.GIntervalTree.ranges(x),
    .valid.GenomicRanges.strand(x),
    .valid.GenomicRanges.mcols(x),
    valid.GenomicRanges.seqinfo(x))
}

setValidity2("GIntervalTree", .valid.GIntervalTree)

#' seqnames accessor 
#' 
#' 
#' @rdname GIntervalTree-class
#' @family GIntervalTree
#' @export
#' @importMethodsFrom GenomicRanges seqnames
setMethod("seqnames", "GIntervalTree", function(x) (x@ranges@partition))

#' ranges accessor
#' 
#' @rdname GIntervalTree-class
#' @family GIntervalTree
#' @export
#' @importMethodsFrom GenomicRanges ranges
setMethod("ranges", "GIntervalTree", function(x) as(x@ranges, "IRanges"))

#' strand accessor
#' 
#' @rdname GIntervalTree-class
#' @family GIntervalTree
#' @export
#' @importMethodsFrom GenomicRanges strand
setMethod("strand", "GIntervalTree", function(x) x@strand)

#' seqinfo accessor
#' 
#' @rdname GIntervalTree-class
#' @family GIntervalTree
#' @export
#' @importMethodsFrom GenomicRanges seqinfo
setMethod("seqinfo", "GIntervalTree", function(x) x@seqinfo)

setMethod("start", "GIntervalTree", function(x, ...) start(x@ranges))
setMethod("end", "GIntervalTree", function(x, ...) end(x@ranges))
setMethod("width", "GIntervalTree", function(x) width(x@ranges))

#' length accessor
#' 
#' 
#' @rdname GIntervalTree-class
#' @family GIntervalTree
#' @export
setMethod("length", "GIntervalTree", function(x) length(x@ranges))

#' construct from GRanges object via coercion
#' 
#' @name as
#' @family GIntervalTree
#' @importClassesFrom GenomicRanges GRanges
setAs("GRanges", "GIntervalTree",
      function(from) {
        if (any(isCircular(from), na.rm=TRUE))
          stop("'GIntervalTree' objects not supported for circular sequences")

        out=new2("GIntervalTree",
                strand=strand(from),
                elementMetadata=mcols(from),
                seqinfo=seqinfo(from),
                ranges=IntervalForest(ranges(from), seqnames(from)),
                check=FALSE)
        out
      }
)

#' constructor function using GRanges object
#' 
#' @family GIntervalTree
#' @export
GIntervalTree <- function(x) {
  as(x, "GIntervalTree")
}

#' coercion from GIntervalTree to GRanges object
#' 
#' @family GIntervalTree
#' @name as
#' @importClassesFrom GenomicRanges GRanges
setAs("GIntervalTree", "GRanges",
      function(from) {
        out=new("GRanges",
                seqnames=(from@ranges@partition),
                strand=strand(from),
                elementMetadata=mcols(from),
                seqinfo=seqinfo(from),
                ranges=ranges(from))
      }
)

#' subsetting
#' 
#' @family GIntervalTree
#' @rdname GIntervalTree-class
#' @export
setMethod("[", "GIntervalTree",
          function(x, i, j, ...) {
            gr <- as(x, "GRanges")[i]
            as(gr, "GIntervalTree")
          })

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
           rngidx="IRanges",
           strand="Rle",
           elementMetadata="DataFrame",
           seqinfo="Seqinfo"),
         prototype(
           strand=Rle(strand()))
)

.valid.GIntervalTree.rngidx <- function(x) {
  x_len <- length(x)
  if (sum(width(x@rngidx)) != x_len)
    return("ivalid 'rngidx' slot")
  NULL
}

.valid.GIntervalTree.ranges <- function(x) {
  if (class(x@ranges) != "IntervalForest")
    return("'x@ranges' must be an IntervalForest instance")
  if (nobj(x@ranges@partitioning) != length(x))
    return("'nobj(x@ranges@partitioning)' must equal 'length(x)'")
  NULL
}

.valid.GIntervalTree <- function(x) {
  c(.valid.GIntervalTree.ranges(x),
    .valid.GIntervalTree.rngidx(x))
}

setValidity2("GIntervalTree", .valid.GIntervalTree)

.GT_getIndex <- function(from) {
  fvals <- as.integer(runValue(seqnames(from)))
  if (IRanges:::isNotSorted(fvals)) {
    flens <- runLength(seqnames(from))
    idx <- IRanges:::orderInteger(fvals)
    rngidx <- successiveIRanges(flens)[idx]
    idx <- IRanges:::orderInteger(start(rngidx))
    return(IRanges(end=cumsum(width(rngidx))[idx], width=width(rngidx)[idx]))
  }
  IRanges()
}


.GT_reorderValue <- function(obj, val, rngidx=NULL)
{
  if (!is(obj, "GIntervalTree") && is.null(rngidx))
    stop("obj must be 'GIntervalTree' object")
  if (is.null(rngidx))
    rngidx <- obj@rngidx  
  if (length(rngidx))
    val <- IRanges:::extractROWS(val, rngidx)
  val
}

#' seqnames accessor 
#' 
#' 
#' @rdname GIntervalTree-class
#' @family GIntervalTree
#' @export
#' @importMethodsFrom GenomicRanges seqnames
setMethod("seqnames", "GIntervalTree",
          function(x) Rle(.GT_reorderValue(x, space(x@ranges))))

#' ranges accessor
#' 
#' @rdname GIntervalTree-class
#' @family GIntervalTree
#' @export
#' @importMethodsFrom GenomicRanges ranges
setMethod("ranges", "GIntervalTree",
          function(x) .GT_reorderValue(x, as(x@ranges, "IRanges")))

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

setMethod("start", "GIntervalTree",
          function(x, ...) .GT_reorderValue(x,start(x@ranges, ...)@unlistData))
setMethod("end", "GIntervalTree",
          function(x, ...) .GT_reorderValue(x, end(x@ranges, ...)@unlistData))
setMethod("width", "GIntervalTree",
          function(x) .GT_reorderValue(x, width(x@ranges)@unlistData))

#' construct from GRanges object via coercion
#' 
#' @name as
#' @family GIntervalTree
#' @importClassesFrom GenomicRanges GRanges
setAs("GRanges", "GIntervalTree",
      function(from) {
        if (any(isCircular(from), na.rm=TRUE))
          stop("'GIntervalTree' objects not supported for circular sequences")

        rl <- split(unname(ranges(from)), seqnames(from))

        new2("GIntervalTree",
              strand=strand(from),
              elementMetadata=mcols(from),
              seqinfo=seqinfo(from),
              ranges=IntervalForest(rl),
              rngidx=.GT_getIndex(from),
              check=FALSE)
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
                seqnames=seqnames(from),
                strand=strand(from),
                elementMetadata=mcols(from),
                seqinfo=seqinfo(from),
                ranges=ranges(from))
        out
      }
)

#' subsetting
#' 
#' @family GIntervalTree
#' @rdname GIntervalTree-class
#' @export
setMethod("[", "GIntervalTree",
          function(x, i, j, ...) {
            gr <- callGeneric(as(x, "GRanges"),i=i, ...)
            as(gr, "GIntervalTree")
          })

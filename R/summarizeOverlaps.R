## -------------------------------------------------------------------------
## summarizeOverlaps methods and related utilities
## -------------------------------------------------------------------------


setGeneric("summarizeOverlaps", signature=c("features", "reads"),
    function(features, reads, mode=Union, ignore.strand=FALSE, ...)
{
    standardGeneric("summarizeOverlaps")
})

## -------------------------------------------------------------------------
## BamFiles and BamViews methods are in Rsamtools
##

## -------------------------------------------------------------------------
## GappedAlignments methods
##

setMethod("summarizeOverlaps", c("GRangesList", "GappedAlignments"),
    function(features, reads, mode, ignore.strand=FALSE, ...)
{
    mode <- match.fun(mode)
    counts <- .dispatchOverlaps(reads, features, mode, ignore.strand)
    SummarizedExperiment(assays=SimpleList(counts=as.matrix(counts)),
        rowData=features, colData=.readsMetadata(reads))
})

setMethod("summarizeOverlaps", c("GRanges", "GappedAlignments"),
    function(features, reads, mode, ignore.strand=FALSE, ...)
{
    mode <- match.fun(mode)
    counts <- .dispatchOverlaps(reads, features, mode, ignore.strand)
    SummarizedExperiment(assays=SimpleList(counts=as.matrix(counts)),
        rowData=features, colData=.readsMetadata(reads))
})

## -------------------------------------------------------------------------
## GappedAlignmentPairs methods
##

setMethod("summarizeOverlaps", c("GRangesList", "GappedAlignmentPairs"),
    function(features, reads, mode, ignore.strand=FALSE, ...)
{
    mode <- match.fun(mode)
    counts <- .dispatchOverlaps(grglist(reads), features, mode, ignore.strand)
    SummarizedExperiment(assays=SimpleList(counts=as.matrix(counts)),
        rowData=features, colData=.readsMetadata(reads))
})

setMethod("summarizeOverlaps", c("GRanges", "GappedAlignmentPairs"),
    function(features, reads, mode, ignore.strand=FALSE, ...)
{
    mode <- match.fun(mode)
    counts <- .dispatchOverlaps(grglist(reads), features, mode, ignore.strand)
    SummarizedExperiment(assays=SimpleList(counts=as.matrix(counts)),
        rowData=features, colData=.readsMetadata(reads))
})

## -------------------------------------------------------------------------
## 'mode' functions 
##

Union <- function(reads, features, ignore.strand=FALSE, ...)
{
    co <- countOverlaps(reads, features, ignore.strand=ignore.strand)
    idx <- co == 1
    if (sum(co == 1) == 0)
        return(integer(length(features)))

    counts <- countOverlaps(features, reads[idx], ignore.strand=ignore.strand)
    names(counts) <- names(features)
    counts 
}

IntersectionStrict <- function(reads, features, ignore.strand = FALSE, ...)
{
    queryseq <- seqlevels(reads)
    circular <- isCircular(features)
    circNames <- intersect(queryseq, names(circular)[circular])
    if (0L != length(circNames)) {
        warning("circular sequence(s) in reads '",
                paste(circNames, sep="' '"), "' ignored")
        if (any(keep <- !seqlevels(reads) %in% circNames))
            reads <- keepSeqlevels(reads, seqlevels(reads)[keep])
        else
            return(integer(length(features)))
    }
 
    fo <- findOverlaps(reads, features, type="within",
        ignore.strand=ignore.strand)
    ## omit reads that hit >1 feature
    co <- tabulate(queryHits(fo), length(reads))
    if (!any(co == 1L))
        return(integer(length(features)))
    unqfo <- fo[queryHits(fo) %in% seq_along(reads)[co == 1L]]
    rle <- Rle(sort(subjectHits(unqfo)))
    counts <- rep(0, length(features))
    counts[runValue(rle)] <- runLength(rle)
    names(counts) <- names(features)
    counts 
}

IntersectionNotEmpty <-  function(reads, features, ignore.strand = FALSE, ...)
{
    co <- countOverlaps(reads, features, ignore.strand=ignore.strand)
    if (sum(co) == 0)
        return(integer(length(features)))

    counts <- .IntersectionNotEmpty(reads, features, ignore.strand)
    names(counts) <- names(features)
    counts 
}

## -------------------------------------------------------------------------
## non-exported helpers 
##

.readsMetadata <- function(reads)
{
    if (length(metadata(reads)) > 0)
        DataFrame(metaData = metadata(reads))
    else
        DataFrame(metaData = character(1))
}

.dispatchOverlaps <-
    function(reads, features, mode, ignore.strand, ...)
{
    if (ignore.strand) {
        if (class(reads) == "GRangesList") {
            r <- unlist(reads)
            strand(r) <- "*"
            reads <- split(r, togroup(reads))
        } else {
            strand(reads) <- "*"
        } 
    }
    mode(reads, features, ignore.strand)
}

.IntersectionNotEmpty <- function(reads, features, ignore.strand=FALSE)
{
    ## unique disjoint regions
    if (class(features) == "GRangesList")
        d <- disjoin(features@unlistData)
    else
        d <- disjoin(features)
    coUnq <- countOverlaps(d, features, ignore.strand=ignore.strand)
    ud <- d[coUnq == 1]

    ## count read if :
    ## (i)  the read or any one of the read fragments (gapped reads) olaps 
    ##      a ud region
    ## (ii) the read or > 1 of the read fragments (gapped reads) olaps >1
    ##      ud region AND the ud regions are from the same feature

    foUD <- findOverlaps(reads, ud, ignore.strand=ignore.strand)
    ## map ud regions back to original features
    foFeatures <- findOverlaps(features, ud, ignore.strand=ignore.strand)
    featuresMap <-
        as.matrix(foFeatures)[order(subjectHits(foFeatures)), 1L]
    backMapFeatures <- featuresMap[subjectHits(foUD)]
    mm <- data.frame(query=queryHits(foUD), subject=backMapFeatures)

    queryRle <- Rle(mm$query)
    qsingle <- runValue(queryRle)[runLength(queryRle) == 1]
    singlehits <- mm$subject[mm$query %in% qsingle]
    qmulti <- runValue(queryRle)[runLength(queryRle) > 1]
    multi <- mm[mm$query %in% qmulti, ]
    lst <- split(multi$subject, multi$query)
    unq <- lapply(lst, function(x) {
             if (length(unique(x)) == 1)
                 unique(x)
             else
                 NA}
           )

    multihits <- do.call(c, unq)
    regions <- c(singlehits, multihits) 

    counts <- rep(0, length(features))
    if (length(regions) == 0)
        return(counts)
    countsRle <- Rle(sort(regions))
    counts[runValue(countsRle)] <- runLength(countsRle) 
    counts
}

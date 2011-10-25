setGeneric("summarizeOverlaps", signature = c("features", "reads"),
    function(features, reads, 
             mode = Union,
             ignore.strand = FALSE, ..., param = ScanBamParam())
{
    standardGeneric("summarizeOverlaps")
})


## summarizeOverlaps methods for BamFiles are found upstairs next to home
## furnishings.  Actually that's not quite right, they are in Rsamtools.


## methods for GappedAlignments
setMethod("summarizeOverlaps", c("GRangesList", "GappedAlignments"),
    function(features, reads, 
             mode, 
             ignore.strand = FALSE, ...)
{
    mode <- match.fun(mode)
    counts <- .dispatch(reads, features, mode, ignore.strand)
    if (length(metadata(reads)) > 0)
        colData <- DataFrame(metaData = metadata(reads))
    else
        colData <- DataFrame(metaData = character(1))
    SummarizedExperiment(assays=SimpleList(counts=as.matrix(counts)),
                         rowData=features, colData=colData)
})

setMethod("summarizeOverlaps", c("GRanges", "GappedAlignments"),
    function(features, reads, 
             mode, 
             ignore.strand = FALSE, ...)
{
    mode <- match.fun(mode)
    counts <- .dispatch(reads, features, mode, ignore.strand)
    if (length(metadata(reads)) > 0)
        colData <- DataFrame(metaData = metadata(reads))
    else
        colData <- DataFrame(metaData = character(1))
    SummarizedExperiment(assays=SimpleList(counts=as.matrix(counts)),
                         rowData=features, colData=colData)
})


.dispatch <- function(reads, features, mode, ignore.strand, ...)
{
    if (ignore.strand)
           strand(reads) <- "*"
    mode(reads, features, ignore.strand)
}

Union <- function(reads, features, ignore.strand=FALSE, ...)
{
    co <- countOverlaps(reads, features, ignore.strand=ignore.strand)
    idx <- co == 1
    if (sum(co == 1) == 0)
        return(integer(length(features)))

    ## gapped reads are treated the same as simple reads for Union 
    simplect <- countOverlaps(features, reads[idx], ignore.strand=ignore.strand)
    counts <- simplect
    names(counts) <- names(features)
    counts 
}

IntersectionStrict <- function(reads, features, ignore.strand = FALSE, ...)
{
    co <- countOverlaps(reads, features, type="within",
        ignore.strand=ignore.strand)
    idx <- co == 1
    if (sum(co == 1) == 0)
        return(integer(length(features)))
    ngaps <- ngap(reads)
    if (sum(ngaps) == 0) {
        gapct <- integer(length(features))
    } else {
        gaps <- reads[ngaps != 0]
        gapct <- .gappedIntersectionStrict(gaps, features,
            ignore.strand=ignore.strand) 
        reads <- reads[ngaps == 0]
        idx <- idx[ngaps == 0]
    }

    fo <- findOverlaps(reads[idx], features, type="within",
        ignore.strand=ignore.strand)
    rle <- Rle(sort(subjectHits(fo)))
    simplect <- rep(0, length(features))
    simplect[runValue(rle)] <- runLength(rle)
    counts <- simplect + gapct
    names(counts) <- names(features)
    counts
}

IntersectionNotEmpty <-  function(reads, features, ignore.strand = FALSE, ...)
{
    co <- countOverlaps(reads, features, ignore.strand=ignore.strand)
    if (sum(co) == 0)
        return(integer(length(features)))
    ngaps <- ngap(reads)
    if (sum(ngaps) == 0) {
       gapct <- integer(length(features))
    } else {
        gaps <- as(reads[ngaps != 0], "GRangesList")
        gapct <- .IntersectionNotEmpty(gaps, features, ignore.strand) 
        reads <- reads[ngaps == 0]
    }

    simplect <- .IntersectionNotEmpty(reads, features, ignore.strand)
    counts <- simplect + gapct
    names(counts) <- names(features)
    counts 
}

.IntersectionNotEmpty <- function(reads, features, ignore.strand=FALSE)
{
    ## disjoint regions
    if (class(features) == "GRangesList")
        d <- disjoin(features@unlistData)
    else
        d <- disjoin(features)

    ## unique disjoint regions (ie, remove regions shared by multiple features)
    coUnq <- countOverlaps(d, features, ignore.strand=ignore.strand)
    ud <- d[coUnq == 1]

    ## count read if :
    ## (i)  the read or exactly 1 of the read fragments (gapped reads) olaps 
    ##      a ud region
    ## (ii) the read or > 1 of the read fragments (gapped reads) olaps >1
    ##      ud region AND the ud regions are from the same feature

    ## overlap between reads and ud regions
    ## create mapping of read fragments back to original reads
    if (class(reads) == "GRangesList") {
        readMap <- rep(seq_len(length(reads)), elementLengths(reads))
        foUD <- findOverlaps(reads@unlistData, ud, ignore.strand=ignore.strand)
    } else {
        foUD <- findOverlaps(reads, ud, ignore.strand=ignore.strand)
        readMap <- seq_len(length(reads))
    }

    ## create mapping of ud regions back to original features
    foFeatures <- findOverlaps(features, ud, ignore.strand=ignore.strand)
    featuresMap <-
        matchMatrix(foFeatures)[order(subjectHits(foFeatures)),"query"]

    ## map read fragments and ud regions back to original
    backMapReads <- readMap[queryHits(foUD)]
    backMapFeatures <- featuresMap[subjectHits(foUD)]

    ## ud regions hit by read fragments must be from same feature
    olaplst <- split(backMapFeatures, backMapReads)
    counts <- rep(0, length(features))

    unq <- lapply(olaplst, unique)
    len <- lapply(unq, length)
    regions <- unq[len == 1]
    if (length(regions) == 0)
        return(counts)

    ## counts
    rle <- Rle(sort(unlist(regions, use.names=FALSE)))
    counts[runValue(rle)] <- runLength(rle) 
    counts
}

.gappedIntersectionStrict <- function(reads, features, ignore.strand=FALSE)
{
    grl <- as(reads, "GRangesList")
    co <- countOverlaps(grl, features, type="within",
        ignore.strand=ignore.strand)

    ## if all parts of gapped read fall within == 1 feature, keep 
    keep <- co == 1
    clean <- grl[keep]
    counts <- rep(0, length(features))
    if (length(grl) > 0) {
        fo <- findOverlaps(clean, features, type="within",
            ignore.strand=ignore.strand)
        shits <- sort(subjectHits(fo))
        rle <- Rle(shits)
        counts[runValue(rle)] <- runLength(rle) 
    } 
    counts
}

.isEqual <- function(x)
{
    diff(range(x)) < .Machine$double.eps ^ 0.5
}

#.gappedUnion <- function(reads, features, ignore.strand=FALSE)
#{
#    grl <- as(reads, "GRangesList")
#    co <- countOverlaps(grl@unlistData, features, ignore.strand=ignore.strand)
#    ## if any portion of the read hits > 1 feature, remove entire read
#    idx <- co != 1
#    map <- rep(seq_len(length(grl)), elementLengths(grl))
#    if (any(idx))
#        clean <- grl[-map[idx]]
#    else
#        clean <- grl
#    countOverlaps(features, clean, ignore.strand=ignore.strand)
#}

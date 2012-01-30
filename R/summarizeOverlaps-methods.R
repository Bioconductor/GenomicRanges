setGeneric("summarizeOverlaps", signature = c("features", "reads"),
    function(features, reads, mode = Union, ignore.strand = FALSE, 
             ..., param = ScanBamParam())
{
    standardGeneric("summarizeOverlaps")
})

## methods for BamFiles and BamViews are in Rsamtools

## methods for GappedAlignments
setMethod("summarizeOverlaps", c("GRangesList", "GappedAlignments"),
    function(features, reads, mode, ignore.strand = FALSE, ...)
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
    function(features, reads, mode, ignore.strand = FALSE, ...)
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

    counts <- countOverlaps(features, reads[idx], ignore.strand=ignore.strand)
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

    fo <- findOverlaps(reads[idx], features, type="within",
        ignore.strand=ignore.strand)
    rle <- Rle(sort(subjectHits(fo)))
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
        matchMatrix(foFeatures)[order(subjectHits(foFeatures)),"query"]
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

.isEqual <- function(x)
{
    diff(range(x)) < .Machine$double.eps ^ 0.5
}


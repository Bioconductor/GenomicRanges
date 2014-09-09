library(digest)

.singleDispatch <-
    c("coverage", "disjointBins", "duplicated", "end", "end<-",
      "flank", "granges", "isDisjoint", "narrow", "ranges",
      "resize", "restrict", "seqinfo", "seqnames", "shift",
      "start", "start<-", "strand", "width", "width<-",
      "strand", "ranges", "seqinfo<-")

.twoDispatch <- c("nearest", "precede", "follow", "distance",
           "distanceToNearest",
           ## FIXME: "findOverlaps", "countOverlaps",
           "subsetByOverlaps")

.otherFuns <- c("order", "rank", "sort")

m <- matrix(1, 5, 3, dimnames=list(NULL, NULL))
mlst <- matrix(1, 3, 3, dimnames=list(NULL, NULL))
mList <- list(m, mlst)
assaysList <- list(gr=SimpleList(m=m), grl=SimpleList(m=mlst))
rowDataList <- 
    list(gr=GRanges("chr1", IRanges(1:5, 10)), 
         grl=split(GRanges("chr1", IRanges(1:5, 10)), c(1,1,2,2,3)))
names(rowDataList[["grl"]]) <- NULL
colData <- DataFrame(x=letters[1:3])

## a list of one SE with GRanges and one with GRangesList
ssetList <- 
    list(SummarizedExperiment(
           assays=assaysList[["gr"]], 
           rowData=rowDataList[["gr"]], 
           colData=colData),
         SummarizedExperiment(
           assays=assaysList[["grl"]], 
           rowData=rowDataList[["grl"]], 
           colData=colData))

test_SummarizedExperiment_GRanges_API <- function() {
    ## are we targetting the correct API? signature for
    ## SummarizedExperiment method should match signature for
    ## GenomicRanges or similar, as in each test below

    for (.fun in .singleDispatch) {
        generic <- getGeneric(.fun)
        method <- getMethod(.fun, "SummarizedExperiment")
        checkIdentical("x", generic@signature)
        checkIdentical(formals(generic@.Data), formals(method@.Data))
    }

    .sig <- c("SummarizedExperiment", "SummarizedExperiment")
    .targ <- c("GenomicRanges", "GenomicRanges")
    for (.fun in .twoDispatch) {
        .siglocal <- body(getMethod(.fun, .sig)@.Data)[[2]][[3]]
        .targlocal <- body(getMethod(.fun, .targ)@.Data)[[2]][[3]]
        checkIdentical(formals(.targlocal), formals(.siglocal))
    }

    ## FIXME: compare, Compare

    .sig <- "SummarizedExperiment"
    for (.fun in .otherFuns) {
        generic <- getGeneric(.fun)
        method <- getMethod(.fun, "SummarizedExperiment")
        checkIdentical(formals(generic@.Data), formals(method@.Data))
    }        
}

test_SummarizedExperiment_GRanges_values <- function()
{
    x <- ssetList[[1]]
    isAssign <- grep("<-$", .singleDispatch, value=TRUE)
    needArgs <- c("flank", "resize")
    isEndomorphism <- c("narrow", "restrict", "shift")
    .funs <- setdiff(.singleDispatch,
                     c(isAssign, needArgs, isEndomorphism))
    ## 'exp' created after manual inspection of results
    exp <- setNames(c("1f7c0", "35e2c", "02dde", "80339", "49a3f",
                      "72f53", "86757", "77198", "ec53a", "35e2c",
                      "625d9", "3c90a"), .funs)
    obs <- sapply(.funs, function(.fun) {
        substr(digest(getGeneric(.fun)(x)), 1, 5)
    })
    checkIdentical(exp, obs)

    .funs <- isAssign
    .gets <- sub("<-$", "", isAssign)
    for (i in seq_along(isAssign)) {
        ## self-assignment isomorphism
        value <- getGeneric(.gets[[i]])(x)
        x1 <- do.call(isAssign[[i]], list(x, value=value))
        checkIdentical(x, x1)
    }

    for (.fun in needArgs) {
        ## all needArgs operate on rowData
        generic <- getGeneric(.fun)
        x1 <- x; rowData(x1) <- generic(rowData(x1), 5)
        checkIdentical(x1, generic(x, 5))
    }
    ## isEndomorphism
    for (.fun in isEndomorphism) {
        generic <- getGeneric(.fun)
        obs <- generic(x)
        checkIdentical(generic(rowData(x)), rowData(obs))
        checkIdentical(assays(x), assays(obs))
    }

    .funs <- c(.twoDispatch[.twoDispatch != "subsetByOverlaps"],
               "findOverlaps", "countOverlaps")
    x1 <- shift(x, seq_len(nrow(x)) * 5)
    for (.fun in .funs) {
        generic <- getGeneric(.fun)
        exp <- generic(rowData(x1), rowData(x1))
        obs <- generic(x1, x1)
        checkIdentical(obs, exp)
    }
    # nearest,SummarizedExperiment,missing-method
    checkIdentical(nearest(rowData(x1)), nearest(x1)) 
    checkIdentical(subsetByOverlaps(rowData(x1), rowData(x1)[3]),
                   rowData(subsetByOverlaps(x1, x1[3])))
}

test_SummarizedExperiment_split <- function() {
    gr <- GRanges(Rle(c("A", "B"), c(2, 3)), IRanges(1:5, 10))
    se <- SummarizedExperiment(m, rowData=gr, colData=colData)
    ## FIXME: unname should not be necessary
    obs <- split(se, seqnames(se))
    exp <- SimpleList(A=se[1:2], B=se[3:5])
    checkEquals(obs, exp)
}

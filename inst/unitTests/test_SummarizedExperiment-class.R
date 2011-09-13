m <- matrix(1, 5, 3, dimnames=list(NULL, NULL))
mlst <- matrix(1, 3, 3, dimnames=list(NULL, NULL))
#assays <- SimpleList(m=m)
#assayslst <- SimpleList(m=mlst)
colData <- DataFrame(x=letters[1:3])
mList <- list(m, mlst)
assaysList <- list(gr=SimpleList(m=m), grl=SimpleList(m=mlst))
rowDataList <- 
    list(gr=GRanges("chr1", IRanges(1:5, 10)), 
         grl=split(GRanges("chr1", IRanges(1:5, 10)), c(1,1,2,2,3)))
names(rowDataList[["grl"]]) <- NULL
ssetList <- 
    list(SummarizedExperiment(assays=assaysList[["gr"]], 
         rowData=rowDataList[["gr"]], colData),
         SummarizedExperiment(assays=assaysList[["grl"]], 
         rowData=rowDataList[["grl"]], colData))

.SubsetGRListAtRangesLevel <- function(grl, idx, ...)
{
    if (is.character(idx)) {
        orig <- idx
        idx <- match(idx, names(grl@unlistData))
    }
    grReduced <-
        GenomicRanges:::deconstructGRLintoGR(grl)[idx,,drop=FALSE]
    grlReduced <-
        GenomicRanges:::reconstructGRLfromGR(grReduced, grl)
    grlReduced[elementLengths(grlReduced) != 0]
}


test_SummarizedExperiment_construction <- function() {
    ## empty-ish
    m1 <- matrix(0, 0, 0)
    checkTrue(validObject(new("SummarizedExperiment")))
    checkTrue(validObject(SummarizedExperiment()))
    checkTrue(validObject(
        SummarizedExperiment(assays=SimpleList(m1))))
    checkException(
        SummarizedExperiment(assays=SimpleList(matrix())),
        "assays dim mismatch", TRUE)
    checkException(
        SummarizedExperiment(assays=SimpleList(m1, matrix())),
        "assays dim mismatch", TRUE)
    checkException(
        SummarizedExperiment(assays=SimpleList(character())),
        "assays class", TRUE)

    ## substance
    for (i in length(ssetList)) {
        sset <- ssetList[[i]] 
        checkTrue(validObject(sset))
        checkIdentical(SimpleList(m=mList[[i]]), assays(sset))
        checkIdentical(rowDataList[[i]], rowData(sset))
        checkIdentical(DataFrame(x=letters[1:3]), colData(sset))
    }
}

test_SummarizedExperiment_getters <- function() {
    for (i in length(ssetList)) {
        sset <- ssetList[[i]] 
        rowData <- rowDataList[[i]] 
        ## dim, dimnames
        checkIdentical(c(length(rowData), nrow(colData)), dim(sset))
        checkIdentical(list(NULL, NULL), dimnames(sset))

        ## row / col / exptData
        checkIdentical(rowData, rowData(sset))
        checkIdentical(colData, colData(sset))
        checkIdentical(SimpleList(), exptData(sset))
    }

    ## assays
    m0 <- matrix(0L, 0, 0, dimnames=list(NULL, NULL))
    m1 <- matrix(0, 0, 0, dimnames=list(NULL, NULL))
    a <- SimpleList(a=m0, b=m1)
    checkIdentical(a, assays(SummarizedExperiment(assays=a)))
    ## assay
    checkException(
        assay(SummarizedExperiment()), "0-length assay", TRUE)
    checkIdentical(m0,
        assay(SummarizedExperiment(assays=a)), "default assay")
    checkIdentical(m1,
        assay(SummarizedExperiment(assays=a), 2),
        "assay, numeric index")
    checkException(
        assay(SummarizedExperiment(assays=a), 3),
        "invalid assay index", TRUE)
    checkIdentical(m1,
        assay(SummarizedExperiment(assays=a), "b"),
        "assay, character index")
    checkException(
        assay(SummarizedExperiment(assays=a), "c"),
        "invalid assay name", TRUE)
}

test_SummarizedExperiment_setters <- function()
{
    for (i in length(ssetList)) {
        sset <- ssetList[[i]] 
        rowData <- rowDataList[[i]] 
        ## row / col / exptData<-
        ss1 <- sset
        revData <- rowData[rev(seq_len(length(rowData))),,drop=FALSE]
        rowData(ss1) <- revData
        checkIdentical(revData, rowData(ss1))
        checkException(rowData(ss1) <- rowData(sset)[1:2,,drop=FALSE],
                       "incorrect row dimensions", TRUE)
        revData <- colData[rev(seq_len(nrow(colData))),,drop=FALSE]
        colData(ss1) <- revData
        checkIdentical(revData, colData(ss1))
        checkException(colData(ss1) <- colData(sset)[1:2,,drop=FALSE],
                       "incorrect col dimensions", TRUE)
        lst <- SimpleList("foo", "bar")
        exptData(ss1) <- lst
        checkIdentical(lst, exptData(ss1))

        ## assay / assays
        ss1 <- sset
        assay(ss1) <- assay(ss1)+1
        checkIdentical(assay(sset)+1, assay(ss1))
        ss1 <- sset
        assay(ss1, 1) <- assay(ss1, 1) + 1
        checkIdentical(assay(sset, "m") + 1, assay(ss1, "m"))
        ss1 <- sset
        assay(ss1, "m") <- assay(ss1, "m") + 1
        checkIdentical(assay(sset, "m")+1, assay(ss1, "m"))

        ## dimnames<-
        ss1 <- sset
        dimnames <- list(letters[seq_len(nrow(ss1))],
                         LETTERS[seq_len(ncol(ss1))])
        rownames(ss1) <- dimnames[[1]]
        colnames(ss1) <- dimnames[[2]]
        checkIdentical(dimnames, dimnames(ss1))
        rowData1 <- rowData
        names(rowData1) <- dimnames[[1]]
        checkIdentical(rowData1, rowData(ss1))
        colData1 <- colData
        row.names(colData1) <- dimnames[[2]]
        checkIdentical(colData1, colData(ss1))
        ss1 <- sset
        dimnames(ss1) <- dimnames
        checkIdentical(dimnames, dimnames(ss1))
        dimnames(ss1) <- NULL
        checkIdentical(list(NULL, NULL), dimnames(ss1))
    }
}

test_SummarizedExperiment_subset <- function()
{
    for (i in length(ssetList)) {
        sset <- ssetList[[i]] 
        rowData <- rowDataList[[i]] 
        ## numeric
        ss1 <- sset[2:3,]
        checkIdentical(c(2L, ncol(sset)), dim(ss1))
        checkIdentical(rowData(ss1), rowData(sset)[2:3,])
        checkIdentical(colData(ss1), colData(sset))
        ss1 <- sset[,2:3]
        checkIdentical(c(nrow(sset), 2L), dim(ss1))
        checkIdentical(rowData(ss1), rowData(sset))
        checkIdentical(colData(ss1), colData(sset)[2:3,,drop=FALSE])
        ss1 <- sset[2:3, 2:3]
        checkIdentical(c(2L, 2L), dim(ss1))
        checkIdentical(rowData(ss1), rowData(sset)[2:3,,drop=FALSE])
        checkIdentical(colData(ss1), colData(sset)[2:3,,drop=FALSE])

        ## character
        ss1 <- sset
        dimnames(ss1) <- list(LETTERS[seq_len(nrow(ss1))],
                               letters[seq_len(ncol(ss1))])
        ridx <- c("B", "C")
        checkIdentical(rowData(ss1[ridx,]), rowData(ss1)[ridx,])
        checkIdentical(rowData(ss1["C",]), rowData(ss1)["C",,drop=FALSE])
        checkException(ss1[LETTERS,], "i-index out of bounds", TRUE)
        cidx <- c("b", "c")
        checkIdentical(colData(ss1[,cidx]), colData(ss1)[cidx,,drop=FALSE])
        checkIdentical(colData(ss1[,"a"]), colData(ss1)["a",,drop=FALSE])
        checkException(ss1[,letters], "j-index out of bounds", TRUE)

        ## logical
        ss1 <- sset
        dimnames(ss1) <- list(LETTERS[seq_len(nrow(ss1))],
                               letters[seq_len(ncol(ss1))])
        checkIdentical(ss1, ss1[TRUE,])
        checkIdentical(c(0L, ncol(ss1)), dim(ss1[FALSE,]))
        checkIdentical(ss1, ss1[,TRUE])
        checkIdentical(c(nrow(ss1), 0L), dim(ss1[,FALSE]))
        idx <- c(TRUE, FALSE)               # recycling
        ss2 <- ss1[idx,]
        checkIdentical(rowData(ss1)[idx,,drop=FALSE], rowData(ss2))
        ss2 <- ss1[,idx]
        checkIdentical(colData(ss1)[idx,,drop=FALSE], colData(ss2))
    }
}

test_SummarizedExperiment_subsetassign <- function()
{
    for (i in length(ssetList)) {
        sset <- ssetList[[i]] 
        ss1 <- sset
        ss1[1:2,] <- ss1[2:1,]
        checkIdentical(rowData(sset)[2:1,], rowData(ss1)[1:2,])
        checkIdentical(rowData(sset[-(1:2),]), rowData(ss1)[-(1:2),])
        checkIdentical(colData(sset), colData(ss1))
        checkIdentical(c(exptData(sset), exptData(sset)), exptData(ss1))

        ss1 <- sset
        ss1[,1:2] <- ss1[,2:1,drop=FALSE]
        checkIdentical(colData(sset)[2:1,,drop=FALSE],
                       colData(ss1)[1:2,,drop=FALSE])
        checkIdentical(colData(sset)[-(1:2),,drop=FALSE],
                       colData(ss1)[-(1:2),,drop=FALSE])
        checkIdentical(rowData(sset), rowData(ss1))
        checkIdentical(c(exptData(sset), exptData(sset)), exptData(ss1))
    }
}

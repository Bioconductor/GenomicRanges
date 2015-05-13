m <- matrix(1, 5, 3, dimnames=list(NULL, NULL))
mlst <- matrix(1, 3, 3, dimnames=list(NULL, NULL))
mList <- list(m, mlst)
assaysList <- list(gr=SimpleList(m=m), grl=SimpleList(m=mlst))
rowRangesList <- 
    list(gr=GRanges("chr1", IRanges(1:5, 10)), 
         grl=split(GRanges("chr1", IRanges(1:5, 10)), c(1,1,2,2,3)))
names(rowRangesList[["grl"]]) <- NULL
colData <- DataFrame(x=letters[1:3])

## a list of one SE with GRanges and one with GRangesList
ssetList <- 
    list(SummarizedExperiment(
           assays=assaysList[["gr"]], 
           rowRanges=rowRangesList[["gr"]], 
           colData=colData),
         SummarizedExperiment(
           assays=assaysList[["grl"]], 
           rowRanges=rowRangesList[["grl"]], 
           colData=colData))


test_SummarizedExperiment_GenomicRanges_coercion <- function()
{
    if (requireNamespace("Biobase", quietly = TRUE)) {
        eset1 <- Biobase::ExpressionSet()

        checkTrue(validObject(eset1))

        se1 <- as(eset1, "SummarizedExperiment")

        checkTrue(validObject(se1))

        data("sample.ExpressionSet", package = "Biobase")

        eset2 <- sample.ExpressionSet
        checkTrue(validObject(eset2))

        se2 <- as(eset2, "SummarizedExperiment")

        checkTrue(validObject(se2))

        checkIdentical(Biobase::experimentData(eset2),
                       metadata(se2)$experimentData)

        checkIdentical(Biobase::annotation(eset2),
                       metadata(se2)$annotation)

        checkIdentical(Biobase::protocolData(eset2),
                       metadata(se2)$protocolData)

        eset2Assays <- SimpleList(as.list(Biobase::assayData(eset2)))
        se2Assays <- assays(se2)
        checkIdentical(eset2Assays$exprs, se2Assays$exprs)
        checkIdentical(eset2Assays$se.exprs, se2Assays$se.exprs)

        checkIdentical(Biobase::featureNames(eset2),
                       rownames(se2))

        checkIdentical(Biobase::sampleNames(eset2),
                       colnames(se2))
    }
}

test_GenomicRanges_SummarizedExperiment_coercion <- function()
{
    if (requireNamespace("Biobase", quietly = TRUE)) {
        assayData <- Biobase::assayData
        experimentData <- Biobase::experimentData
        annotation <- Biobase::annotation
        protocolData <- Biobase::protocolData
        featureNames <- Biobase::featureNames
        featureData <- Biobase::featureData
        sampleNames <- Biobase::sampleNames
        pData <- Biobase::pData

        ## empty SE
        simpleSE <- SummarizedExperiment()

        eset1 <- as(simpleSE, "ExpressionSet")

        checkTrue(validObject(eset1))

        ## Back and forth empty ES
        simpleES <- Biobase::ExpressionSet()

        simpleES2 <- as(as(simpleES, "SummarizedExperiment"), "ExpressionSet")

        checkTrue(validObject(simpleES2))

        checkEquals(as.list(assayData(simpleES)),
                    as.list(assayData(simpleES2)))

        ## Simple SE
        eset2 <- as(ssetList[[1]], "ExpressionSet")
        checkTrue(validObject(eset2))

        ## The ExpressionSet features should have the data from the
        ## SummarizedExperiment rows if they are from GRanges.
        checkIdentical(pData(featureData(eset2)),
                       as.data.frame(rowRanges(ssetList[[1]])))

        # the rowRanges are retained if the object has them to begin with.
        se2_2 <- as(eset2, "SummarizedExperiment")
        rr_se2_2 <- unname(rowRanges(se2_2))
        rr_eset2 <- rowRanges(ssetList[[1]])
        checkEquals(rr_se2_2, rr_eset2)

        eset3 <- as(ssetList[[2]], "ExpressionSet")
        checkTrue(validObject(eset3))

        ## The ExpressionSet features should not have the data from the
        ## SummarizedExperiment rows if they are from GRangesList, but they
        ## should be empty and the same length as the number of ranges.
        checkEquals(unname(NROW(featureData(eset3))),
                       unname(length(rowRanges(ssetList[[2]]))))

        data("sample.ExpressionSet", package = "Biobase")
        eset4 <- sample.ExpressionSet

        eset5 <- as(as(eset4, "SummarizedExperiment"), "ExpressionSet")

        checkTrue(validObject(eset5))

        ## this is necessary because the order in environments is undefined.
        compareLists <- function(x, y) {
            nmsX <- names(x)
            nmsY <- names(y)

            reorderY <- match(nmsY, nmsX)

            checkIdentical(x, y[reorderY])
        }

        compareLists(as.list(assayData(eset4)),
                       as.list(assayData(eset5)))

        checkIdentical(experimentData(eset4),
                       experimentData(eset5))

        checkIdentical(annotation(eset4),
                       annotation(eset5))

        checkIdentical(protocolData(eset4),
                       protocolData(eset5))

        checkIdentical(featureNames(eset4),
                       featureNames(eset5))

        checkIdentical(sampleNames(eset4),
                       sampleNames(eset5))
    }
}

test_GenomicRanges_SummarizedExperiment_coercion_mappingFunctions <- function()
{
    ExpressionSet <- Biobase::ExpressionSet

    ## naiveRangeMapper
    ## valid object from empty object
    checkTrue(validObject(makeSummarizedExperimentFromExpressionSet(ExpressionSet())))

    ## valid object from sample ExpressionSet
    data("sample.ExpressionSet", package = "Biobase")
    eset1 <- sample.ExpressionSet
    checkTrue(validObject(makeSummarizedExperimentFromExpressionSet(eset1)))

    ## makeSummarizedExperimentFromExpressionSet should be the same as `as`
    ## with default args
    checkEquals(makeSummarizedExperimentFromExpressionSet(eset1),
                as(eset1, "SummarizedExperiment"))

    ## probeRangeMapper
    ## valid object from empty object
    checkTrue(validObject(
            makeSummarizedExperimentFromExpressionSet(ExpressionSet(),
                probeRangeMapper)))

    ## valid object from sample ExpressionSet
    se1 <- makeSummarizedExperimentFromExpressionSet(eset1, probeRangeMapper)
    checkTrue(validObject(se1))

    ## Granges returned have rownames that were from the featureNames
    checkTrue(all(rownames(rowRanges(se1)) %in% Biobase::featureNames(eset1)))

    ## geneRangeMapper
    ## valid object from empty object
    checkTrue(validObject(
            makeSummarizedExperimentFromExpressionSet(ExpressionSet(),
                geneRangeMapper(NULL))))

    ## valid object from sample ExpressionSet
    se2 <- makeSummarizedExperimentFromExpressionSet(eset1,
        geneRangeMapper("TxDb.Hsapiens.UCSC.hg19.knownGene"))
    checkTrue(validObject(se2))

    ## Granges returned have rownames that were from the featureNames
    checkTrue(all(rownames(rowRanges(se2)) %in% Biobase::featureNames(eset1)))
}


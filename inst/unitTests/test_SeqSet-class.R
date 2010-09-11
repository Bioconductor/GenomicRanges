## FIXME: not needed when constructor, generics exported
assays <- GenomicRanges:::assays
assay <- GenomicRanges:::assay
SeqSet <- function(...) new("SeqSet", ...)

test_SeqSet_construction <- function() {
    m1 <- matrix(0, 0, 0)
    checkTrue(validObject(new("SeqSet")))
    checkTrue(validObject(SeqSet()))
    checkTrue(validObject(SeqSet(assays=SimpleList(m1))))
    checkException(SeqSet(assays=SimpleList(matrix())),
                   "assays dim mismatch", TRUE)
    checkException(SeqSet(assays=SimpleList(m1, matrix())),
                   "assays dim mismatch", TRUE)
    checkException(SeqSet(assays=SimpleList(character())),
                   "assays class", TRUE)
}

test_SeqSet_accessors <- function() {
    m0 <- matrix(0L, 0, 0)
    m1 <- matrix(0, 0, 0)
    a <- SimpleList(a=m0, b=m1)

    ## assays
    checkIdentical(a, assays(SeqSet(assays=a)))
    ## assay
    checkException(assay(SeqSet()), "0-length assay", TRUE)
    checkIdentical(m0, assay(SeqSet(assays=a)), "default assay")
    checkIdentical(m1, assay(SeqSet(assays=a), 2), "assay, numeric index")
    checkException(assay(SeqSet(assays=a), 3), "invalid assay index", TRUE)
    checkIdentical(m1, assay(SeqSet(assays=a), "b"), "assay, character index")
    checkException(assay(SeqSet(assays=a), "c"), "invalid assay name", TRUE)
}

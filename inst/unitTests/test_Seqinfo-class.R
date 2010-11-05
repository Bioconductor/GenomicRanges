###

test_Seqinfo.merge <- function()
{
    x <- Seqinfo(seqnames=c("chr1", "chr2", "chr3", "chrM"),
                 seqlengths=c(100, 200, NA, 15),
                 isCircular=c(NA, FALSE, FALSE, TRUE))
    checkIdentical(merge(x), x)
    checkIdentical(merge(x, NULL), x)
    checkIdentical(merge(NULL, x), x)
    checkIdentical(merge(x, Seqinfo()), x)
    checkIdentical(merge(Seqinfo(), x), x)
    checkIdentical(merge(x, x), x)
    checkIdentical(merge(x, Seqinfo(rev(names(seqlengths)))), x)

    y <- Seqinfo(seqnames=c("chrM", "chr4", "chr3"),
                 seqlengths=c(15, NA, 300))
    got <- merge(x, y)
    want <- Seqinfo(seqnames=c("chr1", "chr2", "chr3", "chrM", "chr4"),
                    seqlengths=c(100, 200, 300, 15, NA),
                    isCircular=c(NA, FALSE, FALSE, TRUE, NA))
    checkIdentical(got, want)

    ## This contradicts what 'x' says about circularity of chr3 and chrM:
    isCircular(y)[c("chr3", "chrM")] <- c(TRUE, FALSE)
    checkException(merge(x, y), silent=TRUE)
}


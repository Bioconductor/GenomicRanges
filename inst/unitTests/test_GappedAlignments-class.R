test_GappedAlignments_constructor <- function()
{
    checkTrue(validObject(GappedAlignments()))
    checkTrue(validObject(GappedAlignments(seqnames=factor("A"),
                                           pos=1L, cigar="1M",
                                           strand=strand("-"))))
}

test_GappedAlignments_concat <- function() {
    galn <- GappedAlignments(seqnames=factor("A"),
                             pos=1L, cigar="1M",
                             strand=strand("-"))
    galn_c <- GappedAlignments(seqnames=rep(factor("A"), 2),
                               pos=rep(1L, 2), cigar=rep("1M", 2),
                               strand=rep(strand("-"), 2))
    checkIdentical(galn_c, c(galn, galn))
}

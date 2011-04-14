test_GappedAlignments_constructor <- function()
{
    checkTrue(validObject(GappedAlignments()))
    checkTrue(validObject(GappedAlignments(rname=factor("A"),
                                           pos=1L, cigar="1M",
                                           strand=strand("-"))))
}

test_GappedAlignments_concat <- function() {
    galn <- GappedAlignments(rname=factor("A"),
                             pos=1L, cigar="1M",
                             strand=strand("-"))
    checkIdentical("GappedAlignments",
                   class(c(galn, galn)))
}

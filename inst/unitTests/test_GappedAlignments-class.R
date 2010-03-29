test_GappedAlignments_constructor <- function()
{
    checkTrue(validObject(GappedAlignments()))
    checkTrue(validObject(GappedAlignments(rname=factor("A"),
                                           pos=1L, cigar="1M",
                                           strand=strand("-"))))
}


test_GAlignments_constructor <- function()
{
    checkTrue(validObject(GAlignments()))
    checkTrue(validObject(GAlignments(seqnames=factor("A"),
                                      pos=1L, cigar="1M",
                                      strand=strand("-"))))
}

test_GAlignments_concat <- function() 
{
    galn <- GAlignments(seqnames=factor("A"),
                        pos=1L, cigar="1M",
                        strand=strand("-"))
    galn_c <- GAlignments(seqnames=rep(factor("A"), 2),
                          pos=rep(1L, 2), cigar=rep("1M", 2),
                          strand=rep(strand("-"), 2))
    checkIdentical(galn_c, c(galn, galn))
}


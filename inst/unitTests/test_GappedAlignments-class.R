test_Alignments0_constructor <- function()
{
    checkTrue(validObject(Alignments0()))
    checkTrue(validObject(Alignments0(rname=factor("A"),
                                      strand=strand("-"),
                                      pos=1L, cigar="1M")))
}

test_Alignments1_constructor <- function()
{
    checkTrue(validObject(Alignments1()))
    checkTrue(validObject(Alignments1(rname=factor("A"),
                                      strand=strand("-"),
                                      pos=1L, cigar="1M")))
}

test_Alignments2_constructor <- function()
{
    checkTrue(validObject(Alignments2()))
    checkTrue(validObject(Alignments2(rname=factor("A"),
                                      strand=strand("-"),
                                      pos=1L, cigar="1M")))
}


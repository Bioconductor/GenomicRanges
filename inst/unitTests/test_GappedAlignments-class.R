test_GappedAlignments_constructor <- function()
{
    checkTrue(validObject(GappedAlignments()))
    checkTrue(validObject(GappedAlignments(seqnames=factor("A"),
                                           pos=1L, cigar="1M",
                                           strand=strand("-"))))
}

test_GappedAlignments_concat <- function() 
{
    galn <- GappedAlignments(seqnames=factor("A"),
                             pos=1L, cigar="1M",
                             strand=strand("-"))
    galn_c <- GappedAlignments(seqnames=rep(factor("A"), 2),
                               pos=rep(1L, 2), cigar=rep("1M", 2),
                               strand=rep(strand("-"), 2))
    checkIdentical(galn_c, c(galn, galn))
}

test_GappedAlignments_qnarrow <- function() 
{
    gal <- GappedAlignments(seqnames=rep(factor("A"), 8),
                            pos=10:17,
                            cigar=c("5M", "5X", "3M2I3M", "3M2D3M",
                                    "3M2N3M", "3M2S3M", "3M2H3M", "3M2P3M"),
                            strand=Rle(strand(rep("+", 8))))
    n1 <- narrow(gal, start=3)
    q1 <- qnarrow(gal, start=3)
    checkIdentical(qwidth(n1), qwidth(q1)) 
    checkIdentical(width(n1), width(q1))
 
    n2 <- narrow(gal, start=4)
    q2 <- qnarrow(gal, start=4)
    checkIdentical(width(n2), width(q2))
    ## M and X 
    checkIdentical(qwidth(n2[1:2]), qwidth(q2[1:2]))
    ## I 
    checkIdentical(qwidth(q2[3]), width(q2[3]) + 2L)
    ## D, N and P
    checkIdentical(qwidth(q2[c(4,5,8)]), width(q2[c(4,5,8)]))
    ## S and H
    checkIdentical(qwidth(q2[6]), width(q2[6]) + 2L)
    checkIdentical(qwidth(q2[7]), width(q2[7]))
}

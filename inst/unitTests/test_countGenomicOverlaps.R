rng <- function(s, w)
GRanges(seq="chr1", IRanges(s, width=w), strand="+")

rng2 <- function(s, w)
GRanges(seq="chr2", IRanges(s, width=w), strand="+")

make_subject <- function() 
{
    GRangesList(A=rng(1000, 500),
        B=rng2(2000, 900),
        C=rng(c(3000, 3600), c(500, 300)),
        D=rng2(c(7000, 7500), c(600, 300)),
        E1=rng(4000, 500), E2=rng(c(4300, 4500), c(400, 400)),
        F=rng2(3000, 500),
        G=rng(c(5000, 5600), c(500, 300)),
        H1=rng(6000, 500), H2=rng(6600, 400))
}

make_GRquery <- function() 
{
    GRangesList(a=rng(1400, 500),
        b=rng2(2700, 100),
        c=rng(3400, 300),
        d=rng2(7100, 600),
        e=rng(4200, 500),
        f=rng2(c(3100, 3300), 50),
        g=rng(c(5400, 5600), 50),
        h=rng(c(6400, 6600), 50))
}

make_GAquery <- function() 
{
    GappedAlignments(
        rname=(c("chr1", "chr2", "chr1", "chr2", "chr1", "chr2", 
               "chr1", "chr1")),
        pos=as.integer(c(1400, 2700, 3400, 7100, 4200, 3100, 5400, 6400)),
        cigar=c("500M", "100M", "300M", "600M", "500M", "50M150N50M", 
                "50M150N50M", "50M150N50M"),
        strand=strand(rep.int("+", 8L))) 
}

.getCounts <- function(res)
{
    as.vector(values(unlist(res))[["hits"]])
}

test_input_type <- function()
{
    subject <- make_subject()
    GRquery <- make_GRquery()
    GAquery <- make_GAquery() 
    GRq_GRLs <- countGenomicOverlaps(GRquery, subject)
    GAq_GRLs <- countGenomicOverlaps(GAquery, subject)
    GRq_GRs <- countGenomicOverlaps(GAquery, unlist(subject))
    GAq_GRs <- countGenomicOverlaps(GAquery, unlist(subject))

    checkIdentical(GRq_GRLs, GAq_GRLs) 
    checkIdentical(GRq_GRs, GAq_GRs) 
}

test_typeAny <- function(type="any", ...)
{
    .tests <- function(ignore.strand=FALSE)
    {
        subject <- make_subject()
        query <- make_GRquery()

        anyNone <- countGenomicOverlaps(query, subject, type, resolution="none")
        anyDivide <- countGenomicOverlaps(query, subject, type, resolution="divide")
        anyUD <- countGenomicOverlaps(query, subject, type, resolution="uniqueDisjoint")
 
        checkIdentical(.getCounts(anyNone), c(1, 1, rep.int(0L, 12)))
        checkIdentical(round(.getCounts(anyDivide), 3), c(1, 1, rep.int(0.5, 4), 
                       rep.int(0.333, 3), 1, rep.int(0.5, 4)))
        checkIdentical(round(.getCounts(anyUD), 3), c(1, 1, rep.int(0, 4), 1, 
                       0, 0, 1, rep.int(0.5, 4)))
    }
    .tests(ignore.strand=TRUE)
    .tests(ignore.strand=FALSE)
} 


test_typeWithin <- function(type="within", ...)
{
    .tests <- function(ignore.strand=FALSE)
    {
        subject <- make_subject()
        query <- make_GRquery()

        withinNone <- countGenomicOverlaps(query, subject, type=type, resolution="none")
        withinDivide <- countGenomicOverlaps(query, subject, type=type, resolution="divide")
 
        checkIdentical(.getCounts(withinNone), c(0, 1, rep.int(0L, 12)))
        checkIdentical(.getCounts(withinDivide), 
                       c(0, 1, rep.int(0L, 7), 1, 0.5, 0.5, 0.5, 0.5))
    }
    .tests(ignore.strand=TRUE)
    .tests(ignore.strand=FALSE)
} 





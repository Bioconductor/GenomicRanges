test_GenomicData_construction <- function() {
  range1 <- IRanges(start=c(1,2,3), end=c(5,2,8))

  gr <- GenomicData(range1)
  checkTrue(validObject(gr))
  checkIdentical(gr, RangedData(range1))
  gr <- GenomicData(range1, genome = "hg18")
  checkTrue(validObject(gr))
  checkIdentical(genome(gr), "hg18")

  filter <- c(1L, 0L, 1L)
  score <- c(10L, 2L, NA)
  strand <- factor(c("+", NA, "-"), levels = levels(strand()))
  
  gr <- GenomicData(range1, score, genome = "hg18")
  checkTrue(validObject(gr))
  checkIdentical(gr[["score"]], score)
  checkIdentical(strand(gr), factor(rep(NA, 3), levels = levels(strand())))
  gr <- GenomicData(range1, score, filt = filter, strand = strand)
  checkTrue(validObject(gr))
  checkIdentical(gr[["score"]], score)
  checkIdentical(gr[["filt"]], filter)
  checkIdentical(strand(gr), strand)

  range2 <- IRanges(start=c(15,45,20,1), end=c(15,100,80,5))
  ranges <- c(range1, range2)
  score <- c(score, c(0L, 3L, NA, 22L)) 
  chrom <- paste("chr", rep(c(1,2), c(length(range1), length(range2))), sep="")
  
  gr <- GenomicData(ranges, score, genome = "hg18")
  checkTrue(validObject(gr))
  gr <- GenomicData(ranges, score, chrom = chrom, genome = "hg18")
  checkTrue(validObject(gr))
  checkIdentical(gr[["score"]], score)
  checkIdentical(score(gr), score)
  checkIdentical(chrom(gr), factor(chrom))
  checkIdentical(gr[1][["score"]], score[1:3])
  gr <- GenomicData(ranges, foo = score, chrom = chrom, genome = "hg18")
  checkIdentical(score(gr), score) ## if no score, gets first col
  
  checkException(GenomicData(), silent = TRUE)
  checkException(GenomicData(range1, genome = c("hg18", "mm9")), silent = TRUE)
  checkException(GenomicData(range1, genome = 1), silent = TRUE)
  checkException(GenomicData(ranges, chrom = chrom[1:3]), silent = TRUE)
  checkException(GenomicData(range1, strand = strand[1:2]), silent = TRUE)
  checkException(GenomicData(range1, strand = c("+", "-", ".")), silent = TRUE)
}

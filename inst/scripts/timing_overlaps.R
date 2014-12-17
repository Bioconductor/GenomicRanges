### =========================================================================
### Timing findOverlaps() and summarizeOverlaps()
### -------------------------------------------------------------------------

library(RNAseqData.HNRNPC.bam.chr14)
bamfile <- RNAseqData.HNRNPC.bam.chr14_BAMFILES[1]
library(GenomicAlignments)
reads <- readGAlignments(bamfile)
reads_grl <- grglist(reads)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exbygene <- exonsBy(txdb, by="gene")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Timing findOverlaps,GRangesList,GRangesList method
###

system.time(hits1 <- findOverlaps(exbygene, reads_grl))
#digest(hits1)  # "42f770fcf1a2f27f2ceeb18249244677"

# With BioC 3.0 / R 3.1
# ---------------------
# First time:
#   user  system elapsed 
#  1.999   0.020   2.022 
# Subsequent times:
#   user  system elapsed 
#  0.852   0.000   0.854 

# With BioC 3.1/ R 3.2
# --------------------
# First time:
#   user  system elapsed 
#  0.797   0.008   0.806 
# Subsequent times:
#   user  system elapsed 
#  0.141   0.000   0.141 

system.time(hits2 <- findOverlaps(exbygene, reads_grl, ignore.strand=TRUE))
#digest(hits2)  # "3f1a62b338d431ef2705602ce6dfbf9a"

# With BioC 3.0 / R 3.1
# ---------------------
# First time:
#   user  system elapsed 
#  3.605   0.020   3.629 
# Subsequent times:
#   user  system elapsed 
#  1.152   0.000   1.154 

# With BioC 3.1/ R 3.2
# --------------------
# First time:
#   user  system elapsed 
#  1.244   0.008   1.254 
# Subsequent times:
#   user  system elapsed 
#  0.170   0.000   0.171 

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Timing findOverlaps,GRangesList,GAlignments method
###

system.time(hits3 <- findOverlaps(exbygene, reads))
#digest(hits3)  # "42f770fcf1a2f27f2ceeb18249244677"

# With BioC 3.0 / R 3.1
# ---------------------
# First time:
#   user  system elapsed 
#  3.844   0.016   3.864
# Subsequent times:
#   user  system elapsed 
#  1.093   0.000   1.095 

# With BioC 3.1/ R 3.2
# --------------------
# First time:
#   user  system elapsed 
#  1.792   0.008   1.802 
# Subsequent times:
#   user  system elapsed 
#  0.366   0.000   0.367 

system.time(hits4 <- findOverlaps(exbygene, reads, ignore.strand=TRUE))
#digest(hits4)  # "3f1a62b338d431ef2705602ce6dfbf9a"

# With BioC 3.0 / R 3.1
# ---------------------
# First time:
#   user  system elapsed 
#  5.075   0.004   5.084 
# Subsequent times:
#   user  system elapsed 
#  2.642   0.012   2.658 

# With BioC 3.1/ R 3.2
# --------------------
# First time:
#   user  system elapsed 
#  1.766   0.004   1.771 
# Subsequent times:
#   user  system elapsed 
#  0.396   0.000   0.397 

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Timing summarizeOverlaps,GRangesList,GAlignments method
###

system.time(se1 <- summarizeOverlaps(exbygene, reads))
#digest(assay(se1))  # "7e58dc4eb0c56fbf75ff86c09dbafb90"

# With BioC 3.0 / R 3.1
# ---------------------
# First time:
#   user  system elapsed 
#  3.815   0.012   3.831 
# Subsequent times:
#   user  system elapsed 
#  1.361   0.000   1.363 

# With BioC 3.1/ R 3.2
# --------------------
# First time:
#   user  system elapsed 
#  1.783   0.004   1.789 
# Subsequent times:
#   user  system elapsed 
#  0.436   0.004   0.441 

system.time(se2 <- summarizeOverlaps(exbygene, reads, ignore.strand=TRUE))
#digest(assay(se2))  # "7d473582932ee5704acbddd1ffc5c146"

# With BioC 3.0 / R 3.1
# ---------------------
# First time:
#   user  system elapsed 
#  5.545   0.020   5.570 
# Subsequent times:
#   user  system elapsed 
#  2.614   0.000   2.617 

# With BioC 3.1/ R 3.2
# --------------------
# First time:
#   user  system elapsed 
#  2.489   0.004   2.495 
# Subsequent times:
#   user  system elapsed 
#  0.490   0.008   0.498 


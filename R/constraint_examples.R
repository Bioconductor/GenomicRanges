### =========================================================================
### Some examples of concrete Constraint subclasses
### -------------------------------------------------------------------------
###
### Like the Constraint virtual class itself, concrete Constraint subclasses
### cannot have slots.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### EXAMPLE 1: The HasRangeTypeCol constraint.
###
### The HasRangeTypeCol constraint checks that the object it's attached to
### has a unique "rangeType" column in its elementMetadata part and that this
### column is a 'factor' Rle with no NAs and with the following levels (in
### this order): gene, transcript, exon, cds, 5utr, 3utr.

setClass("HasRangeTypeCol", contains="Constraint")

### Like validity methods, "checkConstraint" methods must return NULL or
### a character vector describing the problems found. They should never fail
### i.e. they should never raise an error.
setMethod("checkConstraint", c("GenomicRanges", "HasRangeTypeCol"),
    function(x, constraint, verbose=FALSE)
    {
        emd <- elementMetadata(x)
        idx <- match("rangeType", colnames(emd))
        if (length(idx) != 1L || is.na(idx)) {
            msg <- c("'elementMetadata(x)' must have exactly 1 column ",
                     "named \"rangeType\"")
            return(paste(msg, collapse=""))
        }
        rangeType <- emd[[idx]]
        .LEVELS <- c("gene", "transcript", "exon", "cds", "5utr", "3utr")
        if (!is(rangeType, "Rle") ||
            IRanges:::anyMissing(runValue(rangeType)) ||
            !identical(levels(rangeType), .LEVELS))
        {
            msg <- c("'elementMetadata(x)$rangeType' must be a 'factor' Rle ",
                     "with no NAs and with levels: ",
                     paste(.LEVELS, collapse=", "))
            return(paste(msg, collapse=""))
        }
        NULL
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### EXAMPLE 2: The GeneRanges constraint.
###
### The GeneRanges constraint is defined on top of the HasRangeTypeCol
### constraint. It checks that all the ranges in the object are of type
### "gene".

setClass("GeneRanges", contains="HasRangeTypeCol")

### The checkConstraint() generic will check the HasRangeTypeCol constraint
### first, and, if it's statisfied, it will then check the GeneRanges
### constraint.
setMethod("checkConstraint", c("GenomicRanges", "GeneRanges"),
    function(x, constraint, verbose=FALSE)
    {
        rangeType <- elementMetadata(x)$rangeType
        if (!all(rangeType == "gene")) {
            msg <- c("all elements in 'elementMetadata(x)$rangeType' must be ",
                     "equal to \"gene\"")
            return(paste(msg, collapse=""))
        }
        NULL
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### EXAMPLE 3: The HasGCCol constraint.
###
### The HasGCCol constraint checks that the object it's attached to
### has a unique "GC" column in its elementMetadata part, that this column
### is of type numeric, with no NAs, and that all the values in that column
### are >= 0 and <= 1.

setClass("HasGCCol", contains="Constraint")

setMethod("checkConstraint", c("GenomicRanges", "HasGCCol"),
    function(x, constraint, verbose=FALSE)
    {
        emd <- elementMetadata(x)
        idx <- match("GC", colnames(emd))
        if (length(idx) != 1L || is.na(idx)) {
            msg <- c("'elementMetadata(x)' must have exactly one column ",
                     "named \"GC\"")
            return(paste(msg, collapse=""))
        }
        GC <- emd[[idx]]
        if (!is.numeric(GC) ||
            IRanges:::anyMissing(GC) ||
            any(GC < 0) || any(GC > 1))
        {
            msg <- c("'elementMetadata(x)$GC' must be a numeric vector ",
                     "with no NAs and with values between 0 and 1")
            return(paste(msg, collapse=""))
        }
        NULL
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### EXAMPLE 4: The HighGCRanges constraint.
###
### The HighGCRanges constraint is defined on top of the HasGCCol
### constraint. It checks that all the ranges in the object have a GC content
### >= 0.5.

setClass("HighGCRanges", contains="HasGCCol")

### The checkConstraint() generic will check the HasGCCol constraint
### first, and, if it's statisfied, it will then check the HighGCRanges
### constraint.
setMethod("checkConstraint", c("GenomicRanges", "HighGCRanges"),
    function(x, constraint, verbose=FALSE)
    {
        GC <- elementMetadata(x)$GC
        if (!all(GC >= 0.5)) {
            msg <- "all elements in 'elementMetadata(x)$GC' must be >= 0.5"
            return(msg)
        }
        NULL
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### EXAMPLE 5: The HighGCGeneRanges constraint.
###
### The HighGCGeneRanges constraint is the combination (AND) of the GeneRanges
### and HighGCRanges constraints.

setClass("HighGCGeneRanges", contains=c("GeneRanges", "HighGCRanges"))

### No need to define a method for this constraint: the checkConstraint()
### generic will automatically check the GeneRanges and HighGCRanges
### constraints.


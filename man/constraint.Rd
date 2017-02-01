\name{Constraints}
\alias{Constraints}

\alias{class:Constraint}
\alias{Constraint-class}
\alias{Constraint}
\alias{class:Constraint_OR_NULL}
\alias{Constraint_OR_NULL-class}
\alias{Constraint_OR_NULL}

\alias{constraint}
\alias{constraint<-}
\alias{checkConstraint}


\title{Enforcing constraints thru Constraint objects}

\description{
  Attaching a Constraint object to an object of class A (the "constrained"
  object) is meant to be a convenient/reusable/extensible way to enforce
  a particular set of constraints on particular instances of A.

  THIS IS AN EXPERIMENTAL FEATURE AND STILL VERY MUCH A WORK-IN-PROGRESS!
}

\details{
  For the developer, using constraints is an alternative to the more
  traditional approach that consists in creating subclasses of A and
  implementing specific validity methods for each of them.
  However, using constraints offers the following advantages over the
  traditional approach:
  \itemize{
    \item The traditional approach often tends to lead to a proliferation
          of subclasses of A.
    \item Constraints can easily be re-used across different classes without
          the need to create any new class.
    \item Constraints can easily be combined.
  }

  All constraints are implemented as concrete subclasses of the Constraint
  class, which is a virtual class with no slots.
  Like the Constraint virtual class itself, concrete Constraint subclasses
  cannot have slots.

  Here are the 7 steps typically involved in the process of putting
  constraints on objects of class A:
  \enumerate{
    \item Add a slot named \code{constraint} to the definition of class A.
          The type of this slot must be Constraint_OR_NULL. Note that
          any subclass of A will inherit this slot.
    \item Implements the \code{constraint()} accessors (getter and setter)
          for objects of class A. This is done by implementing the
          \code{"constraint"} method (getter) and replacement method (setter)
          for objects of class A (see the examples below). As a convenience
          to the user, the setter should also accept the name of a constraint
          (i.e. the name of its class) in addition to an instance of that
          class. Note that those accessors will work on instances of any
          subclass of A.
    \item Modify the validity method for class A so it also returns the
          result of \code{checkConstraint(x, constraint(x))} (append this
          result to the result returned by the validity method).
    \item Testing: Create \code{x}, an instance of class A (or subclass of A).
          By default there is no constraint on it (\code{constraint(x)} is
          \code{NULL}). \code{validObject(x)} should return \code{TRUE}.
    \item Create a new constraint (MyConstraint) by extending the
          Constraint class, typically with
          \code{setClass("MyConstraint", contains="Constraint")}.
          This constraint is not enforcing anything yet so you could put
          it on \code{x} (with \code{constraint(x) <- "MyConstraint"}),
          but not much would happen. In order to actually enforce something,
          a \code{"checkConstraint"} method for signature
          \code{c(x="A", constraint="MyConstraint")} needs to be implemented.
    \item Implement a \code{"checkConstraint"} method for signature
          \code{c(x="A", constraint="MyConstraint")}. Like validity methods,
          \code{"checkConstraint"} methods must return \code{NULL} or
          a character vector describing the problems found.
          Like validity methods, they should never fail (i.e. they should
          never raise an error).
          Note that, alternatively, an existing constraint (e.g.
          SomeConstraint) can be adapted to work on objects of class A
          by just defining a new \code{"checkConstraint"} method for
          signature \code{c(x="A", constraint="SomeConstraint")}.
          Also, stricter constraints can be built on top of existing
          constraints by extending one or more existing constraints
          (see the examples below).
    \item Testing: Try \code{constraint(x) <- "MyConstraint"}. It will or
          will not work depending on whether \code{x} satisfies the
          constraint or not. In the former case, trying to modify \code{x}
          in a way that breaks the constraint on it will also raise an
          error.
  }
}

\note{
  WARNING: This note is not true anymore as the \code{constraint} slot has
  been temporarily removed from \link{GenomicRanges} objects (starting with
  package GenomicRanges >= 1.7.9).

  Currently, only \link{GenomicRanges} objects can be constrained, that is:
  \itemize{
    \item they have a \code{constraint} slot;
    \item they have \code{constraint()} accessors (getter and setter)
          for this slot;
    \item their validity method has been modified so it also returns the
          result of \code{checkConstraint(x, constraint(x))}.
  }
  More classes in the GenomicRanges and IRanges packages will support
  constraints in the near future.
}

\author{H. Pagès}

\seealso{
  \code{\link{setClass}},
  \code{\link{is}},
  \code{\link{setMethod}},
  \code{\link{showMethods}},
  \code{\link{validObject}},
  \link{GenomicRanges-class}
}

\examples{
## The examples below show how to define and set constraints on
## GenomicRanges objects. Note that this is how the constraint()
## setter is defined for GenomicRanges objects:
#setReplaceMethod("constraint", "GenomicRanges",
#    function(x, value)
#    {
#        if (isSingleString(value))
#            value <- new(value)
#        if (!is(value, "Constraint_OR_NULL"))
#            stop("the supplied 'constraint' must be a ",
#                 "Constraint object, a single string, or NULL")
#        x@constraint <- value
#        validObject(x)
#        x
#    }
#)

#selectMethod("constraint", "GenomicRanges")  # the getter
#selectMethod("constraint<-", "GenomicRanges")  # the setter

## We'll use the GRanges instance 'gr' created in the GRanges examples
## to test our constraints:
example(GRanges, echo=FALSE)
gr
#constraint(gr)

## ---------------------------------------------------------------------
## EXAMPLE 1: The HasRangeTypeCol constraint.
## ---------------------------------------------------------------------
## The HasRangeTypeCol constraint checks that the constrained object
## has a unique "rangeType" metadata column and that this column
## is a 'factor' Rle with no NAs and with the following levels
## (in this order): gene, transcript, exon, cds, 5utr, 3utr.

setClass("HasRangeTypeCol", contains="Constraint")

## Like validity methods, "checkConstraint" methods must return NULL or
## a character vector describing the problems found. They should never
## fail i.e. they should never raise an error.
setMethod("checkConstraint", c("GenomicRanges", "HasRangeTypeCol"),
    function(x, constraint, verbose=FALSE)
    {
        x_mcols <- mcols(x)
        idx <- match("rangeType", colnames(x_mcols))
        if (length(idx) != 1L || is.na(idx)) {
            msg <- c("'mcols(x)' must have exactly 1 column ",
                     "named \"rangeType\"")
            return(paste(msg, collapse=""))
        }
        rangeType <- x_mcols[[idx]]
        .LEVELS <- c("gene", "transcript", "exon", "cds", "5utr", "3utr")
        if (!is(rangeType, "Rle") ||
            S4Vectors:::anyMissing(runValue(rangeType)) ||
            !identical(levels(rangeType), .LEVELS))
        {
            msg <- c("'mcols(x)$rangeType' must be a ",
                     "'factor' Rle with no NAs and with levels: ",
                     paste(.LEVELS, collapse=", "))
            return(paste(msg, collapse=""))
        }
        NULL
    }
)

#\dontrun{
#constraint(gr) <- "HasRangeTypeCol"  # will fail
#}
checkConstraint(gr, new("HasRangeTypeCol"))  # with GenomicRanges >= 1.7.9

levels <- c("gene", "transcript", "exon", "cds", "5utr", "3utr")
rangeType <- Rle(factor(c("cds", "gene"), levels=levels), c(8, 2))
mcols(gr)$rangeType <- rangeType
#constraint(gr) <- "HasRangeTypeCol"  # OK
checkConstraint(gr, new("HasRangeTypeCol"))  # with GenomicRanges >= 1.7.9

## Use is() to check whether the object has a given constraint or not:
#is(constraint(gr), "HasRangeTypeCol")  # TRUE
#\dontrun{
#mcols(gr)$rangeType[3] <- NA  # will fail
#}
mcols(gr)$rangeType[3] <- NA
checkConstraint(gr, new("HasRangeTypeCol"))  # with GenomicRanges >= 1.7.9

## ---------------------------------------------------------------------
## EXAMPLE 2: The GeneRanges constraint.
## ---------------------------------------------------------------------
## The GeneRanges constraint is defined on top of the HasRangeTypeCol
## constraint. It checks that all the ranges in the object are of type
## "gene".

setClass("GeneRanges", contains="HasRangeTypeCol")

## The checkConstraint() generic will check the HasRangeTypeCol constraint
## first, and, only if it's statisfied, it will then check the GeneRanges
## constraint.
setMethod("checkConstraint", c("GenomicRanges", "GeneRanges"),
    function(x, constraint, verbose=FALSE)
    {
        rangeType <- mcols(x)$rangeType
        if (!all(rangeType == "gene")) {
            msg <- c("all elements in 'mcols(x)$rangeType' ",
                     "must be equal to \"gene\"")
            return(paste(msg, collapse=""))
        }
        NULL
    }
)

#\dontrun{
#constraint(gr) <- "GeneRanges"  # will fail
#}
checkConstraint(gr, new("GeneRanges"))  # with GenomicRanges >= 1.7.9

mcols(gr)$rangeType[] <- "gene"
## This replace the previous constraint (HasRangeTypeCol):
#constraint(gr) <- "GeneRanges"  # OK
checkConstraint(gr, new("GeneRanges"))  # with GenomicRanges >= 1.7.9

#is(constraint(gr), "GeneRanges")  # TRUE
## However, 'gr' still indirectly has the HasRangeTypeCol constraint
## (because the GeneRanges constraint extends the HasRangeTypeCol
## constraint):
#is(constraint(gr), "HasRangeTypeCol")  # TRUE
#\dontrun{
#mcols(gr)$rangeType[] <- "exon"  # will fail
#}
mcols(gr)$rangeType[] <- "exon"
checkConstraint(gr, new("GeneRanges"))  # with GenomicRanges >= 1.7.9

## ---------------------------------------------------------------------
## EXAMPLE 3: The HasGCCol constraint.
## ---------------------------------------------------------------------
## The HasGCCol constraint checks that the constrained object has a
## unique "GC" metadata column, that this column is of type numeric,
## with no NAs, and that all the values in that column are >= 0 and <= 1.

setClass("HasGCCol", contains="Constraint")

setMethod("checkConstraint", c("GenomicRanges", "HasGCCol"),
    function(x, constraint, verbose=FALSE)
    {
        x_mcols <- mcols(x)
        idx <- match("GC", colnames(x_mcols))
        if (length(idx) != 1L || is.na(idx)) {
            msg <- c("'mcols(x)' must have exactly ",
                     "one column named \"GC\"")
            return(paste(msg, collapse=""))
        }
        GC <- x_mcols[[idx]]
        if (!is.numeric(GC) ||
            S4Vectors:::anyMissing(GC) ||
            any(GC < 0) || any(GC > 1))
        {
            msg <- c("'mcols(x)$GC' must be a numeric vector ",
                     "with no NAs and with values between 0 and 1")
            return(paste(msg, collapse=""))
        }
        NULL
    }
)

## This replace the previous constraint (GeneRanges):
#constraint(gr) <- "HasGCCol"  # OK
checkConstraint(gr, new("HasGCCol"))  # with GenomicRanges >= 1.7.9

#is(constraint(gr), "HasGCCol")  # TRUE
#is(constraint(gr), "GeneRanges")  # FALSE
#is(constraint(gr), "HasRangeTypeCol")  # FALSE

## ---------------------------------------------------------------------
## EXAMPLE 4: The HighGCRanges constraint.
## ---------------------------------------------------------------------
## The HighGCRanges constraint is defined on top of the HasGCCol
## constraint. It checks that all the ranges in the object have a GC
## content >= 0.5.

setClass("HighGCRanges", contains="HasGCCol")

## The checkConstraint() generic will check the HasGCCol constraint
## first, and, if it's statisfied, it will then check the HighGCRanges
## constraint.
setMethod("checkConstraint", c("GenomicRanges", "HighGCRanges"),
    function(x, constraint, verbose=FALSE)
    {
        GC <- mcols(x)$GC
        if (!all(GC >= 0.5)) {
            msg <- c("all elements in 'mcols(x)$GC' ",
                     "must be >= 0.5")
            return(paste(msg, collapse=""))
        }
        NULL
    }
)

#\dontrun{
#constraint(gr) <- "HighGCRanges"  # will fail
#}
checkConstraint(gr, new("HighGCRanges"))  # with GenomicRanges >= 1.7.9
mcols(gr)$GC[6:10] <- 0.5
#constraint(gr) <- "HighGCRanges"  # OK
checkConstraint(gr, new("HighGCRanges"))  # with GenomicRanges >= 1.7.9

#is(constraint(gr), "HighGCRanges")  # TRUE
#is(constraint(gr), "HasGCCol")  # TRUE

## ---------------------------------------------------------------------
## EXAMPLE 5: The HighGCGeneRanges constraint.
## ---------------------------------------------------------------------
## The HighGCGeneRanges constraint is the combination (AND) of the
## GeneRanges and HighGCRanges constraints.

setClass("HighGCGeneRanges", contains=c("GeneRanges", "HighGCRanges"))

## No need to define a method for this constraint: the checkConstraint()
## generic will automatically check the GeneRanges and HighGCRanges
## constraints.

#constraint(gr) <- "HighGCGeneRanges"  # OK
checkConstraint(gr, new("HighGCGeneRanges"))  # with GenomicRanges >= 1.7.9

#is(constraint(gr), "HighGCGeneRanges")  # TRUE
#is(constraint(gr), "HighGCRanges")  # TRUE
#is(constraint(gr), "HasGCCol")  # TRUE
#is(constraint(gr), "GeneRanges")  # TRUE
#is(constraint(gr), "HasRangeTypeCol")  # TRUE

## See how all the individual constraints are checked (from less
## specific to more specific constraints):
#checkConstraint(gr, constraint(gr), verbose=TRUE)
checkConstraint(gr, new("HighGCGeneRanges"), verbose=TRUE)  # with
                                                            # GenomicRanges
                                                            # >= 1.7.9

## See all the "checkConstraint" methods:
showMethods("checkConstraint")
}

\keyword{methods}
\keyword{classes}

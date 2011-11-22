\name{Constraints}
\alias{Constraints}

\alias{class:Constraint}
\alias{Constraint-class}
\alias{Constraint}
\alias{class:ConstraintORNULL}
\alias{ConstraintORNULL-class}
\alias{ConstraintORNULL}

\alias{constraint}
\alias{constraint<-}
\alias{checkConstraint}


\title{Enforcing constraints thru Constraint objects}

\description{
  Attaching a Constraint object to an object of class A (the "constrained"
  object) is meant to be a convenient/reusable/extensible way to enforce
  a particular set of constraints on particular instances of A.

  THIS IS A WORK-IN-PROGRESS!
}

\details{
  For the developper, using constraints is an alternative to the more
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
  virtual class, which is a virtual class with no slots.
  Like the Constraint virtual class itself, concrete Constraint subclasses
  cannot have slots.

  Here are the 7 steps typically involved in the process of putting
  constraints on objects of class A (note that the first 4 steps are
  already done for \link{GenomicRanges} objects):
  \enumerate{
    \item Add a slot named \code{constraint} to the definition of class A.
          The type of this slot must be ConstraintORNULL. Note that
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

\author{H. Pages}

\seealso{
  \code{\link{validObject}},
  \link{GenomicRanges-class}
}

\examples{
## The examples below show how to define and set constraints on
## GenomicRanges objects.
## We'll use the GRanges instance 'gr' created in the GRanges examples
## to test our constraints:
example(GRanges, echo=FALSE)
gr
constraint(gr)

## ---------------------------------------------------------------------
## EXAMPLE 1: The HasRangeTypeCol constraint.
## ---------------------------------------------------------------------
## The HasRangeTypeCol constraint checks that the constrained object
## has a unique "rangeType" column in its elementMetadata part and that
## this column is a 'factor' Rle with no NAs and with the following
## levels (in this order): gene, transcript, exon, cds, 5utr, 3utr.

setClass("HasRangeTypeCol", contains="Constraint")

## Like validity methods, "checkConstraint" methods must return NULL or
## a character vector describing the problems found. They should never
## fail i.e. they should never raise an error.
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
            msg <- c("'elementMetadata(x)$rangeType' must be a ",
                     "'factor' Rle with no NAs and with levels: ",
                     paste(.LEVELS, collapse=", "))
            return(paste(msg, collapse=""))
        }
        NULL
    }
)

\dontrun{
constraint(gr) <- "HasRangeTypeCol"  # will fail
}
levels <- c("gene", "transcript", "exon", "cds", "5utr", "3utr")
rangeType <- Rle(factor(c("cds", "gene"), levels=levels), c(8, 2))
elementMetadata(gr)$rangeType <- rangeType
constraint(gr) <- "HasRangeTypeCol"  # OK
\dontrun{
elementMetadata(gr)$rangeType[3] <- NA  # will fail
}

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
        rangeType <- elementMetadata(x)$rangeType
        if (!all(rangeType == "gene")) {
            msg <- c("all elements in 'elementMetadata(x)$rangeType' ",
                     "must be equal to \"gene\"")
            return(paste(msg, collapse=""))
        }
        NULL
    }
)

\dontrun{
constraint(gr) <- "GeneRanges"  # will fail
}
elementMetadata(gr)$rangeType[] <- "gene"
constraint(gr) <- "GeneRanges"  # OK
\dontrun{
elementMetadata(gr)$rangeType[] <- "exon"  # will fail
}

## ---------------------------------------------------------------------
## EXAMPLE 3: The HasGCCol constraint.
## ---------------------------------------------------------------------
## The HasGCCol constraint checks that the constrained object has a
## unique "GC" column in its elementMetadata part, that this column is
## of type numeric, with no NAs, and that all the values in that column
## are >= 0 and <= 1.

setClass("HasGCCol", contains="Constraint")

setMethod("checkConstraint", c("GenomicRanges", "HasGCCol"),
    function(x, constraint, verbose=FALSE)
    {
        emd <- elementMetadata(x)
        idx <- match("GC", colnames(emd))
        if (length(idx) != 1L || is.na(idx)) {
            msg <- c("'elementMetadata(x)' must have exactly ",
                     "one column named \"GC\"")
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

## This replace the previous constraint (GeneRanges):
constraint(gr) <- "HasGCCol"  # OK

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
        GC <- elementMetadata(x)$GC
        if (!all(GC >= 0.5)) {
            msg <- c("all elements in 'elementMetadata(x)$GC' ",
                     "must be >= 0.5")
            return(paste(msg, collapse=""))
        }
        NULL
    }
)

\dontrun{
constraint(gr) <- "HighGCRanges"  # will fail
}
elementMetadata(gr)$GC[6:10] <- 0.5
constraint(gr) <- "HighGCRanges"  # OK

## ---------------------------------------------------------------------
## EXAMPLE 5: The HighGCGeneRanges constraint.
## ---------------------------------------------------------------------
## The HighGCGeneRanges constraint is the combination (AND) of the
## GeneRanges and HighGCRanges constraints.

setClass("HighGCGeneRanges", contains=c("GeneRanges", "HighGCRanges"))

## No need to define a method for this constraint: the checkConstraint()
## generic will automatically check the GeneRanges and HighGCRanges
## constraints.

constraint(gr) <- "HighGCGeneRanges"  # OK

## See how all the individual constraints are checked (from less
## specific to more specific constraints):
checkConstraint(gr, constraint(gr), verbose=TRUE)
}

\keyword{methods}
\keyword{classes}

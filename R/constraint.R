### =========================================================================
### Enforcing constraints thru Constraint objects
### -------------------------------------------------------------------------
###
### Attaching a Constraint object to an object of class A (the "constrained"
### object) is meant to be a convenient/reusable/extensible way to enforce
### a particular set of constraints on particular instances of A. It's an
### alternative to the more traditional approach that consists in creating
### subclasses of A and implementing specific validity methods for each of
### them. However, using constraints offers the following advantages over the
### traditional approach:
###   (a) The traditional approach often tends to lead to a proliferation
###       of subclasses of A.
###   (b) Constraints can easily be re-used across different classes without
###       the need to create any new class.
###   (c) Constraints can easily be combined.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constraint objects.
###

### Virtual class with no slots.
setClass("Constraint", representation("VIRTUAL"))

### Like the Constraint virtual class itself, concrete constraint subclasses
### cannot have slots. A Constraint object doesn't need to contain anything
### anyway so by enforcing this we make sure that "combining" constraints
### (i.e. creating a Constraint subclass that extends more than one concrete
### Constraint subclass) always work smoothly (no slot clashes).
setValidity2("Constraint",
    function(x)
    {
        if (length(slotNames(x)) != 0L)
            return("Constraint objects cannot have slots")
        NULL
    }
)

### Constrained objects can store a Constraint object (or NULL if no
### constraints) in their 'constraint' slot.
setClassUnion("Constraint_OR_NULL", c("Constraint", "NULL"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### constraint() accessor.
###
### Get or set the Constraint object stored in the constrained object.
### Classes that want to support constraints need to implement methods for
### those 2 generics.

setGeneric("constraint", function(x) standardGeneric("constraint"))

setGeneric("constraint<-", signature="x",
    function(x, value) standardGeneric("constraint<-")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The checkConstraint() generic.
###

### Returns a matrix with 1 signature per row. Less specific signatures are
### guaranteed to be *before* (i.e. lower row index) the more specific ones.
.nextCheckConstraintSignatures <- function(xClass, constraintClass)
{
    ## Loosely inspired by validObject().
    xClassDef <- getClassDef(xClass)
    xAncestors <- sapply(xClassDef@contains,
                         slot, "superClass")
    xAncestors <- c("ANY", rev(xAncestors), xClass)
    constraintClassDef <- getClassDef(constraintClass)
    constraintAncestors <- sapply(constraintClassDef@contains,
                                  slot, "superClass")
    constraintAncestors <- c("ANY", rev(constraintAncestors), constraintClass)
    signatures <- NULL
    for (class1 in xAncestors) {
        for (class2 in constraintAncestors) {
            methodDef <- selectMethod("checkConstraint",
                                      c(class1, class2),
                                      optional=TRUE,
                                      doCache=TRUE)
            if (!is.null(methodDef))
                signatures <- rbind(signatures, methodDef@defined)
        }
    }
    unique(signatures)
}

### The checkConstraint() generic function implements its own dispatch
### algorithm. Like validity methods, "checkConstraint" methods must return
### NULL or a character vector describing the problems found.
suppressWarnings(
  setGeneric("checkConstraint", signature=c("x", "constraint"),
    function(x, constraint, verbose=FALSE)
    {
        errors <- NULL
        signatures <- .nextCheckConstraintSignatures(class(x),
                                                     class(constraint))
        if (is.null(signatures))
            return(errors)
        ## We check from less specific to more specific constraints.
        for (i in seq_len(nrow(signatures))) {
            sig <- signatures[i, ]
            sigString <- paste(names(sig),
                               paste0('"', sig, '"'),
                               sep="=", collapse=", ")
            if (verbose)
                message("Calling \"checkConstraint\" method for\n",
                        "    ", sigString)
            method <- getMethod("checkConstraint", sig)
            errors <- method(x, constraint)
            if (length(errors) != 0L) {
                errors <- paste0("from \"checkConstraint\" method for c(",
                                 sigString, "): ", errors)
                ## If a constraint is not satisfied, we don't check the
                ## remaining constraints (so when implementing a constraint
                ## a developer can assume that the less specific constraints
                ## are satisfied).
                break
            }
        }
        errors
    }
  )
)


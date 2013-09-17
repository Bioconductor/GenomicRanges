\name{seqinfo}

\alias{seqinfo}
\alias{seqinfo<-}
\alias{seqnames}
\alias{seqnames<-}
\alias{seqlevels}
\alias{seqlevels,ANY-method}
\alias{seqlevels<-}
\alias{seqlevels<-,ANY-method}
\alias{sortSeqlevels}
\alias{sortSeqlevels,character-method}
\alias{sortSeqlevels,ANY-method}
\alias{seqlevelsInUse}
\alias{seqlevelsInUse,Vector-method}
\alias{seqlevelsInUse,CompressedList-method}
\alias{seqlevels0}
\alias{seqlengths}
\alias{seqlengths,ANY-method}
\alias{seqlengths<-}
\alias{seqlengths<-,ANY-method}
\alias{isCircular}
\alias{isCircular,ANY-method}
\alias{isCircular<-}
\alias{isCircular<-,ANY-method}
\alias{genome}
\alias{genome,ANY-method}
\alias{genome<-}
\alias{genome<-,ANY-method}
\alias{seqnameStyle}
\alias{seqnameStyle,ANY-method}
\alias{seqnameStyle<-}
\alias{seqnameStyle<-,ANY-method}

\alias{seqinfo,List-method}
\alias{seqinfo<-,List-method}
\alias{seqinfo,RangedData-method}
\alias{seqinfo<-,RangedData-method}
\alias{seqinfo,RangesList-method}
\alias{seqnames,RangesList-method}
\alias{seqnames,RangedData-method}

\title{Accessing/modifying sequence information}

\description{
  A set of generic functions for getting/setting/modifying the sequence
  information stored in an object.
}

\usage{
seqinfo(x)
seqinfo(x, new2old=NULL, force=FALSE) <- value

seqnames(x)
seqnames(x) <- value

seqlevels(x)
seqlevels(x, force=FALSE) <- value
sortSeqlevels(x, X.is.sexchrom=NA)
seqlevelsInUse(x)
seqlevels0(x)

seqlengths(x)
seqlengths(x) <- value

isCircular(x)
isCircular(x) <- value

genome(x)
genome(x) <- value

seqnameStyle(x)
seqnameStyle(x) <- value
}

\arguments{
  \item{x}{
    The object from/on which to get/set the sequence information.
  }
  \item{new2old}{
    The \code{new2old} argument allows the user to rename, drop, add and/or
    reorder the "sequence levels" in \code{x}.

    \code{new2old} can be \code{NULL} or an integer vector with one element
    per row in \link{Seqinfo} object \code{value} (i.e. \code{new2old} and
    \code{value} must have the same length) describing how the "new" sequence
    levels should be mapped to the "old" sequence levels, that is, how the
    rows in \code{value} should be mapped to the rows in \code{seqinfo(x)}.
    The values in \code{new2old} must be >= 1 and <= \code{length(seqinfo(x))}.
    \code{NA}s are allowed and indicate sequence levels that are being added.
    Old sequence levels that are not represented in \code{new2old} will be
    dropped, but this will fail if those levels are in use (e.g. if \code{x}
    is a \link{GRanges} object with ranges defined on those sequence levels)
    unless \code{force=TRUE} is used (see below).

    If \code{new2old=NULL}, then sequence levels can only be added to the
    existing ones, that is, \code{value} must have at least as many rows
    as \code{seqinfo(x)} (i.e. \code{length(values) >= length(seqinfo(x))})
    and also \code{seqlevels(values)[seq_len(length(seqlevels(x)))]} must be
    identical to \code{seqlevels(x)}.
  }
  \item{force}{
    Force dropping sequence levels currently in use. This is achieved by
    dropping the elements in \code{x} where those levels are used (hence
    typically reducing the length of \code{x}).
  }
  \item{value}{
    Typically a \link{Seqinfo} object for the \code{seqinfo} setter.

    Either a named or unnamed character vector for the \code{seqlevels}
    setter.

    A vector containing the sequence information to store for the other
    setters.
  }
}

\details{
  The \link{Seqinfo} class plays a central role for the functions described
  in this man page because:

  \itemize{
  \item All these functions (except \code{seqinfo}, \code{seqlevelsInUse},
        and \code{seqlevels0}) work on a \link{Seqinfo} object.

  \item For classes that implement it, the \code{seqinfo} getter should
        return a \link{Seqinfo} object.

  \item Default \code{seqlevels}, \code{seqlengths}, \code{isCircular},
        \code{genome}, and \code{seqnameStyle} getters and setters are
        provided.
        By default, \code{seqlevels(x)} does \code{seqlevels(seqinfo(x))},
        \code{seqlengths(x)} does \code{seqlengths(seqinfo(x))},
        \code{isCircular(x)} does \code{isCircular(seqinfo(x))},
        \code{genome(x)} does \code{genome(seqinfo(x))},
        and \code{seqnameStyle(x)} does \code{seqnameStyle(seqinfo(x))}.
        So any class with a \code{seqinfo} getter will have all the above
        getters work out-of-the-box. If, in addition, the class defines
        a \code{seqinfo} setter, then all the corresponding setters will
        also work out-of-the-box.

        See the \link{GRanges}, \link{GRangesList}, \link{GAlignments},
        and \link{GAlignmentPairs} classes for examples of classes that
        define the \code{seqinfo} getter and setter (those 4 classes are
        defined in the GenomicRanges package).

        See the \link[GenomicFeatures]{TranscriptDb} class (defined in the
        GenomicFeatures package) for an example of a class that defines only
        the \code{seqinfo} getter (no setter).
      }

      The GenomicRanges package defines \code{seqinfo} and
      \code{seqinfo<-} methods for these low-level IRanges data
      structures: \code{List}, \code{RangesList} and
      \code{RangedData}. Those objects do not have the means to formally
      store sequence information. Thus, the wrappers simply store the
      \code{Seqinfo} object within \code{metadata(x)}. Initially, the
      metadata is empty, so there is some effort to generate a
      reasonable default \code{Seqinfo}. The names of any \code{List}
      are taken as the \code{seqnames}, and the \code{universe} of
      \code{RangesList} or \code{RangedData} is taken as the
      \code{genome}.
}

\note{
  The full list of methods defined for a given generic can
  be seen with e.g. \code{showMethods("seqinfo")} or
  \code{showMethods("seqnames")} (for the getters),
  and \code{showMethods("seqinfo<-")} or \code{showMethods("seqnames<-")}
  (for the setters aka \emph{replacement methods}).
  Please be aware that this shows only methods defined in packages
  that are currently attached.
}

\author{H. Pages}

\seealso{
  \itemize{
    \item \link{Seqinfo-class}.

    \item \link{GRanges}, \link{GRangesList},
          \link{GAlignments}, \link{GAlignmentPairs},
          \link{GAlignmentsList}, \link{SummarizedExperiment},
          and \link[GenomicFeatures]{TranscriptDb}, for examples
          of objects on which the functions described in this man page
          work.

    \item \link{seqlevels-utils} for convenience wrappers to the
          \code{seqlevels} getter and setter.

    \item \code{\link{makeSeqnameIds}}, on which \code{sortSeqlevels} is
          based.
  }
}

\examples{
## ---------------------------------------------------------------------
## Finding methods.
## ---------------------------------------------------------------------

showMethods("seqinfo")
showMethods("seqinfo<-")

showMethods("seqnames")
showMethods("seqnames<-")

showMethods("seqlevels")
showMethods("seqlevels<-")

if (interactive())
  ?`GRanges-class`

## ---------------------------------------------------------------------
## Modify seqlevels of an object.
## ---------------------------------------------------------------------

## Overlap and matching operations between objects require matching
## seqlevels. Often the seqlevels in one must be modified to match 
## the other. The seqlevels() function can rename, drop, add and reorder 
## seqlevels of an object. Examples below are shown on TranscriptDb 
## and GRanges but the approach is the same for all objects that have
## a 'Seqinfo' class.

library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
seqlevels(txdb)

## Rename:
seqlevels(txdb) <- sub("chr", "", seqlevels(txdb))
seqlevels(txdb)

seqlevels(txdb) <- paste0("CH", seqlevels(txdb))
seqlevels(txdb)

seqlevels(txdb)[seqlevels(txdb) == "CHM"] <- "M"
seqlevels(txdb)

## Add:
gr <- GRanges(c("chr1", "chr2", "chr3"), IRanges(1:3, 5))
seqlevels(gr)
seqlevels(gr) <- c(seqlevels(gr), "chr4")
seqlevels(gr)

## Reorder:
seqlevels(gr) <- c("chr2", "chr3", "chr1", "chr4")
seqlevels(gr)

## Drop:
seqlevels(gr, force=TRUE) <- c("chr2", "chr1", "chr4")
seqlevels(gr)

## Rename/Add/Reorder:
seqlevels(gr) <- c(chr1="1", chr2="2", "chr4", "chr5")
seqlevels(gr)


## ---------------------------------------------------------------------
## Sort seqlevels in "natural" order
## ---------------------------------------------------------------------

sortSeqlevels(c("11", "Y", "1", "10", "9", "M", "2"))

seqlevels <- c("chrXI", "chrY", "chrI", "chrX", "chrIX", "chrM", "chrII")
sortSeqlevels(seqlevels)
sortSeqlevels(seqlevels, X.is.sexchrom=TRUE)
sortSeqlevels(seqlevels, X.is.sexchrom=FALSE)

seqlevels <- c("chr2RHet", "chr4", "chrUextra", "chrYHet",
               "chrM", "chrXHet", "chr2LHet", "chrU",
               "chr3L", "chr3R", "chr2R", "chrX")
sortSeqlevels(seqlevels)

gr <- GRanges()
seqlevels(gr) <- seqlevels
sortSeqlevels(gr)

## ---------------------------------------------------------------------
## Subset objects by seqlevels. 
## ---------------------------------------------------------------------

tx <- transcripts(txdb)
seqlevels(tx)

## Drop 'M', keep all others.
seqlevels(tx, force=TRUE) <- seqlevels(tx)[seqlevels(tx) != "M"]
seqlevels(tx)

## Drop all except 'ch3L' and 'ch3R'.
seqlevels(tx, force=TRUE) <- c("ch3L", "ch3R")
seqlevels(tx)

## ---------------------------------------------------------------------
## Restore original seqlevels. 
## ---------------------------------------------------------------------

## Applicable to TranscriptDb objects only.
\dontrun{
seqlevels0(txdb)
seqlevels(txdb)
}
}

\keyword{methods}

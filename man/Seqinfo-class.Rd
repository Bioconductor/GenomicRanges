\name{Seqinfo-class}
\docType{class}

% Classes:
\alias{class:Seqinfo}
\alias{Seqinfo-class}
\alias{Seqinfo}

% Constructors:
\alias{Seqinfo}
\alias{updateObject,Seqinfo-method}

% Methods:
\alias{length,Seqinfo-method}
\alias{seqnames,Seqinfo-method}
\alias{seqnames<-,Seqinfo-method}
\alias{names,Seqinfo-method}
\alias{names<-,Seqinfo-method}
\alias{seqlevels,Seqinfo-method}
\alias{seqlevels<-,Seqinfo-method}
\alias{seqlengths,Seqinfo-method}
\alias{seqlengths<-,Seqinfo-method}
\alias{isCircular,Seqinfo-method}
\alias{isCircular<-,Seqinfo-method}
\alias{genome,Seqinfo-method}
\alias{genome<-,Seqinfo-method}
\alias{seqnameStyle}
\alias{seqnameStyle,ANY-method}
\alias{seqnameStyle,Seqinfo-method}
\alias{[,Seqinfo-method}
\alias{coerce,Seqinfo,GenomicRanges-method}
\alias{coerce,Seqinfo,RangesList-method}
\alias{as.data.frame,Seqinfo-method}
\alias{show,Seqinfo-method}
\alias{merge,Seqinfo,missing-method}
\alias{merge,missing,Seqinfo-method}
\alias{merge,Seqinfo,NULL-method}
\alias{merge,NULL,Seqinfo-method}
\alias{merge,Seqinfo,Seqinfo-method}


\title{Seqinfo objects}

\description{
  A Seqinfo object is a table-like object that contains basic information
  about a set of genomic sequences. The table has 1 row per sequence and
  1 column per sequence attribute.  Currently the only attributes are the
  length, circularity flagi, and genome provenance (e.g. hg19) of the
  sequence, but more attributes might be added in the future as the need
  arises.
}

\details{
  Typically Seqinfo objects are not used directly but are part of
  higher level objects. Those higher level objects will generally
  provide a \code{seqinfo} accessor for getting/setting their
  Seqinfo component.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{Seqinfo(seqnames, seqlengths=NA, isCircular=NA, genome=NA)}:
      Creates a Seqinfo object.
    }
  }
}

\section{Accessor methods}{
  In the code snippets below, \code{x} is a Seqinfo object.

  \describe{
    \item{}{
      \code{length(x)}:
      Gets the number of sequences in \code{x}.
    }
    \item{}{
      \code{seqnames(x)}, \code{seqnames(x) <- value}:
      Gets/sets the names of the sequences in \code{x}.
      Those names must be non-NA, non-empty and unique.
      They are also called the \emph{sequence levels} or the \emph{keys}
      of the Seqinfo object.

      Note that, in general, the end-user should not try to alter the
      sequence levels with \code{seqnames(x) <- value}. The recommended way
      to do this is with \code{seqlevels(x) <- value} as described below.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
      Same as \code{seqnames(x)} and \code{seqnames(x) <- value}.
    }
    \item{}{
      \code{seqlevels(x)}:
      Same as \code{seqnames(x)}.
    }
    \item{}{
      \code{seqlevels(x) <- value}:
      Can be used to rename, drop, add and/or reorder the sequence levels.
      \code{value} must be either a named or unnamed character vector.
      When \code{value} has names, the names only serve the purpose of
      mapping the new sequence levels to the old ones.
      Otherwise (i.e. when \code{value} is unnamed) this mapping is
      implicitly inferred from the following rules:

      (1) If the number of new and old levels are the same, and if the
          positional mapping between the new and old levels shows that
          some or all of the levels are being renamed, and if the levels
          that are being renamed are renamed with levels that didn't exist
          before (i.e. are not present in the old levels), then
          \code{seqlevels(x) <- value} will just rename the sequence levels.
          Note that in that case the result is the same as with
          \code{seqnames(x) <- value} but it's still recommended to use
          \code{seqlevels(x) <- value} as it is safer.

      (2) Otherwise (i.e. if the conditions for (1) are not satisfied)
          \code{seqlevels(x) <- value} will consider that the sequence
          levels are not being renamed and will just perform
          \code{x <- x[value]}.

      See below for some examples.
    }
    \item{}{
      \code{seqlengths(x)}, \code{seqlengths(x) <- value}:
      Gets/sets the length for each sequence in \code{x}.
    }
    \item{}{
      \code{isCircular(x)}, \code{isCircular(x) <- value}:
      Gets/sets the circularity flag for each sequence in \code{x}.
    }
    \item{}{
      \code{genome(x)}, \code{genome(x) <- value}:
      Gets/sets the genome identifier or assembly name for each sequence
      in \code{x}.
    }
    \item{}{
      \code{seqnameStyle(x)}:
      Lists the matching seqname styles for \code{x}.
    }
  }
}

\section{Subsetting}{
  In the code snippets below, \code{x} is a Seqinfo object.

  \describe{
    \item{}{
      \code{x[i]}:
      A Seqinfo object can be subsetted only by name i.e. \code{i}
      must be a character vector.
      This is a convenient way to drop/add/reorder the rows (aka the
      sequence levels) of a Seqinfo object.

      See below for some examples.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a Seqinfo object.

  \describe{
    \item{}{
      \code{as.data.frame(x)}:
      Turns \code{x} into a data frame.
    }
    \item{}{
      \code{as(x, "GenomicRanges")}, \code{as(x, "RangesList")}:
      Turns \code{x} (with no \code{NA} lengths) into a GRanges or
      RangesList.
    }
  }
}

\section{Combining Seqinfo objects}{
  There are no \code{c} or \code{rbind} method for Seqinfo objects.
  Both would be expected to just append the rows in \code{y} to the rows
  in \code{x} resulting in an object of length \code{length(x) + length(y)}.
  But that would tend to break the constraint that the seqnames of a Seqinfo
  object must be unique keys.

  So instead, a \code{merge} method is provided.

  In the code snippet below, \code{x} and \code{y} are Seqinfo objects.

  \describe{
    \item{}{
      \code{merge(x, y)}:
      Merge \code{x} and \code{y} into a single Seqinfo object where the
      keys (aka the seqnames) are \code{union(seqnames(x), seqnames(y))}.
      If a row in \code{y} has the same key as a row in \code{x}, and if
      the 2 rows contain compatible information (NA values are compatible
      with anything), then they are merged into a single row in the result.
      If they cannot be merged (because they contain different seqlengths,
      and/or circularity flags, and/or genome identifiers), then an error
      is raised.
      In addition to check for incompatible sequence information,
      \code{merge(x, y)} also compares \code{seqnames(x)} with
      \code{seqnames(y)} and issues a warning if each of them has names not
      in the other. The purpose of these checks is to try to detect situations
      where the user might be combining or comparing objects based on
      different reference genomes.
    }
  }
}

\author{H. Pages}

\seealso{
  \code{\link{seqinfo}}
}

\examples{
  ## Note that all the arguments (except 'genome') must have the
  ## same length. 'genome' can be of length 1, whatever the lengths
  ## of the other arguments are.
  x <- Seqinfo(seqnames=c("chr1", "chr2", "chr3", "chrM"),
               seqlengths=c(100, 200, NA, 15),
               isCircular=c(NA, FALSE, FALSE, TRUE),
               genome="toy")
  x

  x[c("chrY", "chr3", "chr1")]  # subset by names

  ## Rename, drop, add and/or reorder the sequence levels:
  xx <- x
  seqlevels(xx) <- sub("chr", "ch", seqlevels(xx))  # rename
  xx
  seqlevels(xx) <- rev(seqlevels(xx))  # reorder
  xx
  seqlevels(xx) <- c("ch1", "ch2", "chY")  # drop/add/reorder
  xx
  seqlevels(xx) <- c(chY="Y", ch1="1", "22")  # rename/reorder/drop/add
  xx

  y <- Seqinfo(seqnames=c("chr3", "chr4", "chrM"),
               seqlengths=c(300, NA, 15))
  y
  merge(x, y)  # rows for chr3 and chrM are merged
  suppressWarnings(merge(x, y))

  ## Note that, strictly speaking, merging 2 Seqinfo objects is not
  ## a commutative operation, i.e., in general 'z1 <- merge(x, y)'
  ## is not identical to 'z2 <- merge(y, x)'. However 'z1' and 'z2'
  ## are guaranteed to contain the same information (i.e. the same
  ## rows, but typically not in the same order):
  suppressWarnings(merge(y, x))

  ## This contradicts what 'x' says about circularity of chr3 and chrM:
  isCircular(y)[c("chr3", "chrM")] <- c(TRUE, FALSE)
  y
  if (interactive()) {
    merge(x, y)  # raises an error
  }
}

\keyword{methods}
\keyword{classes}

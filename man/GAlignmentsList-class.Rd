\name{GAlignmentsList-class}
\docType{class}

% Class
\alias{class:GAlignmentsList}
\alias{GAlignmentsList-class}
\alias{GAlignmentsList}

% Constructors:
\alias{GAlignmentsList}
\alias{makeGAlignmentsListFromFeatureFragments}
\alias{updateObject,GAlignmentsList-method}

% Accessors:
\alias{names,GAlignmentsList-method}
\alias{names<-,GAlignmentsList-method}
\alias{seqnames,GAlignmentsList-method}
\alias{seqnames<-,GAlignmentsList-method}
\alias{rname,GAlignmentsList-method}
\alias{rname<-,GAlignmentsList-method}
\alias{strand,GAlignmentsList-method}
\alias{strand<-,GAlignmentsList-method}
\alias{cigar,GAlignmentsList-method}
\alias{qwidth,GAlignmentsList-method}
\alias{ngap,GAlignmentsList-method}
\alias{elementMetadata,GAlignmentsList-method}
\alias{elementMetadata<-,GAlignmentsList-method}
\alias{seqinfo,GAlignmentsList-method}
\alias{seqinfo<-,GAlignmentsList-method}
\alias{start,GAlignmentsList-method}
\alias{end,GAlignmentsList-method}
\alias{width,GAlignmentsList-method}

% Coercion:
\alias{as.data.frame,GAlignmentsList-method}
\alias{coerce,GAlignmentsList,CompressedIRangesList-method}
\alias{coerce,GAlignmentsList,IRangesList-method}
\alias{coerce,GAlignmentsList,GRangesList-method}
\alias{coerce,GAlignmentsList,RangesList-method}
\alias{granges,GAlignmentsList-method}
\alias{ranges,GAlignmentsList-method}

% Subsetting:
\alias{[,GAlignmentsList-method}
\alias{[<-,GAlignmentsList,ANY,ANY,ANY-method}
\alias{[[<-,GAlignmentsList-method}
\alias{[[<-,GAlignmentsList,ANY,ANY-method}

% Combining:
\alias{c,GAlignmentsList-method}

% show:
\alias{show,GAlignmentsList-method}

% Other methods:

\title{GAlignmentsList objects}

\description{
  The GAlignmentsList class is a container for storing a collection of
  GappedAlignments objects.
}

\details{
  A GAlignmentsList object contains a list of GappedAlignments. The
  majority of operations on this page are described in greater detail
  on the GappedAlignments man page (?\code{GappedAlignments}).
} 

\section{Constructors}{
  \describe{
    \item{}{
      \code{GAlignmentsList(...)}:
      Creates a GAlignmentsList from a list of GappedAlignments objects.
    }
    \item{}{
      \code{makeGAlignmentsListFromFeatureFragments(seqnames=Rle(factor()), 
                                                    fragmentPos=list(),
                                                    fragmentCigar=list(),
                                                    strand=character(0),
                                                    sep=",")}:
      Constructs a GAlignmentsList from a list of fragmented features.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{x} is a GAlignmentsList object.

  \describe{
    \item{}{
      \code{length(x)}:
      Return the number of elements in \code{x}.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
      Get or set the names of the elements of \code{x}.
    }
    \item{}{
      \code{seqnames(x)}, \code{seqnames(x) <- value}:
      Get or set the name of the reference sequences of the
      alignments in each element of \code{x}.
    }
    \item{}{
      \code{rname(x)}, \code{rname(x) <- value}:
      Same as \code{seqnames(x)} and \code{seqnames(x) <- value}.
    }
    \item{}{
      \code{strand(x)}, \code{strand(x) <- value}:
      Get or set the strand of the alignments in each element 
      of \code{x}.
    }
    \item{}{
      \code{cigar(x)}:
      Returns a character list of length \code{length(x)}
      containing the CIGAR string for the alignments in
      each element of \code{x}.
    }
    \item{}{
      \code{qwidth(x)}:
      Returns an integer list of length \code{length(x)}
      containing the length of the alignments in each element of
      \code{x} *after* hard clipping (i.e. the length of the 
      query sequence that is stored in the corresponding SAM/BAM record).
    }
    \item{}{
      \code{start(x)}, \code{end(x)}:
      Returns an integer list of length \code{length(x)}
      containing the "start" and "end" (respectively) of the 
      alignments in each element of \code{x}. 
    }
    \item{}{
      \code{width(x)}:
      Returns an integer list of length \code{length(x)} containing
      the "width" of the alignments in each element of \code{x}.
    }
    \item{}{
      \code{ngap(x)}:
      Returns an integer list of length \code{x} containing the number
      of gaps (i.e. N operations in the CIGAR) for the alignments
      in each element of \code{x}.
    }
    \item{}{
      \code{seqinfo(x)}, \code{seqinfo(x) <- value}:
      Get or set the information about the underlying sequences in each
      element of \code{x}. \code{value} must be a list of \link{Seqinfo} 
      objects.
    }
    \item{}{
      \code{seqlevels(x)}, \code{seqlevels(x) <- value}:
      Get or set the sequence levels of the alignments in each element
      of \code{x}.
    }
    \item{}{
      \code{seqlengths(x)}, \code{seqlengths(x) <- value}:
      Get or set the sequence lengths for each element of \code{x}.
      \code{seqlengths(x)} is equivalent to \code{seqlengths(seqinfo(x))}.
      \code{value} can be a named non-negative integer or numeric vector
      eventually with NAs.
    }
    \item{}{
      \code{isCircular(x)}, \code{isCircular(x) <- value}:
      Get or set the circularity flags for the alignments in each
      element in \code{x}. \code{value} must be a named logical list 
      eventually with NAs.
    }
    \item{}{
      \code{genome(x)}, \code{genome(x) <- value}:
      Get or set the genome identifier or assembly name for the alignments 
      in each element of \code{x}. \code{value} must be a named character 
      list eventually with NAs.
    }
    \item{}{
      \code{seqnameStyle(x)}:
      Get or set the seqname style for alignments in each element of \code{x}.
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a GAlignmentsList object.

  \describe{
    \item{}{
      \code{granges(x)}, \code{ranges(x)}:
      Return either a \link{GRangesList} or a \link[IRanges]{IRangesList}
      object of length \code{length(x)} where each element is a 
      \link{GRanges} or \link{IRanges} object.

      Note that the gap information (cigar) is lost when a GappedAlignments
      is coerced to GRanges or IRanges. 
    }
    \item{}{
      \code{as(x, "GRangesList")}, \code{as(x, "GRanges")},
      \code{as(x, "IRangesList")}, \code{as(x, "Ranges")}:
      Alternate ways of doing \code{granges(x)}, \code{ranges(x)},
      respectively.
    }
    \item{}{
      \code{as(x, "CompressedIRangesList")}:
      Coerce a GAlignmentsList to CompressedIRangesList.
    }
  }
}

\section{Subsetting and related operations}{
  In the code snippets below, \code{x} is a GAlignmentsList object.

  \describe{
    \item{}{
      \code{x[i]}, \code{x[i] <- value}:
      Get or set list elements \code{i}. \code{i} can be a numeric 
      or logical vector. \code{value} must be a GappedAlignments.
    }
    \item{}{
      \code{x[[i]]}, \code{x[[i]] <- value}:
      Same as \code{x[i]}, \code{x[i] <- value}. 
    }
    \item{}{
      \code{x[i, j]}, \code{x[i, j] <- value}:
      Get or set list elements \code{i} with optional metadata columns
      \code{j}. \code{i} can be a numeric, logical or missing. 
      \code{value} must be a GappedAlignments.
    }
  }
}

\section{Combining}{
  \describe{
    \item{}{
      \code{c(...)}:
      Concatenates the GAlignmentList objects in \code{...}.
    }
  }
}

\references{
  \url{http://samtools.sourceforge.net/}
}

\author{Valerie Obenchain <vobencha@fhcrc.org}

\seealso{
  \itemize{
    \item \link{GappedAlignmentPairs-class}.
    \item \code{\link[Rsamtools]{readBamGappedAlignments}}.
    \item \link{GRangesList-class}.
    \item \link{GRanges-class}.
    \item \link{findOverlaps-methods}.
    \item \link{coverage-methods}.
    \item \code{\link{seqinfo}}.
    \item \link[IRanges]{CompressedIRangesList-class}.
    \item \link{setops-methods}.
  }
}

\examples{
gal1 <- GappedAlignments(
    seqnames=Rle(factor(c("chr1", "chr2", "chr1", "chr3")),
        c(1, 3, 2, 4)),
    pos=1:10, cigar=paste0(10:1, "M"),
    strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
    names=head(letters, 10), score=1:10)

gal2 <- GappedAlignments(
    seqnames=Rle(factor(c("chr2", "chr4")), c(3, 4)), pos=1:7,
    cigar=c("5M", "3M2N3M2N3M", "5M", "10M", "5M1N4M", "8M2N1M", "5M"),
    strand=Rle(strand(c("-", "+")), c(4, 3)),
    names=tail(letters, 7), score=1:7)

gaList <- GAlignmentsList(noGaps=gal1, Gaps=gal2)


## ---------------------------------------------------------------------
## A. BASIC MANIPULATION
## ---------------------------------------------------------------------
length(gaList)
names(gaList)
seqnames(gaList)
strand(gaList)
head(cigar(gaList))
head(qwidth(gaList))
head(start(gaList))
head(end(gaList))
head(width(gaList))
head(ngap(gaList))
seqlevels(gaList)

## Rename the reference sequences:
seqlevels(gaList) <- sub("chr", "seq", seqlevels(gaList))
seqlevels(gaList)

granges(gaList)  # a GRangesList object
ranges(gaList)   # an IRanges object

## ---------------------------------------------------------------------
## B. SUBSETTING
## ---------------------------------------------------------------------

gaList[strand(gaList) == "-"]
gaps <- sapply(gaList, function(x) any(grepl("N", cigar(x), fixed=TRUE)))
gaList[gaps]

## Different ways to subset:
gaList[2]             # a GappedAlignments object of length 1
gaList[[2]]           # a GappedAlignments object of length 1
granges(gaList)[[2]]  # a GRanges object of length 1
ranges(gaList)[[2]]   # a NormalIRanges object of length 1

## ---------------------------------------------------------------------
## C. mcols()/elementMetadata()
## ---------------------------------------------------------------------

## Metadata can exist for both the individual GappedAlignments objects 
## and the overall GAlignmentsList object. By default, 'level=between' 
## which extracts metadata set on the GAlignmentsList. Using 'level=within' 
## will extract the metadata on the GappedAliginments.

mcols(gaList) ## no metadata on the GAlignmentsList object
mcols(gaList, level="within")
}

\keyword{methods}
\keyword{classes}

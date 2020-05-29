\name{DelegatingGenomicRanges-class}
\docType{class}

% Class:
\alias{class:DelegatingGenomicRanges}
\alias{DelegatingGenomicRanges-class}
\alias{class:DelegatingGenomicRanges}
\alias{DelegatingGenomicRanges-class}

% Accessors:
\alias{seqnames,DelegatingGenomicRanges-method}
\alias{ranges,DelegatingGenomicRanges-method}
\alias{strand,DelegatingGenomicRanges-method}
\alias{seqinfo,DelegatingGenomicRanges-method}

% Updating:
\alias{update,DelegatingGenomicRanges-method}

% Coercion:
\title{DelegatingGenomicRanges objects}

\description{
  The \code{DelegatingGenomicRanges} class implements the virtual
  \code{GenomicRanges} class using a delegate \code{GenomicRanges}. This
  enables developers to create \code{GenomicRanges} subclasses that add
  specialized columns or other structure, while remaining agnostic to
  how the data are actually stored. See the Extending GenomicRanges
  vignette.
}

\author{M. Lawrence}

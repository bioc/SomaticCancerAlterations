
\name{SomaticCancerAlterations-helper-functions}
     
\alias{SomaticCancerAlterations-helper-functions}

\alias{ncbi2hg}
\alias{hg2ncbi}
\alias{seqchr}

\title{SomaticCancerAlterations helper functions}

\description{
  Functions to help with common analysis tasks.
}

\usage{
ncbi2hg(x)
 
hg2ncbi(x)
 
seqchr(x)
}

\arguments{
  \item{x}{A GRanges object}
}

\value{
  \item{ncbi2hg, hg2ncbi}{A GRanges object.}
  \item{seqchr}{A character vector with the 'seqnames'.}
}
 
\details{
  These functions provide useful wrappers for converting between genomic
  coordinate systems and classes.
}

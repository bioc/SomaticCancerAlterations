
\name{SomaticCancerAlterations-functions}
    
\alias{SomaticCancerAlterations-functions}

\alias{scaListDatasets}
\alias{scaLoadDatasets}
\alias{scaMetadata}
\alias{mutationDensity}

\title{SomaticCancerAlterations functions}

\description{
  Functions to interact with the 'SomaticCancerAlterations' datasets.
}

\usage{
scaMetadata()

scaListDatasets()

scaLoadDatasets(names, merge = FALSE)

mutationDensity(x, binSize = 1e6, chrs = seqnames(x))
}

\arguments{
  \item{names}{Character vector with data set name(s) to load. The names
    are the same as returned from 'scaListDatasets()'.}
  \item{merge}{Logical.  Should the results be a merged 'GRanges' object
  with data of all studies (TRUE), or a 'GenomicRangesList' with one list
  element per dataset (FALSE [default]).}
  \item{x}{GRanges object}
  \item{binSize}{Integer with the bin size used for computing the
    density of variants.}
  \item{chrs}{Character vector with the chrosomes/seqlevels to include
    (by default, all chromosomes are used.)}
}

\value{
  \item{scaMetadata}{A data frame summarizing the available data
    sets, with rows corresponding to data sets.}
  \item{scaListDatasets}{A character vector with all available study
    names.}
  \item{sca_load_dataset}{A GenomicRangesList or GRanges objecs, depending
    on the 'merge' argument.}
  \item{mutationDensity}{A GRanges object, with ranges corresponding to
    bins along the genome.  The columns 'counts' and 'density' refer to
    the number and the fraction (given per 1Mbp) of alterations,
    respectively.}
}

\details{
  The 'sca_load_dataset' function makes it easy to load the data of a study
  directly into a variable.

  The 'mutationDensity' function computes the density of events in a
  binned-manner.
}
  
\references{
  \url{https://tcga-data.nci.nih.gov/tcga/}
  
  \url{https://wiki.nci.nih.gov/display/TCGA/TCGA+Data+Primer}
  
  \url{https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/}
}

\seealso{
  \code{\link[SomaticCancerAlterations]{SomaticCancerAlterations-data}}
}

\examples{
meta_data = scaMetadata()
head(meta_data)

all_datasets = scaListDatasets()
head(all_datasets)

grl = scaLoadDatasets(all_datasets[1], merge = FALSE)

gr = scaLoadDatasets(all_datasets[1], merge = TRUE)
}

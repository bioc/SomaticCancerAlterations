
\name{SomaticCancerAlterations-functions}
    
\alias{SomaticCancerAlterations-functions}

\alias{sca_list_datasets}
\alias{sca_load_datasets}
\alias{sca_metadata}
\alias{mutation_density}

\title{SomaticCancerAlterations functions}

\description{
  Functions to interact with the 'SomaticCancerAlterations' datasets.
}

\usage{
sca_metadata()

sca_list_datasets()

sca_load_datasets(names, merge = FALSE)

mutation_density(x, bin_size = 1e6, chrs = seqnames(x))
}

\arguments{
  \item{names}{Character vector with data set name(s) to load. The names
    are the same as returned from 'sca_list_datasets()'.}
  \item{merge}{Logical.  Should the results be a merged 'GRanges' object
  with data of all studies (TRUE), or a 'GenomicRangesList' with one list
  element per dataset (FALSE [default]).}
  \item{x}{GRanges object}
  \item{bin_size}{Integer with the bin size used for computing the
    density of variants.}
  \item{chrs}{Character vector with the chrosomes/seqlevels to include
    (by default, all chromosomes are used.)}
}

\value{
  \item{sca_metadata}{A data frame summarizing the available data
    sets, with rows corresponding to data sets.}
  \item{sca_list_datasets}{A character vector with all available study
    names.}
  \item{sca_load_dataset}{A GenomicRangesList or GRanges objecs, depending
    on the 'merge' argument.}
  \item{mutation_density}{A GRanges object, with ranges corresponding to
    bins along the genome.  The columns 'counts' and 'density' refer to
    the number and the fraction (given per 1Mbp) of alterations,
    respectively.}
}

\details{
  The 'sca_load_dataset' function makes it easy to load the data of a study
  directly into a variable.

  The 'mutation_density' function computes the density of events in a
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
meta_data = sca_metadata()
head(meta_data)

all_datasets = sca_list_datasets()
head(all_datasets)

grl = sca_load_datasets(all_datasets[1], merge = FALSE)

gr = sca_load_datasets(all_datasets[1], merge = TRUE)
}

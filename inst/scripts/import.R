
library(SomaticCancerAlterations)
library(GenomicRanges)
library(IRanges)
library(stringr)
#source("code/package/SomaticCancerAlterations/R/import-functions.R")

### parameters ###
open_access_cancer_types = c("GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "SKCM", "THCA")
#open_access_cancer_types = c("GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "SKCM", "THCA", "UCEC")
import_cancer_types = open_access_cancer_types[!grepl("UCEC", open_access_cancer_types)] ## ignore 'UCEC' at the moment

data_dir = "/ebi/research/huber/users/jgehring/projects/28-cbioportal/code/package/SomaticCancerAlterations/data"
dir.create(data_dir, showWarnings = FALSE)

base_regexp = "(.*)_([A-Z]{2,7}).(.*)_(.*).Level_(\\d+.\\d+.\\d+.\\d+)"
base_names = c("institute", "cancer", "sequencer", "source", "version")


### find the maf files ###

maf_dir = "/ebi/research/huber/users/jgehring/projects/28-cbioportal/data/tcga/mutations/maf/new"
maf_select_pattern = ".*broad.mit.edu.*.maf$" ## only import 'broad' studies
maf_paths = dir(path = maf_dir, pattern = ".*.maf", recursive = TRUE, full.names = TRUE)

seq_info_path = "/ebi/research/huber/users/jgehring/projects/28-cbioportal/data/tcga/hs37d5_seq_info.Rda"

maf_select_paths = grep(maf_select_pattern, maf_paths, value = TRUE)
maf_select_paths = sapply(import_cancer_types, grep, maf_select_paths, value = TRUE)

### meta data from path names ###
study_base = basename(dirname(maf_select_paths))
studies = str_match(study_base, base_regexp)[ ,-1] ## ignore full match
colnames(studies) = base_names

#stopifnot(!anyDuplicated())

names(maf_select_paths) = studies[ ,"cancer"]

load(seq_info_path)

meta_list = list()
for(i in seq_along(maf_select_paths)) {
    cancer = names(maf_select_paths)[i]
    print(cancer)
    var_name = sprintf("%s_tcga", tolower(cancer))
    save_path = file.path(data_dir, sprintf("%s.rda", var_name))
    maf = .read_maf(maf_select_paths[i])
    stopifnot(.validate_maf(maf, seq_info_hs37d5))
    maf = .reconstruct_ids(maf)
    meta_list[[var_name]] = .meta_maf(maf, cancer)
    gr = .maf2gr(maf, seq_info_hs37d5)
    metadata(gr) = meta_list[[var_name]]
    assign(var_name, gr)
    save(file = save_path,
         list = c(var_name))
}

meta = as.data.frame(t(sapply(meta_list, c)))
names(meta) = names(meta)

save_path = file.path(data_dir, "meta.rda")
save(file = save_path, "meta")

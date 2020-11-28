
scaListDatasets <- function() {
    x = data(package = "SomaticCancerAlterations")$results[ ,"Item"]
    x = x[(x != "meta")]
    return(x)
}

scaLoadDatasets <- function(names, merge = FALSE) {
    all_datasets = scaListDatasets()
    if(missing(names))
        names = all_datasets
    if(!all(idx <- names %in% all_datasets)) {
        msg = sprintf("Data set %s not found.",
            paste(names[!idx], collapse = ", "))
        stop(msg)
    }
    
    x = sapply(names, .load_dataset, simplify = FALSE, USE.NAMES = TRUE)
    res = GRangesList(unlist(x), compress=FALSE)
    if(merge) {
        datasets = rep(factor(rep(names(res), lengths(res))))
        res = unlist(res)
        res$Dataset = datasets
    }

    return(res)
}

scaMetadata <- function() {
    res = .load_dataset("meta")
    return(res)
}
    
  
.load_dataset <- function(name, package = "SomaticCancerAlterations") {
    tmp_env = new.env()
    data(list = name, package = package, envir = tmp_env)
    res = get(name, envir = tmp_env)
    return(res)
}


ncbi2hg <- function(x) {
    seqnameStyle(x) = "ucsc"
    genome(x) = NA
    return(x)
}


hg2ncbi <- function(x) {
    seqnameStyle(x) = "ncbi"
    genome(x) = NA
    return(x)
}


seqchr <- function(x) {
    y = as.character(seqnames(x))
    return(y)
}

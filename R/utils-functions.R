
sca_list_datasets <- function() {
    x = data(package = "SomaticCancerAlterations")$results[ ,"Item"]
    x = x[(x != "meta")]
    return(x)
}

sca_load_datasets <- function(names, merge = FALSE) {
    all_datasets = sca_list_datasets()
    if(missing(names))
        names = all_datasets
    if(!all(idx <- names %in% all_datasets)) {
        msg = sprintf("Data set %s not found.",
            paste(names[!idx], collapse = ", "))
        stop(msg)
    }
    
    x = sapply(names, .load_dataset, simplify = FALSE, USE.NAMES = TRUE)
    res = GenomicRangesList(unlist(x))
    if(merge) {
        datasets = rep(factor(rep(names(res), elementLengths(res))))
        res = unlist(res)
        res$Dataset = datasets
    }

    return(res)
}

sca_metadata <- function() {
    res = .load_dataset("meta")
    return(res)
}
    
  
.load_dataset <- function(name, package = "SomaticCancerAlterations") {
    tmp_env = new.env()
    data(list = name, package = package, envir = tmp_env)
    res = get(name, envir = tmp_env)
    return(res)
}

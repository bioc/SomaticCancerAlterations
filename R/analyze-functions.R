
### SNV density ###
mutation_density <- function(x, bin_size = 1e6,
    chrs = seqnames(x)) {

    if(!missing(chrs))
        x = keepSeqlevels(x, chrs)
    
    seq_gr = as(seqinfo(x), "GRanges")
    bins = subdivideGRanges(seq_gr, bin_size)

    bins$counts = countOverlaps(bins, x)
    bins$density = bins$counts / width(bins) * 1e6
    
    return(bins)
}

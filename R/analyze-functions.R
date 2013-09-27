
### SNV density ###
mutationDensity <- function(x, binSize = 1e6,
    chrs = seqnames(x)) {

    if(!missing(chrs))
        x = keepSeqlevels(x, chrs)
    
    seq_gr = as(seqinfo(x), "GRanges")
    bins = subdivideGRanges(seq_gr, binSize)

    bins$counts = countOverlaps(bins, x)
    bins$density = bins$counts / width(bins) * 1e6
    
    return(bins)
}

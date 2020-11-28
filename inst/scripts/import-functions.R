
### read_maf ###
.read_maf <- function(path) {
    
    maf = read.delim(path, comment.char = "#", sep = "\t", header = TRUE, stringsAsFactors = TRUE)    

    return(maf)
}

### reconstruct_ids ###
.reconstruct_ids <- function(maf) {

    ## Extract the sample ID
    pattern = "(TCGA-.{2}-.{4})-.{3}-.{3}-.{4}-.{2}"
    #patient_ids = strtrim(maf$Tumor_Sample_Barcode, 12)
    maf$Sample_ID = factor(maf$Tumor_Sample_Barcode)
    maf$Patient_ID = factor(gsub(pattern, "\\1", maf$Tumor_Sample_Barcode))
  
    return(maf)
}


### maf2gr ###
.maf2gr <- function(maf, seq_info,
                    keep_columns = c("Hugo_Symbol", "Entrez_Gene_Id", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Verification_Status", "Validation_Status", "Mutation_Status", "Patient_ID", "Sample_ID") ) {
    
    ## The 'MT' is sometimes coded as 'M', rename this
    ## convert to NCBI notation
    maf$Chromosome = sub("^M$", "MT", maf$Chromosome)

    ## contruct GRanges
    maf$Strand = factor(as.character(maf$Strand), levels = c("+", "-", ".", "*"))
    gr = with(maf,
        GRanges(Chromosome, IRanges(Start_position, End_position), seqinfo = seq_info)
        )
    values(gr) = subset(maf, select = keep_columns)

    gr$index = 1:nrow(maf)

    gr = sort(gr) ## no 'unique', will ignore the data columns
    
    return(gr)
}


.meta_maf = function(maf, cancer = NA) {
    cancer_name = c(GBM = "Glioblastoma multiforme", HNSC = "Head and Neck squamous cell carcinoma", KIRC = "Kidney Chromophobe", LUAD = "Lung adenocarcinoma", LUSC = "Lung squamous cell carcinoma", OV = "Ovarian serous cystadenocarcinoma", SKCM = "Skin Cutaneous Melanoma", THCA = "Thyroid carcinoma ")
    strip_columns = c("Center", "NCBI_Build", "Sequence_Source", "Sequencing_Phase", "Sequencer")
    ## extract the meta data that is the same for all columns
    meta = sapply(strip_columns, function(column, data) {
        x = as.character(unique(data[[column]])) ## TODO: ugly, otherwise sapply will return the index of the factor
        if(length(x) != 1)
            warning(column)
        return(x)
        }, maf)
    meta$Number_Samples = nlevels(maf$Sample_ID)
    meta$Number_Patients = nlevels(maf$Patient_ID)
    meta = c(Cancer_Type = cancer, Cancer_Name = cancer_name[cancer], meta)

    return(meta)
}


### validate_maf ###
.validate_maf <- function(maf, seq_info, drop = FALSE) {

    ## a validator function per column
    validators = c(
        Chromosome = function(x) x %in% seqnames(seq_info),
        Strand = function(x) x %in% "+",
        Variant_Classification = function(x) x %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Silent", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation", "3'UTR", "3'Flank", "5'UTR", "5'Flank", "IGR", "Intron", "RNA", "Targeted_Region"),
        Variant_Type = function(x) x %in% c("SNP", "DNP", "TNP", "ONP", "INS", "DEL", "Consolidated"),
        Reference_Allele = function(x) grepl("^[ACGT-]+$", x),
        Tumor_Seq_Allele1 = function(x) grepl("^[ACGT-]+$", x),
        Tumor_Seq_Allele2 = function(x) grepl("^[ACGT-]+$", x),
        Verification_Status = function(x) x %in% c("Verified", "Unknown", NA),
        Validation_Status = function(x) x %in% c("Untested", "Inconclusive", "Valid", "Invalid", NA),
        Sequencer = function(x) x %in% c("Illumina GAIIx", "Illumina HiSeq", "SOLID", "454", "ABI 3730xl", "Ion Torrent PGM", "Ion Torrent Proton", "PacBio RS", "Illumina MiSeq", "Illumina HiSeq 2500", "454 GS FLX Titanium", "AB SOLiD 4 System")
        )

    ## check each column
    idx_valid = logical(nrow(maf))
    idx_valid[] = TRUE
    for(n in names(validators)) {
        status = validators[[n]](maf[[n]])
        if(!all(status)) {
            idx_valid[status] = FALSE
            msg = sprintf("%s: %d invalid entries: '%s'", n, sum(!status), paste(head(as.character(maf[[n]][!status]), 3), collapse = ", "))
            warning(msg)
        }
    }

    ## return value
    if(drop)
        res = maf[idx_valid, ]
    else
        res = all(status)

    return(res)
}

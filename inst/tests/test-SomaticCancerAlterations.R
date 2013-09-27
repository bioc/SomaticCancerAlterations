
library(testthat)
library(SomaticCancerAlterations)

context("Accessor functions")
  
test_that("'scaListDatasets' works", {
      
    expect_is( class(scaListDatasets()), "character" )
    expect_error( scaListDatasets(1) )
    
})

test_that("'scaLoadDatasets' works", {

    all_datasets = scaListDatasets()
    ds = sample(all_datasets, 1)
    ds2 = sample(all_datasets, 2)

    expect_is( scaLoadDatasets(ds), "GenomicRangesList" )
    expect_is( scaLoadDatasets(ds, merge = FALSE), "GenomicRangesList" )
    expect_is( scaLoadDatasets(ds, merge = TRUE), "GRanges" )
    expect_is( scaLoadDatasets(ds2), "GenomicRangesList" )
    expect_is( scaLoadDatasets(ds2, merge = FALSE), "GenomicRangesList" )
    expect_is( scaLoadDatasets(ds2, merge = TRUE), "GRanges" )
    expect_identical( length(scaLoadDatasets()), length(all_datasets) )
      
    expect_error( scaLoadDatasets("not_present") )
    expect_error( scaLoadDatasets(c(all_datasets, "not_present")) )

    x = scaLoadDatasets(all_datasets[1:2], merge = FALSE)
    expect_identical( names(x), all_datasets[1:2] )

    for(ds in all_datasets) {
        x = scaLoadDatasets(ds, merge = TRUE)
        expect_true( length(x) > 0 )
        expect_is( x$Sample_ID, "factor" )
        expect_is( x$Patient_ID, "factor" )
    }

})

test_that("'scaMetadata' works", {
    
    expect_is( scaMetadata(), "data.frame" )
    expect_error( scaMetadata(1) )
    expect_equal( sort(scaListDatasets()), sort(rownames(scaMetadata())) )
    
})

test_that("'mutationDensity' works", {
    
    ds = sample(scaListDatasets(), 1)  
    x = scaLoadDatasets(ds, merge = TRUE)

    expect_is( mutationDensity(x, 1e8), "GRanges" )
    expect_is( mutationDensity(x, binSize = 1e8), "GRanges" )
    
})

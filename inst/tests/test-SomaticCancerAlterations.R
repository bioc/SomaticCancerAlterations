
library(testthat)
library(SomaticCancerAlterations)

context("Accessor functions")
  
test_that("'sca_list_datasets' works", {
      
    expect_is( class(sca_list_datasets()), "character" )
    expect_error( sca_list_datasets(1) )
    
})

test_that("'sca_load_datasets' works", {

    all_datasets = sca_list_datasets()
    ds = sample(all_datasets, 1)
    ds2 = sample(all_datasets, 2)

    expect_is( sca_load_datasets(ds), "GenomicRangesList" )
    expect_is( sca_load_datasets(ds, merge = FALSE), "GenomicRangesList" )
    expect_is( sca_load_datasets(ds, merge = TRUE), "GRanges" )
    expect_is( sca_load_datasets(ds2), "GenomicRangesList" )
    expect_is( sca_load_datasets(ds2, merge = FALSE), "GenomicRangesList" )
    expect_is( sca_load_datasets(ds2, merge = TRUE), "GRanges" )
    expect_identical( length(sca_load_datasets()), length(all_datasets) )
      
    expect_error( sca_load_datasets("not_present") )
    expect_error( sca_load_datasets(c(all_datasets, "not_present")) )

    x = sca_load_datasets(all_datasets[1:2], merge = FALSE)
    expect_identical( names(x), all_datasets[1:2] )

    for(ds in all_datasets) {
        x = sca_load_datasets(ds, merge = TRUE)
        expect_true( length(x) > 0 )
        expect_is( x$Sample_ID, "factor" )
        expect_is( x$Patient_ID, "factor" )
    }

})

test_that("'sca_metadata' works", {
    
    expect_is( sca_metadata(), "data.frame" )
    expect_error( sca_metadata(1) )
    expect_equal( sort(sca_list_datasets()), sort(rownames(sca_metadata())) )
    
})

test_that("'mutation_density' works", {
    
    ds = sample(sca_list_datasets(), 1)  
    x = sca_load_datasets(ds, merge = TRUE)

    expect_is( mutation_density(x, 1e8), "GRanges" )

})

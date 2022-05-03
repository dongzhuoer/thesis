testthat::context('Testing ftp')
setwd(here::here(''))  # workspace is reset per file


testthat::test_that('geo_ftp()', {
    testthat::expect_identical(geo_ftp('foo2018', 'bar'), 'ftp://ftp.ncbi.nlm.nih.gov/geo/bar/foo2nnn/foo2018')
})


testthat::test_that('*_ftp()', {
    testthat::expect_identical(gds_ftp('GDS479'), 'ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDSnnn/GDS479');
    testthat::expect_identical(gse_ftp('GSE3'), 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE3');
    testthat::expect_identical(gpl_ftp('GPL2570'), 'ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL2nnn/GPL2570');
})


testthat::test_that('*_*_ftp()', {
    testthat::expect_identical(
        gpl_soft_ftp('GPL6947'), 
        'ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6947/soft/GPL6947_family.soft.gz'
    );
    testthat::expect_identical(
        gse_soft_ftp('GSE19161'), 
        'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19161/soft/GSE19161_family.soft.gz'
    );
    testthat::expect_identical(
        gse_matrix_ftp('GSE51280'), 
        'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51280/matrix/GSE51280_series_matrix.txt.gz'
    );
})




testthat::context('Testing read_geo.R')
setwd(here::here(''))  # workspace is reset per file


# GEO datasets summary txt file --------------

testthat::test_that('read_summary()', {
	geo_dataset <- read_summary('inst/extdata/gds_result-cut.txt')

	testthat::expect_true(tibble::is_tibble(geo_dataset));
	testthat::expect_identical(dim(geo_dataset), c(4L, 8L));
	testthat::expect_identical(sapply(geo_dataset, class) %T>% {names(.) <- NULL}, rep('character', 8));
	testthat::expect_identical(
	    colnames(geo_dataset), 
	    c("accession", "type", "species", "ftp", "platform", "data_type", "title", "description")
	);
});





# GPL information html file -----------------

testthat::test_that('read_gpl_html()', {
	GPL98 <- read_gpl_html('inst/extdata/GPL98.html');

	testthat::expect_true(tibble::is_tibble(GPL98$info));
	testthat::expect_identical(dim(GPL98$info), c(16L, 2L));
	testthat::expect_identical(colnames(GPL98$info), c("name", "description"));

	testthat::expect_true(tibble::is_tibble(GPL98$sample));
	testthat::expect_identical(dim(GPL98$sample), c(20L, 16L));


	testthat::expect_identical(read_gpl_html('inst/extdata/GPL10400.html'), list(info = NULL, sample = NULL))
});




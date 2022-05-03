context("Testing data")
setwd(here::here(''))  # workspace is reset per file


test_that('elements in map only contain a single value', {
    expect_true(T)
	cat('\tsee lem4\'s test `make_platform_type() sep_pattern`')
})



test_that('map should be umambiguous (the same id cann\'t correspond map to different values)', {
	#" modify it to see whether the test works (e.g. change `3` to `2`)
	testthat::expect_identical(tibble::tibble(c(1, 2, 3), c('A', 'B', 'C'))[[1]] %>% anyDuplicated, 0L);

	testthat::expect_identical(entrez2symbol[[1]] %>% anyDuplicated, 0L);
	testthat::expect_identical(genbank2symbol[[1]] %>% anyDuplicated, 0L);
	testthat::expect_identical(entrez_or_symbol2symbol[[1]] %>% anyDuplicated, 0L);
	testthat::expect_identical(unigene2entrez[[1]] %>% anyDuplicated, 0L);
});



testthat::test_that('map should be effective (no NA)', {
	#" modify it to see whether the test works (e.g. change `2` to `NA`)
	testthat::expect_true(tibble::tibble(c(1, 2), c('A', 'B')) %>% dplyr::filter_all(dplyr::any_vars(is.na(.))) %>% libzhuoer::print_or_T());

	testthat::expect_true(entrez2symbol %>% dplyr::filter_all(dplyr::any_vars(is.na(.))) %>% libzhuoer::print_or_T());
	testthat::expect_true(genbank2symbol %>% dplyr::filter_all(dplyr::any_vars(is.na(.))) %>% libzhuoer::print_or_T());
	testthat::expect_true(entrez_or_symbol2symbol %>% dplyr::filter_all(dplyr::any_vars(is.na(.))) %>% libzhuoer::print_or_T());
	testthat::expect_true(unigene2entrez %>% dplyr::filter_all(dplyr::any_vars(is.na(.))) %>% libzhuoer::print_or_T());
});



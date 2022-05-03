testthat::context("Testing as-symbol.R")
setwd(here::here(''))  # workspace is reset per file

testthat::test_that('as_symbol_from_symbol()', {
	testthat::expect_identical(as_symbol_from_symbol(c('ANKLE2', 'A1BG', '-1')), c('ANKLE2', 'A1BG', NA));
});

testthat::test_that('as_symbol_from_entrez()', {
	testthat::expect_identical(as_symbol_from_entrez(c('503538', '144568', '-1')), c('A1BG-AS1', 'A2ML1', NA));
});

testthat::test_that('as_symbol_from_ensembl()', {
	testthat::expect_identical(as_symbol_from_ensembl(c('ENSG00000148584', 'ENSG00000184389', '-1')), c('A1CF', 'A3GALT2', NA));
});


testthat::test_that('as_symbol_from_entrez_or_symbol()', {
	testthat::expect_identical(as_symbol_from_entrez_or_symbol(c('503538', 'ANKLE2', 'A1BG', '144568', '-1')), c('A1BG-AS1', 'ANKLE2', 'A1BG', 'A2ML1', NA));
});

testthat::test_that('as_symbol_from_genbank()', {
	testthat::expect_identical(as_symbol_from_genbank(c('BX647329', 'AY005822', '-1')), c('A2M', 'ACOT2', NA));
});

testthat::test_that('as_symbol_from_genbank() ignores transcript id', {
	testthat::expect_identical(as_symbol_from_genbank('AK056168.1'), 'SPATA33');
});

testthat::test_that('as_symbol_from_unigene()', {
	testthat::expect_identical(as_symbol_from_unigene(c('41', 'Hs.96', '-1')), c('CEACAM8', 'PMAIP1', NA));
});

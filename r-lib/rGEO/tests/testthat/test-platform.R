testthat::context('Testing platform.R')
setwd(here::here(''))  # workspace is reset per file


# helper function ------------------

fake_platform_type <- . %>% fake_platform(rGEO.data::gpl_metas) %>% guess_platform_type

#" except for sequence, I directly test guess_platform_type() since it both tests filter_*() and the correct order in guess_platform_type()

testthat::test_that('guess_platform_type() NULL', {
	testthat::expect_identical(guess_platform_type(list(accession = 'doubi', info = NULL)), NULL)
    testthat::expect_identical(fake_platform_type('-1'), NULL)
})



## entrez_id ------------------

testthat::test_that('guess_platform_type() entrez_id', {
	testthat::expect_identical(fake_platform_type('GPL1008')$as_symbol_from, 'entrez')
	#1
	testthat::expect_false(identical(fake_platform_type('GPL3639')$as_symbol_from, 'entrez'))
	#2
	testthat::expect_false(identical(fake_platform_type('GPL8152')$as_symbol_from, 'entrez'))
	testthat::expect_false(identical(fake_platform_type('GPL8153')$as_symbol_from, 'entrez'))
	testthat::expect_false(identical(fake_platform_type('GPL9144')$as_symbol_from, 'entrez'))
	#6
    testthat::expect_false(identical(fake_platform_type('GPL18390')$as_symbol_from, 'entrez'))
	#7
    testthat::expect_false(identical(fake_platform_type('GPL3784')$as_symbol_from, 'entrez'))
    testthat::expect_false(identical(fake_platform_type('GPL18060')$as_symbol_from, 'entrez'))
	#8
    testthat::expect_false(identical(fake_platform_type('GPL6978')$as_symbol_from, 'entrez'))
	#9
    testthat::expect_false(identical(fake_platform_type('GPL5677')$as_symbol_from, 'entrez'))
	#10
    testthat::expect_false(identical(fake_platform_type('GPL10327')$as_symbol_from, 'entrez'))
	#11
    testthat::expect_false(identical(fake_platform_type('GPL1228')$as_symbol_from, 'entrez'))
});


# entrez ---------------
testthat::test_that('guess_platform_type() entrez', {
	testthat::expect_identical(fake_platform_type('GPL13667')$as_symbol_from, 'entrez_or_symbol')
	#1
	testthat::expect_false(identical(fake_platform_type('GPL5055')$as_symbol_from, 'entrez_or_symbol'))
	#4
	testthat::expect_false(identical(fake_platform_type('GPL6978')$as_symbol_from, 'entrez_or_symbol'))
	#5
    testthat::expect_false(identical(fake_platform_type('GPL5677')$as_symbol_from, 'entrez_or_symbol'))
	#6
    testthat::expect_false(identical(fake_platform_type('GPL5694')$as_symbol_from, 'entrez_or_symbol'))
    #7
    testthat::expect_false(identical(fake_platform_type('GPL1228')$as_symbol_from, 'entrez_or_symbol'))
});


# ORF ---------------
testthat::test_that('guess_platform_type() ORF', {
	testthat::expect_identical(fake_platform_type('GPL10184')$as_symbol_from, 'entrez_or_symbol')
	#1
	testthat::expect_false(identical(fake_platform_type('GPL6353')$as_symbol_from, 'entrez_or_symbol'))
	#4
	testthat::expect_false(identical(fake_platform_type('GPL23432')$as_symbol_from, 'entrez_or_symbol'))
	testthat::expect_false(identical(fake_platform_type('GPL6254')$as_symbol_from, 'entrez_or_symbol'))
});


# symbol ---------------
testthat::test_that('guess_platform_type() symbol', {
	testthat::expect_identical(fake_platform_type('GPL10030')$as_symbol_from, 'symbol')
	#1
	testthat::expect_false(identical(fake_platform_type('GPL6152')$as_symbol_from, 'symbol'))
	#2
	testthat::expect_false(identical(fake_platform_type('GPL544')$as_symbol_from, 'symbol'))
	testthat::expect_false(identical(fake_platform_type('GPL545')$as_symbol_from, 'symbol'))
	#4
	testthat::expect_false(identical(fake_platform_type('GPL8798')$measure, 'PlateformeSymbole'))
});


# ensembl ---------------
testthat::test_that('guess_platform_type() ensembl', {
	testthat::expect_identical(fake_platform_type('GPL13916')$as_symbol_from, 'ensembl')
	#3
	testthat::expect_false(identical(fake_platform_type('GPL5678')$as_symbol_from, 'ensembl'))
	#4
	testthat::expect_false(identical(fake_platform_type('GPL4253')$as_symbol_from, 'ensembl'))
});


# refseq ---------------
testthat::test_that('guess_platform_type() refseq', {
	testthat::expect_identical(fake_platform_type('GPL10097')$as_symbol_from, 'genbank')
	#3
	testthat::expect_false(identical(fake_platform_type('GPL19407')$measure, 'RefSeq Single exon'))
	#4
	testthat::expect_false(identical(fake_platform_type('GPL13645')$as_symbol_from, 'genbank'))
});


# genbank ---------------
testthat::test_that('guess_platform_type() genbank', {
	testthat::expect_identical(fake_platform_type('GPL10001')$as_symbol_from, 'genbank')
});


# unigene ---------------
testthat::test_that('guess_platform_type() unigene', {
	testthat::expect_identical(fake_platform_type('GPL1052')$as_symbol_from, 'unigene')
});


# sequence ---------------
testthat::test_that('guess_platform_type() sequence', {
	info10118 <- 'GPL10118' %>% {tibble::add_column(rGEO.data::gpl_metas[[.]]$info, accession = .)}
	testthat::expect_identical(info10118 %>% filter_sequence %>% nrow, 1L)
});



# make_as_symbol ---------------

testthat::test_that('make_as_symbol()', {
	testthat::expect_identical(make_as_symbol('symbol'), hgnc::as_symbol_from_symbol);
	testthat::expect_identical(make_as_symbol('entrez'), hgnc::as_symbol_from_entrez);
	testthat::expect_identical(make_as_symbol('ensembl'), hgnc::as_symbol_from_ensembl);
	testthat::expect_identical(make_as_symbol('genbank'), hgnc::as_symbol_from_genbank);
	testthat::expect_identical(make_as_symbol('entrez_or_symbol'), hgnc::as_symbol_from_entrez_or_symbol);
	testthat::expect_identical(make_as_symbol('unigene'), hgnc::as_symbol_from_unigene);
});

# make_platform_type ---------------

testthat::test_that('make_platform_type()', {
	testthat::expect_identical(
		make_platform_type('Gene Symbol', 'symbol'),
		list(measure = 'Gene Symbol',  sep_pattern = '[^-@\\w]+', as_symbol_from = 'symbol')
	);
	testthat::expect_identical(
		make_platform_type('ENTREZ_GENE_ID', 'entrez'),
		list(measure = 'ENTREZ_GENE_ID',  sep_pattern = '[^\\d]+', as_symbol_from = 'entrez')
	);
	testthat::expect_identical(
		make_platform_type('Ensembl Gene ID', 'ensembl'),
		list(measure = 'Ensembl Gene ID',  sep_pattern = '[^ENSG\\d]+', as_symbol_from = 'ensembl')
	);
	testthat::expect_identical(
		make_platform_type('GB_ACC', 'genbank'),
		list(measure = 'GB_ACC',  sep_pattern = '[^\\w\\.]+', as_symbol_from = 'genbank')
	);
	testthat::expect_identical(
		make_platform_type('Entrez Gene', 'entrez_or_symbol'),
		list(measure = 'Entrez Gene',  sep_pattern = '[^-@\\w]+', as_symbol_from = 'entrez_or_symbol')
	);
	testthat::expect_identical(
		make_platform_type('Unigene ID', 'unigene'),
		list(measure = 'Unigene ID',  sep_pattern = '[^Hs\\.\\d]+', as_symbol_from = 'unigene')
	);
});


testthat::test_that('make_platform_type() sep_pattern', {
	#" demonstrate the effectiveness of pattern `[^-@\w]+`, omitting `[A-Z]` & `\d` since too obvious
	testthat::expect_output(hgnc::hugo_symbol %>% stringr::str_subset('[^-\\w]+') %>% print, '@');
	testthat::expect_output(hgnc::hugo_symbol %>% stringr::str_subset('[^@\\w]+') %>% print, '-');
	testthat::expect_output(hgnc::hugo_symbol %>% stringr::str_subset('[^@\\dA-Z]+') %>% print, '[a-z]');
	testthat::expect_output(hgnc::hugo_symbol %>% stringr::str_subset('[^@a-zA-Z\\d]+') %>% print, '_');


	testthat::expect_identical(stringr::str_subset(hgnc::hugo_symbol, '[^-@\\w]+'), character(0));
	testthat::expect_identical(stringr::str_subset(hgnc::entrez2symbol[[1]], '[^\\d]+'), character(0));
	testthat::expect_identical(stringr::str_subset(hgnc::ensembl2symbol[[1]], '[^ENSG\\d]+'), character(0));
	#" input may contain trailing `.\d`, check Caution in hgnc's readme for details
	testthat::expect_identical(stringr::str_subset(hgnc::genbank2symbol[[1]], '[^\\w\\.]+'), character(0));
	testthat::expect_identical(stringr::str_subset(hgnc::entrez_or_symbol2symbol[[1]], '[^-@\\w]+'), character(0));
	testthat::expect_identical(stringr::str_subset(hgnc::unigene2entrez[[1]], '[^Hs\\.\\d]+'), character(0));
});



testthat::context('Testing read_geo.R')
setwd(here::here(''))  # workspace is reset per file



# matrix --------------------
testthat::test_that('read_gse_matrix() non existent file', {
	testthat::expect_identical(read_gse_matrix('non/existent/file'), NULL)
})

testthat::test_that('read_gse_matrix() GSE51280', {
	matrix <- read_gse_matrix('inst/extdata/GSE51280_series_matrix.txt.gz')

	testthat::expect_true(tibble::is_tibble(matrix))
	testthat::expect_identical(sapply(matrix, class) %>% {names(.) <- NULL; .}, c('character', rep('numeric', 24)))
 	testthat::expect_identical(colnames(matrix)[1], 'ID_REF')
 	testthat::expect_identical(dim(matrix), c(123L, 25L))
})


testthat::test_that('read_gse_matrix() GSE19161', {
	# some duplicated ones are omitted
	matrix <- read_gse_matrix('inst/extdata/GSE19161_series_matrix.txt.gz')

	testthat::expect_identical(sapply(matrix, class) %>% {names(.) <- NULL; .}, c('character', rep('numeric', 61)))
 	testthat::expect_identical(dim(matrix), c(658L, 62L))
})



# soft ----------------------


testthat::test_that('parse_gse_soft() platform with duplicated column name', {
	# use first one
	GPL4133 <- parse_gse_soft('inst/extdata/GSE98737_truncated_family.soft.gz', F)
	testthat::expect_false(GPL4133$info$name %>% duplicated %>% any)
	testthat::expect_identical(GPL4133$table$SPOT_ID[[1]], 'GE_BrightCorner')

	# special case use second one
	GPL7025 <- parse_gse_soft('inst/extdata/GSE12080_family.soft.gz', F)
	testthat::expect_false(GPL7025$info$name %>% duplicated %>% any)
	testthat::expect_identical(GPL7025$table$GB_ACC[[2]], 'NM_006923')

})


testthat::test_that('parse_gse_soft() GSE51280', {
	testthat::expect_output(parse_gse_soft('inst/extdata/GSE51280_family.soft.gz'))

	platform <- parse_gse_soft('inst/extdata/GSE51280_family.soft.gz', F)
	testthat::expect_identical(platform$accession, 'GPL17590')
	testthat::expect_true(tibble::is_tibble(platform$info))
	testthat::expect_identical(colnames(platform$info), c('name', 'description'))
	testthat::expect_identical(platform$info$name, c('ID', 'GENE_SYMBOL', 'GB_ACC', 'Class Name', 'SPOT_ID'))
})


testthat::test_that('parse_gse_soft() GSE19161', {
	# some duplicated ones are omitted

	platform <- parse_gse_soft('inst/extdata/GSE19161_family.soft.gz', F)
	testthat::expect_identical(platform$accession, 'GPL9717')
	testthat::expect_identical(platform$info$name, c('ID', 'Species Scientific Name', 'Sequence Type', 'Sequence Source', 'Transcript ID(Array Design)', 'UniGene ID', 'Gene Title', 'Gene Symbol', 'Ensembl', 'ORF', 'SPOT_ID'))
})


testthat::test_that('read_gse_soft()', {
	# non existent file
	testthat::expect_identical(read_gse_soft('non/existent/file'), NULL)

	# verbose
	testthat::expect_output(read_gse_soft('inst/extdata/GSE51280_family.soft.gz', T))

	# use "GENE_SYMBOL" as symbol
	gse51280_chip <- read_gse_soft('inst/extdata/GSE51280_family.soft.gz')
	testthat::expect_true(tibble::is_tibble(gse51280_chip))
	testthat::expect_identical(colnames(gse51280_chip), c('ID_REF', 'symbol'))
	testthat::expect_true(is.element(gse51280_chip$symbol, hgnc::hugo_symbol) %>% all)

	# use "ORF" as entrez_or_symbol
	gse19161_chip <- read_gse_soft('inst/extdata/GSE19161_family.soft.gz')
	testthat::expect_true(is.element(gse19161_chip$symbol, hgnc::hugo_symbol) %>% all)
})








# platform type ---------------------------------






# testthat::test_that('read_gse() using symbol', {
# 	if (Sys.getenv('TRAVIS') != 'true') skip('skip long test (4s) when not on Travis')
#
# 	gse16571 <- read_gse('inst/extdata/GSE16571_family.soft.gz', 'inst/extdata/GSE16571_series_matrix.txt.gz')
#
# 	testthat::expect_identical(sapply(gse16571, class) %>% {names(.) <- NULL; .}, c(rep('character', 2), rep('numeric', 6)))
# 	testthat::expect_identical(colnames(gse16571)[1:2], c('ID_REF', 'symbol'))
# 	testthat::expect_identical(dim(gse16571), c(20448L, 8L))
# 	testthat::expect_true(is.element(gse16571$symbol, hgnc::hugo_symbol) %>% all)
# });

# entrez GPL6947 (such as GSE16571)




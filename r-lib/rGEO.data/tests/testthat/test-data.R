testthat::context('Testing data')
setwd(here::here(''))  # workspace is reset per file

testthat::test_that('dataset', {
    testthat::expect_true(tibble::is_tibble(dataset));
    testthat::expect_identical(ncol(dataset), 8L);

});

testthat::test_that('platform', {
    testthat::expect_true(tibble::is_tibble(platform));
    testthat::expect_identical(ncol(platform), 9L);
});

testthat::test_that('GPL', {
    testthat::expect_true(
        dataset %>% dplyr::filter(species == 'Homo sapiens', stringr::str_detect(accession, 'GSE')) %>% 
            dplyr::filter(stringr::str_count(platform, 'GPL') == 1L) %>% 
            {stringr::str_extract(.$platform, 'GPL\\d+')} %>% 
            unique %>% setdiff(platform$Accession) %>% libzhuoer::print_or_T(), 
        info = 'not all GPL in dataset are included in platform'
    );
});

testthat::test_that('gpl_metas', {
    testthat::expect_true(is.list(gpl_metas));

    gpl_meta <- sample(gpl_metas, 1) %>% {.[[1]]};
    testthat::expect_true(is.list(gpl_meta));
    testthat::expect_identical(length(gpl_meta), 2L);
    testthat::expect_true(tibble::is_tibble(gpl_meta[[1]]) || is.null(gpl_meta[[1]]));
    testthat::expect_true(tibble::is_tibble(gpl_meta[[2]]) || is.null(gpl_meta[[2]]));
});

context("Testing mapping.R")
setwd(here::here(''))  # workspace is reset per file

# test data --------------------

wide_map <- tibble::tribble(
    ~symbol, ~id,
    'A', '1|2',
    'B', '2',
    'C', '3'
)

long_map <- tibble::tribble(
    ~symbol, ~id,
    'A', '1',
    'A', '2',
    'B', '2',
    'C', '3'
)

square_map <- tibble::tibble(symbol = letters, id = as.character(1:26))


# melt_map -----------------

test_that("melt_map()", {
	expect_identical(melt_map(wide_map, 'symbol', 'id', '\\|'), long_map);
});

test_that("melt_map() doesn\'t transfigure square_map", {
	expect_identical(melt_map(square_map, 'symbol', 'id'), square_map);
});



test_that("melt_map() removes NA", {
	expect_identical(
		melt_map(tibble::add_row(wide_map, 'id' = '4', .before = 2), 'symbol', 'id', '\\|'),
		long_map
	);
	expect_identical(
		melt_map(tibble::add_row(wide_map, 'symbol' = 'D', .before = 3), 'symbol', 'id', '\\|'),
		long_map
	);
	expect_identical(
		melt_map(tibble::add_row(wide_map, .before = 4), 'symbol', 'id', '\\|'),
		long_map
	);
});

test_that("melt_map() removes NA for square_map", {
	expect_identical(
		melt_map(tibble::add_row(square_map, 'id' = '4', .before = 10), 'symbol', 'id'),
		square_map
	);
	expect_identical(
		melt_map(tibble::add_row(square_map, 'symbol' = 'D', .before = 15), 'symbol', 'id'),
		square_map
	);
	expect_identical(
		melt_map(tibble::add_row(square_map, .before = 20), 'symbol', 'id'),
		square_map
	);
});



test_that("melt_map() drops extra variables", {
	expect_identical(
		melt_map(tibble::add_column(wide_map, foo = 'bar'), 'symbol', 'id', '\\|'),
		long_map
	)
});

test_that("melt_map() drops extra variables for square_map", {
	expect_identical(
		melt_map(tibble::add_column(square_map, foo = 'bar'), 'symbol', 'id'),
		square_map
	)
});





# cast_map -----------------

test_that("cast_map()", {
	expect_identical(cast_map(long_map, 'symbol', 'id', '|'), wide_map)
});

test_that("cast_map() doesn\'t transfigure square_map", {
	expect_identical(cast_map(square_map, 'symbol', 'id'), square_map);
});



test_that("cast_map() removes NA", {
	expect_identical(
		cast_map(tibble::add_row(long_map, 'id' = '4', .before = 2), 'symbol', 'id', '|'),
		wide_map
	);
	expect_identical(
		cast_map(tibble::add_row(long_map, 'symbol' = 'D', .before = 3), 'symbol', 'id', '|'),
		wide_map
	);
	expect_identical(
		cast_map(tibble::add_row(long_map, .before = 4), 'symbol', 'id', '|'),
		wide_map
	);
});

test_that("cast_map() removes NA for square_map", {
	expect_identical(
		cast_map(tibble::add_row(square_map, 'id' = '4', .before = 2), 'symbol', 'id'),
		square_map
	);
	expect_identical(
		cast_map(tibble::add_row(square_map, 'symbol' = 'D', .before = 3), 'symbol', 'id'),
		square_map
	);
	expect_identical(
		cast_map(tibble::add_row(square_map, .before = 4), 'symbol', 'id'),
		square_map
	);
});



test_that("cast_map() drops extra variables", {
	expect_identical(
		cast_map(tibble::add_column(long_map, foo = 'bar'), 'symbol', 'id', '|'),
		wide_map
	)
});

test_that("cast_map() drops extra variables for square_map", {
	expect_identical(
		cast_map(tibble::add_column(square_map, foo = 'bar'), 'symbol', 'id'),
		square_map
	)
});






# qualify_map -------

test_that("qualify_map() drops any NA", {
	expect_identical(
		tibble::tribble(
			~from, ~to,
			NA, 1,
			'B', NA,
			NA, NA,
			'C', 3
		) %>% qualify_map,
		tibble::tibble(from = 'C', to = 3)
	)
});



test_that("qualify_map() remove ambiguous mapping", {
	expect_identical(
		tibble::tribble(
			~from, ~to,
			'A', 0.5,
			'A', 1.5,
			'B', 2
		) %>% qualify_map,
		tibble::tibble(from = 'B', to = 2)
	)
});














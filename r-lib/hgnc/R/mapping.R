#' @importFrom magrittr %<>%
#'
NULL


#' @title convert a "wide" map into a "long" one
#'
#' @param df tibble. rows containing NA in any of following variables are removed.
#' @param id string. name of a character variable in `df`, its elements must be unique and can only contain _one_ value each
#' @param measure string. name of a character variable in `df`, its elements needn't to be unique and contain multiple values (we use [stringr::str_split()] to sperate them by `sep_pattern`)
#' @param sep_pattern string. see `measure`, note that it is a _pattern_, for example, '|' should be '\\|'
#'
#' @return tibble. column names are `id` and `measure`.
#' @export
#'
#' @family mapping
#'
#' @examples
#' tibble::tribble(
#'     ~symbol, ~id, ~foobar,
#'     'A|B', '1', 'foo',
#'     'C', '2', 'bar',
#'     'D', '2', 'foobar'
#' ) %>% melt_map('id', 'symbol', '\\|')
#'
melt_map <- function(df, id, measure, sep_pattern = NA) {
	df %<>% dplyr::select(!!id, !!measure) %>% dplyr::filter_all(dplyr::all_vars(!is.na(.)));
	one2many <- sep_pattern %>% {if (is.na(.)) rep(F, nrow(df)) else stringr::str_detect(df[[measure]], .)} ;
    if (all(!one2many)) return(df)

	dplyr::bind_rows(
		dplyr::filter(df, one2many) %>% plyr::dlply(id, . %>% {
			measures <- stringr::str_split(.[[measure]], sep_pattern)[[1]];
			dplyr::tibble(!!id := .[[id]], !!measure := measures)
		}) %>% dplyr::bind_rows(),
		dplyr::filter(df, !one2many)
	)
}


#' @title convert a "long" map into a "wide" one
#'
#' @param df tibble. rows containing NA in any of following variables are removed.
#' @param id string. name of a character variable in `df`, its elements can only contain _one_ value each
#' @param measure string. name of a character variable in `df`, its elements can only contain _one_ value each
#' @param collapse string. used as a separator when collapsing multiple values of `measure` which map to a single value of `id`
#'
#' @return tibble. column names are `id` and `measure`.
#'
#' @family mapping
#'
#' @export
cast_map <- function(df, id, measure, collapse = ' /// ') {
	df %>% dplyr::filter_all(dplyr::all_vars(!is.na(.))) %>%
		dplyr::group_by(!!rlang::sym(id)) %>%
		dplyr::summarise(
			!!measure := paste0(!!rlang::sym(measure), collapse = collapse)
		) %>% dplyr::select(!!id, !!measure)
}





#' @title qualify a map to be "square"
#' 
#' @description qualify a map to be umambiguous and effective (no NA)
#'
#' @param df tibble.
#'
#' @return tibble.
#'
#' @examples
#' tibble::tribble(
#'     ~from, ~to,
#'     NA, 1,
#'     'B', NA,
#'     NA, NA,
#'     'C', 3
#' ) %>% qualify_map
#'
#' tibble::tribble(
#'     ~from, ~to,
#'     'A', 0.5,
#'     'A', 1.5,
#'     'B', 2
#' ) %>% qualify_map
#'
#' @export
qualify_map <- function(df) {
	df %>% dplyr::filter_all(dplyr::all_vars(!is.na(.))) %>%
		dplyr::filter_at(1, dplyr::any_vars(!(. %in% .[duplicated(.)])))
}




#' @title remove problematic rows after you [readr::read_*()][readr::read_delim()] a file
#'
#' @param data tibble
#'
#' @return tibble
#'
#' @examples
#' problematic_df <- readr::read_tsv('a\tb\tc\n1\n2\t3\n1\t2\t3')
#'
#' problematic_df
#'
#' rm_problematic_row(problematic_df)
#'
#' @export
rm_problematic_row <- function(data) {
	problematic_row <- unique(readr::problems(data)$row);
	if (length(problematic_row) > 0L) data %<>% dplyr::slice(-problematic_row);
	data %T>% {attr(., c('problems')) <- NULL}
}



# matrix --------------------

#' @title read GSE series matrix file (contains expression data)
#'
#' @details For now, we assume that the only special value in input file is '' (empty, you may search `\\t\\t`), i.e., no `Inf`, `NaN`, `NA`, etc. And we don't collect sample meta data.
#'
#' @param path string. path to the matrix file
#'
#' @return tibble or NULL. the first variable is `ID_REF` (probe ID), others are gene expression value of each sample
#'
#' @examples
#' read_gse_matrix(system.file('extdata/GSE51280_series_matrix.txt.gz', package = 'rGEO'))
#'
#' @family read raw data
#' 
#' @export

# matrix_file <- 'inst/extdata/GSE51280_series_matrix.txt.gz'
read_gse_matrix <- function(matrix_file) {
	if (!file.exists(matrix_file)) return(NULL)
	matrix_raw <- readr::read_lines(matrix_file);

    matrix_boundary <- stringr::str_which(matrix_raw, 'series_matrix_table');
	if (diff(matrix_boundary) == 2L) return(NULL); # only a header line, no data to read
    
	matrix_raw[seq(matrix_boundary[1] + 1, matrix_boundary[2] - 1)] %>%
    	I() %>% readr::read_tsv(T, readr::cols(.default = 'c')) %>% 
        rm_problematic_row() %>% 
        {suppressWarnings(dplyr::mutate_at(., -1, as.double))}
}
#" maybe used to parse meta data of samples
# meta <- series %>% stringr::str_subset('^!Sample_') %>% stringr::str_replace('^!', '') %>%
# 	stringr::str_split('\t') %>% plyr::laply(as.character) %>% t %>%
# 	plyr::aaply(1, . %>% paste0(collapse = '\t')) %>% {suppressWarnings(read_char_tsv(.))}




# SOFT --------------------

#' @title parse GSE SOFT file (contains platform annotation)
#'
#' @details for now, we drop probes which haven't been mapped to a symbool (mapping to multiple symbols is okay).
#'
#' @param path string. path to the SOFT file
#'
#' @return list
#'
#' 1. accession: string. platform accession
#'
#' 1. info: tibble. chip column name information
#'
#' 1. table: tibble. chip annotation
#'
#' @examples
#' parse_gse_soft(system.file('extdata/GSE19161_family.soft.gz', package = 'rGEO'), verbose = F)
#' 
#' parse_gse_soft(system.file('extdata/GSE51280_family.soft.gz', package = 'rGEO'))
#'
#' @family read raw data
#' 
#' @export
#' 
# soft_file <- 'inst/extdata/GSE12080_family.soft.gz'; verbose = T
parse_gse_soft <- function(soft_file, verbose = T) {
	soft <- readr::read_lines(soft_file);
    accession <- stringr::str_subset(soft, '^\\^PLATFORM') %>% stringr::str_extract('GPL\\d+');
    table_boundary <- stringr::str_which(soft, 'platform_table');

    platform_begin <- stringr::str_which(soft, '^\\^PLATFORM'); # platform section begin
    info <- soft[platform_begin:table_boundary[1]] %>% stringr::str_subset('^#') %>% {
    	if (verbose) {cat('\nplatform meta:\n'); print(.); cat('\n\n\n');}
    	tibble::tibble(
			name = stringr::str_extract(., '(?<=^#)[^=]+(?= =)'),
			description = stringr::str_remove(., '^#[^=]+= +')
		)
    };

	duplicated_name <- info$name %>% {.[duplicated(.)]} %>% unique;
	if (length(duplicated_name) > 0) {
        if (accession %in% dup_use_latter) {
		    col_drop <- sapply(duplicated_name, . %>% {which(info$name == .)[1]})
        } else
			col_drop <- info$name %>% duplicated() %>% which
    	col_types <- rep('c', nrow(info)) %T>% {.[col_drop] = '-'} %>% paste0(collapse = '')
    	info %<>% dplyr::slice(-col_drop)
    } else
    	col_types <- readr::cols(.default = readr::col_character())

    table <- soft[seq(table_boundary[1] + 1, table_boundary[2] - 1)] %>% paste0(collapse = '\n') %>%
    	{suppressWarnings(readr::read_tsv(., T, col_types, name_repair = 'minimal'))} %>% rm_problematic_row()
    aberrant_row <- table[[1]] %>% stringr::str_length() %>% {. > mean(.) * 10} %>% which
    if (length(aberrant_row) > 0L) table %<>% dplyr::slice(-aberrant_row);

	if (length(duplicated_name) > 0)
    	for (i in info$name %in% duplicated_name %>% which) colnames(table)[i] %<>% stringr::str_remove('_\\d$')
    if (!identical(colnames(table), info$name)) stop('parsing platform info failed')
    
    list(accession = accession, info = info, table = table);
}



#' @title read GSE SOFT file (contains platform annotation)
#'
#' @details for now, we drop probes which haven't been mapped to a symbool (mapping multiple symbols is okay)
#'
#' @param path string. path to the SOFT file
#'
#' @return tibble or NULL. the first variable is `ID_REF` (probe ID) and the second one is [HUGO](https://www.genenames.org/) gene symbol
#'
#' @examples
#' read_gse_soft(system.file('extdata/GSE19161_family.soft.gz', package = 'rGEO'))
#' 
#' read_gse_soft(system.file('extdata/GSE51280_family.soft.gz', package = 'rGEO'), verbose = T)
#'
#' @family read raw data
#' 
#' @export

# soft_file <- 'inst/extdata/GSE19161_family.soft.gz'; verbose = T
read_gse_soft <- function(soft_file, verbose = F) {
	if (!file.exists(soft_file)) return(NULL);
	platform <- parse_gse_soft(soft_file, verbose);
    if (verbose) {cat('platform:\n\n'); print(platform)}

	platform_type <- guess_platform_type(platform);
	if (is.null(platform_type)) return(NULL);

	id <- names(platform$table)[1];
	measure = platform_type$measure;
	sep_pattern = platform_type$sep_pattern;
	as_symbol = make_as_symbol(platform_type$as_symbol_from)
	if (verbose) cat('\n\n', platform$accession, ': use "', platform_type$measure, '" as ', platform_type$as_symbol_from, '\n\n\n\n', sep = '')

	platform$table %>%
        hgnc::melt_map(id, measure, sep_pattern) %>%
        dplyr::mutate('symbol' = as_symbol(!!rlang::sym(measure))) %>%
        hgnc::cast_map(id, 'symbol') %>%
		dplyr::rename('ID_REF' = !!rlang::sym(id))
}




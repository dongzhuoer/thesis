
#' @title read GEO datasets summary .txt file
#' 
#' @description refer to [here](https://github.com/dongzhuoer/thesis/r-lib/rGEO.data/blob/master/R-raw/data.Rmd#dataset) for how to get the input file
#' 
#' @param path string. path to the summary .txt file
#' 
#' @return tibble
#' 
#' 1. accession
#' 
#' 1. type
#' 
#' 1. species
#' 
#' 1. ftp
#' 
#' 1. platform
#' 
#' 1. data_type
#' 
#' 1. title
#' 
#' 1. description
#' 
#' @examples 
#' read_summary(system.file('extdata/gds_result-cut.txt', package = 'rGEO.data'))
#'
#' @export 
read_summary <- function(path) {
	content_l <- readr::read_file(path) %>% stringr::str_replace('^\n', '') %>% 
		{stringr::str_split(., '\n\n')[[1]]}  %>% lapply(readr::read_lines)
	
	title <- content_l %>% sapply(. %>% {.[1]}) %>% stringr::str_replace('\\d+\\. ', '')
	description <- content_l %>% sapply(. %>% {.[2]})
	species <- content_l %>% sapply(. %>% stringr::str_subset('^Organism')) %>% stringr::str_replace('Organism:\t', '')
	data_type <- content_l %>% sapply(. %>% stringr::str_subset('^Type')) %>% stringr::str_replace('Type:\t', '')
	platform <- content_l %>% sapply(. %>% stringr::str_subset('Platform[^\n]+Sample') %>% {if (length(.) == 0L) '' else .})
	ftp <- content_l %>% sapply(. %>% stringr::str_subset('^FTP download')) %>% stringr::str_extract('ftp[^\n]+')
	type <- content_l %>% sapply(. %>% stringr::str_subset('Accession: ')) %>% stringr::str_extract('^\\w+')
	accession <- content_l %>% sapply(. %>% stringr::str_subset('Accession: ')) %>% stringr::str_extract('G\\w+')
	
	tibble::tibble(accession, type, species, ftp, platform, data_type, title, description)
}



#' @title read GPL information webpage
#'
#' @param webpage passed on to [xml2::read_html()].  
#'
#' @return list
#'
#' 1. info: tibble (`name`, `description`)
#'
#' 1. sample: tibble, part of the platform table
#' @export
#'
#' @examples 
#' read_gpl_html(system.file('extdata/GPL98.html', package = 'rGEO.data'))
#' read_gpl_html(system.file('extdata/GPL10400.html', package = 'rGEO.data'))
#' 
#' \dontrun{
#'     # feel free to run it, I don't do so to avoid accessing the Internet
#'     read_gpl_html('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16570')
#' }

# read_gpl_html('inst/extdata/GPL98.html')
# xml2::read_html('inst/extdata/GPL98.html') %>% rvest::html_node('table[cellspacing="2"]')
read_gpl_html <- function(webpage) {
	html <- xml2::read_html(webpage);

	info_node <- html %>% rvest::html_node('table[cellspacing="3"]');
	if (rvest::html_text(rvest::html_node(info_node, 'tr strong')) == 'Data table header descriptions') {
		info <- info_node %>% rvest::html_table() %>% tibble::as_tibble() %>% dplyr::slice(-1) %>%
			dplyr::select(name = 1, description = 2);
	} else {
		info <- NULL
	}
	sample_node <- html %>% rvest::html_node('table[cellspacing="2"]')
	if (class(sample_node) == 'xml_missing') {
		sample <- NULL
	} else {
	sample <- sample_node %>% rvest::html_table(fill = T) %>% tibble::as_tibble() %T>%
		{colnames(.) = .[2, , drop = T]} %>% dplyr::slice(-1:-2)
	}

	list(info = info, sample = sample)
}

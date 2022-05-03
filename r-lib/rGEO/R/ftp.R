# helper fucntion -----------

geo_ftp <- function(accession, base_dir) {
	accession %>% {paste0(stringr::str_replace(., '\\d{1,3}$', 'nnn/'), .)} %>%
		paste('ftp://ftp.ncbi.nlm.nih.gov/geo', base_dir, ., sep = '/');
}



# ftp directory -----------

#' @title ftp path to the directory of GEO accession (GDS, GSE or GPL)
#'
#' @param accession character. like 'GDS479', 'GSE3', 'GPL2570'
#'
#' @return character. ftp path to directory,  without trailing `/`.
#' 
#' @examples
#' gds_ftp('GDS479')
#' 
#' gse_ftp('GSE3')
#' 
#' gpl_ftp('GPL2570')
#' 
#' @export

# gds_ftp('GDS479') %>% browseURL()
gds_ftp <- function(accession) {
	geo_ftp(accession, 'datasets');
}



#' @rdname gds_ftp
#' @export

# gse_ftp('GSE3') %>% browseURL()
gse_ftp <- function(accession) {
	geo_ftp(accession, 'series');
}



#' @rdname gds_ftp
#' @export

# gpl_ftp('GPL2570') %>% browseURL()
gpl_ftp <- function(accession) {
	geo_ftp(accession, 'platforms');
}



# ftp file -----------


#' @title ftp path to GEO raw data file
#' 
#' @description `gpl_soft_ftp()` for GPL SOFT file, `gse_soft_ftp()` for GSE SOFT file, `gse_matrix_ftp()` for GSE matrix file
#' 
#' @param accession character. like 'GPL6947', 'GSE19161'
#'
#' @return character. ftp path to file
#' 
#' @examples
#' gpl_soft_ftp('GPL6947')
#' 
#' gse_soft_ftp('GSE19161')
#' 
#' gse_matrix_ftp('GSE19161')
#'
#' @export
gpl_soft_ftp <- function(accession) {
	paste0(gpl_ftp(accession), '/soft/', accession, '_family.soft.gz')
}

#' @rdname gpl_soft_ftp
#' @export
gse_soft_ftp <- function(accession) {
	paste0(gse_ftp(accession), '/soft/', accession, '_family.soft.gz')
}

#' @rdname gpl_soft_ftp
#' @export
gse_matrix_ftp <- function(accession) {
	paste0(gse_ftp(accession), '/matrix/', accession, '_series_matrix.txt.gz')
}


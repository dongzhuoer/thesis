

map_1_to_2 <- function(id, map) {
	map[[2]][match(id, map[[1]])]
}


#' @name as_symbol
#' 
#' @title convert common gene ID (accession) to HUGO gene symbol
#' 
#' @description haha  
#' 
#' @param symbol character. each element is a valid HUGO gene symbol, like `'ANKLE2'`
#' @param entrez character. each element is a valid entrez gene ID, like `'503538'`
#' @param ensembl character. each element is a valid ensembl gene ID, like `'ENSG00000148584'`
#' @param genbank character. each element is a valid Genbank accession, like `'BX647329'`
#' @param entrez_or_symbol character. each element is a valid entrez gene ID or HUGO gene symbol, refer to `vignette('content')` 
#' @param unigene character. each element is a valid unigene ID, with or without `Hs.`, like `'41'`, `'Hs.96'`
#' 
#' @return character. corresponding HUGO gene symbol
NULL

#' @rdname as_symbol
#' 
#' @export   
as_symbol_from_symbol <- function(symbol) {
	hgnc::hugo_symbol[match(symbol, hgnc::hugo_symbol)]
}

#' @rdname as_symbol
#'
#' @export
as_symbol_from_entrez <- function(entrez) {
	map_1_to_2(entrez, hgnc::entrez2symbol)
}

#' @rdname as_symbol
#'
#' @export
as_symbol_from_ensembl <- function(ensembl) {
	map_1_to_2(ensembl, hgnc::ensembl2symbol)
}

#' @section Caution: `as_symbol_from_genbank()` ignores transcript id or version
#'   number, anyhow, the trailing `.\\d`.
#' @rdname as_symbol
#'
#' @export
as_symbol_from_genbank <- function(genbank) {
	genbank %<>% stringr::str_remove('\\.\\d+$');

	map_1_to_2(genbank, hgnc::genbank2symbol)
}

#' @rdname as_symbol
#'
#' @export
as_symbol_from_entrez_or_symbol <- function(entrez_or_symbol) {
	map_1_to_2(entrez_or_symbol, hgnc::entrez_or_symbol2symbol)
}

#' @rdname as_symbol
#'
#' @export
as_symbol_from_unigene <- function(unigene) {
	unigene %<>% stringr::str_replace('^\\d', 'Hs.\\0'); #" add Hs. if necessary

	map_1_to_2(unigene, hgnc::unigene2entrez) %>% as_symbol_from_entrez
}

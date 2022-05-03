#  duplicated column with same name is table
dup_use_latter <- c('GPL13224', 'GPL13915', 'GPL7025', 'GPL7088', 'GPL8335')



# special ---------------------
special <- tibble::tribble(
	~accession, ~name, ~database,
	'GPL9851', 'Primary_Gene_ID', 'entrez',  # there are even numbers in `Gene_Symbol`, for example, 38961 mapped to FlyBase, but corresponding 23157 from `Primary_Gene_ID` is indeed 'septin 6'
	'GPL4757', 'GB_ACC (ref)', 'genbank',  # it fills `symbol` filed with name
	'GPL15684', 'gene_id', 'ensembl'
)

# not dealt with ---------------------
filter_non_human <- function(info_all) {
	info_all %>%
    dplyr::filter(accession %in% c('GPL1883'))  	# Yeast clone vector
}

filter_miRNA <- function(info_all) {
	info_all %>%
    dplyr::filter(stringr::str_detect(full, stringr::regex('miRNA|microRNA', T)) |
    			  accession %in% c('GPL6127', 'GPL9349', 'GPL6513', 'GPL3241', 'GPL17546', 'GPL6743')
    )
}

filter_circRNA <- function(info_all) {
	info_all %>%
    dplyr::filter(stringr::str_detect(full, stringr::fixed('circRNA', T)))
}

# current rules -----------------
filter_entrez_id <- function(info_all) {
	info_all %>%
    dplyr::filter(name != 'GB_ACC') %>%
    dplyr::filter(!(accession %in% c('GPL8152', 'GPL8153', 'GPL9144'))) %>%  # should use symbol
    dplyr::filter(stringr::str_detect(full, stringr::fixed('entrez', T))) %>%
    dplyr::filter(stringr::str_detect(full, stringr::fixed('id', T))) %>%
    dplyr::filter(!duplicated(accession)) %>%
    dplyr::filter(!stringr::str_detect(full, 'IMAGE clone ID')) %>% # should use symbol
    dplyr::filter(!stringr::str_detect(full, stringr::fixed('name', T))) %>%           # ORF means symbol
    dplyr::filter(!stringr::str_detect(name, stringr::fixed('unigene', T))) %>%      # actually means unigene
    dplyr::filter(!stringr::str_detect(full, stringr::fixed('SNP', T))) %>%            # such as 'rs9651273'
	dplyr::filter(!stringr::str_detect(full, stringr::fixed('(LocusLink)/', T))) %>%  # **/**/**/**/**
    dplyr::filter(!stringr::str_detect(full, 'GEO'))  # Column added by GEO staff to ... in Entrez GEO
}

filter_entrez <- function(info_all) {
	info_all %>%
    dplyr::filter(name != 'GB_ACC') %>%
	dplyr::filter(stringr::str_detect(full, stringr::fixed('entrez', T))) %>%
	dplyr::filter(!duplicated(accession)) %>%
	dplyr::filter(!stringr::str_detect(name, stringr::fixed('unigene', T))) %>%        # should use unigene
	dplyr::filter(!stringr::str_detect(full, stringr::fixed('SNP', T))) %>%            # such as 'rs9651273'
	dplyr::filter(!stringr::str_detect(full, stringr::fixed('Protein', T))) %>%        # such as 'BAC87609', 'Q96BN7'
	dplyr::filter(!stringr::str_detect(full, 'GEO'))  # Column added by GEO staff to ... in Entrez GEO
}

filter_ORF <- function(info_all) {
	info_all %>%
    dplyr::filter(!(accession %in% c('GPL6353'))) %>%  # should use ensembl
	dplyr::filter(name %in% c('ORF', 'ORF_LIST')) %>%
	dplyr::filter(!duplicated(accession)) %>%
	dplyr::filter(!stringr::str_detect(full, stringr::regex('ensembl|genbank', T)))
}

filter_symbol <- function(info_all) {
	info_all %>%
    dplyr::filter(!(accession %in% c('GPL6152'))) %>%  	       # RefSeq  is better
	dplyr::filter(!(accession %in% c('GPL544', 'GPL545'))) %>%                 # GB_ACC  is better
	dplyr::filter(stringr::str_detect(name, stringr::fixed('symbol', T))) %>%
	dplyr::filter(name != 'PlateformeSymbole') %>%
	dplyr::filter(!duplicated(accession))
}

filter_ensembl <- function(info_all) {
	info_all %>%
    dplyr::filter(stringr::str_detect(name, stringr::fixed('ensembl', T))) %>%
	dplyr::filter(!duplicated(accession)) %>%
	dplyr::filter(name != 'Ensembl_genes_count') %>%
	dplyr::filter(!stringr::str_detect(full, stringr::fixed('supported', T)))
}

filter_refseq <- function(info_all) {
	info_all %>%
    dplyr::filter(stringr::str_detect(name, stringr::fixed('refseq', T))) %>%
	dplyr::filter(!duplicated(accession)) %>%
	dplyr::filter(name != 'RefSeq Single exon') %>%                     # GPL19407, GB_ACC is better
	dplyr::filter(!stringr::str_detect(full, stringr::fixed('Genome build version', T)))
	# dplyr::filter(stringr::str_detect(full, stringr::fixed('Distance between', T))) %>%  # ... and ... CpG ....
}

filter_genbank <- function(info_all) {
	info_all %>%
    dplyr::filter(name %in% c('GB_ACC' ,'GB_LIST' ,'GENOME_ACC')) %>%
	dplyr::filter(!duplicated(accession))
}

filter_unigene <- function(info_all) {
	info_all %>%
    dplyr::filter(stringr::str_detect(name, stringr::fixed('unigene', T))) %>%
	dplyr::filter(!duplicated(accession))
}

# not supported now -----------------
filter_sequence <- function(info_all) {
	info_all %>%
    dplyr::filter(stringr::str_detect(name, stringr::regex('^sequence$', T))) %>%
	dplyr::filter(!duplicated(accession))
}



# main ------------------
#' @keywords internal
#' @export
guess_platform_type <- function(platform) {
	if (platform$info %>% is.null) return(NULL)

	temp <- special %>% dplyr::filter(accession == platform$accession)
	if (nrow(temp) != 0L) return(make_platform_type(temp$name, temp$database))

	# I add accession column to be compatible with `platform.Rmd`
	info_all <- platform$info %>%
		dplyr::mutate(full = paste0(name, ': ', description)) %>%
		tibble::add_column(accession = platform$accession, .before = 1)

	if (nrow(filter_miRNA(info_all)) != 0L)     return(NULL)
	if (nrow(filter_circRNA(info_all)) != 0L)   return(NULL)
	if (nrow(filter_non_human(info_all)) != 0L) return(NULL)


    temp <- filter_entrez_id(info_all);
    if (nrow(temp) != 0L) return(make_platform_type(temp$name, 'entrez'))

    temp <- filter_entrez(info_all);
    if (nrow(temp) != 0L) return(make_platform_type(temp$name, 'entrez_or_symbol'))

    temp <- filter_ORF(info_all);
    if (nrow(temp) != 0L) return(make_platform_type(temp$name, 'entrez_or_symbol'))

    temp <- filter_symbol(info_all);
    if (nrow(temp) != 0L) return(make_platform_type(temp$name, 'symbol'))

    temp <- filter_ensembl(info_all);
    if (nrow(temp) != 0L) return(make_platform_type(temp$name, 'ensembl'))

    temp <- filter_refseq(info_all);
    if (nrow(temp) != 0L) return(make_platform_type(temp$name, 'genbank'))

    temp <- filter_genbank(info_all);
    if (nrow(temp) != 0L) return(make_platform_type(temp$name, 'genbank'))

    temp <- filter_unigene(info_all);
    if (nrow(temp) != 0L) return(make_platform_type(temp$name, 'unigene'))

    temp <- filter_ORF(info_all);
    if (nrow(temp) != 0L) return(make_platform_type(temp$name, 'entrez_or_symbol'))

    temp <- filter_ORF(info_all);
    if (nrow(temp) != 0L) return(make_platform_type(temp$name, 'entrez_or_symbol'))

    temp <- filter_ORF(info_all);
    if (nrow(temp) != 0L) return(make_platform_type(temp$name, 'entrez_or_symbol'))

	return(NULL)
}

# utilify ------------------

make_as_symbol <- function(as_symbol_from = c('entrez', 'symbol', 'ensembl', 'entrez_or_symbol', 'genbank', 'unigene')) {
	switch(match.arg(as_symbol_from),
		'symbol' = hgnc::as_symbol_from_symbol,
		'entrez' = hgnc::as_symbol_from_entrez,
		'ensembl' = hgnc::as_symbol_from_ensembl,
		'genbank' = hgnc::as_symbol_from_genbank,
		'entrez_or_symbol' = hgnc::as_symbol_from_entrez_or_symbol,
		'unigene' = hgnc::as_symbol_from_unigene
	)
}

make_platform_type <- function(measure, as_symbol_from = c('entrez', 'symbol', 'ensembl', 'entrez_or_symbol', 'genbank', 'unigene')) {
	sep_pattern = switch(match.arg(as_symbol_from),
		'symbol' = '[^-@\\w]+',
		'entrez' = '[^\\d]+',
		'ensembl' = '[^ENSG\\d]+',
		'genbank' = '[^\\w\\.]+',
		'entrez_or_symbol' = '[^-@\\w]+',
		'unigene' = '[^Hs\\.\\d]+'
	)
	list(measure = measure, sep_pattern = sep_pattern, as_symbol_from = as_symbol_from);
}


#' @title make `platform` from GPL accession to test [guess_platform_type()]
#' @keywords internal
fake_platform <- function(accession, gpl_metas) {
	list(accession = accession, info = gpl_metas[[accession]]$info)
}










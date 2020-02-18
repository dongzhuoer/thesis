# Super GSE

GSE may be SuperSeries which contains multiple SubSeries thus multiple GPL. To avoid repetition, these files are no longer needed, following is the code used to parse them:

```r
series_not_in_summary <- dir('data-raw/geo-tsv', 'series', full.names = T) %>% 
	lapply(libzhuoer::read_char_tsv) %>% bind_rows() %>% 
	filter(Taxonomy == 'Homo sapiens', str_detect(`Series Type`, '^Expression profiling by array$')) %>% 
	filter(!(Accession %in% filter(rGEO.data::read_summary('data-raw/gds_result.txt'), str_count(platform, 'GPL') == 1L)$accession))
```

```r
dir('data-raw/GSE-html', full = T) %>% {.[file.size(.) <10]} %>% file.remove

foobar <- function(){
    series_not_in_summary$Accession %>% 
    {setdiff(., dir('data-raw/GSE-html') %>% str_extract('GSE\\d+'))} %>%
    mclapply(. %>% {
        input  <- paste0('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=', .);
        output <- paste0('data-raw/GSE-html/', ., '.html');
        if (!file.exists(output)) download.file(input, output)
    }, mc.cores = 16);
}
foobar()
```

```r
dir('data-raw/GSE-html', full.names = T) %>% mclapply(. %>% {
	tibble(
        series = str_extract(., 'GSE\\d+'), 
        platform = unique(str_extract_all(read_file(.), '(?<=acc=)GPL\\d+')[[1]])
	)
})
```



# GSEA command

I try to run GSEA within R, but finally abort.

```r
#' @title run gsea
#'
#' @param software_dir string. without `/`
#' @param gmt_file string.
#' @param input_file string.
#'
#' @return `NULL`
#'
#' @examples NULL
gsea <- function(bin_file, gmt_file, input_file, output_dir, title, n_permutation, metric, args = '', ...) {
	args %<>% c(
		'-Xmx4096m',
		  paste0('-cp ', bin_file), 'xtools.gsea.Gsea -gui false', paste0('-gmx ', gmt_file),
		  paste0(c('-res', '-cls'), ' ', input_file, c('.txt', '.cls')),
		  paste0('-out ', output_dir), paste0('-rpt_label ', title),
		  paste0('-nperm ', n_permutation), paste0('-metric ', metric), '-collapse false',
		  .
	)

	system2('java', args, ...)
}


# gsea0 < function(software_dir) {gsea(paste0(software_dir, '/gsea-3.0.jar'), paste0(software_dir, '/msigdb_v6.1/msigdb_v6.1_GMTs/msigdb.v6.1.symbols.gmt', 'Pearson'))}
# gsea_command <- . %>%
```



# GDS

GSD is a lot more easier to parse than GSE, but it's only a small fraction of the latter. Finally I decide to try the harder one to utilize more data.

```r
#' @title prepare GSEA input files
#'
#' @description you should run in no collapse mode since `write_gsea_input` would collapse dataset.
#'
#' @param dataset tibble. first column are `id` (hugo gene synbol), the rest are expression data
#' @param output_file string. append extension to be output filename, for example, `'foobar'` will produce files `foobar.txt`, `foobar.cls`.  Make sure `foo/` exist if use `foo/bar`
#' @param interest_gene string. HUGO gene symbol for gene of interest
#' @param collapse_fun function. Select the expression values to use for the single probe that will represent all probe sets for the gene. it should accept a numeric and `na.rm` parameter.
#'
#' @examples
#' \donotrun{
#'     read_gds_dataset('data-raw/GDS817_full.soft.gz') %>% write_gsea_input('data-raw/GDS817', 'ANKLE2')
#' }
#'
#' @section to do: write a separate function: collapse dataset
#'
#' @return logical scalar. `TRUE` if everything is okay
#'
#' @export
write_gsea_input_uncollapse <- function(dataset, output_file, interest_gene, collapse_fun = max) {
	if (!(interest_gene %in% dataset$id)) return(F)

	# dataset <- read_gds_dataset('../GDS/GDS4085_full.soft.gz');
	# interest_gene <- 'ANKLE2';
	# collapse_fun <- max

	max_na <- (ncol(dataset) - 1) * 0.1;

	definite <- dataset %>% dplyr::filter(id %in% lem4::hugo_symbol);
	ambiguous <- dataset %>% dplyr::filter(str_detect(id, '///')) %>%
		dplyr::mutate(id = stringr::str_replace_all(id, ' ', '')) %>%
		dplyr::filter(plyr::laply(stringr::str_split(.data$id, '///'), . %>% is.element(lem4::hugo_symbol) %>% all));

	collapsed <- dplyr::bind_rows(definite, ambiguous) %>% dplyr::mutate_at(-1, as.numeric) %>%
		dplyr::group_by(id) %>% dplyr::summarise_all(collapse_fun, na.rm = T) %>%
		dplyr::mutate_all(. %>% {.[is.infinite(.)] = NA; .}) %>%
	    {dplyr::slice(., dplyr::select(., -1) %>% dplyr::mutate_all(is.na) %>% rowSums %>% {which(. < max_na)})}

	phenotype <- collapsed %>% dplyr::filter(id == interest_gene) %>% dplyr::select(-1);
	if (nrow(phenotype) == 0L) return(F);
	if (as.matrix(phenotype) %>% is.na %>% any) return(F);

	collapsed %>% dplyr::rename('NAME' = 'id') %>%
		tibble::add_column(DESCRIPTION = 'NA', .after = 1) %>%
		readr::write_tsv(paste0(output_file, '.txt'), na = '');
	phenotype %>% readr::format_tsv(col_names = F) %>%
		c('#numeric', paste0('#', interest_gene), .) %>%
		readr::write_lines(paste0(output_file, '.cls'))

	return(T);
}


#' @title prepare GSEA input files
#'
#' @description you should run in no collapse mode since `write_gsea_input` would collapse dataset.
#'
#' @param dataset tibble. first column are `id` (hugo gene synbol), the rest are expression data
#' @param output_file string. append extension to be output filename, for example, `'foobar'` will produce files `foobar.txt`, `foobar.cls`.  Make sure `foo/` exist if use `foo/bar`
#' @param interest_gene string. HUGO gene symbol for gene of interest
#' @param collapse_fun function. Select the expression values to use for the single probe that will represent all probe sets for the gene. it should accept a numeric and `na.rm` parameter.
#'
#' @examples
#' \donotrun{
#'     read_gds_dataset('data-raw/GDS817_full.soft.gz') %>% write_gsea_input('data-raw/GDS817', 'ANKLE2')
#' }
#'
#' @section to do: write a separate function: collapse dataset
#'
#' @return logical scalar. `TRUE` if everything is okay
#'
#' @export
write_gsea_input <- function(dataset, output_file, interest_gene, collapse_fun = max) {
	if (!(interest_gene %in% dataset$id)) return(F)

	# dataset <- read_gds_dataset('../GDS/GDS4085_full.soft.gz');
	# interest_gene <- 'ANKLE2';
	# collapse_fun <- max

	max_na <- (ncol(dataset) - 1) * 0.1;

	definite <- dataset %>% dplyr::filter(id %in% lem4::hugo_symbol);
	ambiguous <- dataset %>% dplyr::filter(str_detect(id, '///')) %>%
		dplyr::mutate(id = stringr::str_replace_all(id, ' ', '')) %>%
		dplyr::filter(plyr::laply(stringr::str_split(.data$id, '///'), . %>% is.element(lem4::hugo_symbol) %>% all));

	collapsed <- dplyr::bind_rows(definite, ambiguous) %>% dplyr::mutate_at(-1, as.numeric) %>%
		dplyr::group_by(id) %>% dplyr::summarise_all(collapse_fun, na.rm = T) %>%
		dplyr::mutate_all(. %>% {.[is.infinite(.)] = NA; .}) %>%
	    {dplyr::slice(., dplyr::select(., -1) %>% dplyr::mutate_all(is.na) %>% rowSums %>% {which(. < max_na)})}

	phenotype <- collapsed %>% dplyr::filter(id == interest_gene) %>% dplyr::select(-1);
	if (nrow(phenotype) == 0L) return(F);
	if (as.matrix(phenotype) %>% is.na %>% any) return(F);

	collapsed %>% dplyr::rename('NAME' = 'id') %>%
		tibble::add_column(DESCRIPTION = 'NA', .after = 1) %>%
		readr::write_tsv(paste0(output_file, '.txt'), na = '');
	phenotype %>% readr::format_tsv(col_names = F) %>%
		c('#numeric', paste0('#', interest_gene), .) %>%
		readr::write_lines(paste0(output_file, '.cls'))

	return(T);
}
```
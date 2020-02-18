

- clean

```bash
lem4 rm output/*/*/edb/{gene_sets.gmt,results.edb}
lem4 "R -e \"dir('output2', recursive = T, '.png$', full = T) %>% file.remove\""
```

- read GSEA output

```r
read_gsea_summary <- function(path) {
	id <- basename(path) %>% str_replace('\\.Gsea\\.\\d+', '');
	# pattern = 'neg_'; correlation = 'negative'
	impl <- function(pattern, correlation) {
		#" empty file just contains a header line, after collapsed into a string, readr will treat it as a file, so I add \n to the end
		report_files <- dir(path, pattern, full.names = T) %>% stringr::str_subset('report');
		if (length(report_files) == 0L) return(tibble::tibble())

		stringr::str_subset(report_files, 'xls$') %>% readr::read_lines() %>% 
            c('\n') %>% stringr::str_replace('\t$', '') %>% paste0(collapse = '\n') %>%
			readr::read_tsv(T, 'cccnnnnnnnc') %>% tibble::add_column(correlation, .before = 1)
	}

	dplyr::bind_rows(impl('neg_', 'negative'), impl('pos_', 'positive')) %>% add_column(id, .before = 1);
}
read_gsea_summary('/path/to/qGSEA/tests/testthat/output/GSE19161.Gsea.1523948718159')
#> # A tibble: 15 x 13
#>    id    correlation NAME  `GS<br> follow … `GS DETAILS`  SIZE     ES    NES
#>    <chr> <chr>       <chr> <chr>            <chr>        <dbl>  <dbl>  <dbl>
#>  1 GSE1… negative    HALL… HALLMARK_HYPOXIA Details ...     19 -0.375 -1.08 
#>  2 GSE1… negative    HALL… HALLMARK_ESTROG… Details ...     17 -0.370 -1.01
#> ...

read_collection_output <- function(collection) {
    dir(collection, full = T) %>% 
        mclapply(qGSEA::read_gsea_summary) %>% bind_rows %>% 
        add_column(collection = basename(collection), .before = 1)
}

dir('output', full = T) %>% 
    lapply(read_collection_output) %>% bind_rows %>% 
    rename(accession = id) %>% write_rds('data-raw/gsea_output_full.rds')
```


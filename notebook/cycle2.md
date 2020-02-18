```r
library(tidyverse)
```



# 02/02

```r
BiocInstaller::biocLite('GEOquery')
```

基本没用



# 02/03


```bash
# copy raw data from workstation to local
rsync --ignore-existing -vh -e "ssh -p port" user@ip:/path/to/project/raw /path/to/project
```

```r
# check gzip file integrity
libzhuoer::check_gzip('raw')
```



# 02/08

explore clusterProfiler vignette



# 03/01 gene sets in MSigDB

```yaml
MSigDB:
    genesets < 1000:
        c1.all.v6.1.symbols.gmt: 326
        c2.cp.biocarta.v6.1.symbols.gmt: 217
        c2.cp.kegg.v6.1.symbols.gmt: 186
        c2.cp.reactome.v6.1.symbols.gmt: 674
        c3.mir.v6.1.symbols.gmt: 221
        c3.tft.v6.1.symbols.gmt: 615
        c4.cgn.v6.1.symbols.gmt: 427
        c4.cm.v6.1.symbols.gmt: 431
        c5.cc.v6.1.symbols.gmt: 580
        c5.mf.v6.1.symbols.gmt: 901
        c6.all.v6.1.symbols.gmt: 189
        h.all.v6.1.symbols.gmt: 50
    genesets > 1000:
        c2.cgp.v6.1.symbols.gmt: 3409
        c5.bp.v6.1.symbols.gmt: 4436
        c7.all.v6.1.symbols.gmt: 4872
    not atom:
        c2.all.v6.1.symbols.gmt: 4738
        c2.cp.v6.1.symbols.gmt: 1329
        c3.all.v6.1.symbols.gmt: 836
        c4.all.v6.1.symbols.gmt: 858
        c5.all.v6.1.symbols.gmt: 5917
        msigdb.v6.1.symbols.gmt: 5917
```

I remember large collection would cause error (maybe it just waste too much time). According to size and biological significance, finally I decide to use the first 12 (<1000), and `c2.cgp`.


# 03/02

```r
geo_all <- rGEO.data::read_summary(system.file('extdata/gds_result.txt', package = 'rGEO.data'))
gds <- geo_all %>% filter(species == 'Homo sapiens', type == 'DataSet') 
gds_ftp <- gds$accession %>% {paste0(rGEO::gds_ftp(.), '/soft/', ., '_full.soft.gz')}
head(gds_ftp)
#> [1] "ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS6nnn/GDS6177/soft/GDS6177_full.soft.gz"
#> [2] "ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS6nnn/GDS6100/soft/GDS6100_full.soft.gz"
#> [3] "ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS6nnn/GDS6083/soft/GDS6083_full.soft.gz"
#> [4] "ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS6nnn/GDS6082/soft/GDS6082_full.soft.gz"
#> [5] "ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS6nnn/GDS6063/soft/GDS6063_full.soft.gz"
#> [6] "ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS6nnn/GDS6010/soft/GDS6010_full.soft.gz"

write_lines(gds_ftp, 'download.md');
# then use `aria2c` to download links in `download.md` on workstation
```



# 03/03

```r
# read_gds ------------------
#' @title read GDS files
#'
#' @param path string. path to GDS DataSet full SOFT file, such as https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3958
#'
#' @return tibble
NULL 

#' @rdname  read_gds
#' @export
read_gds_dataset <- function(path) {
	content <- readr::read_lines(path);
	
    n_sample <- content %>% stringr::str_subset('dataset_sample_count') %>% 
    	stringr::str_extract('\\d+') %>% as.integer;

	content %>% {.[stringr::str_which(., 'dataset_table') %>% {seq(.[1] + 1, .[2] - 1)}]} %>% 
		paste0(collapse = '\n') %>% libzhuoer::read_char_tsv() %>% 
		select('id' = 'Gene symbol', seq(3, 2 + n_sample))
}

#' @rdname  read_gds
#' @export
read_gds_annotation <- function(path) {
	# path <- 'data-raw/GDS817_full.soft.gz'
	content <- readr::read_lines(path) %>% stringr::str_subset('^[!#\\^]')
	
	get_value <- function(name) {stringr::str_subset(content, name) %>% 
			stringr::str_replace(name, '') %>% {if (length(.) == 0L) '' else .}};
	
	tibble::tibble(
		accession = stringr::str_extract(path, 'GDS\\d+'),
		pubmed = get_value('!dataset_pubmed_id = '),
		subsets = get_value('!subset_description = ') %>% paste0(collapse = ', ')
	)
}
```



# 03/04 run GSEA on GDS

For GDS, almost all "Gene symbol" is same as HUGO symbol (MGC12488 is [deprecated](https://www.ncbi.nlm.nih.gov/gene/84786)).

```r
lem4::read_gds_dataset('path/to/GDS1381_full.soft.gz') %>% lem4::write_gsea_input('gsea/GDS1381', 'ANKLE2')
# write *cls
```


```bash
java -Xmx512m -cp /path/to/gsea-3.0.jar xtools.gsea.Gsea -gmx /path/to/msigdb_v6.1/msigdb_v6.1_GMTs/c5.mf.v6.1.symbols.gmt -res gsea/GDS1381.txt -cls gsea/GDS1381.cls#ANKLE2  -chip gsea/GDS1381.chip -out output -collapse false -nperm 1000 -rpt_label my_analysis -metric Pearson -gui false

java -Xmx512m -cp /path/to/gsea-3.0.jar xtools.gsea.Gsea -gmx /path/to/msigdb_v6.1/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt -res gsea/GDS4090.txt -cls gsea/GDS4090.cls#ANKLE2  -chip gsea/GDS4090.chip -out output -collapse true -nperm 10 -permute phenotype -rpt_label my_analysis -metric Pearson -gui false

# java -Xmx512m -cp /path/to/gsea-3.0.jar xtools.gsea.LeadingEdgeTool -dir output/my_analysis.Gsea.1520140414563 -gsets set_names_comma_delimited
```



# 03/07

```r
left_join(
	dir('path/to/GDS', full.names = T) %>% mclapply(lem4::read_gds_annotation) %>% bind_rows,
	rGEO.data::read_summary(system.file('extdata/gds_result.txt', package = 'rGEO.data'))
) -> gds

gds %>% mutate(sample = str_extract(platform, '\\d+(?= Sample)')) %>% select(1:2, 11, 9, 3, 10) %>% 
	filter(!(accession %in% read_lines('data-raw/non_breast-ovarian_cancer_GDS.md'))) %>% 
	write_tsv('output/GDS_info.xlsx')
```



# 03/08

```r
libzhuoer::read_char_csv('data-raw/platform.csv') %>% 
    filter(Taxonomy == 'Homo sapiens') %>% {.$Accession} %>% rGEO::gpl_soft_ftp() %>% 
    write_lines('download.md')
# then use `aria2c` to download links in `download.md` on workstation
```



```r
library(tidyverse)
```



# 01/28

```r
system.file('extdata/gds_result.txt', package = 'rGEO.data') %>% read_lines %>% 
    str_extract('ftp://[\\w\\W]+GDS\\d+') %>% na.omit %>% as.character
#> ...
#> [1772] "ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDSnnn/GDS12"
```

only 1772 data has GDS file, others has only GSE, such as https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77386



# 01/29

```r
system.file('extdata/gds_result.txt', package = 'rGEO.data') %>% read_lines %>% 
    str_extract('Accession: \\w+') %>% na.omit %>% 
	str_replace('Accession: ', '') %>% str_replace('\\d+', '') %>% unique
#> [1] "GDS" "GSE"
```

There are two type of accession: GDS & GSE.

```r
contents <- system.file('extdata/gds_result.txt', package = 'rGEO.data') %>% 
    read_file %>% str_replace('^\n', '') %>% 
    {str_split(., '\n\n')[[1]]} %>% lapply(read_lines)

contents %>% sapply(length) %>% table()
#> .
#>     6     7     8     9 
#>    16 23154    38     2
```

I try various ways to explore the format of `gds_result.txt` and finally write `rGEO.data::read_summary`.

```r
geo_all <- rGEO.data::read_summary(system.file('extdata/gds_result.txt', package = 'rGEO.data'))

geo_all$accession %>% str_replace('\\d+', '') %>% table
#> .
#>   GDS   GSE 
#>  1772 21438 

# SuperSeries would contain more than one platform
geo_all %>% mutate(n_gpl = str_count(platform, 'GPL')) %>% 
    group_by(n_gpl) %>% summarise(n = n())
#> # A tibble: 4 x 2
#>   n_gpl     n
#>   <int> <int>
#> 1     0   329
#> 2     1 20210
#> 3     2  2208
#> 4     3   463

# check
geo_all$ftp %>% str_detect('^ftp://') %>% all
#> [1] TRUE
geo_all$data_type %>% str_detect('Expression profiling by array') %>% all
#> [1] TRUE
```

- explore relationship between GDS & GSE

```r
setdiff(
    geo_all %>% filter(type == 'Series') %>% {.$platform} %>% 
        str_extract_all('GDS\\d{1,5}') %>% unlist,
    geo_all %>% filter(type == 'DataSet') %>% {.$accession}
)
#> [1] "GDS4998" "GDS4234" "GDS3953" "GDS3958" "GDS3960" "GDS1649" "GDS592" 

# inspect them in the browser
.Last.value %>% paste0('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=', .) %>% 
    libzhuoer::browse_url()
```

For example, GSE48200 contains GDS4998-GPL6246-Mus_musculus and GDS4997-GPL6244-Homo_sapiens.

For more thant one species, 

|                | title | description |
| -------------- | ----- | ----------- |
| breast cancer  | 18    | 16          |
| ovarian cancer | 4     | 0           |

Thus I use 'Homo sapiens' only.

```r
geo <- geo_all %>% filter(!str_detect(data_type, 'Third-party reanalysis')) %>% 
	filter(species == 'Homo sapiens', type == 'Series') 
```

|                | title | description |
| -------------- | ----- | ----------- |
| breast cancer  | 1090  | 1195        |
| ovarian cancer | 273   | 237         |



# 01/30

- Prof Zhu wants to see the correlation between lem4 and many individual genes, so I search NCBI for their Entrez ID.

- [lem-4](https://www.ncbi.nlm.nih.gov/gene/23141)
- [E2F1](https://www.ncbi.nlm.nih.gov/gene/1869)
- [CCND1](https://www.ncbi.nlm.nih.gov/gene/595)
- [CCND2](https://www.ncbi.nlm.nih.gov/gene/894)
- [CCNE1](https://www.ncbi.nlm.nih.gov/gene/898)
- [MYC](https://www.ncbi.nlm.nih.gov/gene/4609)
- [CDK1](https://www.ncbi.nlm.nih.gov/gene/983)
- [CDK4](https://www.ncbi.nlm.nih.gov/gene/1019)
- [CDKN1A](https://www.ncbi.nlm.nih.gov/gene/1026)
- [CDKN2A](https://www.ncbi.nlm.nih.gov/gene/1029)
- [AURKA](https://www.ncbi.nlm.nih.gov/gene/6790)
- [ERα](https://www.ncbi.nlm.nih.gov/gene/2099)



# 01/31

- single-gene-correlation

```r
zhu <- tribble(
	~name, ~symbol, ~id,
	'lem-4', 'lem4', '23141',
	'E2F1', 'E2F1', '1869',
	'CCND1', 'CCND1', '595',
	'CCND2', 'CCND2', '894',
	'CCNE1', 'CCNE1', '898',
	'MYC', 'MYC', '4609',
	'CDK1', 'CDK1', '983',
	'CDK4', 'CDK4', '1019',
	'p21', 'CDKN1A', '1026',
	'p16', 'CDKN2A', '1029',
	'AURKA', 'AURKA', '6790',
	'ERα', 'ESR1', '2099'
)
```

```r
calc_expression <- function(gse, Gene) {
	Probe <- gse$chip %>% filter(str_detect(gene, Gene)) %>% `$`(probe);
	gse$mat %>% filter(probe %in% Probe) %>% {.[, -1]} %>% 
		{if (nrow(.) > 0) colMeans(., na.rm = T) else rep(NA, ncol(.))};
}
# calc_expression(gse69078, '-1')

extract_submat <- function(gse) {
	plyr::laply(zhu$id, . %>% calc_expression(gse, .)) %>% t %>% as_tibble() %>%
		{colnames(.) <- zhu$symbol; .} %>%
		add_column(Sample = colnames(gse$mat)[-1], .before = 1)
}
# submat69078 <- extract_submat(gse69078)	
```

On workstation, I extract expression matrix and save them under `gse/`, like `gse/GSE69031.rds`.

```r
x <- dir('gse', full.names = T) %>% 
    mclapply(. %>% {tryCatch(read_rds(.) %>% extract_submat, error = function(e){message(.)})} )
names(x) = dir('gse', full.names = T)
write_rds(x, 'submats.rds')
```

Then copy `submats.rds` to local.

```r
submats <- read_rds('data-raw/submats.rds')[-211]

r_t_p <- function(x, y) {
	na_nan <- . %>% {any(is.na(.) | is.nan(.))};
	if (na_nan(x) || na_nan(y)) return(NA);
	
	n = length(x);
	if (n <= 2) return(NA)
	
	r <- cor(x, y);
	abs(r)*sqrt((n - 2)/(1 - r^2)) %>% pt(n-2)
}

plyr::ldply(
	seq_along(submats),
	function(i) {
		sapply(
			seq_len(nrow(zhu))[-1],
			function(j){
				submats[[i]] %>% {r_t_p(.[[2]], .[[j + 1]])}
			}
		) %>% tibble(p = ., gene = zhu$name[-1], sample = names(submats)[i])	
	}	
) %>% as_tibble -> sing_cor_df

sing_cor_df %>% filter(!is.na(p)) %>% ggplot() + 
	geom_histogram(aes(p, fill = p < 0.05), binwidth = 0.01) + facet_wrap(~gene)
```


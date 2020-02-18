```r
library(tidyverse)
```



# 03/20

```r
hgnc <- libzhuoer::read_char_tsv('path/to/hgnc_complete_set.txt.gz') %>% 
	filter(status == 'Approved')

hgnc_sample <- hgnc %>% sample_n(100)

read_rds('data-raw/platform_column_name.rds') %>% 
    {.[!str_detect(sapply(., class), 'error')]} %>% 
    yaml::write_yaml('output/platform_header.yaml')
```



# 03/24 try to convert various ID in platform file to HUGO symbol

- HGNC's ena actually means Genbank accessions which is not included in RefSeq.

```r
hgnc %>% select(symbol, entrez_id:refseq_accession) %>% filter(refseq_accession == ena)
intersect(hgnc$ena, hgnc$refseq_accession)
```

- GSE 只能说没得救了

For example, in GSE93601 platform table, on line 1129 (search 'TC0101128'), there is a _line break_ that split one row, you can't even read tsv.

- we shouldn't believe genebank

```r
filter(gse78958_raw$platform, str_detect(ENTREZ_GENE_ID, '23141')) %>% print(width = Inf)
# A tibble: 3 x 16
  ID          GB_ACC   SPOT_ID `Species Scientific Name` `Annotation Date`
  <chr>       <chr>    <chr>   <chr>                     <chr>
1 212200_at   AW274877 NA      Homo sapiens              Oct 6, 2014
2 212201_at   AW274877 NA      Homo sapiens              Oct 6, 2014
3 213962_s_at AI924382 NA      Homo sapiens              Oct 6, 2014
  `Sequence Type`    `Sequence Source`
  <chr>              <chr>
1 Consensus sequence GenBank
2 Consensus sequence GenBank
3 Consensus sequence GenBank
  `Target Description`
  <chr>
1 Consensus includes gb:AK025933.1 /DEF=Homo sapiens cDNA: FLJ22280 fis, clone…
2 gb:AW274877 /DB_XREF=gi:6661907 /DB_XREF=xm62a09.x1 /CLONE=IMAGE:2688760 /FE…
3 gb:AI924382 /DB_XREF=gi:5660346 /DB_XREF=wn60d01.x1 /CLONE=IMAGE:2449825 /FE…
  `Representative Public ID` `Gene Title`
  <chr>                      <chr>
1 AW274877                   ankyrin repeat and LEM domain containing 2
2 AW274877                   ankyrin repeat and LEM domain containing 2
3 AI924382                   ankyrin repeat and LEM domain containing 2
  `Gene Symbol` ENTREZ_GENE_ID
  <chr>         <chr>
1 ANKLE2        23141
2 ANKLE2        23141
3 ANKLE2        23141
  `RefSeq Transcript ID`
  <chr>
1 NM_015114 /// XM_005266159 /// XM_005266160 /// XM_005266161 /// XM_006719735
2 NM_015114 /// XM_005266159 /// XM_005266160 /// XM_005266161 /// XM_006719735
3 NM_015114 /// XM_005266159 /// XM_005266160 /// XM_005266161 /// XM_006719735
  `Gene Ontology Biological Process`
  <chr>
1 0000278 // mitotic cell cycle // traceable author statement /// 0007049 // c…
2 0000278 // mitotic cell cycle // traceable author statement /// 0007049 // c…
3 0000278 // mitotic cell cycle // traceable author statement /// 0007067 // m…
  `Gene Ontology Cellular Component`
  <chr>
1 0005737 // cytoplasm // inferred from direct assay /// 0005783 // endoplasmi…
2 0005737 // cytoplasm // inferred from direct assay /// 0005783 // endoplasmi…
3 0005789 // endoplasmic reticulum membrane // traceable author statement /// …
  `Gene Ontology Molecular Function`
  <chr>
1 0005515 // protein binding // inferred from physical interaction /// 0008601…
2 0005515 // protein binding // inferred from physical interaction /// 0008601…
3 0005515 // protein binding // inferred from physical interaction /// 0008601…
```



# 03/25 

```r
breast_cancer <- rGEO.data::dataset %>% 
    filter(species == 'Homo sapiens') %>% 
	filter(str_count(platform, 'GPL') == 1L) %>% 
 	filter(str_detect(title, '[Bb]reast|[Oo]varian')) %>%  
	filter(str_detect(title, '[Cc]ancer|[Aa]denocarcinoma|[Cc]arcinoma|[Tt]umor')) %>% 
	mutate(description = str_replace(description, '\\(Submitter supplied\\) ', '')) %>% 
    mutate(sample = as.integer(str_extract(platform, '\\d+(?= Sample)'))) %>%
    select(accession, sample, type:description) %T>% write_rds('data/breast_cancer.rds')
#> # A tibble: 1,332 x 9
#>    accession sample type   species  ftp    platform  data_type title description
#>    <chr>      <int> <chr>  <chr>    <chr>  <chr>     <chr>     <chr> <chr>      
#>  1 GDS5819       13 DataS… Homo sa… ftp:/… Platform… "\tExpre… Meta… Analysis o…
#>  2 GDS5816        9 DataS… Homo sa… ftp:/… Platform… "\tExpre… Meth… Analysis o…
#> ...
```

```r
geo_breast <- breast_cancer %>% 
    filter(str_detect(title, '[Bb]reast'), str_count(platform, 'GPL') == 1L) %T>% 
    {cowplot::plot_grid( 
        ggplot(., aes(sample)) + geom_density(adjust = 0.25), 
        ggplot(filter(., sample > 100), aes(sample)) + geom_density(adjust = 0.2) + 
            scale_x_continuous(breaks = seq(0,1000,100)) + geom_vline(xintercept = 232, color = 'blue'),
        ncol = 2
    ) %>% print} %>% 
    filter(., sample > 225) %>% arrange(desc(sample)) %T>% print(n = 50) 
#> # A tibble: 37 x 9
#>    accession sample type   species  ftp   platform data_type title  description 
#>    <chr>      <int> <chr>  <chr>    <chr> <chr>    <chr>     <chr>  <chr>       
#>  1 GSE93601    1110 Series Homo sa… ftp:… Platfor… "\tExpre… Alcoh… Investigati…
#> ...
#> 37 GSE50811     241 Series Homo sa… ftp:… Platfor… "\tExpre… Breas… Eribulin me…
#> ...

geo_ovary <- breast_cancer %>% 
    filter(str_detect(title, '[Oo]varian'), str_count(platform, 'GPL') == 1L) %T>% 
    {cowplot::plot_grid(
        ggplot(., aes(sample)) + geom_density(adjust = 0.25), 
        ggplot(filter(., sample > 10), aes(sample)) + geom_density(adjust = 0.2) + 
            scale_x_continuous(breaks = seq(0,1000,50)) + geom_vline(xintercept = 50, color = 'blue'),
        ncol = 2
    ) %>% print} %>% 
    filter(., sample > 50) %>% arrange(desc(sample)) %T>% print(n = 100)
#> # A tibble: 52 x 9
#>    accession sample type   species  ftp    platform  data_type title description
#>    <chr>      <int> <chr>  <chr>    <chr>  <chr>     <chr>     <chr> <chr>      
#>  1 GSE74357     529 Series Homo sa… ftp:/… Platform… "\tExpre… Gene… The goal o…
#> ...
#> 52 GSE8057       51 Series Homo sa… ftp:/… Platform… "\tExpre… Expr… Time-cours…
#> ...
```

- select small data for unittest

```r
platform_type_df <- breast_cancer %>% 
    filter(str_count(platform, 'GPL') == 1L) %>% 
	mutate(platform = str_extract(platform, 'GPL\\d+')) %>% 
	mutate(database = sapply(
        platform, 
        . %>% guess_platform_type %>% {if (is.null(.)) NA else formalArgs(.$as_symbol)}
    )) %>% arrange(sample) %>% select(accession, sample, database, platform)

dir('raw', recursive = T, full = T) %>% file.info %>% 
    as_tibble(rownames = 'file') %>% write_rds('raw.rds')
file_size_df <- read_rds('data-raw/raw.rds') %>% 
    mutate(accession = str_extract(file, 'GSE\\d+')) %>% 
	filter(!is.na(accession)) %>% group_by(accession) %>% 
	summarise(size = prod(size), both = n() == 2L) %>% filter(both)

find_small_gse <- inner_join(platform_type_df, file_size_df) %>% arrange(size)
```



# 03/26

- make GSEA input from GSE on workstation

```r
mclapply(
    str_extract(dir('raw/GSE'), 'GSE\\d+'),
    . %>% {
        data <- lem4::read_gse(
            paste0('raw/GSE/', ., '_family.soft.gz'), 
            paste0('raw/GSE/', ., '_series_matrix.txt.gz'), 
            T
        )
        lem4::make_gsea_input(data, 'ANKLE2', paste0('gsea/', .))
    }, mc.preschedule = F
)
```


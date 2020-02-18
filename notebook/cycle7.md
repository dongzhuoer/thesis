```r
library(tidyverse)
```



# 04/29 synthesis GSEA output

- make `gsea_output_full` and `leading_edge`

```r
bad_collection <- c('c3.mir', 'c4.cgn', 'c7.all')

sample_df <- read_rds('data/breast_cancer.rds') %>% filter(str_detect(accession, 'GSE')) %>%
    select(accession, n_sample = sample) %T>% print
#> # A tibble: 1,193 x 2
#>    accession n_sample
#>    <chr>        <int>
#>  1 GSE114359        6
#>  2 GSE106549        4
#> ...

gsea_output_full <- 'data/gsea_output_full.rds' %>% read_rds %>%
    rename(gene_set = NAME, p = 'NOM p-val', FDR = 'FDR q-val') %>%
    select(collection, correlation, gene_set, p, FDR, accession) %>%
    inner_join(sample_df) %>% filter(!(collection %in% bad_collection)) %T>% print
#> # A tibble: 8,000,434 x 7
#>    collection correlation gene_set      p   FDR accession n_sample
#>    <chr>      <chr>       <chr>     <dbl> <dbl> <chr>        <int>
#>  1 c1.all     negative    CHRXQ12  0.0975     1 GSE10046         6
#>  2 c1.all     negative    CHR12Q12 0.0966     1 GSE10046         6
#> ...

leading_edge <- read_rds('data/leading_edge.rds') %>% 
	inner_join(sample_df) %>% filter(!(collection %in% bad_collection)) %T>% print
#> # A tibble: 8,000,434 x 7
#>    collection correlation gene_set      p   FDR accession n_sample
#>    <chr>      <chr>       <chr>     <dbl> <dbl> <chr>        <int>
#>  1 c1.all     negative    CHRXQ12  0.0975     1 GSE10046         6
#>  2 c1.all     negative    CHR12Q12 0.0966     1 GSE10046         6
#> ...
```

- drop aberrant dataset (two many significant result)

```r
gsea_output_aberrant <- gsea_output_full %>% group_by(collection, accession) %>% 
    summarise(aberrant = sum(FDR < 0.25) / n() >= 0.25) %>% 
    filter(aberrant) %>% select(-aberrant) %T>% print
#> # A tibble: 532 x 2
#> # Groups:   collection [10]
#>    collection accession
#>    <chr>      <chr>    
#>  1 c1.all     GSE107251
#>  2 c1.all     GSE10879 
#> ...

gsea_output_full %<>% anti_join(gsea_output_aberrant) %T>% print
#> # A tibble: 532 x 2
#> # Groups:   collection [10]
#>    collection accession
#>    <chr>      <chr>    
#>  1 c1.all     GSE107251
#>  2 c1.all     GSE10879 
#> ...

leading_edge %<>% anti_join(gsea_output_aberrant) %T>% print
#> # A tibble: 358,180 x 6
#>    PROBE    `RANK METRIC SCORE` gene_set collection accession n_sample
#>    <chr>                  <dbl> <chr>    <chr>      <chr>        <int>
#>  1 ANKLE2                 1     CHR12Q24 c1.all     GSE102484      683
#>  2 POLE                   0.528 CHR12Q24 c1.all     GSE102484      683
#> ...
```

- summary result (filter small sample)

```r
msigdb <- xml2::read_xml('/path/to/GSEA/msigdb_v6.1/msigdb_v6.1.xml') %>% 
    xml2::as_list() %>% {lapply(.[[1]], attributes)} %T>% 
    {names(.) <- sapply(., . %>% {.$STANDARD_NAME}) %>% toupper}
str(head(msigdb))
#> List of 6
#>  $ AAANWWTGC_UNKNOWN:List of 25
#> ...
#>   ..$ ORGANISM            : chr "Homo sapiens"
#> ...
#>   ..$ CATEGORY_CODE       : chr "C3"
#>   ..$ SUB_CATEGORY_CODE   : chr "TFT"
#> ...
#>   ..$ MEMBERS_SYMBOLIZED  : chr "MEF2C,ATP1B1,RORA,CITED2,APP,| __truncated__
#>   ..$ MEMBERS_EZID        : chr "4208,481,6095,10370,351,4216,| __truncated__
#> ...
#>  $ AAAYRNCTG_UNKNOWN:List of 25
#> ...
min_n_sample <- 100L;

gsea_summary <- gsea_output_full %>% 
	filter(n_sample >= min_n_sample) %>% group_by(collection, correlation, gene_set) %>% 
	summarise(N = n(), n = sum(FDR < 0.25) - sum(FDR*(FDR < 0.25))) %>% ungroup() %>% 
	mutate(
		title = msigdb[gene_set] %>% sapply(. %>% {.$DESCRIPTION_BRIEF}), 
		url = msigdb[gene_set] %>% sapply(. %>% {.$EXTERNAL_DETAILS_URL}), 
		description = msigdb[gene_set] %>% sapply(. %>% {.$DESCRIPTION_FULL})
	) %>% group_by(collection) %>% arrange(desc(n)) %T>% 
    print %>% write_rds('data/gsea_summary.rds', 'bz2')
#> # A tibble: 19,610 x 8
#> # Groups:   collection [10]
#>    collection correlation gene_set       N     n title     url     description  
#>    <chr>      <chr>       <chr>      <int> <dbl> <chr>     <chr>   <chr>        
#>  1 c1.all     positive    CHR12Q24     102  35.8 Genes in… http:/… Genes in cyt…
#>  2 c2.cgp     positive    HASLINGER…    89  33.4 Genes ch… "\"NON… PURPOSE: Gen…
#> ...

leading_edge_per_data <- leading_edge %>% 
    filter(n_sample >= min_n_sample) %>% 
    left_join(gsea_output_full) %>% 
    group_by(collection, accession, PROBE) %>% summarise(score = sum(1 - FDR)) %>% 
    group_by(collection, accession) %>% arrange(desc(score)) %>% slice(1:10) %T>% 
    print() %>% write_rds('data/leading_edge_per_data.rds', 'bz2')
#> # A tibble: 5,163 x 4
#> # Groups:   collection, accession [517]
#>    collection accession PROBE    score
#>    <chr>      <chr>     <chr>    <dbl>
#>  1 c1.all     GSE102484 SART3     1.99
#>  2 c1.all     GSE102484 ABCB9     1 
#> ...
```



# # 04/30 select final result

```r
gsea_summary <- read_rds('data/gsea_summary.rds')

core <- bind_rows(
	gsea_summary %>% filter(collection == 'h.all', n >= 15),
	gsea_summary %>% filter(collection == 'c1.all') %>% slice(1),
	gsea_summary %>% filter(collection == 'c2.cgp', n >= 29),
	gsea_summary %>% filter(collection == 'c2.cp.biocarta', n >= 20),
	gsea_summary %>% filter(collection == 'c2.cp.kegg') %>% slice(1),
	gsea_summary %>% filter(collection == 'c2.cp.reactome') %>% slice(1),
	gsea_summary %>% filter(collection == 'c3.tft') %>% slice(1),
	gsea_summary %>% filter(collection == 'c4.cm', n >= 24),
	gsea_summary %>% filter(collection == 'c5.all', n >= 27),
	gsea_summary %>% filter(collection == 'c6.all') %>% slice(1)
) %>% select(1,3,4,5) %>% mutate_at('n', round) %T>% print
#> # A tibble: 20 x 4
#> # Groups:   collection [10]
#>    collection     gene_set                                               N     n
#>    <chr>          <chr>                                              <int> <dbl>
#>  1 h.all          HALLMARK_MITOTIC_SPINDLE                              71    25
#>  2 h.all          HALLMARK_G2M_CHECKPOINT                               76    24
#>  3 h.all          HALLMARK_E2F_TARGETS                                  71    23
#>  4 h.all          HALLMARK_MYC_TARGETS_V1                               67    22
#>  5 h.all          HALLMARK_MYC_TARGETS_V2                               73    19
#>  6 c1.all         CHR12Q24                                             102    36
#>  7 c2.cgp         HASLINGER_B_CLL_WITH_CHROMOSOME_12_TRISOMY            89    33
#>  8 c2.cgp         LI_WILMS_TUMOR_VS_FETAL_KIDNEY_1_DN                   85    30
#>  9 c2.cgp         TOYOTA_TARGETS_OF_MIR34B_AND_MIR34C                   86    30
#> 10 c2.cp.biocarta BIOCARTA_G1_PATHWAY                                   79    23
#> 11 c2.cp.biocarta BIOCARTA_G2_PATHWAY                                   79    22
#> 12 c2.cp.biocarta BIOCARTA_MCM_PATHWAY                                  79    20
#> 13 c2.cp.kegg     KEGG_CELL_CYCLE                                       87    24
#> 14 c2.cp.reactome REACTOME_TRANSPORT_OF_MATURE_TRANSCRIPT_TO_CYTOPL…    72    15
#> 15 c3.tft         E2F_Q3_01                                             72    30
#> 16 c4.cm          MODULE_98                                             75    24
#> 17 c4.cm          MODULE_198                                            76    24
#> 18 c5.all         GO_NUCLEAR_ENVELOPE_ORGANIZATION                      86    28
#> 19 c5.all         GO_NUCLEAR_PORE                                       84    28
#> 20 c6.all         RB_P107_DN.V1_UP                                      76    28

extra <- bind_rows(
	gsea_summary %>% filter(collection == 'h.all') %>% slice(6:8),
	gsea_summary %>% filter(collection == 'c1.all', gene_set == 'CHR12Q23', correlation == 'positive'),
	gsea_summary %>% filter(collection == 'c2.cgp', gene_set == 'NIKOLSKY_BREAST_CANCER_12Q24_AMPLICON', correlation == 'positive'),
	gsea_summary %>% filter(collection == 'c2.cp.biocarta', gene_set == 'BIOCARTA_CELLCYCLE_PATHWAY', correlation == 'positive'),
	gsea_summary %>% filter(collection == 'c2.kegg', gene_set == 'KEGG_DNA_REPLICATION', correlation == 'positive'),
	# so many
	gsea_summary %>% filter(collection == 'c2.cp.reactome', gene_set %in% c('REACTOME_MITOTIC_G2_G2_M_PHASES', 'REACTOME_CELL_CYCLE', 'REACTOME_TRANSPORT_OF_MATURE_MRNA_DERIVED_FROM_AN_INTRONLESS_TRANSCRIPT'), correlation == 'positive'),
	gsea_summary %>% filter(collection == 'c4.cm', gene_set == 'MODULE_320', correlation == 'positive'),
	# so many
	gsea_summary %>% filter(collection == 'c5.all', gene_set %in% c('GO_REGULATION_OF_MICROTUBULE_POLYMERIZATION_OR_DEPOLYMERIZATION', 'GO_ATP_DEPENDENT_DNA_HELICASE_ACTIVITY', 'GO_DNA_HELICASE_ACTIVITY', 'GO_RNA_HELICASE_ACTIVITY', 'GO_TRNA_TRANSPORT', 'GO_REGULATION_OF_MICROTUBULE_BASED_PROCESS', 'GO_MITOTIC_SPINDLE'), correlation == 'positive'),
	gsea_summary %>% filter(collection == 'c6.all', gene_set %in% c('E2F1_UP.V1_UP', 'MYC_UP.V1_UP'), correlation == 'positive')
) %>% select(1,3,4,5) %>% mutate_at('n', round) %T>% print
#> # A tibble: 19 x 4
#> # Groups:   collection [8]
#>    collection    gene_set                                                N     n
#>    <chr>         <chr>                                               <int> <dbl>
#>  1 h.all         HALLMARK_MTORC1_SIGNALING                              62    15
#>  2 h.all         HALLMARK_DNA_REPAIR                                    65    14
#>  3 h.all         HALLMARK_UNFOLDED_PROTEIN_RESPONSE                     64    13
#>  4 c1.all        CHR12Q23                                               91    11
#>  5 c2.cgp        NIKOLSKY_BREAST_CANCER_12Q24_AMPLICON                  66    28
#>  6 c2.cp.biocar… BIOCARTA_CELLCYCLE_PATHWAY                             77    18
#>  7 c2.cp.reacto… REACTOME_MITOTIC_G2_G2_M_PHASES                        66    14
#>  8 c2.cp.reacto… REACTOME_CELL_CYCLE                                    64    14
#>  9 c2.cp.reacto… REACTOME_TRANSPORT_OF_MATURE_MRNA_DERIVED_FROM_AN_…    64    14
#> 10 c4.cm         MODULE_320                                             76    22
#> 11 c5.all        GO_REGULATION_OF_MICROTUBULE_POLYMERIZATION_OR_DEP…    85    26
#> 12 c5.all        GO_ATP_DEPENDENT_DNA_HELICASE_ACTIVITY                 89    26
#> 13 c5.all        GO_DNA_HELICASE_ACTIVITY                               91    26
#> 14 c5.all        GO_RNA_HELICASE_ACTIVITY                               85    26
#> 15 c5.all        GO_TRNA_TRANSPORT                                      83    25
#> 16 c5.all        GO_REGULATION_OF_MICROTUBULE_BASED_PROCESS             85    25
#> 17 c5.all        GO_MITOTIC_SPINDLE                                     82    25
#> 18 c6.all        E2F1_UP.V1_UP                                          78    23
#> 19 c6.all        MYC_UP.V1_UP                                           67    14

# analysis of `leading_edge` is omitted
```



# 05/15 store GSEA output

GSEA output directory structure:

- *label*, the label of analysis you specified, I use GSE accession, such as GSE62555
- *hash*, series of numbers to distinguish different runs of same analysis, such as 1524003938698 (user often run same analysis with different parameters, so GSEA add this filed to the end of the output dir)

I once considered the following as essential files,

```
edb/*label*_collapsed_to_symbols.rnk            ranked gene scores
edb/results.edb                                 permutations' ES
index.html,xtools.css                           main report
*label*.Gsea.*hash*.rpt                         parameters
gsea_report_for_ANKLE2_{neg,pos}_*hash*.xls     summary
*label*.Gsea.*hash*/*gene_set*.xls              per gene set summary
!*label*.Gsea.*hash*/gene_set_sizes.xls         useless 
!*label*.Gsea.*hash*/ranked_gene_list_*.*.xls   ranked gene scores and other (big) information
```

and attempted to store them for reproducibility

```r
dir.create('output-minimal')

gsea_output_dirs <- dir('output', full = T) %>% dir(full = T)
gsea_output_dirs %>% str_replace('^output', 'output-minimal') %>% 
    paste0('/edb') %>% parallel::mclapply(dir.create, recursive = T)

parallel::mclapply(
    gsea_output_dirs,
    function(gsea_output_dir) {
        essential_files <- c(
            dir(gsea_output_dir, full = T, recursive = T, 'edb|rpt|rnk$'),
            paste0(gsea_output_dir, c('/index.html', '/xtools.css')),
            dir(gsea_output_dir, full = T, 'xls$') 
        )
        essential_files %>% {file.copy(., str_replace(., '^output', 'output-minimal'), T)}    
    }, mc.preschedule = F
) -> dev.null

# check for progress
gsea_output_dirs %>% str_replace('^output', 'output-minimal') %>% 
    paste0('/edb/results.edb') %>% {.[file.exists(.)]}
```

finally I gave up, and decided to only store `.rpt` file (especially for the 'rnd_seed' field).

```r
dir.create('rpt')

dir('output', recursive = T, '.rpt$', full = T) %>% 
    {file.copy(., str_replace_all(., '/', '_') %>% str_replace('output_', 'rpt/'))}
```

then compress `rpt/` into [rpt.tar.gz](output/rpt.tar.gz).



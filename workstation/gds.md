```r
dir('gsea', full.names = T) %>% file.remove

dir('raw', 'full') %>% str_replace('_full.soft.gz', '') %>% mclapply(
    . %>% {
        input  <- paste0('raw/', ., '_full.soft.gz')
        output <- paste0('gsea/', .)
        lem4::read_gds_dataset(input) %>% lem4::write_gsea_input(output, 'ANKLE2')
    }
)

dir('output', full.names = T) %>% unlink(recursive = T)
system('rm log/*')

gene_sets <- c(
    'c1.all',
    'c2.cp.biocarta',
    'c2.cp.kegg',
    'c2.cp.reactome',
    'c3.mir',
    'c3.tft',
    'c4.cgn',
    'c4.cm',
    'c5.cc',
    'c5.mf',
    'c6.all',
    'h.all',

    'c2.cgp'
)

gds_command <- . %>% paste0(
    'java -Xmx512m -cp /path/to/GSEA/gsea-3.0.jar xtools.gsea.Gsea', 
    ' -gmx /path/to/GSEA/msigdb_v6.1/msigdb_v6.1_GMTs/gene_set.v6.1.symbols.gmt', 
    ' -res gsea/', ., '.txt -cls gsea/', ., '.cls#ANKLE2 -chip gsea/', ., '.chip', 
    ' -out output/gene_set -collapse false -nperm 1000 -permute phenotype', 
    ' -rpt_label ', ., ' -metric Pearson -gui false &> log/gene_set_', ., ' &'
)


bash <- dir('gsea', 'txt') %>% str_replace('.txt', '') %>% gds_command %>% c('wait')
#bash <- dir('../../GDS', full.names = T) %>% str_extract('GDS\\d+') %>% gds_command %>% c('wait')


lapply(gene_sets, . %>% str_replace_all(bash, 'gene_set', .)) %>% unlist %>% 
    c('#!/bin/bash', .) %>% write_lines('temp.sh')





column_name <- . %>% read_lines %>% 
    {.[1:str_which(., '^!platform_table_begin$')]} %>% str_subset('^#');

dir('raw/GPL', full.names = T) %>% 
    {x <- mclapply(., column_name); names(x) <- .; x} %>% 
    write_rds('platform_column_name.rds')


dir('raw/GPL', full.names = T)[str_which(sapply(x, class), 'error')] %>% 
    mclapply(. %>% read_lines %>% str_subset('^#'))
```


```bash
ls raw/*full* | wc -l
ls gsea/*.txt | wc -l
ls -d output/GDS* | wc -l
```


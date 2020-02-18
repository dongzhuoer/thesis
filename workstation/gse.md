```
GSE54002: 
    c2.cgp+2048: T
    c5.all+2048: T
    c5.all+4096: T
    msigdb+4096: F
    msigdb+8192: T
```


```r
cores <- 24;
gse_command <- . %>% paste0(
    'java -Xmx4000m -cp /path/to/GSEA/gsea-3.0.jar xtools.gsea.Gsea', 
    ' -gmx /path/to/GSEA/msigdb_v6.1/msigdb_v6.1_GMTs/gene_set.v6.1.symbols.gmt', 
    ' -res gsea/', ., '.txt -cls gsea/', ., '.cls#ANKLE2 -chip gsea/', ., '.chip', 
    ' -out output/gene_set -collapse true -nperm 1000 -permute phenotype', 
    ' -rpt_label ', ., ' -metric Pearson -gui false &> log/gene_set_', .
)
raw_command <- dir('gsea', 'GSE\\d+.cls') %>% str_extract('GSE\\d+') %>% unique %>% gse_command

gene_set_core <- c('h.all', 'c1.all', 'c2.cp.biocarta', 'c2.cp.kegg', 
    'c2.cp.reactome', 'c3.mir', 'c3.tft', 'c4.cgn', 'c4.cm', 'c6.all'
)
gene_set_extra <- c('c5.all', 'c2.cgp', 'c7.all')

make_shell <- function(total_commands) {
    shell <- vector('list', cores)

    for (i in seq_len(cores)) {
        shell[[i]] = total_commands %>% {.[seq(i, length(.), cores)]}
    }

    shell
}

finished <- dir('output', full.names = T) %>% dir(full.names = T) %>% 
    paste0('/index.html') %>% {.[file.exists(.)]} %>% 
    str_remove('output/') %>% str_remove('.Gsea.\\d+/index.html') %>% 
    str_remove('error_') %>% str_replace('/', '_')

core_command  <- lapply(gene_set_core,  . %>% str_replace_all(raw_command, 'gene_set', .)) %>% 
    unlist %>% {sample(., length(.))} 
extra_command <- lapply(gene_set_extra, . %>% str_replace_all(raw_command, 'gene_set', .)) %>% 
    unlist %>% {sample(., length(.))} 

shells <- list(
    core_command %>% {.[!(str_extract(., '(?<=log/)[\\w\\.]+') %in% finished)]}, 
    extra_command %>% {.[!(str_extract(., '(?<=log/)[\\w\\.]+') %in% finished)]}
) %>% lapply(make_shell)

lapply(
    seq_len(cores), 
    function(i) {
        lapply(
            shells, 
            . %>% {.[[i]]}) %>% unlist %>% c('#!/bin/bash', .) %>% 
                write_lines(., paste0('shell/', i, '.sh')
        )
    }
) -> dev.null
```

- on local


```bash
R --slave -e "paste0('nohup shell/', 1:24, '.sh &> /dev/null &') %>% cat(sep = '\n')" | code -

rm -r output/*/*GSE*; rm log/*GSE*

while true; do echo; lem4 "cd /path/to/project; ls -d output/*/* | wc -l | sed 's/^/normal: /'"; lem4 "cd /path/to/project; ls -d output/*/error* | wc -l | sed 's/^/error: /'"; echo;  sleep 10; done
while true; do echo; echo total ----------------------------------; lem4 "cd /path/to/project; ls -d output/*/* | wc -l | sed 's/^/MYC: /'; ls -d output6/*/* | wc -l | sed 's/^/BRCA1: /'"; echo; echo error -----------------------------------; lem4 "cd /path/to/project; ls -d output5/*/error* | wc -l | sed 's/^/MYC: /'; ls -d output6/*/error* | wc -l | sed 's/^/BRCA1: /';"; echo;  sleep 10; done
```


```bash
setdiff2 $known_error $total_error 
setdiff2 <(lem4 "grep -P 'none' log/* | grep -oP 'log/[\w\.]+\w+'| sed 's/log\///'") <(lem4 "ls output/*/error* -d | grep -oP 'output/[\w\.]+/\w+' | sed 's/output\///' | sed 's/\/error//'")
```




```r
dir()
dir('output', full = T) 
dir('output', full = T) %>% sapply(. %>% dir %>% length)
```
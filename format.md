
# url

- GSE & GPL, `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=` + accession
- GDS, `https://www.ncbi.nglm.nih.gov/sites/GDSbrowser?acc=` + accession



# GDS raw data

1. DataSet full SOFT file: `*_full.soft`, 2 + 5
2. DataSet SOFT file: `*.soft`, raw data
3. Series family SOFT file: `*_family.*`, original submitter-supplied file
4. Series family MINiML file: same as above
5. Annotation SOFT file: `*.annot`, some annotation header, and some annotation column in data



# GSE raw data

- `*series_matrix.txt`, expression data
- `*family.soft`, platform_table contains annotation for probe (includes entrez id)



# GEO structure

1. GDS is a subset of well annotated GSE (~10%)

   some GDS can't be found in `geo` because they are not from human, such as GDS592
   	
   	
2. GSE can be SuperSeries, which is composed of several SubSeries 
   
   Some contains multiple data type, like GSE97484 (chip + ChIP-seq + Methylation); some contain multiple platform, like GSE4198. And it's impossible to know which sample belongs to which platform from meta data

3. Except for SuperSeries, every GSE should have only a GPL (platform)




# `gds_result.txt`

```

record 1

record 2

record 3
```

every record conform the following format,

```
\d+. $Title
$description
Organism:\t$species
Type:\t\t$type
*Platform*Sample*
FTP download: GEO[(***)] $ftp
DataSet|Series\t\tAccession: (GDS|GSE)\d+\tID: \d+
```

1. line number
   - 6, without "Platform" line (mostly `third-party reanalysis`, with one exception)
   - 7, most common case   
   - 8 or 9, optional `Project: ...` or `SRA Run Selector: ` line

2. for "Third-party reanalysis", the source data is the same

3. an examplpe of multi-species dataset

   ```
   14866. Niche modulated versus niche modulating genes in multiple myeloma
   (Submitter supplied) Background. Multiple myeloma (MM) cells depend on the bone marrow (BM) niche for growth and survival. However, the tumor genes regulated by the niche are largely unknown.  Design and Methods. BM aspiration samples were obtained from MM-patients with a high tumor load. Gene expression profile (GEP) was recorded immediately following aspiration and at subsequent time points. Identification of niche-regulated genes relied on spontaneous gene modulation following loss of niche regulation. more...
   Organism:	Homo sapiens
   Type:		Expression profiling by array
   Datasets: GDS4007 GDS4008 Platform: GPL6244 32 Samples
   FTP download: GEO (CEL) ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE36nnn/GSE36036/
   Series		Accession: GSE36036	ID: 200036036
   ```









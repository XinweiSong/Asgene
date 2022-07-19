
# AsgeneDB: A curated orthology arsenic metabolism gene database for metagenome annotation

<br> **Description**: A manually curated arsenic functional gene
database (AsgeneDB) and R package (Asgene package) are developed for
rapid and accurate metagenomic analysis.<br> <br> **Authors**: Xinwei
Song, Yongguan Zhu, Yongming Luo, Jianming Xu, Bin Ma\* <br>

## Overview

Arsenic (As) is a kind of toxic metal-like element widely distributed in
the world. To understand the microbial community of arsenic metabolism
in the environment, we developed a curated arsenic functional gene
database (AsgeneDB) covering five arsenic metabolic pathways (transport,
respiratory, reduction, oxidative and methylation processes), 59 arsenic
biotransformation functional gene families and 414773 representative
sequences. Here, protein sequences for As gene families were recruited
from multiple public databases such as UniProt, NCBI RefSeq, KEGG, COG,
eggNOG, arCOG and KOG. AsgeneDB covers 46 phyla and 1653 genera of
bacterial, archaea and fungi. It can quickly analyze the arsenic
metabolism and transformation function of microbial communities by
integrating multiple lineal homology databases with high specificity,
comprehensiveness, representativeness and accuracy. AsgeneDB and the
associated R Package will greatly promote the study of arsenic
metabolism in microbial communities in various environments.

## AsgeneDB Documentation

### Database files

##### 1. The database files user needs are built into the R package (Asgene). Therefore, when using Asgene R package to analyze metagenomic data, there is no need to download AsgeneDB separately.

##### 2. Example datasets are provided as input and output to help users better understand this package. The example datasets consists of three protein sequence files, stored as testdata1.rda, testdata2.rda, testdata3.rda.

##### 3. In order to facilitate users to use the AsgeneDB individually for personalized analysis, AsgeneDB can be downloaded from <https://data.cyverse.org/dav-anon/iplant/home/xinwei/AsgeneDB/AsgeneDB.zip>.

###### **AsgeneDB.zip details: includes 4 files**

1.  **AsgeneDB.fa**: Fasta format representative sequences obtained by
    clustering curated sequences at 100% sequence identity. This file
    can be used for “BLAST” searching arsenic genes in shotgun
    metagenomes.

2.  **asgene.map**: A mapping file that maps sequence IDs to gene names,
    only sequences belonging to arsenic gene families are included. This
    file is used to generate arsenic gene profiles from BLAST-like
    results against the database.

3.  **id\_gene\_tax\_pathway\_total.csv**: Species table of sequences in
    AsgeneDB. <br>

    **Columns included:**<br>

    1.  Gene name (colnames:gene) <br>
    2.  Corresponding arsenic metabolic pathway (colnames:pathway) <br>
    3.  Taxid (colnames:taxid) <br>
    4.  Protein ID (colnames:protein\_id) <br>
    5.  Kindom classification of species (colnames:kindom) <br>
    6.  Phylum classification of species (colnames:phylum) <br>
    7.  Class classification of species (colnames:class) <br>  
    8.  Order classification of species (colnames:order) <br>
    9.  Family classification of species (colnames:family) <br>
    10. Genus classification of species (colnames:genus) <br>
    11. Species classification of species (colnames:species) <br>

4.  **length.txt**: The file contains the length of amino acid sequences
    in AsgeneDB for standardizing arsenic gene abundance statistics.

## Dependent Tools

1.  **R Studio **

-   Depends: R ≥ 3.4.0
-   Imports: dplyr, seqinr

2.  **database searching tools:**<br>

-   usearch: <https://www.drive5.com/usearch/download.html>
-   diamond: <https://github.com/bbuchfink/diamond/releases>
-   blast: <https://ftp.ncbi.nlm.nih.gov/blast/executables/>

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Asgene Package

### Installation

You can install the development version of Asgene from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("XinweiSong/Asgene")
```

### Usage

**Description**:<br> we provide Asgene Package for metagenomic alignment
(nucleic acid or protein sequence), subsequent gene family abundance
statistics and sample abundance standardization. The database files user
needs are built into the Asgene. Therefore, users only need to choose a
database search tool according to their needs (e.g., USEARCH, BLAST and
DIAMOND) and input three parameters (e.g., working path, search
parameters of tool and filetype) to automatically analyze statistics and
output statistical results. Users can select gene abundance statistics
(Option: abundance) to normalize read counts per kilobase per million
reads (RPKM) to eliminate differences in sequencing depth and reference
sequence length between samples. In addition, if the user selects
functional species statistics (Option: taxonomy), the driveing species
of each arsenic metabolism gene at different classification levels in
the sample can be generated automatically.

### Example

This is a basic example which shows you how to use the package:

``` r
library(Asgene)
#
#Arsenic metabolism gene abundance analysis
Asgene(anlysis = "abundance", workdir = "./", method = "diamond", toolpath = "./", search_parameters = "-e 1e-4 -p 28 --query-cover 80 --id 50",seqtype = "nucl", filetype = "fasta", PE = TRUE , output = "./")
#Arsenic metabolism taxonomy analysis
Asgene(anlysis = "taxonomy", workdir = "./", method = "diamond", toolpath = "./", search_parameters = "-e 1e-4 -p 28 --query-cover 80 --id 50",seqtype = "nucl", filetype = "fasta",PE = TRUE, output = "./")
```

### Output

##### Output 1

<div class="figure">

<img src="/Users/xinweisong/Desktop/results1.png" alt="Output of As metabolic gene abundance analysis" width="20%" height="10%" />
<p class="caption">
Output of As metabolic gene abundance analysis
</p>

</div>

#### **NOTE:**

-   Before you begin, place AsgeneDB files in the current working path.
-   You need to place a tab-separated file in your working path that
    contains the sample name and the number of sequences. Note that file
    extensions should not be included here. For example:
    <https://data.cyverse.org/dav-anon/iplant/home/xinwei/AsgeneDB/sampleinfo.txt>

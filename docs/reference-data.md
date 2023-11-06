---
title: "Reference Data"
layout: single
permalink: /reference-data/
toc: true
toc_label: "On This Page"
sidebar:
  nav: "docs"
---

## Reference Data 

Reference data used by modules.

### Ensemble Data

For applications that require Ensemble data.  (Cache data used by VEP.)
```
/shared/biodata/ngs/Reference/Genomes/Homo_sapiens/Ensembl-110
```

### 10X Genomics
The Hutch has serveral packages from 10X Genomoics. 
 - Cell Ranger
 - Cellranger ATAC 
 - Cell Ranger ARC
 - Space Ranger
Refernce data for 10X Genomics tools is located
```
/shared/biodata/ngs/Reference/10X`
```

### Kraken2
Kraken 2 is the newest version of Kraken, a taxonomic classification system using exact k-mer matches to achieve high accuracy and fast classification speeds. This classifier matches each k-mer within a query sequence to the lowest common ancestor (LCA) of all genomes containing the given k-mer. The k-mer assignments inform the classification algorithm.

```
/shared/biodata/ngs/Reference/Kraken2
```

### STAR2
[STAR](https://github.com/alexdobin/STAR) aligns RNA-seq reads to a reference genome using uncompressed suffix arrays.
Modules: `STAR/2.7.10b-GCC-12.2.0`

```
/shared/biodata/ngs/Reference/Genomes/Homo_sapiens/GRCh38/Sequence/STAR2Index
```

### AlphaFold

AlphaFold is an AI system developed by DeepMind that predicts a protein’s 3D structure from its amino acid sequence.

```
/shared/biodata/ngs/Reference/protein/
```

### CHOPCHOP
CHOPCHOP (version 3) is a tool for selecting target sites for CRISPR/Cas9, CRISPR/Cpf1, CRISPR/Cas13 or NICKASE/TALEN-directed mutagenesis.

```
/shared/biodata/ngs/Reference/chopchop
```

### gatk (Broad Institute)  Funcotator Data
The Genome Analysis Toolkit or [GATK](https://www.broadinstitute.org/gatk/) is a software package developed 
at the Broad Institute to analyse next-generation resequencing data.

```
/shared/biodata/ngs/Reference/Funcotator
```

### Bowtie2 Index data

```
/shared/biodata/ngs/Reference/Genomes/Drosophila_melanogaster/UCSC/BDGP6_dm6/Sequence/Bowtie2Index
```

### Oncotator

```
/shared/biodata/ngs/Reference/Oncotator
```

### PCGR - Personal Cancer Genome Reporter
```
/shared/biodata/ngs/Reference/pcgr/data
```

### RSeQC (Picard-Style) Interval Files
```
/shared/biodata/ngs/Reference/RNASeQC/hg38_rRNA_intervals.list
```

### [SnpEff](http://pcingola.github.io/SnpEff/)

```
/shared/biodata/ngs/Reference/snpEff
```

### tcrdist3

tcrdist3 is an open-source python package that enables a broad array of T cell receptor sequence analyses.
There are many other applications for analyzing TCR data; DeepTCR, clusTCR, [pubtcrs](https://github.com/phbradley/pubtcrs).

```
/shared/biodata/ngs/Reference/tcrdist
```


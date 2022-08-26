---
title: "Vendor Software"
layout: single
permalink: /vendorsoftware/
toc: true
toc_label: "On This Page"
sidebar:
  nav: "docs"
---

### Software from Scientific Instrument Vendors 

Vendor supplied software available as modules.

### Oxford Nanopore
 - Nanopolish - Software package for signal-level analysis of Oxford Nanopore sequencing data.
 - fast5 - A lightweight C++ library for accessing Oxford Nanopore Technologies sequencing data.
 - ont-guppy-cpu  Guppy is a production basecaller provided by Oxford Nanopore
 - Python libraries
    - NanoComp - Compare multiple runs of long read sequencing data and alignments. Creates violin plots or box plots of length, quality and percent identity and creates dynamic, overlaying read length histograms and a cumulative yield plot.
    - NanoFilt - Filtering and trimming of Oxford Nanopore Sequencing data
    - nanoget - Functions to extract information from Oxford Nanopore sequencing data and alignments.
    - nanomath - Provides a few simple math and statistics functions for other scripts processing Oxford Nanopore sequencing data
    - NanoPlot - Plotting suite for Oxford Nanopore sequencing data and alignments
    - nanoplotter - Plotting functions of Oxford Nanopore sequencing data
    - NanoStat - Calculate statistics for Oxford Nanopore sequencing data and alignments
    - pauvre - Tools for plotting Oxford Nanopore and other long-read data.
    - Porechop - Porechop is a tool for finding and removing adapters from Oxford Nanopore reads. Adapters on the ends of reads are trimmed off, and when a read has an adapter in its middle, it is treated as chimeric and chopped into separate reads. Porechop performs thorough alignments to effectively find adapters, even at low sequence identity.


### PacBio
 - bam2fastx - Conversion of PacBio BAM files into gzipped fasta and fastq files, including splitting of barcoded dat.

### Illumina
 - blc2fastq
 - illumina-dump (SRA-Toolkit)
 - BaseSpace Command line tools
 - Strelka 2

### 10X Genomics
Refernce data for 10X Genomics tools is located `/shared/biodata/ngs/Reference/10X`
 - Cell Ranger
 - Cellranger ATAC 
 - Cell Ranger ARC
 - Space Ranger

## Open Source Tools
Open source tools for reading sequence data from 
 - SRA Toolkit - The Sequence Read Archive (SRA Toolkit) stores raw sequence data from "next-generation" sequencing technologies
 - Strelka - Strelka2 is a fast and accurate small variant caller
 optimized for analysis of germline variation in small cohorts and
 somatic variation in tumor/normal sample pairs.
 - Unicycler - Unicycler is an assembly pipeline for bacterial genomes. It can assemble Illumina-only read sets where it functions as a SPAdes-optimiser. It can also assembly long-read-only sets (PacBio or Nanopore) where it runs a miniasm+Racon pipeline. For the best possible assemblies, give it both Illumina reads and long reads, and it will conduct a hybrid assembly.

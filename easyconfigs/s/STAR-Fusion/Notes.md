# STAR-Fusion

### Data Sources
STAR Fusion uses CTAT genome lib for human fusion transcript detection. Or users
can build their own reference library. The CTAT data is
from the [Trinity Cancer Transcriptome Analysis Toolkit](https://github.com/NCIP/Trinity_CTAT/wiki)

The CTAT can be downloaded from here: `https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/`
The reference files are in the format: `CTAT_resource_lib.tar.gz`
The 'plug-n-play' libs are that, just download, unpack it (tar -zxvf filename.tar.gz)
The `source` downloads takes hours to days to build. 

As of Aug, 2023 the CTAT reference data is located here /fh/scratch/app/CTAT.

Star Fusion does not use an environment variable to locate the CTAT data.

Command Line Argument to locate CTAT data
```
--genome_lib /path/to/ctat_genome_lib_build_dir
```



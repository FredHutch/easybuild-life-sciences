#!/bin/bash

# manually load the Kraken2 databases after Kraken2 is installed
set -e

source /app/Lmod/lmod/lmod/init/profile
module load Kraken2

KRAKEN2_DB_PATH=/shared/biodata/microbiome/kraken2
DATABASES='viral bacteria human nr nt'
DATABASES='nt'

for DBNAME in $DATABASES; do
    echo installing $DBNAME
    kraken2-build --download-library archaea --db $DBNAME
done

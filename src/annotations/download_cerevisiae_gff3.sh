#!/bin/bash
DIR=.
wget -nc -O $DIR/Saccharomyces_cerevisiae.R64-1-1.105.gff3.gz http://ftp.ensembl.org/pub/release-105/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.105.gff3.gz
gzip -dfc $DIR/Saccharomyces_cerevisiae.R64-1-1.105.gff3.gz > $DIR/Saccharomyces_cerevisiae.R64-1-1.105.gff3

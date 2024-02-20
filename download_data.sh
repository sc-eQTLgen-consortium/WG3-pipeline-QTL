#!/bin/bash

which wget >/dev/null || { echo "wget: command not found"; exit 1; }
which md5sum >/dev/null || { echo "md5sum: command not found"; exit 1; }

mkdir -p data
cd data || exit

echo "Ribosomal_genes.txt"
wget https://raw.githubusercontent.com/sc-eQTLgen-consortium/WG1-pipeline-QC/master/Demultiplexing/Ribosomal_genes.txt \
  && md5sum -c - <<<"1f232f4815abe534768a32e3203600a3  Ribosomal_genes.txt"

echo "Mitochondrial_genes.txt"
wget https://raw.githubusercontent.com/sc-eQTLgen-consortium/WG1-pipeline-QC/master/Demultiplexing/Mitochondrial_genes.txt \
  && md5sum -c - <<<"c5f65ae3bce28b9340615d605b2bb792  Mitochondrial_genes.txt"
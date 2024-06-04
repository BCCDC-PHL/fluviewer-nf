#!/bin/bash

mkdir -p .github/data/assemblies

curl -o .github/data/assemblies/NC_026423.1.fa "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=NC_026423.1&db=nucleotide&rettype=fasta"
curl -o .github/data/assemblies/NC_026431.1.fa "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=NC_026431.1&db=nucleotide&rettype=fasta"

#!/bin/bash

set -eo pipefail

sed -i 's/cpus = 8/cpus = 4/g' nextflow.config
sed -i "s/memory = '32 GB'/memory = '2 GB'/g" nextflow.config 

nextflow run main.nf \
	 -profile conda \
	 --cache ${HOME}/.conda/envs \
	 --fastq_input .github/data/fastq \
	 --outdir .github/data/test_output

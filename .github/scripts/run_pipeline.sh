#!/bin/bash

set -eo pipefail

sed -i 's/cpus = 8/cpus = 4/g' nextflow.config
sed -i "s/memory = '32 GB'/memory = '2 GB'/g" nextflow.config 

nextflow run main.nf \
	 -profile conda \
	 --cache ${HOME}/.conda/envs \
	 --fastq_input .github/data/fastq \
	 --db .github/data/fluviewer_db/FluViewer_db_v_0_1_8.fa \
	 --outdir .github/data/test_output
#!/bin/bash

set -eo pipefail

source ${HOME}/.bashrc
 
eval "$(conda shell.bash hook)"

conda activate base

# Check for a sign that we're in the GitHub Actions environment.
# Prevents these settings from being applied in other environments.
if [ -n "${GITHUB_ACTIONS}" ]; then 
    sed -i 's/cpus = 8/cpus = 4/g' nextflow.config
    sed -i '/memory/d' modules/FluViewer.nf
fi
nextflow run main.nf \
	 -profile conda \
	 --cache ${HOME}/.conda/envs \
	 --fastq_input .github/data/fastq \
	 --db .github/data/fluviewer_db/FluViewer_db_v_0_1_8.fa \
	 --outdir .github/data/test_output

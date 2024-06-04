#!/bin/bash

mkdir -p apptainer-build
pushd apptainer-build > /dev/null

for env_yaml in ../environments/*.yml; do
    image_name=$(head -n 1 $env_yaml | cut -d ' ' -f 2)
    echo "building image ${image_name} from file ${env_yaml}..."
    cp ${env_yaml} ./environment.yml
    apptainer build ${image_name}.sif ../.github/data/conda-env.def
    ls -lh ${image_name}.sif
done



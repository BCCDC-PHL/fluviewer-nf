#!/bin/bash

for image_file in ./apptainer-build/*.sif; do
    echo "Pushing image ${image_file}..."
    IMAGE_NAME=$(basename ${image_file} .sif)
    echo IMAGE_NAME=$IMAGE_NAME
    IMAGE_ID=ghcr.io/${GITHUB_REPOSITORY_OWNER}/${IMAGE_NAME}
    IMAGE_ID=$(echo $IMAGE_ID | tr '[A-Z]' '[a-z]')
    VERSION=$(echo "${GITHUB_REF}" | sed -e 's,.*/\(.*\),\1,')
    [ "$VERSION" == "main" ] && VERSION=latest
    echo "Pushing $image to $IMAGE_ID:$VERSION..., using username ${GHCR_USERNAME}"
    echo ${GITHUB_TOKEN} | apptainer registry login -u ${GHCR_USERNAME} --password-stdin oras://ghcr.io
    apptainer push ${image_file} oras://${IMAGE_ID}:${VERSION}
done

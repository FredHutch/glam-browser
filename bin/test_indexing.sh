#!/bin/bash

set -e

# Run this from the base level of the glam-browser repository

GENESHOT_VER=v0.6.3
TAR_GZ_URL=https://github.com/Golob-Minot/geneshot/releases/download/${GENESHOT_VER}/geneshot.results.tar.gz

# Get the testing data
if [[ ! -s testing/output/geneshot.results.hdf5 ]]; then

    [[ ! -d testing ]] && \
    mkdir testing

    [[ ! -s testing/geneshot.results.tar.gz ]] && \
    curl -s -L ${TAR_GZ_URL} --output testing/geneshot.results.tar.gz 

    tar -xvzf testing/geneshot.results.tar.gz -C testing

fi

# Index the testing data
NXF_VER=20.04.1 \
nextflow \
    run \
    index.nf \
    -profile testing \
    --input testing/geneshot.results.hdf5 \
    --output_folder testing/ \
    --output_prefix glam \
    -with-docker ubuntu:latest

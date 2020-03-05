# Gene-Level Association of Microbiomes (GLAM) - Browser

The GLAM Browser is intended to make it easier to visualize and understand the
results of a microbiome experiment which has been analyzed by whole-genome
shotgun (WGS) sequencing and analyzed with gene-level metagenomic analysis. The
input files for the GLAM Browser are the output files from the `geneshot` analysis
tool ([link](https://github.org/golob-minot/geneshot)).

The GLAM Browser displays:

  * Table of CAG associations (output from corncob)
  * Volcano plot (estimated coefficient vs. -log10(p-value))
  * Table of genes found in a selected CAG
  * Customizable display of CAG abundance across samples

### Running the GLAM Browser

First, build the Docker image for the GLAM Browser:

```#!/bin/bash
docker build . -t glam
```

To run the app, you need to point it to a directory which contains the needed HDF5
files using the `DATA_DIR` environment variable in the container. This requires two
steps:

  1. Mount your local data folder to the container
  2. Point the app to read data from that mounted folder

You can set this to any path on your system, but for the example we will use `$PWD/data/`.

```#!/bin/bash
docker run -v $PWD/data:/share/data -e DATA_DIR=/share/data -p 3838:7777 glam
```

The app should be accessible on the 3838 port, at `localhost:3838`.

#### Development Mode

If you would like to run a version of the app which you are actively editing, such
as a copy in the `$PWD/app` folder, then use the following command.

```#!/bin/bash
docker run -v $PWD/data:/share/data -e DATA_DIR=/share/data -v $PWD/app:/share/app -e APP_DIR=/share/app -p 3838:7777 glam
```

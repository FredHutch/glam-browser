# Gene-Level Association of Microbiomes (GLAM) - Browser

The GLAM Browser is intended to make it easier to visualize and understand the
results of a microbiome experiment which has been analyzed by whole-genome
shotgun (WGS) sequencing and analyzed with gene-level metagenomic analysis. The
input files for the GLAM Browser are the output files from the `geneshot` analysis
tool ([link](https://github.org/golob-minot/geneshot)).

The GLAM Browser displays:

  * Summary of specimens by the number of genes detected and reads aligned
  * Description of CAG metrics including size, entropy, etc.
  * Summary of user-defined association testing
  * Taxonomic summary for selected CAGs
  * Customizable abundance display for selected CAGs

The input to the Browser is the single results HDF5 file output by `geneshot`.

## Running the GLAM Browser

The GLAM Browser can be run in a few different ways. Depending on your preferred method, you will have to make sure to make the input data accessible in a compatible manner.

### Run Locally with Docker

The most straightforward for a user to run the GLAM Browser is to run a Docker container on your local computer. The two things to think about are that (1) the Docker container needs to be able to access your input data (with a combination of mounting the folder and setting the `HDF5_FP` environment variable inside the container), and (2) it also needs to be able to expose the port (`8050`) which the app is hosted on.

1. Check for the latest tagged Docker image at [Quay](https://quay.io/repository/fhcrc-microbiome/glam?tab=tags). We use `latest` in the example below.
2. Find the absolute path to your input HDF5 on your local filesystem. In this example we will assume that the file is found at `/path/to/folder/containing/data/results.hdf5`
3. Set up [Docker Desktop](https://www.docker.com/products/docker-desktop) on your computer
4. Run the Docker container as follows:

```#!/bin/bash

docker run -it -v /path/to/folder/containing/data/:/share --env HDF5_FP=/share/results.hdf5 -p 8050:8050 quay.io/fhcrc-microbiome/glam:latest
```

You can now access the app at `127.0.0.1:8050` in your browser.

Explanation of flags:

* `-it`: Run in interactive mode (issuing Control+C to close)
* `-v`: Mount the data-containing folder to `/share` in the container
* `--env`: Set an environment variable in the container
* `-p`: Expose the port used by GLAM

### Run Locally from Source

You can also run the app by:

1. Cloning the GitHub repo locally
2. Setting up a Python3 virtual environment
3. Installing the Python requirements `pip3 install -r requirements.txt`
4. Running the app with the following command

```#!/bin/bash
HDF5_FP=/path/to/folder/containing/data/results.hdf5 python3 app.py
```

Because of the challenges of establishing a common working environment, this local run option is less preferred than using Docker.

### Web Deployment

The GLAM Browser is set up as a Dash app, and can be deployed as such following the instructions laid out [in the Dash documentation.](https://dash.plotly.com/deployment)

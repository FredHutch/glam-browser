# Gene-Level Association of Microbiomes (GLAM) - Browser

The GLAM Browser is intended to make it easier to visualize and understand the
results of a microbiome experiment which has been analyzed by whole-genome
shotgun (WGS) sequencing and analyzed with gene-level metagenomic analysis. The
input file for the GLAM Browser is the output file from the `geneshot` analysis
tool ([link](https://github.org/golob-minot/geneshot)).

The GLAM Browser displays:

  * Summary of specimens by the number of genes detected and reads aligned
  * Description of CAG metrics including size, entropy, etc.
  * Summary of user-defined association testing
  * Taxonomic summary for selected CAGs
  * Customizable abundance display for selected CAGs
  * ... and more

In order to view a set of `geneshot` results, you must first generate a GLAM index
file using the `index.nf` utility in this repository. That `*.glam.hdf5` file may
then be viewed in the GLAM Browser by placing it in a single folder and optionally
providing a `manifest.json` file in that folder which adds a name and description
for each dataset.

## Generating the GLAM Index File

To process a single `geneshot` output file to produce the GLAM index file, use the
`index.nf` utility as follows (e.g., on `geneshot.results.hdf5`):

```#!/bin/bash
nextflow run FredHutch/glam-browser \
    --input geneshot.results.hdf5 \
    --output glam.index.hdf5
```

Run `nextflow run FredHutch/glam-browser --help` for a complete list of options.

## Hosted at Fred Hutch Cancer Research Center

Researchers at Fred Hutch have access to a hosted version of the GLAM Browser at
[https://glam-mri.fredhutch.org/](https://glam-mri.fredhutch.org/).
You must be connected to the Fred Hutch VPN to access this site.
Please contact the maintainers of this repository with any questions.

## Running the GLAM Browser

The GLAM Browser can be run in a few different ways.
As part of starting the GLAM Browser, you will be instructed to start a redis server. This database is used to cache data for rapid loading.
To run the GLAM Browser on your personal computer or laptop, follow the instructions to run locally with Docker.
To run on an institutional high-performance computing cluster (HPC), run locally with Singularity.
For other use cases, you may choose to run directly from source.
Built using the Dash framework, the GLAM Browser can be deployed using a variety of approaches which are outlined in community documentation linked below.

### Run Locally with Docker Compose

The most straightforward for a user to run the GLAM Browser is to run a set of Docker containers on your local computer. We will use Docker Compose, which will start GLAM and redis each in their own container and network them together.

1. Clone this GitHub repository to your computer and navigate to that folder.
2. Make a folder called `data/` and place your input data files (and optional `manifest.json`) in that folder.
3. Set up [Docker Desktop](https://www.docker.com/products/docker-desktop) on your computer
4. Run the Docker containers as follows:

```#!/bin/bash

docker-compose -f public-docker-compose.yml up
```

Once the app loads (usually takes <1min) you can access the app at `0.0.0.0:8050` in your browser.

### Run Locally with Singularity

It is common for researchers to not have access to a computing system which provides access to Docker 'containers' (or 'images') as well as a generous amount of memory and storage. Instead, many HPC systems are now providing access to 'containers' using Singularity, which has the strong benefit of not needing to be run as root (which would be a security flaw for a shared HPC).

NOTE: If you are running Singularity at Fred Hutch, do not use the shared `rhino` system, instead use `grabnode` to log into a dedicated worker node.

Before starting the app, make a folder in your working directory to contain all of the geneshot files to display. In this example that folder will be named `$PWD/data`. Next copy the redis configuration file from this repository to a local file called `redis.conf`. For example:

```
maxmemory 4000mb
maxmemory-policy allkeys-lru
maxmemory-samples 5
```

Next, start two containers -- one for redis which runs in the background and another for GLAM.

```#!/bin/bash

# Start the redis server
SINGULARITY_CACHEDIR=$PWD/cache/ singularity run --bind $PWD:/share docker://redis:6 redis-server /share/redis.conf &

SINGULARITY_CACHEDIR=$PWD/cache/ singularity run --bind $PWD/data:/share docker://quay.io/fhcrc-microbiome/glam:latest

```

Explained:

* `SINGULARITY_CACHEDIR`: Sets the directory where Singularity container files are written
* `--bind $PWD:/share`: Mounts the local directory as `/share` inside the container
* `quay.io/fhcrc-microbiome/glam:latest`: Docker container to be run from Singularity. Replace `latest` with any tag of interest if you want to run a specific version of GLAM.

Once GLAM has started, you will be able to access it from a web browser with the name of the compute node with `:8050`. For example, if I am logged into `gizmof13` when running the command above, then I can point my browser to `gizmof13:8050` to access the browser.

### Run Locally from Source

You can also run the app by:

1. Cloning the GitHub repo locally
2. Setting up a Python3 virtual environment
3. Installing the Python requirements `pip3 install -r requirements.txt`
4. Running the app with the following command

```#!/bin/bash
DATA_FOLDER=/path/to/folder/containing/data/ python3 app.py
```

Because of the challenges of establishing a common working environment, this local run option is less preferred than using Docker.

## Naming and Describing Datasets

The GLAM Browser natively displays the results from multiple datasets. To allow the user to describe and name those datasets in a completely flexible manner, the user can add a `manifest.json` file to the data directory (right next to the geneshot output files) which contains some descriptions they might find useful. The format of that manifest is as follows:

```json
{
    "page_title": "Title at the top of the page",
    "page_description": "Text with more room to describe the collection of data which is found in this folder.",
    "contents": [
        {
            "fp": "name-of-file-within-folder.hdf5",
            "name": "Short name for dataset",
            "description": "Longer description of this dataset."
        },
        {
            "fp": "name-of-second-file-within-folder.hdf5",
            "name": "Short name for second dataset",
            "description": "Longer description of this second dataset."
        }
    ]
}
```

### Web Deployment

The GLAM Browser is set up as a Dash app, and can be deployed as such following the instructions laid out [in the Dash documentation.](https://dash.plotly.com/deployment)

## GLAM Screenshots

The number of genes detected for each sample is displayed as a scatter plot or histogram, as well as the proportion of reads which align uniquely to a single gene. Individual specimens can be masked from this plot (as well as any plot showing specimens) using the interactive manifest table at the bottom of GLAM.

![Gene Detection](https://github.com/FredHutch/glam-browser/blob/master/assets/richness_example.png?raw=true)

The similarity of samples is displayed as PCA or t-SNE plots, with options to overlay user-defined metadata.

![Ordination](https://github.com/FredHutch/glam-browser/blob/master/assets/ordination_example.png?raw=true)

The characteristics of CAGs, including size, entropy, etc., are displayed as a scatter plot or histogram. Clicking on any point in this scatter plot will 'select' a CAG in multiple plots in GLAM.

![CAG Summary](https://github.com/FredHutch/glam-browser/blob/master/assets/cag_summary_example.png?raw=true)

If the user provided a `--formula` when running `geneshot`, the results of that association analysis will be presented as a volcano plot. Clicking on any point in this scatter plot will 'select' a CAG in multiple plots in GLAM. Note that masking specimens in the manifest table will _not_ change the data in this plot.

![Volcano](https://github.com/FredHutch/glam-browser/blob/master/assets/volcano_example.png?raw=true)

The distribution of taxonomic annotations across genes is shown for the selected CAG.

![Taxonomy](https://github.com/FredHutch/glam-browser/blob/master/assets/taxonomy_example.png?raw=true)

The abundance of a single CAG across samples can be displayed in many different ways, including overlaid metadata provided by the user in the manifest.

![Single CAG](https://github.com/FredHutch/glam-browser/blob/master/assets/single_cag_example.png?raw=true)

The complete list of specimens in the manifest is provided as an interactive table. The user can filter and sort the table in a fairly complete manner. Crucially, by deselecting the checkbox on any row, that specimen is temporarily masked from all of the other displays in GLAM.

![Manifest](https://github.com/FredHutch/glam-browser/blob/master/assets/manifest_example.png?raw=true)

### Details

## Details on GLAM indexing

The `index.nf` utility creates an HDF5 file by (1) running the `index.py` script on the `geneshot` results file to create a new HDF5 file and then (2) repacking that HDF5 to reduce storage size. The tables in the GLAM index are:

* `/manifest`: A subset of the original manifest file, but with FASTQ input paths removed (`R1`, `R2`, `I1`, and `I2`) explicitly requiring that every `specimen` only has a single entry (row) in the table.
* `/experiment_metrics`: A direct copy of `/summary/experiment` in the input, containing information on the parameters used to run the experiment.
* `/specimen_metrics`: A direct copy of `/summary/all` table in the input, but with the `prop_reads` column added (`aligned_reads / n_reads`). Contains information on the number of genes aligned, number of genes assembled, number of reads, etc. for each specimen in the analysis.
* `/analysis_metrics`: A small summary of the analysis, including flags for whether there was taxonomic analysis, functional analysis, corncob analysis, and betta enrichment analysis included in the outputs. This makes it easier to execute logic to only display those results which are available.
* `/cag_annotations`: A direct copy of `/annot/cag/all` in the input, but with the `size_log10` column added (`log10(size)`).
* `/cag_abundances`: A direct copy of (`/abund/cag/wide`) in the input.
* `/distances/*`: A direct copy of all of the distance matrices in the input.
* `/cag_associations/results/{parameter}`: For each of the `parameter` values in the input (from `/stats/cag/corncob`), a table with those corncob results.
* `/cag_associations/manifest`: A small table with the list of parameters which were computed in the corncob analysis.
* `/enrichments/results/{parameter}/{annotation_group}`: For each of the `parameter` values and `annotation_group` values in the input (from `/stats/enrichment/betta`), a table with those functional enrichment results.
* `/enrichments/manifest`: A small table with the list of parameters and annotation groups which were computed in the analysis.
* `/gene_annotations/detailed`: A direct copy of `/annot/gene/all` in the input. This is a very large table, and so we will make sure to index by the `CAG` column. GLAM will only read in subsets of this table for a particular CAG at a time.
* `/gene_annotations/summary/parameter/{}`: For each `parameter` in the analysis, the table will include a summary of the number of genes which had a given annotation (taxonomic and functional) for each CAG. Notably, this will only include those CAGs which pass the FDR threshold of 0.01 for that parameter. By limiting the number of CAGs in this table, we will be able to read the complete table into GLAM and keep it in memory to render figures from.

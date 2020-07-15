#!/usr/bin/env nextflow

/*
  GLAM Browser for the Gene-Level Analysis of Metagenomes

  This utility can be used to process a set of geneshot results
  for visualization in the GLAM Browser.

  By preprocessing the geneshot outputs, we remove the need for
  server-side computation during callbacks, except as needed to
  render individual plots.
*/

// Using DSL-2
nextflow.preview.dsl=2

// Default values for flags
// If these are not set by the user, then they will be set to the values below
// This is useful for the if/then control syntax below
params.input = false
params.output = false
params.help = false

// Set the containers to user
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/glam-browser <ARGUMENTS>
    
    Arguments:
      --input               Geneshot results (".hdf5") file to process
      --output_folder       Folder for output GLAM index
      --output_prefix       Prefix used to name the GLAM index (".index.hdf5" will be appended)

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.input == false || params.output_folder == false || params.output_prefix == false){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Make CAGs for each set of samples, with the subset of genes for this shard
process indexGeneshotResults {
    container "${container__pandas}"
    label "mem_veryhigh"

    input:
    path input_hdf

    output:
    file "${params.output_prefix}.index.hdf5"

    """#!/bin/bash

index.py ${input_hdf} ${params.output_prefix}.index.hdf5

    """
}

// Repack an HDF5 file
process repackHDF {

    container "${container__pandas}"
    label "mem_medium"
    publishDir "${params.output_folder}", mode: "copy", overwrite: true
    
    input:
    file output_hdf5
        
    output:
    file "${output_hdf5}"

    """
#!/bin/bash

set -e

[ -s ${output_hdf5} ]

h5repack -f GZIP=5 ${output_hdf5} TEMP && mv TEMP ${output_hdf5}
    """
}


workflow {

  // Index the input file
  indexGeneshotResults(
    Channel.fromPath(params.input)
  )

  // Repack the output file
  repackHDF(
    indexGeneshotResults.out
  )

}

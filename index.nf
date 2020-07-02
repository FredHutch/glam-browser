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

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/glam-browser <ARGUMENTS>
    
    Arguments:
      --input               Geneshot results (".hdf5") file to process
      --output              Location for output GLAM index (".hdf5")

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.input == false || params.output == false){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}


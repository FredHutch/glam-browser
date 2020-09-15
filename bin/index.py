#!/usr/bin/env python3

# Script used to process a set of geneshot results for visualization in the GLAM Browser
import argparse
from collections import defaultdict
from functools import lru_cache
import h5py
import logging
import os
import pandas as pd
import numpy as np


class Taxonomy:

    def __init__(self, store):
        """Read in the taxonomy table."""

        # Read the taxonomy table
        logging.info("Reading in /ref/taxonomy")

        self.taxonomy_df = pd.read_hdf(
            store,
            "/ref/taxonomy"
        ).apply(
            lambda c: c.fillna(0).apply(float).apply(
                int) if c.name in ["parent", "tax_id"] else c,
        ).set_index(
            "tax_id"
        )

        self.all_taxids = set(self.taxonomy_df.index.values)

    @lru_cache(maxsize=None)
    def path_to_root(self, tax_id, max_steps=100):
        """Parse the taxonomy to yield a list with all of the taxa above this one."""

        visited = []

        for _ in range(max_steps):

            # Skip taxa which are missing
            if tax_id not in self.all_taxids:
                break

            # Add to the path we have visited
            visited.append(tax_id)

            # Get the parent of this taxon
            parent_id = self.taxonomy_df.loc[tax_id, "parent"]

            # If the chain has ended, stop
            if parent_id in visited or parent_id == 0:
                break

            # Otherwise, keep walking up
            tax_id = parent_id

        return visited

    @lru_cache(maxsize=None)
    def anc_at_rank(self, tax_id, rank):
        for anc_tax_id in self.path_to_root(tax_id):
            if self.taxonomy_df.loc[anc_tax_id, "rank"] == rank:
                return anc_tax_id

    def name(self, tax_id):
        if tax_id in self.all_taxids:
            return self.taxonomy_df.loc[tax_id, "name"]

    def parent(self, tax_id):
        if tax_id in self.all_taxids:
            return self.taxonomy_df.loc[tax_id, "parent"]

    def make_cag_tax_df(
        self,
        taxa_vc, 
    ):
        """Return a nicely formatted taxonomy table from a list of tax IDs and the number of assignments for each."""

        # We will construct a table with all of the taxa in the tree, containing
        # The ID of that taxon
        # The name of that taxon
        # The number of genes assigned to that taxon (or its children)
        # The rank of that taxon
        # The ID of the parent of that taxon

        # The number of genes found at that taxon or in its decendents
        counts = defaultdict(int)

        # Keep track of the total number of genes with a valid tax ID
        total_genes_assigned = 0

        # Iterate over each terminal leaf
        for tax_id, n_genes in taxa_vc.items():

            # Skip taxa which aren't in the taxonomy
            if tax_id not in self.taxonomy_df.index.values:
                continue

            # Count all genes part of this analysis
            total_genes_assigned += n_genes

            # Walk up the tree from the leaf to the root
            for anc_tax_id in self.path_to_root(tax_id):

                # Add to the sum for every node we visit along the way
                counts[anc_tax_id] += n_genes

        if len(counts) == 0:
            return pd.DataFrame([{
                "tax_id": None
            }])

        # Make a DataFrame
        df = pd.DataFrame({
            "count": counts,
        })

        # Add the name, parent, rank
        df = df.assign(
            tax_id=df.index.values,
            parent=self.taxonomy_df["parent"],
            rank=self.taxonomy_df["rank"],
            name=self.taxonomy_df["name"],
            total=total_genes_assigned,
        )

        return df


def parse_manifest(store):
    """Read in the manifest and filter columns for visualization."""

    # Read the whole manifest
    logging.info("Reading in /manifest")
    manifest = pd.read_hdf(store, "/manifest")
    logging.info("Read in data from {:,} rows and {:,} columns".format(
        manifest.shape[0],
        manifest.shape[1],
    ))

    # Make sure that we have a specimen column
    assert "specimen" in manifest.columns.values

    # Drop any columns with paths to reads
    for k in ["R1", "R2", "I1", "I2"]:
        if k in manifest.columns.values:
            logging.info("Removing column {}".format(k))
            manifest = manifest.drop(columns=k)

    # Drop any columns which do not have unique values for each specimen
    for k in manifest.columns.values:
        if k == "specimen":
            continue
        # Look at just the unique values for this column
        d = manifest.reindex(columns=["specimen", k]).drop_duplicates()
        # Make sure that every specimen has only one value in this column
        if d["specimen"].value_counts().max() > 1:
            logging.info("Removing column with duplicate values: {}".format(k))
            manifest = manifest.drop(columns=k)

    # Now drop duplicates and make sure the specimen column is still unique
    manifest = manifest.drop_duplicates()
    assert manifest["specimen"].value_counts().max() == 1

    # Return the filtered manifest
    return manifest


def parse_specimen_metrics(store):
    """Read in the specimen metrics"""
    key_name = "/summary/all"

    logging.info("Reading in {}".format(key_name))

    df = pd.read_hdf(store, key_name)

    # Compute `prop_reads`
    return df.assign(
        prop_reads = df["aligned_reads"] / df["n_reads"]
    )


def parse_cag_annotations(store, max_n_cags=250000):
    """Read in the CAG annotations."""
    key_name = "/annot/cag/all"

    logging.info("Reading in {}".format(key_name))

    df = pd.read_hdf(store, key_name)

    # Store information for no more than `max_n_cags`
    if max_n_cags is not None:
        df = df.head(
            max_n_cags
        )

    # Compute `prop_reads`
    return df.assign(
        size_log10 = df["size"].apply(np.log10)
    ).drop(
        columns=["std_abundance"]
    )


def parse_gene_annotations(store, tax, max_n_cags=250000):
    """Make a summary of the gene-level annotations for this subset of CAGs."""
    key_name = "/annot/gene/all"

    logging.info("Reading in {}".format(key_name))

    df = pd.read_hdf(store, key_name)

    logging.info("Read in {:,} annotations for {:,} CAGs".format(
        df.shape[0],
        df["CAG"].unique().shape[0],
    ))

    # Store information for no more than `max_n_cags`
    if max_n_cags is not None:
        logging.info("Limiting annotations to {:,} CAGs".format(max_n_cags))
        df = df.query(
            "CAG < {}".format(max_n_cags)
        )
        logging.info("Filtered down to {:,} annotations for {:,} CAGs".format(
            df.shape[0],
            df["CAG"].unique().shape[0],
        ))

    # Trim the `eggNOG_desc` to 100 characters, if present
    df = df.apply(
        lambda c: c.apply(lambda n: n[:100] if isinstance(n, str) and len(n) > 100 else n) if c.name == "eggNOG_desc" else c
    )

    # Summarize the number of genes with each functional annotation, if available
    if "eggNOG_desc" in df.columns.values:
        logging.info("Summarizing functional annotations")
        functional_df = summarize_annotations(df, "eggNOG_desc")
    else:
        functional_df = None

    # Summarize the number of genes with each taxonomic annotation, if available
    # This function also does some taxonomy parsing to count the assignments to higher levels
    if "tax_id" in df.columns.values and tax is not None:
        logging.info("Summarizing taxonomic annotations")
        counts_df, rank_summaries = summarize_taxonomic_annotations(df, tax)
    else:
        counts_df, rank_summaries = None, None

    return functional_df, counts_df, rank_summaries


def summarize_taxonomic_annotations(df, tax):
    # Make a DataFrame with the number of taxonomic annotations per CAG
    counts_df = pd.concat([
        tax.make_cag_tax_df(
            cag_df["tax_id"].dropna().apply(int).value_counts()
        ).assign(
            CAG = cag_id
        )
        for cag_id, cag_df in df.reindex(
            columns=["CAG", "tax_id"]
        ).groupby(
            "CAG"
        )
    ]).reset_index(
        drop=True
    )

    logging.info("Made tax ID count table")

    rank_summaries = {
        rank: counts_df.query(
            "rank == '{}'".format(rank)
        ).drop(
            columns=["rank", "parent"]
        )
        for rank in ["phylum", "class", "order", "family", "genus", "species"]
    }

    return counts_df, rank_summaries


def summarize_annotations(df, col_name):
    """Count up the unique annotations for a given CAG."""
    assert col_name in df.columns.values, (col_name, df.columns.values)
    assert "CAG" in df.columns.values, ("CAG", df.columns.values)

    return pd.DataFrame([
        {
            "label": value,
            "count": count,
            "CAG": cag_id
        }
        for cag_id, cag_df in df.groupby("CAG")
        for value, count in cag_df[col_name].dropna().value_counts().items()
    ])

def parse_cag_abundances(store, max_n_cags=250000):
    """Read in the CAG abundances."""
    key_name = "/abund/cag/wide"

    logging.info("Reading in {}".format(key_name))

    df = pd.read_hdf(store, key_name)

    # Store information for no more than `max_n_cags`
    if max_n_cags is not None:
        df = df.head(
            max_n_cags
        )

    return df


def parse_genome_containment(store, max_n_cags=250000, min_cag_prop=0.25):
    """Read in a summary of CAGs aligned against genomes."""
    key_name = "/genomes/cags/containment"

    if key_name in store:

        logging.info("Reading in {}".format(key_name))

        df = pd.read_hdf(store, key_name)

        # Store information for no more than `max_n_cags`
        if max_n_cags is not None:
            logging.info("Subsetting to {:,} CAGs".format(max_n_cags))
            df = df.query(
                "CAG < {}".format(max_n_cags)
            )
            logging.info("Retained {:,} alignments for {:,} CAGs".format(
                df.shape[0],
                df["CAG"].unique().shape[0]
            ))

        # Filter alignments by `min_cag_prop`
        if min_cag_prop is not None:
            logging.info("Filtering to cag_prop >= {}".format(min_cag_prop))
            df = df.query("cag_prop >= {}".format(min_cag_prop))
            logging.info("Retained {:,} alignments for {:,} CAGs".format(
                df.shape[0],
                df["CAG"].unique().shape[0]
            ))

        # Yield the subset for each CAG
        for group_ix, group_df in df.reindex(
            columns=[
                "CAG",
                "genome",
                "n_genes",
                "cag_prop",
            ]
        ).assign(
            group=df["CAG"].apply(lambda v: v % 1000)
        ).groupby("group"):
            yield group_ix, group_df.drop(columns="group")

    else:

        logging.info("No genome alignments found")

        yield None, None


def parse_distance_matrices(store, all_keys):
    """Read in each of the distance matrices in the store."""

    for k in all_keys:
        if k.startswith("/distances/"):
            logging.info("Reading in {}".format(k))
            yield k.replace("/distances/", ""), pd.read_hdf(store, k)


def parse_corncob_results(store, all_keys, max_n_cags=250000):
    """Read in and parse the corncob results from the store."""

    key_name = "/stats/cag/corncob"
    if key_name in all_keys:

        logging.info("Reading in {}".format(key_name))

        # Read in the complete set of results
        df = pd.read_hdf(
            store,
            key_name
        )

        # Store information for no more than `max_n_cags`
        if max_n_cags is not None:
            logging.info("Limiting annotations to {:,} CAGs".format(max_n_cags))
            df = df.query(
                "CAG < {}".format(max_n_cags)
            )
            logging.info("Filtered down to {:,} associations for {:,} CAGs".format(
                df.shape[0],
                df["CAG"].unique().shape[0],
            ))

        # Compute the log10 p_value and q_value
        logging.info("Corncob: Calculating -log10 p-values and q-values")
        df = add_neg_log10_values(df)

        for parameter, parameter_df in df.groupby(
            "parameter"
        ):
            if parameter != "(Intercept)":
                yield parameter, parameter_df.drop(
                    columns="parameter"
                ).reset_index(
                    drop=True
                )


def parse_enrichment_results(store, all_keys):
    """Read in and parse the betta results for annotation-level associations."""

    key_name = "/stats/enrichment/betta"

    if key_name in all_keys:

        logging.info("Reading in {}".format(key_name))

        # Read in the complete set of results
        df = pd.read_hdf(
            store,
            key_name
        )

        # Compute the log10 p_value and q_value
        logging.info("Betta: Calculating -log10 p-values and q-values")
        df = add_neg_log10_values(df)

        for parameter, parameter_df in df.groupby(
            "parameter"
        ):

            if parameter != "(Intercept)":

                for annotation, annotation_df in parameter_df.groupby(
                    "annotation"
                ):
                    yield parameter, annotation, annotation_df.drop(
                        columns=["parameter", "annotation"]
                    ).reset_index(
                        drop=True
                    )


def add_neg_log10_values(df):
    for k in ["p", "q"]:

        old_col = "{}_value".format(k)
        new_col = "neg_log10_{}value".format(k)

        df = df.assign(
            NEW_COL=df[
                old_col
            ].clip(
                lower=df[old_col][df[old_col] > 0].min()
            ).apply(
                np.log10
            ) * -1
        ).rename(
            columns={
                "NEW_COL": new_col
            }
        )

    # Also add the abs(wald) as abs_wald
    if "wald" in df.columns.values:
        df = df.assign(
            abs_wald = df["wald"].abs()
        )

    return df


def index_geneshot_results(input_fp, output_fp, skip_enrichments=["eggNOG_desc"]):

    # Keep all of the data in a dict linking the key to the table
    dat = {}

    # Keep a summary of the analysis results which are present
    analysis_features = []

    # Some data will be copied directly from input to output
    # Those tables which need to be transformed in any way will be copied
    # into memory and written out at the end

    # First go through and copy some data elements directly
    input_store = h5py.File(input_fp, "r")
    output_store = h5py.File(output_fp, "w")

    # Copy a set of input tables to their destination
    for source_key, dest_key in [
        ("summary/experiment", "experiment_metrics"),
        ("genomes/manifest", "genome_manifest"),
        ("genomes/summary", "genome_summary"),
        ("genomes/annotations", "genome_annotations"),
        ("genomes/detail", "genome_details"),
    ]:
        if source_key in input_store:
            logging.info("Copying %s to %s" % (source_key, dest_key))
            input_store.copy(
                source_key,
                output_store,
                name=dest_key
            )
        else:
            logging.info("Table(s) in %s not found" % source_key)

    # Close the connections
    input_store.close()
    output_store.close()

    # Open a connection to the input HDF5 for writing processed tables
    with pd.HDFStore(input_fp, "r") as store:
        all_keys = store.keys()

        # Read in the manifest
        dat["/manifest"] = parse_manifest(store)

        # Read in the specimen metrics
        dat["/specimen_metrics"] = parse_specimen_metrics(store)

        # Read in the CAG annotations
        dat["/cag_annotations"] = parse_cag_annotations(store)

        # Read in the CAG abundances
        dat["/cag_abundances"] = parse_cag_abundances(store)

        # Read in the genome containments
        for group_ix, group_df in parse_genome_containment(store):
            if group_ix is not None:
                dat["/genome_containment/{}".format(group_ix)] = group_df

        # Read in the distance matrices
        for metric_name, metric_df in parse_distance_matrices(store, all_keys):

            # Record which distance matrices are present
            analysis_features.append({
                "group": "distances",
                "key": "metric",
                "value": metric_name
            })

            # Store the actual distance matrix
            dat["/distances/{}".format(metric_name)] = metric_df

        # Read in the corncob results
        for parameter, parameter_df in parse_corncob_results(store, all_keys):

            # Record which corncob results are present
            analysis_features.append({
                "group": "cag_associations",
                "key": "parameter",
                "value": parameter
            })

            # Store the actual corncob results
            dat["/cag_associations/{}".format(parameter)] = parameter_df

        # Read in the betta results (for enrichment of corncob associations by annotation)
        for parameter, annotation, df in parse_enrichment_results(store, all_keys):

            # Skip a user-defined set of annotations
            # By default this is intended to include eggNOG_desc
            if annotation in skip_enrichments:
                continue

            # NOTE: The eggNOG_desc entrichments will be skipped by default. However,
            # if the user provides skip_enrichments=[], then this code block will help
            # to format that long text string.
            # Trim the eggNOG_desc labels to 100 character
            if annotation == "eggNOG_desc":
                # Trim the `eggNOG_desc` to 100 characters, if present
                df = df.apply(
                    lambda c: c.apply(lambda n: n[:100] if isinstance(n, str) and len(n) > 100 else n) if c.name == "label" else c
                )

            # Record which enrichment results are present
            analysis_features.append({
                "group": "enrichments",
                "key": "annotation",
                "value": annotation
            })

            # Store the actual betta results
            dat["/enrichments/{}/{}".format(parameter, annotation)] = df

        # Assemble the `analysis_features` table
        dat["/analysis_features"] = pd.DataFrame(
            analysis_features
        ).drop_duplicates()

        # Read in the taxonomy, if present
        if "/ref/taxonomy" in all_keys:
            tax = Taxonomy(store)
        else:
            tax = None

        # Read in the gene annotations for just those CAGs
        functional_annot_df, counts_df, rank_summaries = parse_gene_annotations(store, tax)

    # Store the summary annotation tables if the annotations are available
    if functional_annot_df is not None:
        for group_ix, group_df in functional_annot_df.assign(
            group = functional_annot_df["CAG"].apply(lambda v: v % 1000)
        ).groupby("group"):
            key_name = "/gene_annotations/functional/{}".format(group_ix)
            dat[key_name] = group_df.drop(columns="group")

    if counts_df is not None:
        for group_ix, group_df in counts_df.assign(
            group = counts_df["CAG"].apply(lambda v: v % 1000)
        ).groupby("group"):
            key_name = "/gene_annotations/taxonomic/all/{}".format(group_ix)
            dat[key_name] = group_df.drop(columns="group")

        for rank, rank_df in rank_summaries.items():
            for group_ix, group_df in rank_df.assign(
                group = rank_df["CAG"].apply(lambda v: v % 1000)
            ).groupby("group"):
                key_name = "/gene_annotations/taxonomic/{}/{}".format(
                    rank, 
                    group_ix
                )
                dat[key_name] = group_df.drop(columns="group")

    # Write out all of the tables to HDF5
    with pd.HDFStore(output_fp, "a") as store:
        logging.info("Writing to {}".format(output_fp))
        for key_name, df in dat.items():

            if df is None:
                logging.info(
                    "No information found for {}".format(key_name)
                )

            else:
                logging.info("Writing a table with {:,} rows and {:,} columns to {}".format(
                    df.shape[0], df.shape[1], key_name
                ))

                df.to_hdf(
                    store,
                    key_name
                )


if __name__ == "__main__":

    log_formatter = logging.Formatter(
        "%(asctime)s %(levelname)-8s [GLAM Index] %(message)s"
    )
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    # Write logs to STDOUT
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_formatter)
    root_logger.addHandler(console_handler)

    parser = argparse.ArgumentParser(
        description="""
        Index a set of geneshot results for visualization with the GLAM Browser.

        Example Usage:

        index.py <INPUT_HDF_FP> <OUTPUT_HDF_FP>

        """
    )

    parser.add_argument(
        "input",
        type=str,
        help="Path to results HDF5 file generated by geneshot"
    )

    parser.add_argument(
        "output",
        type=str,
        help="Path to write out index HDF5 file which can be visualized with the GLAM Browser"
    )

    # Parse the arguments
    args = parser.parse_args()

    # Make sure the input file exists
    assert os.path.exists(args.input), "Cannot find {}".format(args.input)

    # Make sure that the output file does not exist
    assert os.path.exists(args.output) is False, "{} already exists".format(args.output)

    index_geneshot_results(
        args.input,
        args.output
    )

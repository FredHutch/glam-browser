#!/usr/bin/env python3

# Script used to process a set of geneshot results for visualization in the GLAM Browser
import argparse
import logging
import os
import pandas as pd
import numpy as np


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


def parse_taxonomy(store):
    """Read in the taxonomy table."""

    # Read the taxonomy table
    logging.info("Reading in /ref/taxonomy")
    return pd.read_hdf(store, "/ref/taxonomy")


def parse_experiment_metrics(store):
    """Read in the experiment metrics"""
    key_name = "/summary/experiment"

    logging.info("Reading in {}".format(key_name))

    return pd.read_hdf(store, key_name)


def parse_specimen_metrics(store):
    """Read in the specimen metrics"""
    key_name = "/summary/all"

    logging.info("Reading in {}".format(key_name))

    df = pd.read_hdf(store, key_name)

    # Compute `prop_reads`
    return df.assign(
        prop_reads = df["aligned_reads"] / df["n_reads"]
    )


def parse_cag_annotations(store):
    """Read in the CAG annotations."""
    key_name = "/annot/cag/all"

    logging.info("Reading in {}".format(key_name))

    df = pd.read_hdf(store, key_name)

    # Compute `prop_reads`
    return df.assign(
        size_log10 = df["size"].apply(np.log10)
    ).drop(
        columns=["std_abundance"]
    )


def parse_gene_annotations(store):
    """Read in the gene-level annotations."""
    key_name = "/annot/gene/all"

    logging.info("Reading in {}".format(key_name))

    df = pd.read_hdf(store, key_name)

    # Trim the `eggNOG_desc` to 100 characters, if present
    df = df.apply(
        lambda c: c.apply(lambda n: n[:100] if isinstance(n, str) and len(n) > 100 else n) if c.name == "eggNOG_desc" else c
    )

    # Yield a table for each individual CAG
    for cag_id, cag_df in df.groupby("CAG"):

        # Just save the number of unique values for tax ID and eggNOG annotation
        yield cag_id, gene_annotation_value_counts(cag_df)


def gene_annotation_value_counts(cag_df):
    """Count up the unique annotations for a given CAG."""
    return pd.DataFrame([
        {
            "annotation": col_name,
            "value": value,
            "count": count
        }
        for col_name in ["tax_id", "eggNOG_desc"]
        for value, count in cag_df[col_name].dropna().value_counts().items()
    ])

def parse_cag_abundances(store):
    """Read in the CAG abundances."""
    key_name = "/abund/cag/wide"

    logging.info("Reading in {}".format(key_name))

    return pd.read_hdf(store, key_name)


def parse_distance_matrices(store, all_keys):
    """Read in each of the distance matrices in the store."""

    for k in all_keys:
        if k.startswith("/distances/"):
            logging.info("Reading in {}".format(k))
            yield k.replace("/distances/", ""), pd.read_hdf(store, k)


def parse_corncob_results(store, all_keys):
    """Read in and parse the corncob results from the store."""

    key_name = "/stats/cag/corncob"
    if key_name in all_keys:

        logging.info("Reading in {}".format(key_name))

        # Read in the complete set of results
        df = pd.read_hdf(
            store,
            key_name
        )

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


def parse_gene_annotation_summaries(store, dat, annotation_columns=["eggNOG_desc", "tax_id"]):
    """If we have corncob results, make a summary table for CAGs passing FDR for each."""
    for key_name, df in dat.items():
        if key_name.startswith("/cag_associations/"):
            parameter_name = key_name.replace("/cag_associations/", "")

            # Get the list of CAGs which pass the FDR threshold
            cag_id_list = df.query(
                "q_value <= 0.01"
            )["CAG"].tolist()

            if len(cag_id_list) > 0:

                # Get the annotations for all of the genes in these CAGs
                summary_df = [
                    dat["/gene_annotations/CAG/CAG{}".format(cag_id)]
                    for cag_id in cag_id_list
                ]
                if len(summary_df) > 0:
                    summary_df = pd.concat(summary_df)

                    # If this group of CAGs has any annotations, save it
                    if any([
                        k in summary_df.columns.values and summary_df[k].dropna().shape[0] > 0 
                        for k in annotation_columns
                    ]):

                        # Compute summary metrics by CAG
                        yield parameter_name, pd.DataFrame([
                            {
                                "CAG": cag_id,
                                "annotation": col_name,
                                "key": label,
                                "count": label_df.shape[0]
                            }
                            for cag_id, cag_df in summary_df.groupby("CAG")
                            for col_name in annotation_columns
                            if col_name in cag_df.columns.values
                            for label, label_df in cag_df.groupby(col_name)
                        ])


def get_cags_to_include(dat, top_n=10000):
    """Limit the number of CAGs for which we will save information."""
    cags_to_include = set([])

    for key_name, df in dat.items():

        # Table with CAG annotations
        if key_name == "/cag_annotations":

            # Keep the top_n CAGs by size and mean abundance
            for col_name in ["size", "mean_abundance"]:
                
                logging.info("Selecting the top {:,} CAGs by {}".format(
                    top_n,
                    col_name
                ))

                cags_to_include.update(set(
                    df.sort_values(
                        by=col_name,
                        ascending=False
                    ).head(
                        top_n
                    )["CAG"].tolist()
                ))
                logging.info("Total number of CAGs to display: {:,}".format(
                    len(cags_to_include)
                ))

        # Table with the abundances of each CAG
        elif key_name == "/cag_abundances":

            logging.info("Selecting the top {:,} CAGs by maximum abundance".format(
                top_n
            ))

            # Compute the maximum relative abundance of every CAG
            max_abund = df.set_index("CAG").max(axis=1)

            # Keep the top_n CAGs by maximum relative abundances
            cags_to_include.update(set(
                [
                    cag_id
                    for cag_id in max_abund.sort_values(
                        ascending=False
                    ).head(top_n).index.values
                ]
            ))
            logging.info("Total number of CAGs to display: {:,}".format(
                len(cags_to_include)
            ))

        # Table with CAG association metrics
        elif key_name.startswith("/cag_associations"):

            logging.info("Selecting the top {:,} CAGs by absolute Wald ({})".format(
                top_n,
                key_name
            ))

            # Keep the top_n CAGs by absolute Wald metric
            cags_to_include.update(set(
                df.sort_values(
                    by="abs_wald",
                    ascending=False
                ).head(
                    top_n
                )["CAG"].tolist()
            ))
            logging.info("Total number of CAGs to display: {:,}".format(
                len(cags_to_include)
            ))

    return cags_to_include


def filter_data_to_selected_cags(dat, cags_to_include):
    """Reduce the size of the data that will be saved by filtering to these CAGs."""

    # ABUNDANCES #
    # Only keep abundances for the CAGs in the list
    logging.info("Retaining abundances for {:,} / {:,} CAGs".format(
        len(cags_to_include), dat["/cag_abundances"].shape[0]
    ))
    # Subset the table 
    dat["/cag_abundances"] = dat["/cag_abundances"].loc[
        dat["/cag_abundances"]["CAG"].isin(cags_to_include)
    ]
    # Every CAG should have an annotation
    assert dat["/cag_abundances"].shape[0] == len(cags_to_include)

    # Iterate over the tables to access programmatically named tables
    tables_to_delete = []
    for table_name in dat:
        # ASSOCIATIONS #
        # Consider all of the CAG association tables
        if table_name.startswith("/cag_associations/"):
            # Subset the table
            dat[table_name] = dat[table_name].loc[
                dat[table_name]["CAG"].isin(cags_to_include)
            ]
            # Every set of associations should have some CAGs included
            assert dat[table_name].shape[0] > 0

        # DETAILED ANNOTATIONS #
        elif table_name.startswith("/gene_annotations/CAG/CAG"):
            cag_id = int(table_name.replace("/gene_annotations/CAG/CAG", ""))

            # Delete tables for CAGs that aren't in this list
            if cag_id not in cags_to_include:
                tables_to_delete.append(table_name)

    for key_name in tables_to_delete:
        del dat[key_name]


def index_geneshot_results(input_fp, output_fp):

    # Keep all of the data in a dict linking the key to the table
    dat = {}

    # Keep a summary of the analysis results which are present
    analysis_features = []

    # Open a connection to the input HDF5
    with pd.HDFStore(input_fp, "r") as store:
        all_keys = store.keys()

        # Read in the manifest
        dat["/manifest"] = parse_manifest(store)

        # Read in the experiment metrics
        dat["/experiment_metrics"] = parse_experiment_metrics(store)

        # Read in the specimen metrics
        dat["/specimen_metrics"] = parse_specimen_metrics(store)

        # Read in the CAG annotations
        dat["/cag_annotations"] = parse_cag_annotations(store)

        # Read in the CAG abundances
        dat["/cag_abundances"] = parse_cag_abundances(store)

        # Read in the taxonomy
        if "/ref/taxonomy" in all_keys:
            dat["/taxonomy"] = parse_taxonomy(store)

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

        # If we have corncob results, make aggregate tables for each parameter
        # which include all CAGs with FDR alpha=0.01 and summarize the number
        # of genes from each CAG which have a given annotation
        items_to_add = {
            "/gene_annotations/parameter/{}".format(parameter): parameter_df
            for parameter, parameter_df in parse_gene_annotation_summaries(store, dat)
        }
        # Add those items to the larger `dat` object as a second isolated step
        for k, v in items_to_add.items():
            dat[k] = v

    # Assemble the `analysis_features` table
    dat["/analysis_features"] = pd.DataFrame(
        analysis_features
    ).drop_duplicates()

    # Limit the number of CAGs for which we will save information
    cags_to_include = get_cags_to_include(dat)
    logging.info("Indexing details for {:,} CAGs out of {:,} total".format(
        len(cags_to_include),
        dat["/cag_abundances"].shape[0]
    ))

    # Add the gene annotations for those selected CAGs
    for cag_id, cag_annotations in parse_gene_annotations(store):
        if cag_id in cags_to_include:
            dat["/gene_annotations/CAG/CAG{}".format(cag_id)] = cag_annotations

    # Now actually subset the information to this set of CAGs
    filter_data_to_selected_cags(dat, cags_to_include)

    # Write out all of the tables to HDF5
    with pd.HDFStore(output_fp, "w") as store:
        logging.info("Writing to {}".format(output_fp))
        for key_name, df in dat.items():
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

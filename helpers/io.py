#!/usr/bin/env python3

from collections import defaultdict
import json
import os
import numpy as np
import pandas as pd
from filelock import Timeout, FileLock
from time import sleep

# All of the data in the indicated folder is loaded as a
# list of Geneshot objects
def parse_directory(folder_path):
    # Set up the data object
    page_data = {
        "page_title": "GLAM Browser",
        "page_description": "Gene-level analysis of microbiome datasets",
    }

    # Read the files in the directory
    file_list = os.listdir(folder_path)

    # If there is a file named 'manifest.json', read its information
    if "manifest.json" in file_list:
        # Load the manifest as JSON format
        manifest = json.load(open(
            os.path.join(
                folder_path,
                "manifest.json"
            )
        ))

        # Iterate over the elements to import from the manifest
        for k in ["page_title", "page_description"]:
            # Skip any keys which are not found
            if k in manifest:
                # Add to the data structure
                page_data[k] = str(manifest[k])

        # If there is a list of "contents" in the manifest, use that
        if "contents" in manifest:
            # Create the list, which will be filled with Geneshot objects
            page_data["contents"] = []
            # Make sure that each element has the expected contents
            for i in manifest["contents"]:
                assert "fp" in i
                i["fp"] = os.path.join(folder_path, i["fp"])
                assert os.path.exists(i["fp"]), "File not found: {}".format(i["fp"])
                if "name" not in i:
                    i["name"] = i["fp"].split("/")[-1]
                page_data["contents"].append(i)

    # If the manifest did not have a "contents" list, 
    # iterate over the files in this folder
    if "contents" not in page_data:
        page_data["contents"] = [
            {
                "fp": os.path.join(folder_path, fp),
                "name": fp
            }
            for fp in file_list
            if fp.endswith(".hdf5")
        ]

    # Return the data
    return page_data


def hdf5_get_item(
    fp, 
    key_path, 
    index_col=None, 
    f=None, 
    where=None, 
    columns=None,
    timeout=5, 
    retry=5,
):
    """Read data from the HDF5 store."""

    # Set up a file lock to prevent multiple concurrent access attempts
    lock = FileLock("{}.lock".format(fp), timeout=timeout)

    # Read in the table
    print("Reading in {} from {}".format(key_path, fp))

    try:
        with lock:
            with pd.HDFStore(fp, "r") as store:
                try:
                    df = pd.read_hdf(
                        store,
                        key_path,
                        where=where,
                        columns=columns
                    )
                except KeyError:
                    return None
    except Timeout:
        print("Another instance of this application currently holds the lock on {}.".format(fp))

        sleep(retry)
        return hdf5_get_item(
            fp, 
            key_path, 
            index_col=index_col,
            f=f,
            where=where,
            columns=columns,
            timeout=timeout,
            retry=retry,
        )

    # Set the index
    if index_col is not None:
        df.set_index(
            index_col,
            inplace=True
        )
    # Apply a function
    if f is not None:
        df = f(df)

    return df


def hdf5_manifest(fp):
    return hdf5_get_item(fp, "/manifest", index_col="specimen")

def hdf5_richness(fp):
    return hdf5_get_item(
        fp,
        "/summary/all", 
        index_col="specimen", 
        f=hdf5_richness_f
    )

def hdf5_richness_f(df):
    return df.assign(
        prop_reads_aligned=df["aligned_reads"] / df["n_reads"]
    )

def hdf5_cag_summary(fp):
    return hdf5_get_item(fp, "/annot/cag/all", index_col="CAG", f=hdf5_cag_summary_f)

def hdf5_cag_summary_f(df):
    return df.assign(
        size_log10=df["size"].apply(np.log10)
    )

def hdf5_metrics(fp):
    return hdf5_get_item(
        fp,
        "/summary/experiment", 
        index_col="variable"
    )["value"]

def hdf5_distances(fp, metric):
    return hdf5_get_item(
        fp,
        "/distances/{}".format(metric),
        index_col="specimen"
    )

def hdf5_corncob(fp):
    return hdf5_get_item(fp, "/stats/cag/corncob", f=hdf5_corncob_f)
def hdf5_corncob_f(df):
    return df.assign(
        neg_log_pvalue=df["p_value"].apply(np.log10) * -1,
        neg_log_qvalue=df["q_value"].apply(np.log10) * -1,
    )

def hdf5_taxonomy(fp):
    return hdf5_get_item(fp, "/ref/taxonomy", f=hdf5_taxonomy_f)
def hdf5_taxonomy_f(df):
    # Format the parent tax ID as an integer
    return df.apply(
        lambda c: c.fillna(0).apply(float).apply(
            int) if c.name in ["parent", "tax_id"] else c,
    ).set_index("tax_id")

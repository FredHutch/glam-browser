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
    timeout=5, 
    retry=5,
):
    """Read data from the HDF5 store."""

    # Set up a file lock to prevent multiple concurrent access attempts
    lock = FileLock("{}.lock".format(fp), timeout=timeout)

    # Read in the table
    try:
        with lock:
            with pd.HDFStore(fp, "r") as store:
                try:
                    df = pd.read_hdf(
                        store,
                        key_path
                    )
                except KeyError:
                    return None
    except Timeout:

        sleep(retry)
        return hdf5_get_item(
            fp, 
            key_path, 
            timeout=timeout,
            retry=retry,
        )

    return df


def hdf5_taxonomy(fp):
    return hdf5_get_item(
        fp, 
        "/taxonomy"
    ).apply(
        lambda c: c.fillna(0).apply(float).apply(int) if c.name in ["parent", "tax_id"] else c,
    ).set_index(
        "tax_id"
    )

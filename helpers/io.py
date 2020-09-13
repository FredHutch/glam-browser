#!/usr/bin/env python3

from collections import defaultdict
from filelock import Timeout, FileLock
import h5py
import json
import numpy as np
import os
import pandas as pd
from time import sleep

# All of the data in the indicated folder is loaded as a
# list of Geneshot objects
class Manifest:

    def __init__(self, folder_path):

        # Save the folder path
        self.folder_path = folder_path

        # Keep a list of all files in this folder
        self.all_filepaths = []

        # Read in the page data
        self.page_data = self.parse_directory(folder_path)

    def parse_fp(self, selected_dataset, page=None, key=None):
        """Function to return the file path for whichever dataset is selected."""
        if selected_dataset in [[-1], ["-1"], -1, "-1"]:
            # No dataset was selected
            return None

        if isinstance(selected_dataset, list):
            selected_dataset = selected_dataset[0]
        selected_dataset = int(selected_dataset)
        
        # Get the list of datasets
        dataset_list = self.dataset_list(page=page, key=key)

        if len(dataset_list) < (selected_dataset + 1):
            return None

        return os.path.join(
            self.folder_path,
            dataset_list[selected_dataset]["fp"]
        )

    def page_title(self, page=None, key=None):
        return self.get_page_info("page_title", page=page, key=key)

    def page_description(self, page=None, key=None):
        return self.get_page_info("page_description", page=page, key=key)

    def dataset_list(self, page=None, key=None):
        dataset_list = self.get_page_info("contents", page=page, key=key)
        if dataset_list is None:
            return []
        else:
            return dataset_list

    def default(self, selected_dataset, default_field, page=None, key=None):
        if selected_dataset in [[-1], ["-1"], -1, "-1"]:
            # No dataset was selected
            return None

        if isinstance(selected_dataset, list):
            selected_dataset = selected_dataset[0]
        selected_dataset = int(selected_dataset)

        # Get the list of datasets
        dataset_list = self.dataset_list(page=page, key=key)

        if len(dataset_list) < (selected_dataset + 1):
            return None

        if "defaults" not in dataset_list[selected_dataset]:
            return None
        
        return dataset_list[selected_dataset]["defaults"].get(default_field)

    def get_page_info(self, info_name, page=None, key=None):

        # Limit the amount of data available to prevent accessing the key directly
        if isinstance(info_name, str):
            if info_name not in ["page_title", "page_description", "contents"]:
                return None
        elif isinstance(info_name, tuple):
            if info_name[0] != "default":
                return None
        else:
            return None

        # If no page was selected, use 'main'
        if page is None:
            page = 'main'

        # If the page starts with a '/', remove it
        if page.startswith("/"):
            if page == "/":
                page = "main"
            else:
                page = page[1:]

        # Make sure that the indicated page has data
        if page not in self.page_data:
            return None

        # If a key is required, only return the file if the provided key matches
        if "key" in self.page_data[page]:
            if key is not None and key == self.page_data[page]["key"]:
                if isinstance(info_name, str):
                    return self.page_data[page][info_name]
                else:
                    return self.page_data[page]["defaults"][info_name[1]]
            else:
                return None
        
        # No key is required
        else:
            if isinstance(info_name, str):
                return self.page_data[page][info_name]
            else:
                return self.page_data[page]["defaults"][info_name[1]]


    def parse_directory(self, folder_path):

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

            # Make sure that the manifest conforms to the expected schema
            assert isinstance(manifest, dict), "Manifest must be a dict"

            # Make sure that the manifest is EITHER a single page, or multiple valid pages
            all_valid = all([self.is_valid_group(v) for v in manifest.values()])
            assert self.is_valid_group(manifest) or all_valid, "Manifest is not formatted correctly"

            # Either the manifest contains information for a single page, or multiple

            # The manifest is a single page
            if self.is_valid_group(manifest):

                # Validate that all files are present
                self.assert_valid_filepaths(manifest)

                # Return this single page as the manifest
                return {
                    "main": manifest
                }

            else:
                # The manifest must contain multiple pages

                # Validate that all files are present
                for v in manifest.values():
                    self.assert_valid_filepaths(v)

                # Return this data structure
                return manifest

        # If the manifest does not exist, just use all of the files in the folder
        else:
            # Make a single page with just the files in this folder
            return {
                "main": {
                    "page_title": "GLAM Browser",
                    "page_description": "Gene-level analysis of microbiome datasets",
                    "contents": [
                    {
                        "fp": os.path.join(folder_path, fp),
                        "name": fp
                    }
                    for fp in file_list
                    if fp.endswith(".hdf5")
                    ]
                }
            }


    def is_valid_group(self, page_data, verbose=False):
        """Return True if the page data conforms to the schema."""

        if isinstance(page_data, dict) is False:
            if verbose:
                print("Page is not formatted as a dict: {}".format(json.dumps(page_data)))
            return False

        for k in ["page_title", "page_description", "contents"]:
            if k not in page_data:
                if verbose:
                    print("Missing key '{}' in '{}'".format(k, json.dumps(page_data)))
                return False

        # Check for defaults
        if "defaults" in page_data:
            if isinstance(page_data["defaults"], dict) is False:
                if verbose:
                    print("Section 'defaults' must be a dict")
                return False
            for k, v in page_data["defaults"].items():
                for t in [k, v]:
                    if isinstance(t, str) is False:
                        if verbose:
                            print("All defaults must be formatted as strings")
                        return False

        if isinstance(page_data["page_title"], str) is False:
            return False

        if isinstance(page_data["page_description"], str) is False:
            return False

        if isinstance(page_data["contents"], list) is False:
            return False


        for file_data in page_data["contents"]:
            if isinstance(file_data, dict) is False:
                if verbose:
                    print("File data is not a dict: '{}'".format(json.dumps(file_data)))
                return False
            
            for k in ["fp", "name", "description"]:
                if k not in file_data:
                    if verbose:
                        print("Missing key '{}' in '{}'".format(k, json.dumps(file_data)))
                    return False

        return True


    def assert_valid_filepaths(self, page_data):
        """Assert that all of the file paths in the page data point to files that exist."""

        assert "contents" in page_data, \
            "Data does not have 'contents' key"

        for file_data in page_data["contents"]:

            assert "fp" in file_data, \
                "File does not contain 'fp' key {}".format(json.dumps(file_data))

            fp = os.path.join(self.folder_path, file_data["fp"])

            # Make sure that the file exists
            assert os.path.exists(fp), \
                "File does not exist: {}".format(fp)

            # Add the filepath to the list of all filepaths
            self.all_filepaths.append(fp)


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


def hdf5_get_keys(
    fp, 
    group_path, 
    timeout=5, 
    retry=5,
):
    """Read keys from a group in the HDF5 store."""

    # Set up a file lock to prevent multiple concurrent access attempts
    lock = FileLock("{}.lock".format(fp), timeout=timeout)

    # Read in the keys
    try:
        with lock:
            with h5py.File(fp, "r") as f:
                try:
                    key_list = list(f[group_path].keys())
                except:
                    return None
    except Timeout:

        sleep(retry)
        return hdf5_get_keys(
            fp, 
            group_path, 
            timeout=timeout,
            retry=retry,
        )

    return key_list


def hdf5_taxonomy(fp):
    return hdf5_get_item(
        fp, 
        "/taxonomy"
    ).apply(
        lambda c: c.fillna(0).apply(float).apply(int) if c.name in ["parent", "tax_id"] else c,
    ).set_index(
        "tax_id"
    )

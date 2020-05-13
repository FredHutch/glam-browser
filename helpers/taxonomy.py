#!/usr/bin/env python3

from collections import defaultdict
import pandas as pd

def path_to_root(tax_id, taxonomy_df, max_steps=100):
    """Yield a list with all of the taxa above this one."""

    visited = set([])

    for _ in range(max_steps):

        # Yield this taxon
        yield tax_id

        # Add to the path we have visited
        visited.add(tax_id)

        # Get the parent of this taxon
        parent_id = taxonomy_df.loc[tax_id, "parent"]

        # If the chain has ended, stop
        if parent_id in visited or parent_id == 0:
            break

        # Otherwise, keep walking up
        tax_id = parent_id


def make_cag_tax_df(tax_id_list, taxonomy_df, ranks_to_keep=["phylum", "class", "order", "family", "genus", "species"]):
    """Return a nicely formatted taxonomy table from a list of tax IDs."""

    # We will construct a table with all of the taxa in the tree, containing
    # The name of that taxon
    # The rank of that taxon
    # The number of genes found at that taxon or in its decendents
    # The name of the parent of that taxon

    counts = defaultdict(int)

    # Iterate over each terminal leaf
    for tax_id, n_genes in tax_id_list.apply(int).value_counts().items():

        # Skip taxa which aren't in the taxonomy
        if tax_id not in taxonomy_df.index.values:
            continue

        # Walk up the tree from the leaf to the root
        for anc_tax_id in path_to_root(tax_id, taxonomy_df):

            # Add to the sum for every node we visit along the way
            counts[anc_tax_id] += n_genes

    if len(counts) == 0:
        return

    # Make a DataFrame
    df = pd.DataFrame({"count": counts})

    # Add the name, parent, rank
    df = df.assign(
        tax_id=df.index.values,
        parent_tax_id=taxonomy_df["parent"],
        rank=taxonomy_df["rank"],
    )

    # Set the parent of the root as ""
    df.loc[0, "parent"] = ""

    # Remove any taxa which aren't at the right rank (but keep the root)
    df = df.assign(
        to_remove=df.apply(
            lambda r: r["rank"] not in (ranks_to_keep) and r["tax_id"] != 0,
            axis=1
        )
    )

    # Remove all of the taxa which aren't at the right rank
    while df["to_remove"].any():
        ix_to_remove = df.index.values[df["to_remove"]][0]

        # Drop this row from the filtered table
        # Also update any rows which include this taxon as a parent
        df = df.drop(
            index=ix_to_remove
        ).replace(to_replace={
            "parent_tax_id": {
                df.loc[ix_to_remove, "tax_id"]: df.loc[ix_to_remove, "parent_tax_id"]
            }
        })

    # Make a dict linking the tax ID with the taxon name
    tax_names = taxonomy_df["name"].reindex(index=df.index)

    # Count up the frequency of each name
    tax_names_vc = tax_names.value_counts()

    # Identify any repeated names
    repeated_names = tax_names_vc.index.values[tax_names_vc.values > 1]

    # Make each name unique by adding the tax ID (if needed)
    for n in repeated_names:
        for tax_id in tax_names.index.values[tax_names == n]:
            tax_names.loc[tax_id] = "{} ({})".format(n, tax_id)

    # Add this nicely formatted name to the output to replace the tax IDs
    df = df.assign(
        name=df["tax_id"].apply(lambda t: tax_names.get(t, "")),
        parent=df["parent_tax_id"].apply(lambda t: tax_names.get(t, "")),
    )

    return df.reindex(columns=["name", "parent", "count"])

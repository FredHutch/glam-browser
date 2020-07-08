#!/usr/bin/env python3

from collections import defaultdict
import pandas as pd

def path_to_root(tax_id, taxonomy_df, max_steps=100):
    """Yield a list with all of the taxa above this one."""

    visited = set([])

    for _ in range(max_steps):

        # Add to the path we have visited
        visited.add(tax_id)

        # Get the parent of this taxon
        parent_id = taxonomy_df.loc[tax_id, "parent"]

        # If the chain has ended, stop
        if parent_id in visited or parent_id == 0:
            break

        # Otherwise, keep walking up
        tax_id = parent_id

    return visited


def make_cag_tax_df(taxa_vc, taxonomy_df, ranks_to_keep=["phylum", "class", "order", "family", "genus", "species"]):
    """Return a nicely formatted taxonomy table from a list of tax IDs."""

    # We will construct a table with all of the taxa in the tree, containing
    # The name of that taxon
    # The rank of that taxon
    # The name of the parent of that taxon

    # The number of genes found at that taxon or in its decendents
    counts = defaultdict(int)

    # Get the ancestors of every observed taxa
    ancestors = {}

    # To do so, start by iterating over every observed taxon
    for tax_id in taxa_vc.index.values:

        # Skip taxa which aren't in the taxonomy
        if tax_id not in taxonomy_df.index.values:
            continue

        # Make sure to add the ancestors for all taxids at-or-above those which were observed
        for anc_tax_id in path_to_root(tax_id, taxonomy_df):

            # Only add it if we haven't yet
            ancestors[anc_tax_id] = ancestors.get(
                anc_tax_id, 
                path_to_root(anc_tax_id, taxonomy_df)
            )

    # The 'consistent_counts' are the number of hits to taxa which
    # are consistent with this taxa, either above or below it in the taxonomy
    consistent_counts = {
        tax_id: 0
        for tax_id, ancestors in ancestors.items()
    }

    # Keep track of the total number of genes with a valid tax ID
    total_genes_assigned = 0

    # Iterate over each terminal leaf
    for tax_id, n_genes in taxa_vc.items():

        # Skip taxa which aren't in the taxonomy
        if tax_id not in taxonomy_df.index.values:
            continue

        # Count all genes part of this analysis
        total_genes_assigned += n_genes

        # Walk up the tree from the leaf to the root
        for anc_tax_id in ancestors[tax_id]:

            # Add to the sum for every node we visit along the way
            counts[anc_tax_id] += n_genes

        # Iterate over every observed taxon and add to the
        # 'consistent_counts' if it is an ancestor or descendent
        # of this taxon
        for other_tax_id, other_taxon_ancestors in ancestors.items():
            if other_tax_id == tax_id:
                consistent_counts[other_tax_id] += n_genes
            elif tax_id in other_taxon_ancestors:
                consistent_counts[other_tax_id] += n_genes
            elif other_tax_id in ancestors[tax_id]:
                consistent_counts[other_tax_id] += n_genes

    if len(counts) == 0:
        return

    # Make a DataFrame
    df = pd.DataFrame({
        "count": counts,
        "consistent": consistent_counts,
    }).fillna(
        0
    )

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
        total=total_genes_assigned,
    )

    return df.reindex(columns=["name", "parent", "count", "consistent", "rank", "total"])

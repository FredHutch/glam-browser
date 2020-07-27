#!/usr/bin/env python3

import json
import numpy as np
import os
import pandas as pd
import dash_core_components as dcc
import plotly.express as px
import plotly.graph_objs as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import leaves_list
from scipy.stats import zscore
from seaborn import color_palette
from skbio.stats.distance import permanova, DistanceMatrix, anosim
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


def parse_manifest_json(manifest_json, manifest_df):
    """Parse the filtered manifest stored in the browser."""
    if manifest_json is None:
        # Return the full table from the HDF5 if nothing is stored in the browser
        return manifest_df
    else:
        # Otherwise, return the table stored in the browser
        return pd.DataFrame(
            json.loads(manifest_json)
        ).set_index(
            "specimen"
        )

def add_metadata_to_dataframe(plot_df, plot_manifest_df, metadata):
    # Add the specified metadata to the DataFrame
    plot_df = plot_df.assign(
        METADATA=plot_manifest_df.reset_index(
        ).reindex(
            columns=["specimen", metadata]
        ).drop_duplicates(
        ).set_index(
            "specimen"
        )[metadata]
    ).rename(columns={
        "METADATA": metadata
    })

    # Remove specimens which lack the specified metadata
    assert plot_df[metadata].isnull().mean(
    ) < 1.0, "Metadata is missing for all specimens"

    # At least one sample is missing metadata
    if plot_df[metadata].isnull().any():
        # Subset to specimens which contain the metadata
        plot_df = plot_df.dropna()

    # Make sure we didn't drop all of the samples
    assert plot_df.shape[0] > 0

    # Make a numeric transform of the metadata
    if plot_df[metadata].apply(
        lambda n: isinstance(n, float) or isinstance(n, int)
    ).all():

        # The metadata is already numeric
        plot_df = plot_df.assign(
            METADATA_FLOAT=plot_df[metadata]
        )

    # Try to convert to datetime
    elif pd.to_datetime(plot_df[metadata], errors="coerce").isnull().sum() == 0:
        plot_df = plot_df.assign(
            METADATA_FLOAT=pd.to_datetime(
                plot_df[metadata],
                errors="raise"
            ).apply(str)
        )

    # Treat as categorical
    else:
        # Assign a rank order
        rank_order = {
            n: ix
            for ix, n in enumerate(plot_df[metadata].drop_duplicates().sort_values().values)
        }
        plot_df = plot_df.assign(
            METADATA_FLOAT=plot_df[metadata].apply(rank_order.get)
        )

    # Set up a metadata color map depending on how many distinct values
    # there are. For small numbers, use "colorblind", otherwise "coolwarm"
    if plot_df["METADATA_FLOAT"].unique().shape[0] <= 5:
        palette_name = "colorblind"

    else:
        palette_name = "coolwarm"

    # Here is the actual color map
    cmap = dict(zip(
        plot_df["METADATA_FLOAT"].drop_duplicates().sort_values().values,
        color_palette(
            palette_name,
            plot_df["METADATA_FLOAT"].unique().shape[0]
        ).as_hex()
    ))

    # Now add that color to the plot DataFrame
    plot_df = plot_df.assign(
        METADATA_COLOR=plot_df["METADATA_FLOAT"].apply(cmap.get)
    )
    assert plot_df["METADATA_COLOR"].isnull().sum() == 0, (plot_df.head(), cmap)

    return plot_df


##################
# RICHNESS GRAPH #
##################
def update_richness_graph(
    richness_df,
    selected_metric,
    selected_type,
    selected_metadata,
    log_x,
    manifest_json,
    full_manifest_df,
):
    # Make sure the selected parameters are in scope
    assert selected_type in ["hist", "scatter"]

    # Get the filtered manifest from the browser
    plot_manifest_df = parse_manifest_json(manifest_json, full_manifest_df)

    # Subset the richness table based on the filtered manifest
    plot_richness_df = richness_df.reindex(
        index=list(set([n for n in plot_manifest_df.index.values]))
    )
    
    # Calculate the percent of reads aligned
    plot_richness_df = plot_richness_df.assign(
        pct_reads_aligned = (plot_richness_df["prop_reads"] * 100).apply(lambda v: round(v, 2))
    )

    # If metadata was selected, add it to the plot
    if selected_metadata != "none":

        # Add the indicated metadata to the plotting dataframe
        # This will also add a METADATA_COLOR column as appropriate
        plot_richness_df = add_metadata_to_dataframe(
            plot_richness_df,
            plot_manifest_df,
            selected_metadata
        )

        # Make sure that we have the needed metadata
        assert selected_metadata in plot_richness_df.columns.values
        assert "METADATA_COLOR" in plot_richness_df.columns.values

    assert selected_metric in plot_richness_df.columns.values, (selected_metric, richness_df.columns.values)

    metric_names = {
        "pct_reads_aligned": "Pct. Reads Aligned",
        "n_genes_aligned": "Num. Genes Aligned",
        "n_genes_assembled": "Num. Genes Assembled",
    }

    if selected_type == "scatter":
        if selected_metric == "n_genes_aligned":
            if selected_metadata == "none":
                hovertemplate = "Sample: %{id}<br>%{x:,} reads<br>%{y:,} genes detected by alignment<extra></extra>"
            else:
                hovertemplate = "Sample: %{id}<br>%{x:,} reads<br>%{y:,} genes detected by alignment<br>%{text}<extra></extra>"
        elif selected_metric == "n_genes_assembled":
            if selected_metadata == "none":
                hovertemplate = "Sample: %{id}<br>%{x:,} reads<br>%{y:,} genes detected by assembly<extra></extra>"
            else:
                hovertemplate = "Sample: %{id}<br>%{x:,} reads<br>%{y:,} genes detected by assembly<br>%{text}<extra></extra>"
        else:
            if selected_metadata == "none":
                hovertemplate = "Sample: %{id}<br>%{x:,} reads<br>%{y:.2f} percent of reads aligned uniquely<extra></extra>"
            else:
                hovertemplate = "Sample: %{id}<br>%{x:,} reads<br>%{y:.2f} percent of reads aligned uniquely<br>%{text}<extra></extra>"

        if selected_metadata == "none":
            fig = go.Figure(
                data=go.Scatter(
                    x=plot_richness_df["n_reads"],
                    y=plot_richness_df[selected_metric],
                    ids=plot_richness_df.index.values,
                    hovertemplate=hovertemplate,
                    mode="markers",
                )
            )
        else:
            fig = go.Figure(
                data=go.Scatter(
                    x=plot_richness_df["n_reads"],
                    y=plot_richness_df[selected_metric],
                    ids=plot_richness_df.index.values,
                    hovertemplate=hovertemplate,
                    mode="markers",
                    marker_color=plot_richness_df["METADATA_COLOR"],
                    text=plot_richness_df[selected_metadata].apply(
                        lambda n: "{}: {}".format(selected_metadata, n)
                    ),
                )
            )

        fig.update_layout(
            xaxis_title="Number of Reads",
        )
        if log_x == "on":
            fig.update_layout(
                xaxis_type="log",
            )
        else:
            fig.update_layout(
                xaxis_range=[0, plot_richness_df["n_reads"].max() * 1.05],
            )

    else:
        assert selected_type == "hist"

        fig = go.Figure(
            data=[
                go.Histogram(
                    y=plot_richness_df[selected_metric],
                    hovertemplate="Range: %{y}<br>Count: %{x}<extra></extra>",
                )
            ],
        )
        fig.update_layout(
            xaxis_title="Number of Specimens"
        )

    fig.update_layout(
        yaxis_range=[0, plot_richness_df[selected_metric].max() * 1.05],
        yaxis_title=metric_names[selected_metric],
        template="simple_white",
        height=500,
    )

    return fig


#######################
# SINGLE SAMPLE GRAPH #
#######################
def plot_sample_vs_cag_size(
    sample_abund,
    sample_name,
    cag_summary_df,
    display_metric,
):
    # Set up the DataFrame for plotting
    plot_df = cag_summary_df.reindex(
        columns=["size", "size_log10"]
    ).assign(
        prop=sample_abund
    ).query(
        "prop > 0"
    )
    if display_metric == "clr":
        # Calculate the CLR
        clr = plot_df["prop"].apply(np.log10)
        clr = clr - clr.mean()

        # Round to 2 decimals
        clr = clr.apply(lambda v: round(v, 2))

        # Add to the DataFrame
        plot_df = plot_df.assign(
            clr = clr
        )
    else:
        # Round the proportional abundance to 4 decimals
        plot_df = plot_df.apply(
            lambda c: c.apply(lambda v: round(v, 4)) if c.name == "prop" else c
        )

    # Reset the index to put CAG into a columns
    plot_df = plot_df.reset_index()

    # Rename the columns
    column_names = {
        "size": "Number of Genes",
        "clr": "Relative Abundance (Centered Log-Ratio)",
        "prop": "Relative Abundance (Proportion)"
    }
    plot_df = plot_df.rename(
        columns=column_names
    ).reindex(
        columns=["CAG", "Number of Genes", column_names[display_metric]]
    )

    # Make the plot
    fig = px.scatter(
        plot_df,
        x="Number of Genes",
        y=column_names[display_metric],
        hover_data=["CAG", "Number of Genes", column_names[display_metric]]
    )

    # Set the display theme
    fig.update_layout(
        showlegend=False,
        template="simple_white",
        height=500,
        title={
            'text': sample_name,
            'y': 0.9,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top',
        },
    )
    # Log scale the CAG size (horizontal) axis
    fig.update_xaxes(type="log")

    return fig

def calc_clr(v):
    """Calculate the CLR for a vector of abundances."""
    # Index of non-zero values
    ix = (v > 0)

    # Lowest non-zero value
    min_val = v[ix].min()

    # Geometric mean
    gmean_val = v[ix].apply(np.log10).mean()

    # Compute the CLR
    return v.clip(lower=min_val).apply(np.log10) - gmean_val

def plot_samples_pairwise(
    primary_sample_abund_df,
    primary_sample_name,
    secondary_sample_abund_df,
    secondary_sample_name,
    display_metric,
    cag_summary_df,
):
    # Make an abundance DataFrame to plot
    plot_df = pd.DataFrame({
        "primary": primary_sample_abund_df[primary_sample_name],
        "secondary": secondary_sample_abund_df[secondary_sample_name],
    })
    # Mask CAGs that are zero in both samples
    plot_df = plot_df.assign(
        max_val = plot_df.max(axis=1)
    ).query(
        "max_val > 0"
    ).drop(
        columns="max_val"
    )

    # If required, calculate the CLR
    if display_metric == "clr":
        
        plot_df = plot_df.apply(calc_clr)

    # Round the proportional abundance to 2 (clr) or 4 (prop) decimals
    plot_df = plot_df.apply(
        lambda c: c.apply(lambda v: round(v, 2 if display_metric == "clr" else 4))
    )

    # Add the CAG size
    plot_df = plot_df.assign(
        SIZE = cag_summary_df["size"]
    ).rename(columns={
        "SIZE": "Number of Genes",
    })

    # Reset the index to put the CAG label into a column
    plot_df = plot_df.reset_index()

    # Make the plot
    fig = px.scatter(
        plot_df,
        x="primary",
        y="secondary",
        labels=dict(
            primary="Abundance in {}".format(primary_sample_name),
            secondary="Abundance in {}".format(secondary_sample_name),
        ),
        hover_data=plot_df.columns.values,
    )

    # Set the display theme
    fig.update_layout(
        showlegend=False,
        template="simple_white",
        height=500,
    )

    return fig

####################
# ORDINATION GRAPH #
####################
def run_pca(df):
    """Dimensionality reduction with PCA."""

    # Initialize the PCA object
    pca = PCA()

    # Fit to the data
    pca.fit(df)

    # Make an output DataFrame
    return pd.DataFrame(
        pca.transform(
            df
        ),
        index=df.index.values,
        columns=[
            "PC%d (%s%s)" % (
                ix + 1,
                round(100 * r, 1) if r > 0.01 else "%.1E" % (100 * r),
                '%'
            )
            for ix, r in enumerate(
                pca.explained_variance_ratio_
            )
        ]
    )

def run_tsne(df, perplexity=30, n_components=2):
    """Dimensionality reduction with t-SNE."""

    # Initialize the TSNE object
    tsne = TSNE(
        n_components=n_components,
        perplexity=perplexity,
    )

    # Make an output DataFrame with the transformed data
    return pd.DataFrame(
        tsne.fit_transform(
            df
        ),
        index=df.index.values,
        columns=[
            "t-SNE %d" % (
                ix + 1
            )
            for ix in range(n_components)
        ]
    )
def update_ordination_graph(
    distances_df,
    algorithm,
    primary_pc,
    secondary_pc,
    perplexity,
    metadata,
    manifest_json,
    full_manifest_df,
):
    """Perform ordination and make the display plots."""

    # Get the filtered manifest from the browser
    plot_manifest_df = parse_manifest_json(manifest_json, full_manifest_df)

    if algorithm == "pca":
        plot_df = run_pca(
            distances_df.reindex(
                index=plot_manifest_df.index.values
            )
        )

    else:
        assert algorithm == "tsne", "Algorithm not found: %s" % algorithm

        plot_df = run_tsne(
            distances_df.reindex(
                index=plot_manifest_df.index.values
            ),
            perplexity=perplexity
        )

        # Always plot the first and second axes
        primary_pc = 1
        secondary_pc = 2

    # Set up the figure, which will be populated below
    fig = go.Figure()

    # The plot will depend on whether metadata has been selected
    if metadata == "none":

        # No metadata

        # Scatter plot
        fig.add_trace(
            go.Scatter(
                x=plot_df[plot_df.columns.values[primary_pc - 1]],
                y=plot_df[plot_df.columns.values[secondary_pc - 1]],
                ids=plot_df.index.values,
                text=plot_df.index.values,
                hovertemplate="%{id}<extra></extra>",
                mode="markers",
                marker_color="blue"
            )
        )

    else:
        # Add the indicated metadata to the plotting dataframe
        # This will also add a METADATA_COLOR column as appropriate
        plot_df = add_metadata_to_dataframe(
            plot_df, 
            plot_manifest_df, 
            metadata
        )

        # Limited number of unique metadata values
        if plot_df[metadata].unique().shape[0] <= 5:

            # Iterate over each of the metadata groups
            for metadata_label, metadata_plot_df in plot_df.groupby(metadata):
                # Scatter plot
                fig.add_trace(
                    go.Scatter(
                        x=metadata_plot_df[
                            plot_df.columns.values[primary_pc - 1]
                        ],
                        y=metadata_plot_df[
                            plot_df.columns.values[secondary_pc - 1]
                        ],
                        name=metadata_label,
                        ids=plot_df.index.values,
                        text=metadata_plot_df[metadata].apply(
                            lambda n: "{}: {}".format(metadata, n)
                        ),
                        hovertemplate="Sample: %{id}<br>%{text}<extra></extra>",
                        mode="markers",
                        marker_color=metadata_plot_df["METADATA_COLOR"].values[0],
                    )
                )

        else:

            # Scatter plot
            fig.add_trace(
                go.Scatter(
                    x=plot_df[plot_df.columns.values[primary_pc - 1]],
                    y=plot_df[plot_df.columns.values[secondary_pc - 1]],
                    ids=plot_df.index.values,
                    text=plot_df[metadata].apply(
                        lambda n: "{}: {}".format(metadata, n)
                    ),
                    hovertemplate="Sample: %{id}<br>%{text}<extra></extra>",
                    mode="markers",
                    marker_color=plot_df["METADATA_COLOR"],
                )
            )

        fig.update_yaxes(
            title_text=metadata,
        )


    fig.update_xaxes(
        title_text=plot_df.columns.values[primary_pc - 1],
    )
    fig.update_yaxes(
        title_text=plot_df.columns.values[secondary_pc - 1],
    )

    fig.update_layout(
        showlegend=False,
        template="simple_white",
        height=500,
    )

    return fig

def print_anosim(
    distances_df,
    metadata,
    manifest_json,
    full_manifest_df,
    permutations=9999,
):
    """Run anosim and return a Markdown summary."""

    # Get the filtered manifest from the browser
    plot_manifest_df = parse_manifest_json(
        manifest_json, 
        full_manifest_df
    ).reset_index(
    ).reindex(
        columns = ["specimen", metadata]
    ).drop_duplicates(
    ).set_index(
        "specimen"
    )

    # Remove any samples with NaN for this field
    samples_to_analyze = plot_manifest_df[metadata].dropna().index.values

    # Filter down the distance matrix and run permanova
    r = anosim(
        DistanceMatrix(
            distances_df.reindex(
                index=samples_to_analyze,
                columns=samples_to_analyze,
            ).values
        ),
        plot_manifest_df[metadata].reindex(index=samples_to_analyze),
        permutations=permutations
    )

    return dcc.Markdown("""
        _ANOSIM_ ([ref](http://scikit-bio.org/docs/0.2.3/generated/generated/skbio.stats.distance.anosim.html)):

        * R: {:.2} (Range: -1 to 1)
        * p: {:.2E}
        * Permutations: {:,}
        * Sample size: {:,}
        * Number of groups: {:,}
        """.format(
            r["test statistic"],
            r["p-value"],
            r["number of permutations"],
            int(r["sample size"]),
            int(r["number of groups"]),
        ))


####################
# CAG SUMMARY CARD #
####################
def draw_cag_summary_graph_hist(
    cag_summary_df,
    metric_primary,
    nbinsx,
    log_scale,
    metric,
):
    # Apply the filters
    plot_df = cag_summary_df.assign(
        CAGs = 1
    )

    axis_names = {
        "CAG": "CAG ID",
        "size": "Number of Genes (log10)",
        "mean_abundance": "Mean Abundance",
        "std_abundance": "Std. Abundance",
        "prevalence": "Prevalence",
        "entropy": "Entropy",
    }

    # Draw a histogram
    fig = px.histogram(
        plot_df,
        x=metric_primary if metric_primary != "size" else "size_log10",
        y="CAGs" if metric == "cags" else "size",
        histfunc="sum",
        nbins=nbinsx,
    )

    if metric == "cags":
        ylabel = "Total number of CAGs per bin"

    else:
        ylabel = "Total number of genes per bin"

    fig.update_layout(
        xaxis_title=axis_names[metric_primary],
        yaxis_title=ylabel,
        template="simple_white",
        height=400,
        width=600,
    )

    # Apply the log transform
    if log_scale == "on":
        fig.update_yaxes(type="log")

    return fig


###############
# CAG HEATMAP #
###############
def draw_cag_abundance_heatmap(
    cag_abund_df,
    metadata_selected,
    abundance_metric,
    cluster_by,
    taxa_rank,
    manifest_json,
    full_manifest_df,
    cag_taxa_dict,
    cag_estimate_dict,
    figure_width = 800,
    figure_height = 800,
):
    # Get the filtered manifest from the browser
    plot_manifest_df = parse_manifest_json(manifest_json, full_manifest_df)

    # Sort the manifest by the indicated fields
    if len(metadata_selected) > 0:
        plot_manifest_df = plot_manifest_df.sort_values(
            by=metadata_selected,
        ).reindex(
            columns=metadata_selected[::-1],
        )
    
    # Subset the CAG abundances to just those selected samples
    plot_df = cag_abund_df.reindex(
        index=plot_manifest_df.index
    )

    if abundance_metric in ["log10", "zscore"]:
        # Transform to log10 relative abundance

        # First find the lowest non-zero value
        lowest_value = plot_df.apply(
            lambda c: c[c > 0].min()
        ).min()

        # Now transform and set the floor
        plot_df = plot_df.clip(
            lower=lowest_value
        ).applymap(
            np.log10
        ).applymap(
            lambda i: round(i, 1)
        )
    
    if abundance_metric == "zscore":
        # Transform into the Z-score per-sample
        plot_df = (plot_df - plot_df.mean()) / plot_df.std()

    # If selected, cluster the specimens by CAG abundance
    if cluster_by == "cag":
        plot_df = plot_df.reindex(
            index=plot_df.index.values[
                leaves_list(
                    linkage(
                        plot_df,
                        method="ward"
                    )
                )
            ]
        )
        plot_manifest_df = plot_manifest_df.reindex(
            index=plot_df.index
        )

    # Group the CAGs based on their abundance
    plot_df = plot_df.reindex(
        columns=plot_df.columns.values[
            leaves_list(
                linkage(
                    plot_df.T,
                    method="ward"
                )
            )
        ]
    )

    # Depending on whether metadata or taxonomic information has
    # been provided, the plot will be set up in different ways
    # Metadata - & taxonomic annotations - : single plot
    # Metadata + & taxonomic annotations - : two plots, metadata on top of abund
    # Metadata - & taxonomic annotations + : two plots, tax to the right of abund
    # Metadata + & taxonomic annotations + : three plots, combining tax and metadata

    has_metadata = len(metadata_selected) > 0

    # Only plot the taxonomic annotation if we have sufficient taxonomic information
    has_taxonomy = any(
        cag_taxa_dict.get(cag_id) is not None
        for cag_id in plot_df.columns.values
    )

    # Set the mouseover text template
    if abundance_metric == "raw":
        hovertemplate = "%{y}<br>%{x}<br>Rel. Abund.: %{z}<extra></extra>"
    elif abundance_metric == "log10":
        hovertemplate = "%{y}<br>%{x}<br>Rel. Abund. (log10): %{z}<extra></extra>"
    elif abundance_metric == "zscore":
        hovertemplate = "%{y}<br>%{x}<br>Rel. Abund. (z-score): %{z}<extra></extra>"

    # No metadata
    if has_metadata is False:

        # No taxonomic information
        if has_taxonomy is False:

            # No estimated coefficients of association
            if cag_estimate_dict is None:

                # Make a very simple plot
                fig = go.Figure(
                    data=draw_cag_abund_heatmap_panel(
                        plot_df,
                        hovertemplate=hovertemplate
                    ),
                )

            else:
                
                # Only add the estimated coefficient
                fig = draw_cag_abund_heatmap_with_estimate(
                    plot_df, cag_estimate_dict,
                    hovertemplate=hovertemplate
                )

        else:

            # We have taxonomic information, but no estimated coefficients
            if cag_estimate_dict is None:

                fig = draw_cag_abund_heatmap_with_tax(
                    plot_df, cag_taxa_dict, taxa_rank,
                    hovertemplate=hovertemplate
                )

            else:

                # We have taxonomic information and estimated coefficients
                fig = draw_cag_abund_heatmap_with_tax_and_estimates(
                    plot_df, cag_taxa_dict, taxa_rank, cag_estimate_dict,
                    hovertemplate=hovertemplate
                )

    else:
        # We have metadata

        # No taxonomic information
        if has_taxonomy is False:

            # No estimated coefficients of association
            if cag_estimate_dict is None:

                # Make a plot with just metadata
                fig = draw_cag_abund_heatmap_with_metadata(
                    plot_df, plot_manifest_df, metadata_selected,
                    hovertemplate=hovertemplate
                )

            else:
                # We have estimated coefficients
                fig = draw_cag_abund_heatmap_with_metadata_and_estimate(
                    plot_df,
                    plot_manifest_df,
                    metadata_selected,
                    cag_estimate_dict,
                    hovertemplate=hovertemplate,
                )
        
        else:

            # We have taxonomic information

            # No estimated coefficients of association
            if cag_estimate_dict is None:

                # Make a plot with metadata and taxonomic annotations
                fig = draw_cag_abund_heatmap_with_metadata_and_tax(
                    plot_df, 
                    plot_manifest_df, 
                    metadata_selected, 
                    cag_taxa_dict,
                    taxa_rank,
                    hovertemplate=hovertemplate,
                )

            else:

                # We will make a plot with all three: metadata, taxa, and estimates
                fig = draw_cag_abund_heatmap_with_metadata_tax_and_estimates(
                    plot_df, 
                    plot_manifest_df, 
                    metadata_selected, 
                    cag_taxa_dict,
                    taxa_rank,
                    cag_estimate_dict,
                    hovertemplate=hovertemplate,
                )

    fig.update_layout(
        width=figure_width,
        height=figure_height,
        template="simple_white",
    )
    return fig


def draw_cag_abund_heatmap_with_tax(
    cag_abund_df, 
    cag_tax_dict,
    taxa_rank,
    hovertemplate = "Specimen: %{y}<br>CAG: %{x}<br>Rel. Abund.: %{z}<extra></extra>",
):

    # Make a plot with two panels, one over the other, sharing the x-axis

    # The taxa plot will be very small
    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True,
        row_heights=[
            0.05, 0.95
        ],
        vertical_spacing=0.005,
    )

    # Plot the abundances on the left
    fig.add_trace(
        draw_cag_abund_heatmap_panel(cag_abund_df, hovertemplate=hovertemplate), 
        row=2, 
        col=1
    )

    # Plot the taxonomic annotations on the right
    fig.add_trace(
        draw_cag_abund_taxon_panel(cag_tax_dict, taxa_rank, cag_abund_df.columns.values), 
        row=1, 
        col=1
    )
    # Rotate the angle of the x-tick labels
    fig.update_xaxes(tickangle=90)

    return fig

def draw_cag_abund_heatmap_with_estimate(
    cag_abund_df, 
    cag_annot_dict,
    hovertemplate = "Specimen: %{y}<br>CAG: %{x}<br>Rel. Abund.: %{z}<extra></extra>",
):

    # Make a plot with two panels, side-by-side, sharing the y-axis

    # The estimate plot will be smaller
    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True,
        row_heights=[
            0.15, 0.85
        ],
        vertical_spacing=0.005,
    )

    # Plot the abundances on the left
    fig.add_trace(
        draw_cag_abund_heatmap_panel(cag_abund_df, hovertemplate=hovertemplate), 
        row=2, col=1
    )

    # Plot the estimated coefficients on the right
    fig.add_trace(
        draw_cag_estimate_panel(
            cag_annot_dict, 
            cag_abund_df.columns.values,
            orientation="horizontal"
        ), 
        row=1, 
        col=1,
    )
    # Rotate the angle of the x-tick labels
    fig.update_xaxes(tickangle=90)

    return fig

def draw_cag_abund_heatmap_with_metadata(
    cag_abund_df, 
    plot_manifest_df, 
    metadata_selected,
    hovertemplate = "Specimen: %{y}<br>CAG: %{x}<br>Rel. Abund.: %{z}<extra></extra>",
):

    # Make a plot with two panels, side by side, sharing the y-axis

    # The relative height of the subplots will be set dynamically
    metadata_height = max(
        0.1,
        min(
            0.5,
            len(metadata_selected) / 20.
        ),
    )
    fig = make_subplots(
        rows=1, cols=2, shared_yaxes=True,
        column_widths=[
            1 - metadata_height,
            metadata_height,
        ],
        horizontal_spacing=0.01,
    )

    # Plot the metadata on the top
    fig.add_trace(
        draw_metadata_heatmap_panel(plot_manifest_df), row=1, col=2
    )

    # Plot the abundances on the bottom
    fig.add_trace(
        draw_cag_abund_heatmap_panel(cag_abund_df, hovertemplate=hovertemplate), row=1, col=1
    )

    return fig

def draw_cag_abund_heatmap_with_metadata_and_tax(
    cag_abund_df, 
    plot_manifest_df, 
    metadata_selected, 
    cag_tax_dict, 
    taxa_rank,
    hovertemplate="Specimen: %{y}<br>CAG: %{x}<br>Rel. Abund.: %{z}<extra></extra>",
):

    # Make a plot with four panels:
    # taxon - blank
    # cag-abun - metadata

    # The relative height of the subplots will be set dynamically
    metadata_height = max(
        0.1,
        min(
            0.5,
            len(metadata_selected) / 20.
        ),
    )

    data = [
        draw_cag_abund_heatmap_panel(
            cag_abund_df, 
            hovertemplate = hovertemplate,
            xaxis = "x",
            yaxis = "y"
        ),
        draw_cag_abund_taxon_panel(
            cag_tax_dict, 
            taxa_rank, 
            cag_abund_df.columns.values,
            xaxis = "x",
            yaxis = "y2"
        ),
        draw_metadata_heatmap_panel(
            plot_manifest_df,
            xaxis = "x2",
            yaxis = "y"
        )
    ]

    layout = go.Layout(
        yaxis = dict(
            domain = [0., 0.94]
        ),
        yaxis2 = dict(
            domain = [0.95, 1.]
        ),
        xaxis = dict(
            domain = [0., 0.99 - metadata_height]
        ),
        xaxis2 = dict(
            domain = [1.01 - metadata_height, 1.0]
        ),
    )

    fig = go.Figure(
        data=data,
        layout=layout
    )
    # Rotate the angle of the x-tick labels
    fig.update_xaxes(tickangle=90)

    return fig

def draw_cag_abund_heatmap_with_metadata_and_estimate(
    cag_abund_df, 
    plot_manifest_df, 
    metadata_selected, 
    cag_annot_dict, 
    hovertemplate="Specimen: %{y}<br>CAG: %{x}<br>Rel. Abund.: %{z}<extra></extra>",
):

    # Make a plot with four panels:
    # estimate - blank
    # cag-abun - metadata

    # First set up the data which goes into the plots
    data = [
        draw_cag_abund_heatmap_panel(
            cag_abund_df, 
            hovertemplate = hovertemplate,
            xaxis = "x",
            yaxis = "y",
        ),
        draw_metadata_heatmap_panel(
            plot_manifest_df,
            xaxis = "x2",
            yaxis = "y"
        ),
        draw_cag_estimate_panel(
            cag_annot_dict, 
            cag_abund_df.columns.values,
            xaxis = "x",
            yaxis = "y2",
            orientation = "horizontal"
        ),
    ]

    # Dynamically set the amount of plot area taken up by metadata
    metadata_height = max(
        0.1,
        min(
            0.5,
            len(metadata_selected) / 20.
        ),
    )

    # Now set the relative plot area taken up by each axis
    # The proportion taken up by each plot also includes 2% internal padding
    layout = go.Layout(
        yaxis=dict(
            domain=[0, 0.84]
        ),
        xaxis=dict(
            domain=[0, 0.99 - metadata_height]
        ),
        yaxis2=dict(
            domain=[0.86, 1.0]
        ),
        xaxis2=dict(
            domain=[1.01 - metadata_height, 1]
        ),
    )

    # Make the figure
    fig = go.Figure(
        data=data,
        layout=layout
    )

    # Rotate the angle of the x-tick labels
    fig.update_xaxes(tickangle=90)


    return fig

def draw_cag_abund_heatmap_with_tax_and_estimates(
    cag_abund_df, 
    cag_tax_dict,
    taxa_rank,
    cag_annot_dict,
    hovertemplate="Specimen: %{y}<br>CAG: %{x}<br>Rel. Abund.: %{z}<extra></extra>",
):

    # Make a plot with three panels:
    # cag-abun - taxa - estimate

    data = [
        draw_cag_abund_heatmap_panel(
            cag_abund_df, 
            hovertemplate=hovertemplate,
            yaxis="y",
        ),
        draw_cag_abund_taxon_panel(
            cag_tax_dict, 
            taxa_rank, 
            cag_abund_df.columns.values,
            yaxis="y2",
        ),
        draw_cag_estimate_panel(
            cag_annot_dict, 
            cag_abund_df.columns.values,
            orientation="horizontal",
            yaxis="y3",
        )
    ]

    layout = go.Layout(
        yaxis=dict(
            domain=[0, 0.795]
        ),
        yaxis2=dict(
            domain=[0.805, 0.845]
        ),
        yaxis3=dict(
            domain=[0.855, 1.0]
        ),
    )

    fig = go.Figure(
        data=data, layout=layout
    )

    # Rotate the angle of the x-tick labels
    fig.update_xaxes(tickangle=90)

    return fig


def draw_cag_abund_heatmap_with_metadata_tax_and_estimates(
    cag_abund_df, 
    plot_manifest_df, 
    metadata_selected, 
    cag_tax_dict,
    taxa_rank,
    cag_annot_dict,
    hovertemplate="Specimen: %{y}<br>CAG: %{x}<br>Rel. Abund.: %{z}<extra></extra>",
):

    # Make a plot with six panels:
    # estimate - blank
    # taxa - blank
    # cag-abund - metadata

    # The relative height of the subplots will be set dynamically
    metadata_height = max(
        0.1,
        min(
            0.5,
            len(metadata_selected) / 20.
        ),
    )

    data = [
        draw_cag_abund_heatmap_panel(
            cag_abund_df, 
            hovertemplate=hovertemplate,
            xaxis="x",
            yaxis="y",
        ),
        draw_cag_abund_taxon_panel(
            cag_tax_dict, 
            taxa_rank, 
            cag_abund_df.columns.values,
            xaxis="x",
            yaxis="y2",
        ),
        draw_cag_estimate_panel(
            cag_annot_dict, 
            cag_abund_df.columns.values,
            orientation="horizontal",
            xaxis="x",
            yaxis="y3"
        ),
        draw_metadata_heatmap_panel(
            plot_manifest_df,
            xaxis="x2",
            yaxis="y"
        )
    ]

    layout = go.Layout(
        yaxis = dict(domain=[0, 0.8]),
        yaxis2 = dict(domain=[0.81, 0.85]),
        yaxis3 = dict(
            domain=[0.86, 1.0],
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor='Grey',
        ),
        xaxis = dict(domain=[0, 0.99 - metadata_height]),
        xaxis2 = dict(domain=[1 - metadata_height, 1.]),
    )

    fig = go.Figure(
        data=data,
        layout=layout
    )

    # Rotate the angle of the x-tick labels
    fig.update_xaxes(tickangle=90)

    return fig

def draw_cag_abund_taxon_panel(
    cag_tax_dict, 
    taxa_rank,
    cag_order,
    xaxis="x",
    yaxis="y",
):

    # For each CAG, pick out the top hit
    summary_df = pd.DataFrame([
        summarize_cag_taxa(cag_id, cag_tax_df)
        for cag_id, cag_tax_df in cag_tax_dict.items()
    ])

    # Format the taxa name as a scalar for the color value
    summary_df = summary_df.assign(
        name_scalar = summary_df["name"].apply(
            dict(zip(
                summary_df["name"].unique(),
                range(summary_df["name"].unique().shape[0])
            )).get
        )
    )

    # Reorder the display to match the abundances
    summary_df = summary_df.set_index(
        "CAG"
    ).reindex(
        index=cag_order
    ).reset_index()

    return go.Heatmap(
        y=[taxa_rank],
        x=["CAG {}".format(i) for i in summary_df["CAG"].values],
        z=summary_df.reindex(columns=["name_scalar"]).T.values,
        text=summary_df.reindex(columns=["label"]).T.values,
        showscale=False,
        colorscale='Viridis',
        hovertemplate="%{x}<br>%{text}<extra></extra>",
        xaxis=xaxis,
        yaxis=yaxis,
    )


def draw_cag_estimate_panel(
    cag_annot_dict,
    cag_order,
    xaxis = "x",
    yaxis = "y",
    orientation = "vertical"
):
    """Render the subplot with the estimated coefficients for each CAG."""
    # Explicitly order the values for plotting
    plot_values = {
        col_name: [
            cag_annot_dict[int(cag_id)][col_name]
            for cag_id in cag_order
        ]
        for col_name in ["estimate", "std_error"]
    }
    # Render the hover text
    hovertext = [
        "Wald: {:.2}".format(
            cag_annot_dict[int(cag_id)]["wald"]
        )
        for cag_id in cag_order
    ]

    # Set the values that will be plotted
    estimate_values = plot_values["estimate"]
    index_values = [
        "CAG {}".format(cag_id)
        for cag_id in cag_order
    ]
    error_values = dict(
        type='data',
        array=plot_values["std_error"],
        visible=True
    )

    if orientation == "vertical":
        return go.Scatter(
            x = estimate_values,
            y = index_values,
            error_x = error_values,
            ids = cag_order,
            text = hovertext,
            hovertemplate = "CAG %{id}<br>Estimated Coefficient: %{x}<br>%{text}<extra></extra>",
            mode = "markers",
            marker_color = "LightSkyBlue",
            xaxis = xaxis,
            yaxis = yaxis,
            showlegend=False
        )
    else:
        assert orientation == "horizontal"

        return go.Scatter(
            y = estimate_values,
            x = index_values,
            error_y = error_values,
            ids = cag_order,
            text = hovertext,
            hovertemplate = "CAG %{id}<br>Estimated Coefficient: %{y}<br>%{text}<extra></extra>",
            mode = "markers",
            marker_color = "LightSkyBlue",
            xaxis = xaxis,
            yaxis = yaxis,
            showlegend=False
        )


def summarize_cag_taxa(cag_id, cag_tax_df):
    """Helper function to summarize the top hit at a given rank."""

    # If there are no hits at this level, return None
    if cag_tax_df is None:
        return {
            "CAG": cag_id,
            "name": 'none',
            "label": "No genes assigned at this level"
        }

    # Return the top hit
    return {
        "CAG": cag_id,
        "name": cag_tax_df["name"].values[0],
        "label": "{}<br>{:,} genes assigned".format(
            cag_tax_df["name"].values[0],
            int(cag_tax_df["count"].values[0])
        )
    }


def draw_metadata_heatmap_panel(
    plot_manifest_df,
    hovertemplate="Specimen: %{y}<br>Label: %{x}<br>Value: %{text}<extra></extra>",
    xaxis="x",
    yaxis="y",
):

    return go.Heatmap(
        z=plot_manifest_df.apply(
            lambda r: r.apply(dict(zip(
                r.drop_duplicates().sort_values(), 
                np.arange(0, 1, 1 / r.unique().shape[0]))).get
            ),
        ).values,
        text=plot_manifest_df.values,
        x=["{}".format(i) for i in plot_manifest_df.columns.values],
        y=["Specimen: {}".format(i) for i in plot_manifest_df.index.values],
        colorscale='Viridis',
        showscale=False,
        hovertemplate=hovertemplate,
        xaxis=xaxis,
        yaxis=yaxis,
    )

def draw_cag_abund_heatmap_panel(
    cag_abund_df,
    hovertemplate = "%{y}<br>%{x}<br>Rel. Abund.: %{z}<extra></extra>",
    xaxis="x",
    yaxis="y",
):
    return go.Heatmap(
        z=cag_abund_df.values,
        x=["CAG {}".format(i) for i in cag_abund_df.columns.values],
        y=["Specimen: {}".format(i) for i in cag_abund_df.index.values],
        colorbar={"title": "Abundance (log10)"},
        colorscale='blues',
        hovertemplate=hovertemplate,
        xaxis = xaxis,
        yaxis = yaxis,

    )

##########################
# CAG ANNOTATION HEATMAP #
##########################
def draw_cag_annotation_heatmap(
    cag_annot_df,
    annotation_type,
    enrichment_df,
    cag_estimate_dict,
    n_annots,
    annotation_names,
    figure_width=1000,
    figure_height=800,
):
    """Render the heatmap with CAG annotations."""
    if cag_annot_df is None:
        return go.Figure()


    # If the annotation is taxonomic, then we will use a dedicated function to prepare the data for plotting
    if annotation_type == "taxonomic":


        # Three data structures are needed to plot the full taxonomy, the number of counts
        # at each terminal node, the taxonomy linking each terminal node to its parents,
        # and (for convenience) a list of ancestors for each terminal node
        plot_df, tax_df = format_taxonomic_annot_df(
            cag_annot_df,
            enrichment_df,
            n_annots
        )

    # Otherwise, just format the annotation by pivoting to wide format and selecing the top
    # N annotations by either frequency or the enrichment absolute Wald (if provided)
    else:


        # Format the annotation table
        plot_df = format_annot_df(
            cag_annot_df, 
            annotation_type, 
            enrichment_df, 
            n_annots
        )

        # Replace the annotation names, if provided
        if annotation_names is not None:
            plot_df = plot_df.rename(
                columns=annotation_names.get
            )

        # Sort the rows and columns with linkage clustering
        plot_df = cluster_dataframe(plot_df)

    # Lacking functional/taxonomic enrichments or CAG-level estimated coefficients of association
    if enrichment_df is None and (cag_estimate_dict is None or len(cag_estimate_dict) == 0):

        # If the plot type is Taxonomic
        if annotation_type == "taxonomic":

            # Make a plot with the proportion of genes assigned, alongside the taxonomy
            fig = plot_taxonomic_annotations_without_enrichment(
                plot_df,
                tax_df,
            )

        # all other plot types
        else:
            # Make a very simple plot
            fig = go.Figure(
                data=draw_cag_annotation_panel(
                    plot_df
                )
            )

    # We have CAG association metrics, but no label enrichment
    elif enrichment_df is None and cag_estimate_dict is not None and len(cag_estimate_dict) > 0:

        if annotation_type == "taxonomic":

            fig = plot_taxonomic_annotations_with_cag_associations_only(
                plot_df,
                tax_df,
                cag_estimate_dict
            )

        else:

            # Just plot the association metrics for the CAGs
            fig = draw_cag_annot_heatmap_with_cag_estimates(
                plot_df,
                cag_estimate_dict,
            )

    # We have information on parameter association for both CAGs and annotation labels
    else:

        # Plot the full taxonomy with a dedicated function to render the etree
        if annotation_type == "taxonomic":
            fig = plot_taxonomic_annotations_with_enrichment(
                plot_df,
                tax_df,
                cag_estimate_dict,
                enrichment_df
            )

        # All other annotation types
        else:

            # Make a plot with association metrics on both the rows and columns        
            # Four panels, side-by-side, sharing the x-axis and y-axis

            fig = draw_cag_annot_heatmap_with_cag_estimates_and_enrichments(
                plot_df,
                cag_estimate_dict,
                enrichment_df,
            )


    fig.update_layout(
        width=figure_width,
        height=figure_height,
        template="simple_white",
    )

    return fig

def draw_path_to_root_tree(path_to_root_df):

    # Set up the figure
    fig = go.Figure()

    # Keep track of where each taxon has been plotted
    taxon_position = {}

    # Add to the figure by walking over each of the ranks in reverse order
    for rank_ix, rank in enumerate(path_to_root_df.columns.values[::-1]):

        # Iterate over each organism at this level
        for org_name in path_to_root_df[rank].dropna().unique():
            
            # Set the vertical position based on the position in the DataFrame
            y = np.mean([
                ix
                for ix, v in enumerate(path_to_root_df[rank].values)
                if v == org_name
            ])

            # Set the horizontal position based on the rank
            x = path_to_root_df.shape[1] - rank_ix

            # Save the position
            taxon_position["{}-{}".format(rank, org_name)] = (x, y)

            # Draw a point for this organism
            fig.add_trace(
                go.Scatter(
                    x = [x],
                    y = [y],
                    ids = ["{}: {}".format(rank, org_name)],
                    hovertemplate = "%{id}<extra></extra>",
                    showlegend=False,
                    mode="markers",
                    marker=dict(color="Blue", size=6)
                )
            )

    # There will be some taxa in the table which do not have terminal assignments
    # To capture those in the tree, we will add them to the left-most column
    # purely for the purpose of drawing lines to ancestors
    first_rank = path_to_root_df.columns.values[0]
    for ix, v in path_to_root_df[first_rank].items():
        if v is None:
            # Add a value to the table
            path_to_root_df.loc[ix, first_rank] = ix
            # Also add a pseudo position for plotting dashed lines
            taxon_position[
                "{}-{}".format(
                    first_rank, ix
                )
            ] = (
                0, 
                np.where(path_to_root_df.index.values == ix)[0][0]
            )
    
    # For each taxon in the tree, add a line to its parent
    for rank_ix, rank in enumerate(path_to_root_df.columns.values[:-1]):

        # Iterate over each unique organism at this level
        for org_name in path_to_root_df[rank].dropna().unique():

            # Find the parent of this organism in the table
            parent_taxa = None
            for parent_rank, parent_orgs in path_to_root_df.loc[
                path_to_root_df[rank] == org_name
            ].drop(
                columns=path_to_root_df.columns.values[:(rank_ix + 1)]
            ).iteritems():

                # Check to see that there is a unique parent taxon
                if parent_orgs.isnull().sum() == 0 and parent_orgs.unique().shape[0] == 1:
                    parent_taxa = "{}-{}".format(
                        parent_rank,
                        parent_orgs.values[0]
                    )
                    break

            # Get the position for this organism
            org_position = taxon_position["{}-{}".format(rank, org_name)]

            # Get the position for the parent
            if parent_taxa is None:
                # If there is no parent in the table, draw a line back to the root
                parent_position = (path_to_root_df.shape[1], org_position[1])
            else:
                # Look up the position used earlier
                parent_position = taxon_position[parent_taxa]
            
            # Add to the figure
            fig.add_trace(
                go.Scatter(
                    x = [org_position[0], parent_position[0] - 0.5, parent_position[0]],
                    y = [org_position[1], org_position[1], parent_position[1]],
                    showlegend=False,
                    mode="lines",
                    hoverinfo='skip',
                    line=dict(color='grey', width=1)
                )
            )


    # Set the tick values and text to sync with other subplots
    fig.update_layout(
        yaxis = dict(
            tickmode = 'array',
            tickvals = list(range(path_to_root_df.shape[0])),
            ticktext = path_to_root_df.index.values
        )
    )

    return fig


def plot_taxonomic_annotations_with_enrichment(
    plot_df,
    tax_df,
    cag_estimate_dict,
    enrichment_df
):
    """Plot the number of genes assigned to taxa alongside the tree, with association metrics."""

    # Make a DataFrame with the ancestors of each terminal node, which can be used for clustering
    path_to_root_df = format_path_to_root_df(plot_df.index.values, tax_df)
    
    # Make the dendrogram
    fig = draw_path_to_root_tree(
        path_to_root_df
    )

    # Get the order of taxa from the dendrogram
    dendro_leaves = fig['layout']['yaxis']['ticktext']
    dendro_ticks = fig['layout']['yaxis']['tickvals']

    # Add the heatmap panel
    fig.add_trace(
        draw_taxonomic_annotation_heatmap_panel(
            plot_df,
            tax_df,
            dendro_leaves,
            dendro_ticks,
            xaxis="x3",
        )
    )

    # Add the CAG estimates
    fig.add_trace(
        draw_cag_estimate_panel(
            cag_estimate_dict,
            plot_df.columns.values,
            yaxis="y2",
            xaxis="x3",
            orientation="horizontal"
        ),

    )

    fig.add_trace(
        draw_enrichment_estimate_panel(
            enrichment_df.loc[
                enrichment_df["label"].isin(dendro_leaves)
            ].set_index(
                "label"
            ),
            dendro_leaves,
            dendro_ticks,
            yaxis="y",
            xaxis="x2",
            orientation="vertical"
        )
    )

    # Edit yaxis for the dendrogram
    fig.update_layout(
        xaxis={
            'domain': [0.8, 1.0],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'showticklabels': False,
            'zeroline': False,
            'ticks': ""
        }
    )
    # Edit yaxis for the enrichment metrics, style the zeroline
    fig.update_layout(
        xaxis2={
            'domain': [0.7, 0.79],
            'zeroline': True,
            'zerolinewidth': 1,
            'zerolinecolor': 'Grey',
        }
    )
    # Edit yaxis for the heatmap and CAG association metrics
    fig.update_layout(
        xaxis3={
            'domain': [0, 0.69],
            'showline': False,
            'zeroline': False,
            'showspikes': True,
            'spikethickness': 2,
            'spikedash': "dot",
            'spikecolor': "#999999",
            'spikemode': "across",
            'spikesnap': 'cursor',
        }
    )
    # Edit xaxis shared by heatmap, enrichment values and dendrogram
    fig.update_layout(
        yaxis={
            'domain': [0, 0.9],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'ticks': "",
            'anchor': "x3",
            'showspikes': True,
            'spikethickness': 2,
            'spikedash': "dot",
            'spikecolor': "#999999",
            'spikemode': "across",
            'spikesnap': 'cursor',
        }
    )
    # Edit xaxis used for the CAG association metrics
    fig.update_layout(
        yaxis2={
            'domain': [0.91, 1],
            'anchor': "x3",
            'zeroline': True,
            'zerolinewidth': 1,
            'zerolinecolor': 'Grey',
        }
    )

    return fig

def plot_taxonomic_annotations_with_cag_associations_only(
    plot_df,
    tax_df,
    cag_estimate_dict
):
    """Plot the number of genes assigned to taxa alongside the tree, with CAG association metrics only."""

    # Make a DataFrame with the ancestors of each terminal node, which can be used for clustering
    path_to_root_df = format_path_to_root_df(plot_df.index.values, tax_df)

    # Make the dendrogram
    fig = draw_path_to_root_tree(
        path_to_root_df
    )
    # Get the order of taxa from the dendrogram
    dendro_leaves = fig['layout']['yaxis']['ticktext']
    dendro_ticks = fig['layout']['yaxis']['tickvals']

    # Add the heatmap panel
    fig.add_trace(
        draw_taxonomic_annotation_heatmap_panel(
            plot_df,
            tax_df,
            dendro_leaves,
            dendro_ticks,
            xaxis="x2",
        )
    )

    # Add the CAG estimates
    fig.add_trace(
        draw_cag_estimate_panel(
            cag_estimate_dict,
            plot_df.index.values,
            xaxis="x2",
            yaxis="y2",
            orientation="horizontal"
        ),

    )

    # Edit yaxis for the dendrogram
    fig.update_layout(
        yaxis={
            'domain': [0.8, 1.0],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'showticklabels': False,
            'zeroline': False,
            'ticks': ""
        }
    )
    # Edit yaxis for the heatmap and CAG association metrics
    fig.update_layout(
        yaxis2={
            'domain': [0, 0.79],
            'showline': False,
            'zeroline': False,
            'showspikes': True,
            'spikethickness': 2,
            'spikedash': "dot",
            'spikecolor': "#999999",
            'spikemode': "across",
            'spikesnap': 'cursor',
        }
    )
    # Edit xaxis shared by heatmap, enrichment values and dendrogram
    fig.update_layout(
        xaxis={
            'domain': [0, 0.9],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'ticks': "",
            'anchor': "y2",
            'showspikes': True,
            'spikethickness': 2,
            'spikedash': "dot",
            'spikecolor': "#999999",
            'spikemode': "across",
            'spikesnap': 'cursor',
        }
    )
    # Edit xaxis used for the CAG association metrics
    fig.update_layout(
        xaxis2={
            'domain': [0.91, 1],
            'anchor': "y2",
            'zeroline': True,
            'zerolinewidth': 1,
            'zerolinecolor': 'Grey',
        }
    )

    return fig


def plot_taxonomic_annotations_without_enrichment(
    plot_df,
    tax_df,
):
    """Plot the number of genes assigned to taxa alongside the tree."""
    # Make a DataFrame with the ancestors of each terminal node, which can be used for clustering
    path_to_root_df = format_path_to_root_df(plot_df.index.values, tax_df)

    # Make the dendrogram
    fig = draw_path_to_root_tree(
        path_to_root_df
    )

    # Get the order of taxa from the dendrogram
    dendro_leaves = fig['layout']['yaxis']['ticktext']
    dendro_ticks = fig['layout']['yaxis']['tickvals']

    # Add the heatmap panel
    fig.add_trace(
        draw_taxonomic_annotation_heatmap_panel(
            plot_df,
            tax_df,
            dendro_leaves,
            dendro_ticks,
            xaxis="x2",
        )
    )

    # Edit xaxis for the dendrogram
    fig.update_layout(
        xaxis={
            'domain': [0.8, 1.0],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'showticklabels': False,
            'zeroline': False,
            'ticks':""
        }
    )
    # Edit xaxis for the heatmap
    fig.update_layout(
        xaxis2={
            'domain': [0, 0.79],
            'showline': False,
            'zeroline': False,
        }
    )
    # Edit shared yaxis
    fig.update_layout(
        yaxis={
            'domain': [0, 1],
            'mirror': False,
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'ticks':"",
            'anchor': "x2",
            'showspikes': True,
            'spikethickness': 2,
            'spikedash': "dot",
            'spikecolor': "#999999",
            'spikemode': "across",
            'spikesnap': 'cursor',
        }
    )

    return fig


def draw_taxonomic_annotation_heatmap_panel(
    plot_df,
    tax_df,
    dendro_leaves,
    dendro_ticks,
    xaxis="x",
    yaxis="y"
):
    # Cluster the table with gene abundances
    plot_df = cluster_dataframe(plot_df)
    
    # Switch the rows and columns
    plot_df = plot_df.T

    # Rename the tax IDs to organism names
    plot_df = plot_df.rename(
        columns=tax_df["name"].get
    ).reindex(
        columns=dendro_leaves
    )

    # Make the table with text to display
    text_df = plot_df.apply(
        lambda c: pd.Series(
            ["CAG {}<br>Label: {}<br>Genes assigned: {:,}".format(
                c.name,
                org_name,
                int(count),
            ) for org_name, count in c.items()],
            index=c.index
        ),
        axis=1
    )

    # Scale to show as a proportion of the maximum number of genes assigned per CAG
    plot_df = plot_df.T.apply(
        lambda c: 100. * c / c.max() if c.max() > 0 else c
    ).T

    # Create the Heatmap
    return go.Heatmap(
        text=text_df.T.values,
        y=dendro_ticks,
        x=["CAG {}".format(cag_id) for cag_id in plot_df.index.values],
        z=plot_df.T.values,
        colorscale='Blues',
        colorbar={"title": "Percent of gene assignments"},
        hovertemplate="%{text}<extra></extra>",
        zmin=0.,
        zmax=100.,
        xaxis=xaxis,
        yaxis=yaxis,
    )


def format_path_to_root_df(
    tax_id_list, 
    tax_df, 
    rank_list = [
        "species", 
        "genus", 
        "family",
        "order",
        "class",
        "phylum",
    ]
):
    return pd.DataFrame({
        tax_id: {
            rank: anc_at_rank(tax_id, tax_df, rank)
            for rank in rank_list
        }
        for tax_id in tax_id_list
    }).T.rename(
        index=tax_df["name"].get
    ).reindex(
        columns=rank_list
    ).sort_values(
        by=rank_list[::-1]
    )


def draw_cag_annot_heatmap_with_cag_estimates_and_enrichments(
    plot_df,
    cag_estimate_dict,
    enrichment_df,
):

    data = [
        draw_cag_annotation_panel(
            plot_df,
            xaxis="x",
            yaxis="y"
        ),
        draw_cag_estimate_panel(
            cag_estimate_dict,
            plot_df.index.values,
            xaxis="x",
            yaxis="y2",
            orientation = "horizontal"
        ),
        draw_enrichment_estimate_panel(
            enrichment_df,
            plot_df.columns.values,
            plot_df.columns.values,
            xaxis="x2",
            yaxis="y",
            orientation="vertical"
        ),
    ]

    layout = go.Layout(
        yaxis=dict(
            domain=[0, 0.85],
            showspikes= True,
            spikethickness= 2,
            spikedash= "dot",
            spikecolor= "#999999",
            spikemode= "across",
            spikesnap="cursor"
        ),
        xaxis=dict(
            domain=[0, 0.85],
            showspikes=True,
            spikethickness=2,
            spikedash="dot",
            spikecolor="#999999",
            spikemode="across",
            spikesnap="cursor"
        ),
        yaxis2=dict(
            domain=[0.86, 1.0],
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor='Grey',
        ),
        xaxis2=dict(
            domain=[0.86, 1.0],
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor='Grey',
        ),
    )

    fig = go.Figure(data=data, layout=layout)

    return fig


def draw_cag_annot_heatmap_with_cag_estimates(
    plot_df,
    cag_estimate_dict,
):

    fig = make_subplots(
        rows=2, 
        cols=1, 
        shared_xaxes=True,
        row_heights=[
            0.15, 0.85
        ],
        horizontal_spacing=0.005,
    )

    # Plot the abundances on the left
    fig.add_trace(
        draw_cag_annotation_panel(
            plot_df
        ),
        row=2, 
        col=1
    )

    # Plot the estimates on the right
    fig.add_trace(
        draw_cag_estimate_panel(
            cag_estimate_dict,
            plot_df.index.values,
            orientation="horizontal",
        ),
        row=1,
        col=1,
    )
    fig.update_layout(
        xaxis2=dict(
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor='Grey',
        ),
        yaxis=dict(
            showspikes= True,
            spikethickness= 2,
            spikedash= "dot",
            spikecolor= "#999999",
            spikemode= "across",
            spikesnap="cursor"
        ),
    )

    return fig

def draw_enrichment_estimate_panel(
    enrichment_df,
    label_order,
    label_ids,
    xaxis="x",
    yaxis="y",
    orientation="vertical"
):
    """Render the subplot with the estimated coefficients for annotation label."""
    # Explicitly order the values for plotting
    plot_df = enrichment_df.reindex(
        index=label_order
    )
    
    # Render the hover text
    hovertext = plot_df["wald"].apply(
        "Wald: {:.2}".format
    )

    # Pull out the values which will be plotted
    estimate_values = plot_df["estimate"]
    error_values = dict(
        type='data',
        array=plot_df["std_error"],
        visible=True
    )

    # Trim the label names
    label_ids = [
        n[:30] + "..." if isinstance(n, str) and len(n) > 30 else n
        for n in label_ids
    ]

    if orientation == "horizontal":
        return go.Scatter(
            x=label_ids,
            y=estimate_values,
            error_y=error_values,
            ids=plot_df.index.values,
            text=hovertext,
            hovertemplate="%{id}<br>Estimated Coefficient: %{y}<br>%{text}<extra></extra>",
            mode="markers",
            marker_color="LightSkyBlue",
            xaxis=xaxis,
            yaxis=yaxis,
            showlegend=False,
        )
    else:
        assert orientation == "vertical"
        return go.Scatter(
            y=label_ids,
            x=estimate_values,
            error_x=error_values,
            ids=plot_df.index.values,
            text=hovertext,
            hovertemplate="%{id}<br>Estimated Coefficient: %{x}<br>%{text}<extra></extra>",
            mode="markers",
            marker_color="LightSkyBlue",
            xaxis=xaxis,
            yaxis=yaxis,
            showlegend=False,
        )

def format_taxonomic_annot_df(cag_annot_df, enrichment_df, n_annots):
    """Format the table of taxonomic assignment for a set of CAGs."""

    # Make a table with just the taxonomy
    tax_df = cag_annot_df.drop(
        columns=["count", "total", "CAG"]
    ).drop_duplicates(
    ).set_index(
        "tax_id"
    )

    # Calculate the number of counts assigned specifically to each taxon
    cag_annot_df = cag_annot_df.groupby(
        "CAG"
    ).apply(
        find_counts_at_taxon
    )

    # Keep track of how many nodes we're adding to the plot
    terminal_nodes = set([])
    
    # If the enrichment values are provided, add the top n_annots / 2
    if enrichment_df is not None:

        # Iterate over the most consistently associated labels
        for _, r in enrichment_df.sort_values(
            by="abs_wald",
            ascending=False
        ).iterrows():
            # If there is an exact match to this label, add it
            matched_org = tax_df.query(
                "rank == '{}'".format(r["rank"])
            ).query(
                "name == '{}'".format(r["label"])
            )

            if matched_org.shape[0] == 1:
                terminal_nodes.add(
                    matched_org.index.values[0],
                )

            # Stop once (or if) the limit on terminal nodes has been reached
            if len(terminal_nodes) >= int(n_annots / 2.):
                break

    # Now add more taxa according to their frequency in the underlying data
    for tax_id in cag_annot_df.loc[
        cag_annot_df["rank"].isin([
            "species", "genus", "family", "class", "order", "phylum"
        ])
    ].sort_values(
        by="specific_prop",
        ascending=False
    )[
        "tax_id"
    ].values:

        # Add this node to the plot
        terminal_nodes.add(
            tax_id
        )

        # Stop once (or if) the limit on terminal nodes has been reached
        if len(terminal_nodes) >= int(n_annots):
            break

    # For each CAG, we will format the wide table so that each row contains
    # a terminal node, and the value in each cell is the number of genes
    # which are assigned to that taxa or its children
    plot_df = cag_annot_df.loc[
        cag_annot_df["tax_id"].isin(terminal_nodes)
    ].reset_index(
        drop=True
    ).pivot_table(
        columns="CAG",
        index="tax_id",
        values="at_or_below"
    ).fillna(
        0
    )

    return plot_df, tax_df


def find_counts_at_taxon(cag_df):

    # Sum up the number of counts below each taxon
    counts_below = cag_df.query(
        "tax_id != parent"
    ).groupby(
        "parent"
    )[
        "count"
    ].sum(
    )

    # Make a vector with the number of counts below each taxon
    counts_below = cag_df["tax_id"].apply(
        lambda t: counts_below.get(t, 0)
    )
    # Count up the number at each taxon
    specific_count = cag_df["count"] - counts_below
    # Compute the proportion of all assignments
    specific_prop = specific_count / specific_count.sum()

    return cag_df.assign(
        counts_below = counts_below,
        specific_count = specific_count,
        at_or_below = counts_below + specific_count,
        specific_prop = specific_prop,
    )


def add_tax_node(tax_id, terminal_nodes, internal_nodes, tax_df):
    """Add a single tax ID and its ancestors, as appropriate."""

    # This tax ID is already present
    if tax_id in terminal_nodes or tax_id in internal_nodes:
        return

    # Add the tax ID as a terminal node
    terminal_nodes.add(tax_id)

    # Now walk up the taxonomy to the root
    for anc_tax_id in path_to_root(tax_id, tax_df):

        # If the ancestor tax ID is a terminal node, move it to the internal node set
        if anc_tax_id in terminal_nodes:
            terminal_nodes.remove(anc_tax_id)

        # Make sure to add this to the internal node set
        internal_nodes.add(anc_tax_id)


def anc_at_rank(tax_id, tax_df, rank):
    """Return the name of the parent of a taxon at a given rank."""

    # Check to see if we are already at this rank
    if tax_df.loc[tax_id, "rank"] == rank:
        return tax_df.loc[tax_id, "name"]

    # Otherwise, walk up the parents until you find it
    else:
        for anc_tax_id in path_to_root(tax_id, tax_df):
            # Check to see if we have reached the rank of interest
            if tax_df.loc[anc_tax_id, "rank"] == rank:
                return tax_df.loc[anc_tax_id, "name"]

    # If none was found, return None
    return


def path_to_root(tax_id, tax_df, max_iter=1000):
    """Parse the taxonomy to walk up to the root (will not yield the query itself)."""

    # Get the parent tax ID
    anc_tax_id = tax_df.loc[tax_id, "parent"]

    # Keep track of all of the tax IDs that we've visited
    visited_taxa = set([tax_id])

    # Walk up the taxonomy
    for _ in range(max_iter):

        # Stop walking up when we reach the end
        if anc_tax_id in visited_taxa:
            break
        elif anc_tax_id is None:
            break
        elif anc_tax_id == tax_id:
            break

        # Add the ancestor
        visited_taxa.add(anc_tax_id)

        # Yield the ancestor tax ID
        yield anc_tax_id

        # Now find the parent of this one
        if anc_tax_id not in tax_df.index.values:
            break

        # In the next iteration, process the parent tax ID
        tax_id = anc_tax_id
        anc_tax_id = tax_df.loc[tax_id, "parent"]
        

def format_annot_df(cag_annot_df, annotation_type, enrichment_df, n_annots):
    """Format the table of CAG annotations."""

    # If the annotations are functional, we can just pivot across those functions
    if annotation_type == "eggNOG_desc":

        wide_df = cag_annot_df.pivot_table(
            index="CAG",
            columns="label",
            values="count",
            aggfunc=sum,
        ).fillna(
            0
        )

    else:
        # Otherwise, the gene annotations are in tax ID space, while the annotation type is either 'species', 'genus', or 'family'
        
        # If we are including the non-specific taxa, use the 'consistent' number of hits, otherwise use 'count'
        wide_df = cag_annot_df.pivot_table(
            index="CAG",
            columns="name",
            values="count",
            aggfunc=sum,
        ).fillna(
            0
        )

    # If enrichments are provided, keep the top `n_annots` by enrichment (absolute Wald)
    if enrichment_df is not None:
        wide_df = wide_df.reindex(
            columns = enrichment_df.reindex(
                index=list(set(wide_df.columns.values))
            ).dropna(
            )[
                "abs_wald"
            ].sort_values(
                ascending=False
            ).head(
                n_annots
            ).index.values
        )

    # Otherwise keep the most abundant annotations (normalizing to the total number of genes per CAG)
    else:
        wide_df = wide_df.reindex(
            columns = (
                wide_df.T / wide_df.sum(axis=1)
            ).T.sum().sort_values(
                ascending=False
            ).head(
                n_annots
            ).index.values
        )

    return wide_df

def draw_cag_annotation_panel(
    plot_df,
    xaxis = "x",
    yaxis = "y",
):

    # Scale the assignments to the maximum for each CAG
    prop_df = 100 * (plot_df.T / plot_df.max(axis=1))

    # Format the mouseover text
    text_df = plot_df.apply(
        lambda c: [
            "CAG {}<br>{}<br>Genes assigned: {:,}".format(
                cag_id,
                c.name,
                int(ncounts),
            )
            for cag_id, ncounts in c.items()
        ]
    ).T

    return go.Heatmap(
        text=text_df.values,
        z=prop_df.values,
        x=["CAG {}".format(i) for i in plot_df.index.values],
        y=[
            n[:30] + "..." if len(n) > 30 else n
            for n in plot_df.columns.values
        ],
        colorbar={"title": "Percent of gene assignments"},
        colorscale='blues',
        hovertemplate = "%{text}<extra></extra>",
        zmin=0.,
        zmax=1.,
        xaxis=xaxis,
        yaxis=yaxis,
    )

#################
# VOLCANO GRAPH #
#################
def draw_volcano_graph(
    corncob_df,
    cag_summary_df,
    parameter, 
    comparison_parameter,
    cag_size_range,
    neg_log_pvalue_min, 
    fdr_on_off, 
):
    if corncob_df is None or neg_log_pvalue_min is None:
        return go.Figure()

    # Filter the corncob_df by CAG size
    corncob_df = corncob_df.assign(
        size = cag_summary_df["size"]
    ).query(
        "size >= {}".format(10**cag_size_range[0])
    ).query(
        "size <= {}".format(10**cag_size_range[1])
    )
    
    # If a comparison parameter was selected, plot the p-values against each other
    if comparison_parameter != "coef":
        return draw_double_volcano_graph(
            corncob_df,
            parameter,
            comparison_parameter,
            neg_log_pvalue_min,
            fdr_on_off,
        )

    # Subset to the pvalue threshold
    plot_df = corncob_df.query(
        "neg_log10_pvalue >= {}".format(neg_log_pvalue_min)
    )

    if fdr_on_off == "off":
        plot_y = "neg_log10_pvalue"
        hovertemplate = "CAG %{id}<br>Estimate: %{x}<br>p-value (-log10): %{y}<extra></extra>"
        yaxis_title = "p-value (-log10)"
    else:
        plot_y = "neg_log10_qvalue"
        hovertemplate = "CAG %{id}<br>Estimate: %{x}<br>q-value (-log10): %{y}<extra></extra>"
        yaxis_title = "q-value (-log10)"

    fig = go.Figure(
        data=go.Scattergl(
            x=plot_df["estimate"],
            y=plot_df[plot_y],
            ids=plot_df.index.values,
            text=plot_df.index.values,
            hovertemplate=hovertemplate,
            mode="markers",
            opacity=0.5,
        ),
    )

    fig.update_layout(
        xaxis_title="Estimated Coefficient",
        yaxis_title=yaxis_title,
        title={
            'text': "Estimated Associations",
            'y': 0.9,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top',
        },
        template="simple_white",
        height=400,
        width=600
    )

    return fig

def draw_double_volcano_graph(
    corncob_df,
    parameter, 
    comparison_parameter,
    neg_log_pvalue_min, 
    fdr_on_off, 
):

    # Set the metric to plot
    if fdr_on_off == "off":
        plot_y = "neg_log10_pvalue"
        hovertemplate = "CAG %{id}<br>" + comparison_parameter + " p-value (-log10): %{x}<br>" + parameter + " p-value (-log10): %{y}<extra></extra>"
        axis_suffix = "p-value (-log10)"
    else:
        plot_y = "neg_log10_qvalue"
        hovertemplate = "CAG %{id}<br>" + comparison_parameter + " q-value (-log10): %{x}<br>" + parameter + " q-value (-log10): %{y}<extra></extra>"
        axis_suffix = "q-value (-log10)"

    # Subset to these two parameters and pivot to be wide
    plot_df = pd.concat([
        corncob_df.query(
            "parameter == '{}'".format(param_name)
        ).query(
            "neg_log10_pvalue >= {}".format(neg_log_pvalue_min)
        )
        for param_name in [parameter, comparison_parameter]
    ]).pivot_table(
        index="CAG",
        columns="parameter",
        values=plot_y
    )
    
    # Set the minimum value after filtering
    plot_df.fillna(
        plot_df.apply(lambda c: c.dropna().min()).min(),
        inplace=True
    )

    fig = go.Figure(
        data=go.Scattergl(
            x=plot_df[comparison_parameter],
            y=plot_df[parameter],
            ids=plot_df.index.values,
            text=plot_df.index.values,
            hovertemplate=hovertemplate,
            mode="markers",
            opacity=0.5,
        ),
    )

    fig.update_layout(
        xaxis_title="{} {}".format(comparison_parameter, axis_suffix),
        yaxis_title="{} {}".format(parameter, axis_suffix),
        template="simple_white",
        height=400,
        width=600
    )

    return fig


####################
# ENRICHMENT GRAPH #
####################

def draw_enrichment_graph(
    enrichment_df, 
    annotation, 
    parameter,
):

    fig = go.Figure(
        data=go.Scatter(
            x = enrichment_df["estimate"],
            y = list(range(enrichment_df.shape[0])),
            error_x = dict(
                type='data',
                array=enrichment_df["std_error"],
                visible=True
            ),
            ids = enrichment_df.index.values,
            text = enrichment_df.apply(
                lambda r: "FDR-adjusted p-value: {:.2E}".format(r['q_value']),
                axis=1
            ),
            hovertemplate = "%{id}<br>Estimated Coefficient: %{x}<br>%{text}<extra></extra>",
            mode = "markers",
            marker_color = "LightSkyBlue",
        )
    )

    # Add a vertical dashed line at x=0 
    fig.add_shape(
        dict(
            type="line",
            x0=0,
            y0=0,
            x1=0,
            y1=enrichment_df.shape[0],
            line=dict(
                dash="dash",
                width=1,
            )
        )
    )
    # # Adjust the axis limits so that it is visible
    xmin = (enrichment_df["estimate"] - (enrichment_df["std_error"] / 2)).min()
    if xmin > 0:
        xmin = 0
    xmax = (enrichment_df["estimate"] + (enrichment_df["std_error"] / 2)).max()
    if xmax < 0:
        xmax = 0
    xpadding = (xmax - xmin) * 0.05
    fig.update_xaxes(
        range=[
            xmin - xpadding,
            xmax + xpadding,
        ]
    )

    fig.update_layout(
        xaxis_title="Estimated Coefficient of Association",
        template="simple_white",
        height=600,
        width=600,
        yaxis = dict(
            tickmode = "array",
            tickvals = list(range(enrichment_df.shape[0])),
            ticktext=[
                n[:30] + "..." if len(n) > 30 else n
                for n in enrichment_df.index.values
            ],
            automargin = True,
        )
    )

    return fig


##################
# TAXONOMY GRAPH #
##################
def draw_taxonomy_sunburst(
    cag_tax_df, 
    plot_title,
    min_ngenes,
    ranks_to_plot = [
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
):
    # If no assignments were made, just show an empty plot
    if cag_tax_df is None:
        fig = go.Figure(data=[])
        fig.update_layout(
            template="simple_white",
        )
        return fig

    # Set the index on the tax ID
    cag_tax_df = cag_tax_df.set_index("tax_id")

    # Make sure that the parent column is formatted as an integer
    cag_tax_df = cag_tax_df.apply(
        lambda c: c.fillna(0).apply(int) if c.name == "parent" else c
    )

    # Add the tax ID to any duplicated names (to make them unique)
    # This makes a Series with the name assigned to tax ID, 
    # but does not add to the table itself
    name_vc = cag_tax_df["name"].value_counts()
    name_dict = {
        tax_id: r["name"] if name_vc[r["name"]] == 1 else "{} (NCBI ID {})".format(
            r["name"], tax_id
        )
        for tax_id, r in cag_tax_df.iterrows()
    }

    # Remove any rows which are below the threshold, or are at the wrong rank
    taxa_to_remove = [
        tax_id
        for tax_id, r in cag_tax_df.iterrows()
        if r["count"] < min_ngenes or r["rank"] not in ranks_to_plot
    ]

    # Make sure that we have some data left to plot
    if len(taxa_to_remove) == cag_tax_df.shape[0]:
        fig = go.Figure(data=[])
        fig.update_layout(
            template="simple_white",
        )
        return fig

    # Walk through and remove the rows, reassigning the 'parent' for each
    for tax_id in taxa_to_remove:
        # Get the parent of this taxon
        parent_tax_id = cag_tax_df.loc[tax_id, "parent"]

        # For any taxa which have this taxon (to be removed) as the parent,
        # replace that value with the parent of this taxon 
        cag_tax_df = cag_tax_df.replace(
            to_replace={
                "parent": {
                    tax_id: parent_tax_id
                }
            }
        ).drop(
            index=tax_id
        )

    # Now we will fill in the name of the parent
    #  (which is currently encoded as a tax ID)
    # Crucially, set the 'parent_name' to None for any taxon
    #  which is its own parent
    cag_tax_df = cag_tax_df.assign(
        unique_name = pd.Series(name_dict),
        parent_name = cag_tax_df.apply(
            lambda r: name_dict[r["parent"]] if r["parent"] != r.name else None,
            axis=1
        ),
    )

    fig = go.Figure(
        data=go.Sunburst(
            labels=cag_tax_df["unique_name"],
            parents=cag_tax_df["parent_name"],
            values=cag_tax_df["count"],
            branchvalues="total",
        )
    )

    fig.update_layout(
        title={
            'text': plot_title,
            'y': 0.9,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top',
        },
    )

    return fig


####################
# SINGLE CAG GRAPH #
####################
def draw_single_cag_graph(
    plot_df,
    plot_title,
    axis_label,
    xaxis,
    plot_type,
    color,
    facet,
    log_scale
):

    # Make a list of the columns needed for plotting
    columns_for_plotting = [xaxis, "CAG_ABUND"]
    if color != "none":
        columns_for_plotting.append(color)
    if facet != "none":
        columns_for_plotting.append(facet)

    # Drop any samples which are missing the required data
    plot_df = plot_df.reindex(
        columns=columns_for_plotting
    ).dropna()

    # Protect against completely missing data with an empty plot
    empty_fig = go.Figure()
    empty_fig.update_layout(
        template="simple_white",
        yaxis_title=axis_label
    )
    if plot_df.shape[0] == 0 or (plot_df["CAG_ABUND"] > 0).sum() == 0:
        return empty_fig
    if xaxis == color and xaxis != "none":
        return empty_fig
    if xaxis == facet and xaxis != "none":
        return empty_fig
    if color == facet and color != "none":
        return empty_fig

    # For plotting on a log scale, replace zero values with the minimum
    if log_scale == "on":
        min_value = plot_df.query("CAG_ABUND > 0")["CAG_ABUND"].min() / 2
        plot_df.replace(
            to_replace = {"CAG_ABUND": 0},
            value = min_value,
            inplace=True
        )

    if plot_type == "scatter":
        plot_func = px.scatter

    elif plot_type == "boxplot":
        plot_func = px.box

    elif plot_type == "strip":
        plot_func = px.strip

    else:
        assert plot_type == "line"
        plot_func = px.line

    fig = plot_func(
        plot_df.sort_values(by=xaxis),
        x = xaxis,
        y = "CAG_ABUND",
        color = None if color == "none" else color,
        facet_col = None if facet == "none" else facet
    )

    # Apply the log transform
    if log_scale == "on":
        fig.update_yaxes(type="log")

    fig.update_layout(
        template="simple_white",
        yaxis_title=axis_label,
        title={
            'text': plot_title,
            'y': 1.0,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top',
        },
    )
    return fig

##################
# GENOME HEATMAP #
##################
def plot_genome_scatter(genome_df, parameter, genome_manifest):

    # Get the genome names to plot
    genome_names = genome_manifest.set_index(
        "id"
    )[
        "name"
    ].reindex(index=genome_df[
        "genome_id"
    ])

    # Draw a scatter plot
    fig = go.Figure(
        go.Scattergl(
            x=genome_df["prop_pass_fdr"].apply(lambda v: round(v, 2)),
            y=genome_df["mean_est_coef"].apply(lambda v: round(v, 2)),
            ids=genome_df["genome_id"],
            text=genome_names,
            marker_color="blue",
            hovertemplate="Genome %{text}<br>Accession: %{id}<br>Genome Proportion Passing FDR: %{x}<br>Mean Estimated Coefficient: %{y}<extra></extra>",
            mode="markers",
            opacity=0.5,
        ),
        layout=go.Layout(
            clickmode="event+select"
        )
    )

    # Set the style of the entire plot
    fig.update_layout(
        xaxis_title="Proportion of Genome Passing FDR",
        yaxis_title="Mean Estimated Coefficient",
        template="simple_white",
        showlegend=False,
        height=400,
        width=600,
    )

    return fig

def cluster_dataframe(plot_df):
    # Cluster on the basis of the proportion of counts from each column
    prop_df = plot_df / plot_df.sum().clip(lower=1)

    # Only cluster on each axis if there are >3 elements on the axis
    if plot_df.shape[0] > 3:

        # Get the labels on the sorted axis
        label_list = plot_df.index.values[
            leaves_list(linkage(
                prop_df,
                method="ward"
            ))
        ]

        # Apply the new order to the tables
        plot_df = plot_df.reindex(index=label_list)
        prop_df = prop_df.reindex(index=label_list)

    if plot_df.shape[1] > 3:
        # Get the labels on the sorted axis
        label_list = plot_df.columns.values[
            leaves_list(linkage(
                prop_df.T,
                method="ward"
            ))
        ]

        # Apply the new order to the tables
        plot_df = plot_df.reindex(columns=label_list)
        prop_df = prop_df.reindex(columns=label_list)

    return plot_df

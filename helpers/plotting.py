#!/usr/bin/env python3

import json
import numpy as np
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import leaves_list
from scipy.stats import zscore
from seaborn import color_palette
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

##################
# RICHNESS GRAPH #
##################
def update_richness_graph(
    richness_df,
    selected_metric,
    selected_type,
    log_x,
    manifest_json,
    full_manifest_df,
):
    # Make sure the selected parameters are in scope
    assert selected_type in ["hist", "scatter"]
    assert selected_metric in richness_df.columns.values, (selected_metric, richness_df.columns.values)

    # Get the filtered manifest from the browser
    plot_manifest_df = parse_manifest_json(manifest_json, full_manifest_df)

    # Subset the richness table based on the filtered manifest
    plot_richness_df = richness_df.reindex(
        index=plot_manifest_df.index.values
    )

    metric_names = {
        "prop_reads_aligned": "Prop. Reads Aligned",
        "n_genes_aligned": "Num. Genes Aligned",
        "n_genes_assembled": "Num. Genes Assembled",
    }

    if selected_type == "scatter":
        if "genes" in selected_metric:
            hovertemplate = "%{text}: %{x:,} reads - %{y:,} genes<extra></extra>"
        else:
            hovertemplate = "%{text}: %{x:,} reads - %{y:.2f} aligned<extra></extra>"

        fig = go.Figure(
            data=go.Scatter(
                x=plot_richness_df["n_reads"],
                y=plot_richness_df[selected_metric],
                text=plot_richness_df.index.values,
                hovertemplate=hovertemplate,
                mode="markers",
            ),
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

    # Make a plot with two panels, one on top of the other, sharing the x-axis
    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True
    )

    # The plot will depend on whether metadata has been selected
    if metadata == "none":

        # No metadata

        # Histogram on the top panel
        fig.add_trace(
            go.Histogram(
                x=plot_df[plot_df.columns.values[primary_pc - 1]],
                hovertemplate="Range: %{x}<br>Count: %{y}<extra></extra>",
            ),
            row=1, col=1
        )
        fig.update_yaxes(
            title_text="Number of Specimens",
            row=1, col=1
        )
        # Scatter on the bottom panel
        fig.add_trace(
            go.Scatter(
                x=plot_df[plot_df.columns.values[primary_pc - 1]],
                y=plot_df[plot_df.columns.values[secondary_pc - 1]],
                ids=plot_df.index.values,
                text=plot_df.index.values,
                hovertemplate="%{id}<extra></extra>",
                mode="markers",
                marker_color="blue"
            ),
            row=2, col=1
        )

    else:

        # Add the specified metadata to the DataFrame
        plot_df = plot_df.assign(
            METADATA = plot_manifest_df[metadata]
        ).rename(columns={
            "METADATA": metadata
        })

        # Make a numeric transform of the metadata
        if plot_df[metadata].apply(
            lambda n: isinstance(n, float) or isinstance(n, int)
        ).all():

            # The metadata is already numeric
            plot_df = plot_df.assign(
                METADATA_FLOAT = plot_df[metadata]
            )

        # Try to convert to datetime
        elif pd.to_datetime(plot_df[metadata], errors="coerce").isnull().sum() == 0:
            plot_df = plot_df.assign(
                METADATA_FLOAT = pd.to_datetime(
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
                METADATA_FLOAT = plot_df[metadata].apply(rank_order.get)
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
            METADATA_COLOR = plot_df["METADATA_FLOAT"].apply(cmap.get)
        )
        assert plot_df["METADATA_COLOR"].isnull().sum() == 0, (plot_df.head(), cmap)

        # Scatterplot or Boxplot on the top panel
        if plot_df[metadata].unique().shape[0] <= 5:

            # Iterate over each of the metadata groups
            for metadata_label, metadata_plot_df in plot_df.groupby(metadata):
                # Boxplot on the upper panel
                fig.add_trace(
                    go.Box(
                        x=metadata_plot_df[
                            plot_df.columns.values[primary_pc - 1]
                        ],
                        name=metadata_label,
                        marker_color=metadata_plot_df["METADATA_COLOR"].values[0],
                    ),
                    row=1, col=1
                )
                # Scatter on the bottom panel
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
                        text=metadata_plot_df[metadata],
                        hovertemplate="%{id}<br>%{text}<extra></extra>",
                        mode="markers",
                        marker_color=metadata_plot_df["METADATA_COLOR"].values[0],
                    ),
                    row=2, col=1
                )

        else:
            # Scatter on the upper panel
            fig.add_trace(
                go.Scatter(
                    x=plot_df[
                        plot_df.columns.values[primary_pc - 1]
                    ],
                    y=plot_df[metadata],
                    ids=plot_df.index.values,
                    text=plot_df[metadata],
                    hovertemplate="%{id}<br>%{text}<extra></extra>",
                    mode="markers",
                    marker_color=plot_df["METADATA_COLOR"],
                ),
                row=1, col=1
            )

            # Scatter on the bottom panel
            fig.add_trace(
                go.Scatter(
                    x=plot_df[plot_df.columns.values[primary_pc - 1]],
                    y=plot_df[plot_df.columns.values[secondary_pc - 1]],
                    ids=plot_df.index.values,
                    text=plot_df[metadata],
                    hovertemplate="%{id}<br>%{text}<extra></extra>",
                    mode="markers",
                    marker_color=plot_df["METADATA_COLOR"],
                ),
                row=2, col=1
            )

        fig.update_yaxes(
            title_text=metadata,
            row=1, col=1
        )


    fig.update_xaxes(
        title_text=plot_df.columns.values[primary_pc - 1],
        row=1, col=1
    )
    fig.update_xaxes(
        title_text=plot_df.columns.values[primary_pc - 1],
        row=2, col=1
    )
    fig.update_yaxes(
        title_text=plot_df.columns.values[secondary_pc - 1],
        row=2, col=1
    )

    fig.update_layout(
        showlegend=False,
        template="simple_white",
        height=800,
        width=600,
    )

    return fig


####################
# CAG SUMMARY CARD #
####################
def draw_cag_summary_graph_hist(
    cag_summary_df,
    metric_primary,
    size_range,
    entropy_range,
    prevalence_range,
    abundance_range,
    nbinsx,
    log_scale,
):
    # Apply the filters
    plot_df = cag_summary_df.query(
        "size >= {}".format(10**size_range[0])
    ).query(
        "size <= {}".format(10**size_range[1])
    ).query(
        "prevalence >= {}".format(prevalence_range[0])
    ).query(
        "prevalence <= {}".format(prevalence_range[1])
    ).query(
        "mean_abundance >= {}".format(abundance_range[0])
    ).query(
        "mean_abundance <= {}".format(abundance_range[1])
    ).query(
        "entropy >= {}".format(entropy_range[0])
    ).query(
        "entropy <= {}".format(entropy_range[1])
    ).apply(
        lambda c: c.apply(np.log10) if c.name == "size" else c
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
    fig = go.Figure(
        go.Histogram(
            x=plot_df[metric_primary],
            histfunc="count",
            nbinsx=nbinsx,
            hovertemplate="Range: %{x}<br>Count: %{y}<extra></extra>",
        ),
    )
    fig.update_layout(
        xaxis_title=axis_names[metric_primary],
        yaxis_title="Number of CAGs",
        template="simple_white",
        height=400,
        width=600,
    )

    # Apply the log transform
    if log_scale == "on":
        fig.update_yaxes(type="log")

    return fig


def draw_cag_summary_graph_scatter(
    cag_summary_df,
    metric_primary,
    metric_secondary,
    size_range,
    entropy_range,
    prevalence_range,
    abundance_range,
    selected_cag_json,
):
    # Apply the filters
    plot_df = cag_summary_df.query(
        "size >= {}".format(10**size_range[0])
    ).query(
        "size <= {}".format(10**size_range[1])
    ).query(
        "prevalence >= {}".format(prevalence_range[0])
    ).query(
        "prevalence <= {}".format(prevalence_range[1])
    ).query(
        "mean_abundance >= {}".format(abundance_range[0])
    ).query(
        "mean_abundance <= {}".format(abundance_range[1])
    ).query(
        "entropy >= {}".format(entropy_range[0])
    ).query(
        "entropy <= {}".format(entropy_range[1])
    ).apply(
        lambda c: c.apply(np.log10) if c.name == "size" else c
    )

    axis_names = {
        "CAG": "CAG ID",
        "size": "Number of Genes (log10)",
        "mean_abundance": "Mean Abundance",
        "std_abundance": "Std. Abundance",
        "prevalence": "Prevalence",
        "entropy": "Entropy",
    }

    # Parse the selected CAG data
    if selected_cag_json is not None:
        # Set the points which are selected in the scatter plot
        selectedpoints = np.where(
            plot_df.index.values == json.loads(selected_cag_json)["id"]
        )
    else:
        selectedpoints = None

    # Draw a scatter plot
    fig = go.Figure(
        go.Scattergl(
            x=plot_df[metric_primary],
            y=plot_df[metric_secondary],
            ids=plot_df.index.values,
            text=plot_df.index.values,
            marker_color="blue",
            hovertemplate="CAG %{id}<br>X-value: %{x}<br>Y-value: %{y}<extra></extra>",
            mode="markers",
            opacity=0.5,
            selectedpoints=selectedpoints
        )
    )

    # Set the style of the entire plot
    fig.update_layout(
        xaxis_title=axis_names[metric_primary],
        yaxis_title=axis_names[metric_secondary],
        template="simple_white",
        showlegend=False,
        height=400,
        width=600,
    )

    return fig


###############
# CAG HEATMAP #
###############
def draw_cag_heatmap(
    cag_abund_df,
    metadata_selected,
    abundance_metric,
    cluster_by,
    taxa_rank,
    manifest_json,
    full_manifest_df,
    cag_tax_dict,
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
    cag_abund_df = cag_abund_df.reindex(
        index=plot_manifest_df.index
    )

    if abundance_metric in ["log10", "zscore"]:
        # Transform to log10 relative abundance

        # First find the lowest non-zero value
        lowest_value = cag_abund_df.apply(
            lambda c: c[c > 0].min()
        ).min()

        # Now transform and set the floor
        cag_abund_df = cag_abund_df.clip(
            lower=lowest_value
        ).applymap(
            np.log10
        ).applymap(
            lambda i: round(i, 1)
        )
    
    if abundance_metric == "zscore":
        # Transform into the Z-score per-sample
        cag_abund_df = (cag_abund_df - cag_abund_df.mean()) / cag_abund_df.std()

    # If selected, cluster the specimens by CAG abundance
    if cluster_by == "cag":
        cag_abund_df = cag_abund_df.reindex(
            index=cag_abund_df.index.values[
                leaves_list(
                    linkage(
                        cag_abund_df,
                        method="ward"
                    )
                )
            ]
        )
        plot_manifest_df = plot_manifest_df.reindex(
            index=cag_abund_df.index
        )

    # Make the specimens columns and CAGs rows
    cag_abund_df = cag_abund_df.T
    plot_manifest_df = plot_manifest_df.T

    # Set the figure width
    figure_width = 800
    # Set the figure height
    figure_height = 800

    # Depending on whether metadata or taxonomic information has
    # been provided, the plot will be set up in different ways
    # Metadata - & taxonomic annotations - : single plot
    # Metadata + & taxonomic annotations - : two plots, metadata on top of abund
    # Metadata - & taxonomic annotations + : two plots, tax to the right of abund
    # Metadata + & taxonomic annotations + : three plots, combining tax and metadata

    has_metadata = len(metadata_selected) > 0
    has_taxonomy = (taxa_rank != "none") & any(df is not None for _, df in cag_tax_dict.items())

    if has_metadata is False and has_taxonomy is False:

        fig = go.Figure(
            data=draw_cag_abund_heatmap_panel(cag_abund_df),
        )

    elif has_metadata is False and has_taxonomy:

        fig = draw_cag_abund_heatmap_with_tax(
            cag_abund_df, cag_tax_dict, taxa_rank
        )

    elif has_metadata and has_taxonomy is False:

        fig = draw_cag_abund_heatmap_with_metadata(
            cag_abund_df, plot_manifest_df, metadata_selected
        )

    else:

        fig = draw_cag_abund_heatmap_with_metadata_and_tax(
            cag_abund_df, 
            plot_manifest_df, 
            metadata_selected, 
            cag_tax_dict, 
            taxa_rank,
        )

    fig.update_layout(
        width=figure_width,
        height=figure_height,
    )
    return fig


def draw_cag_abund_heatmap_with_tax(
    cag_abund_df, 
    cag_tax_dict,
    taxa_rank
):

    # Make a plot with two panels, side-by-side, sharing the y-axis

    # The taxa plot will be very small
    fig = make_subplots(
        rows=1, cols=2, shared_yaxes=True,
        column_widths=[
            0.95, 0.05
        ],
        horizontal_spacing=0.005,
    )

    # Plot the abundances on the left
    fig.add_trace(
        draw_cag_abund_heatmap_panel(cag_abund_df), row=1, col=1
    )

    # Plot the taxonomic annotations on the right
    fig.add_trace(
        draw_cag_abund_taxon_panel(cag_tax_dict, taxa_rank), row=1, col=2
    )

    return fig

def draw_cag_abund_heatmap_with_metadata(
    cag_abund_df, 
    plot_manifest_df, 
    metadata_selected
):

    # Make a plot with two panels, one on top of the other, sharing the x-axis

    # The relative height of the subplots will be set dynamically
    metadata_height = max(
        0.1,
        min(
            0.5,
            len(metadata_selected) / 20.
        ),
    )
    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True,
        row_heights=[
            metadata_height,
            1 - metadata_height
        ],
        vertical_spacing=0.01,
    )

    # Plot the metadata on the top
    fig.add_trace(
        draw_metadata_heatmap_panel(plot_manifest_df), row=1, col=1
    )

    # Plot the abundances on the bottom
    fig.add_trace(
        draw_cag_abund_heatmap_panel(cag_abund_df), row=2, col=1
    )

    return fig

def draw_cag_abund_heatmap_with_metadata_and_tax(
    cag_abund_df, 
    plot_manifest_df, 
    metadata_selected, 
    cag_tax_dict, 
    taxa_rank,
):

    # Make a plot with four panels:
    # metadata - blank
    # cag-abun - taxon

    # The relative height of the subplots will be set dynamically
    metadata_height = max(
        0.1,
        min(
            0.5,
            len(metadata_selected) / 20.
        ),
    )

    fig = make_subplots(
        rows=2, 
        cols=2, 
        shared_xaxes=True,
        shared_yaxes=True,
        row_heights=[
            metadata_height,
            1 - metadata_height
        ],
        vertical_spacing=0.01,
        column_widths=[
            0.95, 0.05
        ],
        horizontal_spacing=0.005,
    )

    # Plot the abundances on the bottom-left
    fig.add_trace(
        draw_cag_abund_heatmap_panel(cag_abund_df), row=2, col=1
    )

    # Plot the taxonomic annotations on the bottom-right
    fig.add_trace(
        draw_cag_abund_taxon_panel(cag_tax_dict, taxa_rank), row=2, col=2
    )

    # Plot the metadata on the top-left
    fig.add_trace(
        draw_metadata_heatmap_panel(plot_manifest_df), row=1, col=1
    )

    return fig

def draw_cag_abund_taxon_panel(
    cag_tax_dict, 
    taxa_rank
):

    # For each CAG, pick out the top hit
    summary_df = pd.DataFrame([
        summarize_cag_taxa(cag_id, cag_tax_df, taxa_rank)
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

    return go.Heatmap(
        x=[taxa_rank],
        y=["CAG {} -".format(i) for i in summary_df["CAG"].values],
        z=summary_df.reindex(columns=["name_scalar"]).values,
        text=summary_df.reindex(columns=["label"]).values,
        showscale=False,
        colorscale='Viridis',
        hovertemplate="%{y}<br>%{text}<extra></extra>"
    )


def summarize_cag_taxa(cag_id, cag_tax_df, taxa_rank):
    """Helper function to summarize the top hit at a given rank."""

    # If there are no hits at this level, return None
    if cag_tax_df is None or ((cag_tax_df["rank"] == taxa_rank).sum() == 0):
        return {
            "CAG": cag_id,
            "name": 'none',
            "label": "No genes assigned at this level"
        }

    # Filter down to this rank
    df = cag_tax_df.query("rank == '{}'".format(taxa_rank))

    # Sort by 'consistent' hits
    df.sort_values(by="consistent", inplace=True, ascending=False)

    # Return the top hit
    return {
        "CAG": cag_id,
        "name": df["name"].values[0],
        "label": "{}<br>{:,} genes assigned at the {} level or above<br>{:,} genes consistent with {}".format(
            df["name"].values[0],
            int(df["count"].values[0]),
            taxa_rank,
            int(df["consistent"].values[0]),
            df["name"].values[0],
        )
    }


def draw_metadata_heatmap_panel(
    plot_manifest_df,
    hovertemplate="Specimen: %{x}<br>Label: %{y}<br>Value: %{text}<extra></extra>",
):
    return go.Heatmap(
        z=plot_manifest_df.apply(
            lambda r: r.apply(dict(zip(r.unique(), np.arange(0, 1, 1 / r.unique().shape[0]))).get),
            axis=1
        ).values,
        text=plot_manifest_df.values,
        y=["{}".format(i) for i in plot_manifest_df.index.values],
        x=["- {}".format(i) for i in plot_manifest_df.columns.values],
        colorscale='Viridis',
        showscale=False,
        hovertemplate=hovertemplate,
    )

def draw_cag_abund_heatmap_panel(
    cag_abund_df,
    hovertemplate = "Specimen: %{x}<br>CAG: %{y}<br>Rel. Abund. (log10): %{z}<extra></extra>",
):
    return go.Heatmap(
        z=cag_abund_df.values,
        y=["CAG {} -".format(i) for i in cag_abund_df.index.values],
        x=["- {}".format(i) for i in cag_abund_df.columns.values],
        colorbar={"title": "Abundance (log10)"},
        colorscale='Bluered_r',
        hovertemplate=hovertemplate,
    )

#################
# VOLCANO GRAPH #
#################
def draw_volcano_graph(
    corncob_df,
    parameter, 
    neg_log_pvalue_min, 
    fdr_on_off, 
    selected_cag_json
):
    if corncob_df is None or neg_log_pvalue_min is None:
        return go.Figure()

    # Subset to just this parameter
    plot_df = corncob_df.query(
        "parameter == '{}'".format(parameter)
    ).query(
        "neg_log_pvalue >= {}".format(neg_log_pvalue_min)
    )

    # Parse the selected CAG data
    if selected_cag_json is not None:
        # Set the points which are selected in the scatter plot
        selectedpoints = np.where(
            plot_df["CAG"].values == json.loads(selected_cag_json)["id"]
        )
    else:
        selectedpoints = None

    if fdr_on_off == "off":
        plot_y = "neg_log_pvalue"
        hovertemplate = "CAG %{id}<br>Estimate: %{x}<br>p-value (-log10): %{y}<extra></extra>"
        yaxis_title = "p-value (-log10)"
    else:
        plot_y = "neg_log_qvalue"
        hovertemplate = "CAG %{id}<br>Estimate: %{x}<br>q-value (-log10): %{y}<extra></extra>"
        yaxis_title = "q-value (-log10)"

    fig = go.Figure(
        data=go.Scattergl(
            x=plot_df["estimate"],
            y=plot_df[plot_y],
            ids=plot_df["CAG"].values,
            text=plot_df["CAG"].values,
            hovertemplate=hovertemplate,
            mode="markers",
            opacity=0.5,
            selectedpoints=selectedpoints,
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


##################
# TAXONOMY GRAPH #
##################
def draw_taxonomy_sunburst(cag_tax_df, cag_id, min_ngenes):
    # If no assignments were made, just show an empty plot
    if cag_tax_df is None:
        fig = go.Figure(data=[])
        fig.update_layout(
            template="simple_white",
        )
        return fig, 1, {}

    # Filter by the number of genes
    if (cag_tax_df["count"] >= min_ngenes).any():
        cag_tax_df = cag_tax_df.query("count >= {}".format(min_ngenes))

    fig = go.Figure(
        data=go.Sunburst(
            labels=cag_tax_df["name"],
            parents=cag_tax_df["parent"],
            values=cag_tax_df["count"],
            branchvalues="total",
        )
    )

    fig.update_layout(
        title={
            'text': "CAG {}".format(cag_id),
            'y': 0.9,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top',
        },
    )

    # Format the marks for the updated slider
    marks = {
        str(int(n)): str(int(n))
        for n in [
            1,
            int(cag_tax_df["count"].max() / 2),
            cag_tax_df["count"].max(),
        ]
    }

    return fig, cag_tax_df["count"].max(), marks


####################
# SINGLE CAG GRAPH #
####################
def draw_single_cag_graph(
    plot_df,
    cag_id,
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
        yaxis_title="CAG {}".format(cag_id)
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
        yaxis_title="CAG {}".format(cag_id)
    )
    return fig

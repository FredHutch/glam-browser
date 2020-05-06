#!/usr/bin/env python3

import click
from collections import defaultdict
import dash
import dash_table
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import json
import numpy as np
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from time import time

# Read in data from the HDF5 defined in `HDF5_fp`
hdf5_fp = os.getenv("HDF5_FP")
assert os.path.exists(hdf5_fp), "Path does not exist: {}".format(hdf5_fp)
with pd.HDFStore(hdf5_fp, "r") as store:

    # Table with summary of the entire experiment
    exp_metrics = pd.read_hdf(
        store, "/summary/experiment"
    ).set_index(
        "variable"
    )[
        "value"
    ]

    # Table for RICHNESS GRAPH
    richness_df = pd.read_hdf(store, "/summary/all")

    # Table for CAG SIZE GRAPH
    cag_summary_df = pd.read_hdf(store, "/annot/cag/all")

    # Specimen metadata table (manifest provided by user)
    manifest_df = pd.read_hdf(store, "/manifest")

    # Pairwise distance tables
    distance_df_dict = {
        metric: pd.read_hdf(store, "/distances/%s" % metric).set_index("specimen")
        for metric in ["euclidean", "aitchison", "braycurtis", "jaccard"]
    }

    # Corncob results (for volcano plot)
    corncob_df = pd.read_hdf(store, "/stats/cag/corncob")

    # Taxonomy table
    taxonomy_df = pd.read_hdf(store, "/ref/taxonomy")

    # Abundance of all genes
    gene_abundance_df = pd.read_hdf(store, "/abund/gene/wide")
    
    # Abundance of all CAGs
    cag_abundance_df = pd.read_hdf(store, "/abund/cag/wide")

# Precompute some useful metrics

# Precompute the proportion of reads which align
richness_df = richness_df.assign(
    prop_reads_aligned=richness_df["aligned_reads"] / richness_df["n_reads"]
)

# Round to 4 significant figures on the CAG summary metrics
cag_summary_limits = cag_summary_df.describe()
cag_summary_df = cag_summary_df.apply(
    lambda c: c.apply(
        lambda i: '{:g}'.format(float('{:.4g}'.format(i)))
    ) if c.name in [
        "std_abundance", "prevalence", "mean_abundance"
    ] else c
)

# The limits of CAG sizes
cag_size_log = cag_summary_df["size"].apply(np.log10)
cag_size_min = cag_size_log.min()
cag_size_max = cag_size_log.max()
cag_size_range = cag_size_max - cag_size_min

# Set index on manifest
manifest_df.set_index("specimen", inplace=True)

# Metadata fields found in the manifest
metadata_fields = [n for n in manifest_df.columns.values if n not in ["R1", "R2", "I1", "I2"]]

# Calculate the -log10(p_value)
corncob_df = corncob_df.assign(
    neg_log_pvalue = corncob_df["p_value"].apply(np.log10) * -1,
    neg_log_qvalue = corncob_df["q_value"].apply(np.log10) * -1,
)
max_neg_log_pvalue = corncob_df.groupby("parameter")["neg_log_pvalue"].max()

# Format the parent tax ID as an integer
taxonomy_df = taxonomy_df.apply(
    lambda c: c.fillna(0).apply(float).apply(int) if c.name in ["parent", "tax_id"] else c,
).set_index("tax_id")

# external CSS stylesheets
external_stylesheets = [
    {
        'href': 'https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/css/bootstrap.min.css',
        'rel': 'stylesheet',
        'integrity': 'sha384-Vkoo8x4CGsO3+Hhxv8T/Q5PaXtkKtu6ug5TOeNV6gBiFeWPGFN9MuhOf23Q9Ifjh',
        'crossorigin': 'anonymous'
    }
]

##################
# HELPER METHODS #
##################
def path_to_root(tax_id, max_steps=100):
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


def make_cag_tax_df(tax_id_list, ranks_to_keep=["phylum", "class", "order", "family", "genus", "species"]):
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
        for anc_tax_id in path_to_root(tax_id):

            # Add to the sum for every node we visit along the way
            counts[anc_tax_id] += n_genes

    if len(counts) == 0:
        return

    # Make a DataFrame
    df = pd.DataFrame({"count": counts})

    # Add the name, parent, rank
    df = df.assign(
        tax_id = df.index.values,
        parent_tax_id = taxonomy_df["parent"],
        rank = taxonomy_df["rank"],
    )

    # Set the parent of the root as ""
    df.loc[0, "parent"] = ""

    # Remove any taxa which aren't at the right rank (but keep the root)
    df = df.assign(
        to_remove = df.apply(
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
        name = df["tax_id"].apply(lambda t: tax_names.get(t, "")),
        parent = df["parent_tax_id"].apply(lambda t: tax_names.get(t, "")),
    )

    return df.reindex(columns=["name", "parent", "count"])

###############################
# REUSABLE DISPLAY COMPONENTS #
###############################
def card_wrapper(card_title, card_body):
    """Nest a set of display elements within a card, including body and title"""
    return html.Div(
        [
            html.Div(
                [
                    card_title
                ],
                className="card-header"
            ),
            html.Div(
                card_body,
                className="row"
            )

        ],
        className="card"
    )


def exp_table_row(header1, value1, header2, value2, header_bg="#F4F6F6", value_bg="white", spacer_bg="white"):
    return [     # Table Body
        html.Tr([    # Row
            html.Td(
                header1,
                style={"background-color": header_bg}
            ),
            html.Td(
                value1,
                style={"background-color": value_bg}
            ),
            html.Td(
                "",
                style={"background-color": spacer_bg}
            ),
            html.Td(
                header2,
                style={"background-color": header_bg}
            ),
            html.Td(
                value2,
                style={"background-color": value_bg}
            ),
        ]
    )]


def graph_div(anchor_id, graph_id, className="col-sm-8"):
    """Return a div containing a dcc.Graph and anchor, all within a col-sm-8."""
    return html.Div(
        [
            html.A(id=anchor_id),
            dcc.Graph(
                id=graph_id
            )
        ],
        className=className,
    )


def metadata_field_dropdown(
    dropdown_id, 
    label_text='Metadata Label', 
    default_value="none"
):
    return [
        html.Label(label_text),
        dcc.Dropdown(
            id=dropdown_id,
            options=[
                {'label': 'None', 'value': 'none'},
            ] + [
                {'label': f, "value": f}
                for f in metadata_fields
            ],
            value=default_value
        ),
        html.Br(),

    ]


def plot_type_dropdown(
    dropdown_id, 
    label_text='Plot Type',
    options=[
        {'label': 'Points', 'value': 'scatter'},
        {'label': 'Histogram', 'value': 'hist'},
    ],
    default_value="scatter"
):
    return [
        html.Label(label_text),
        dcc.Dropdown(
            id=dropdown_id,
            options=options,
            value=default_value
        ),
        html.Br(),
    ]

def ordination_pc_slider(
    slider_id,
    label_text,
    default_value
):
    return basic_slider(
        slider_id,
        label_text,
        min_value=1,
        max_value=min(manifest_df.shape[0], 10),
        default_value=default_value,
        marks = range(
            1,
            min(
                manifest_df.shape[0], 
                13
            ),
            2
        ),
        included=False
    )

def basic_slider(
    slider_id,
    label_text,
    min_value=1,
    max_value=10,
    step_value=1,
    default_value=1,
    marks=[],
    included=True
):
    return [
        html.Div([
            html.Label(label_text),
            dcc.Slider(
                id=slider_id,
                min=min_value,
                max=max_value,
                step=step_value,
                marks={
                    str(n): str(n)
                    for n in marks
                },
                value=default_value,
                included=included
            ),
            html.Br()
        ], id="%s-div" % slider_id)
    ]

def volcano_parameter_dropdown(
    dropdown_id, 
    label_text='Parameter',
):
    if corncob_df is None:

        return [
            html.Label(label_text),
            dcc.Dropdown(
                id=dropdown_id,
                options=[
                    {'label': 'None', 'value': 'none'},
                ],
                value="none"
            ),
            html.Br(),
        ]

    else:

        parameter_list = corncob_df[
            "parameter"
        ].drop_duplicates(
        ).sort_values(
        ).tolist(
        )

        if len(parameter_list) > 1:
            default_value = parameter_list[1]
        else:
            default_value = parameter_list[0]

        return [
            html.Label(label_text),
            dcc.Dropdown(
                id=dropdown_id,
                options=[
                    {'label': l, 'value': l}
                    for l in parameter_list
                ],
                value=default_value
            ),
            html.Br(),
        ]


def volcano_pvalue_slider(slider_id, label_text='P-Value Filter'):
    """This slider is missing the max and marks, which will be updated by a callback."""
    return [
        html.Label(label_text),
        dcc.Slider(
            id=slider_id,
            min=0,
            step=0.1,
            value=1,
        ),
        html.Br()
    ]
    
def cag_size_slider(slider_id, label_text='CAG Size Filter'):
    return [
        html.Label(label_text),
        dcc.RangeSlider(
            id=slider_id,
            min=0,
            max=cag_size_max,
            step=0.1,
            marks={
                str(n): str(10**n)
                for n in range(int(cag_size_max))
            },
            value=[
                max(cag_size_min, np.log10(5)),
                cag_size_max
            ]
        ),
        html.Br()
    ]
    

def cag_metric_slider(slider_id, metric, label_text):
    min_val = cag_summary_limits.loc["min", metric]
    max_val = cag_summary_limits.loc["max", metric]
    range_val = max_val - min_val
    step_val = range_val / 5
    return [
        html.Label(label_text),
        dcc.RangeSlider(
            id=slider_id,
            min=min_val,
            max=max_val,
            step=step_val / 10.,
            marks={
                str(n): str(round(n, 2))
                for n in np.arange(
                    min_val,
                    max_val,
                    step_val
                )
            },
            value=[
                min_val,
                max_val
            ]
        ),
        html.Br()
    ]
    

def cag_metric_dropdown(slider_id, label_text='Metric', default_value="size"):
    return [
        html.Label(label_text),
        dcc.Dropdown(
            id=slider_id,
            options=[
                {'label': 'CAG ID', 'value': 'CAG'},
                {'label': 'Size', 'value': 'size'},
                {'label': 'Entropy', 'value': 'entropy'},
                {'label': 'Mean Abund.', 'value': 'mean_abundance'},
                {'label': 'Prevalence', 'value': 'prevalence'},
                {'label': 'Std. Abund', 'value': 'std_abundance'},
            ],
            value=default_value
        ),
        html.Br()
    ]

def log_scale_radio_button(id_string, default="off", label_text="Log Scale"):
    return [
        html.Br(),
        html.Label(label_text),
        dcc.RadioItems(
            id=id_string,
            options=[
                {'label': 'On ', 'value': 'on'},
                {'label': 'Off', 'value': 'off'},
            ],
            value=default,
        ),
        html.Br()
    ]

def nbins_slider(id_string):
    return [
        html.Label('Number of Bins'),
        dcc.Slider(
            id=id_string,
            min=5,
            max=100,
            step=1,
            marks={
                "5": "5",
                "20": "20",
                "100": "100",
            },
            value=20
        ),
        html.Br()
    ]

#################################
# \ REUSABLE DISPLAY COMPONENTS #
#################################


# Set up the app
app = dash.Dash(
    __name__, 
    external_stylesheets=external_stylesheets
)
app.title = "GLAM Browser"
app.config.suppress_callback_exceptions = True

app.layout = html.Div(
    children=[
        dbc.NavbarSimple(
            brand="GLAM Browser",
            light=True,
            sticky="top",
            children=[
                dbc.NavItem(
                    dbc.NavLink(
                        os.getenv("HDF5_FP").split("/")[-1],
                        href="#"
                    )
                ),
                dbc.DropdownMenu(
                    children = [
                        dbc.DropdownMenuItem("Genes Detected", href="#richness"),
                        dbc.DropdownMenuItem("CAG Size", href="#cag-size"),
                        dbc.DropdownMenuItem("PCA / t-SNE", href="#ordination"),
                    ],
                    nav=True,
                    in_navbar=True,
                    label="Contents",
                )
            ],
        ),
        ########
        # BODY #
        ########
        ############################
        # EXPERIMENT SUMMARY TABLE #
        ############################
        card_wrapper(
            "Experiment",
            [
                html.Div(
                    [
                        dbc.Table([             # Table
                            html.Tbody(
                                exp_table_row(  # Wrapper for each row
                                    "Total Reads",
                                    "{:,}".format(exp_metrics["total_reads"]),
                                    "Aligned Reads",
                                    "{}%".format(round(
                                        100 * exp_metrics["aligned_reads"] / exp_metrics["total_reads"],
                                        1
                                    ))
                                ) + exp_table_row(
                                    "Genes (#)",
                                    "{:,}".format(exp_metrics["num_genes"]),
                                    "CAGs (#)",
                                    "{:,}".format(exp_metrics["num_cags"])
                                ) + exp_table_row(
                                    "Specimens (#)",
                                    "{:,}".format(exp_metrics["num_samples"]),
                                    "Formula",
                                    "{}".format(exp_metrics["formula"])
                                )
                            )
                        ], bordered=True, hover=True, responsive=True)
                    ],
                    className="col-sm-12 my-auto",
                )
            ]
        ),
        ##############################
        # / EXPERIMENT SUMMARY TABLE #
        ##############################        
        ##################
        # RICHNESS GRAPH #
        ##################
        card_wrapper(
            "Gene Analysis Summary",
            [
                graph_div("richness", 'richness-graph'),
                html.Div(
                    [
                        html.Br(),
                        html.Label('Display Values'),
                        dcc.Dropdown(
                            id="richness-metric-dropdown",
                            options=[
                                {'label': 'Genes Assembled (#)', 'value': 'n_genes_assembled'},
                                {'label': 'Genes Aligned (#)', 'value': 'n_genes_aligned'},
                                {'label': 'Reads Aligned (%)', 'value': 'prop_reads_aligned'},
                            ],
                            value='prop_reads_aligned'
                        ),
                        html.Br(),
                    ] + plot_type_dropdown(
                        "richness-type-dropdown"
                    ),
                    className="col-sm-4 my-auto",
                )
            ]
        ),
        ####################
        # / RICHNESS GRAPH #
        ####################        
        ####################
        # ORDINATION GRAPH #
        ####################
        card_wrapper(
            "Ordination Analysis",
            [
                graph_div("ordination", 'ordination-graph'),
                html.Div(
                    [
                        html.Label('Distance Metric'),
                        dcc.Dropdown(
                            id="ordination-metric",
                            options=[
                                {'label': 'Euclidean', 'value': 'euclidean'},
                                {'label': 'Aitchison', 'value': 'aitchison'},
                                {'label': 'Braycurtis', 'value': 'braycurtis'},
                            ],
                            value='euclidean'
                        ),
                        html.Br(),
                        html.Label('Ordination Method'),
                        dcc.Dropdown(
                            id="ordination-algorithm",
                            options=[
                                {'label': 'PCA', 'value': 'pca'},
                                {'label': 't-SNE', 'value': 'tsne'},
                            ],
                            value='pca'
                        ),
                        html.Br(),
                    ] + ordination_pc_slider(
                        "ordination-primary-pc",
                        'Primary Axis',
                        1
                    ) + ordination_pc_slider(
                        "ordination-secondary-pc",
                        'Secondary Axis',
                        2
                    ) + basic_slider(
                        "ordination-tsne-perplexity",
                        "Perplexity",
                        max_value=100,
                        default_value=30,
                        marks=[0, 10, 20, 30, 50, 70, 100],
                        included=False,
                    ) + metadata_field_dropdown(
                        "ordination-metadata"
                    ),
                    className="col-sm-4 my-auto",
                )
            ],
        ),
        ######################
        # / ORDINATION GRAPH #
        ######################
        #####################
        # CAG SUMMARY GRAPH #
        #####################
        card_wrapper(
            "CAG Summary",
            [
                html.Div(
                    [
                        html.A(id="cag-summary"),
                        dcc.Graph(
                            id='cag-summary-graph-hist'
                        ),
                        dcc.Graph(
                            id='cag-summary-graph-scatter'
                        ),
                    ],
                    className="col-sm-8 my-auto",
                ),
                html.Div(
                    cag_metric_dropdown(
                        "cag-summary-metric-primary",
                        default_value="size",
                        label_text="Primary Metric (x-axis)",
                    ) + cag_metric_dropdown(
                        "cag-summary-metric-secondary",
                        default_value="entropy",
                        label_text="Secondary Metric (y-axis)",
                    ) + cag_size_slider(
                        "cag-summary-size-slider"
                    ) + cag_metric_slider(
                        "cag-summary-entropy-slider",
                        "entropy",
                        "CAG Entropy Filter",
                    ) + cag_metric_slider(
                        "cag-summary-prevalence-slider",
                        "prevalence",
                        "CAG Prevalence Filter"
                    ) + cag_metric_slider(
                        "cag-summary-abundance-slider",
                        "mean_abundance",
                        "CAG Abundance Filter",
                    ) + nbins_slider(
                        "cag-summary-nbinsx-slider"
                    ) + log_scale_radio_button(
                        "cag-summary-log",
                        default="on",
                        label_text="Histogram Log Scale"
                    ),
                    className="col-sm-4 my-auto",
                ),
                html.Div(id='global-selected-cag', style={"display": "none"}),
                html.Div(id='cag-summary-selected-cag', style={"display": "none"}),
            ]
        ),
        #######################
        # / CAG SUMMARY GRAPH #
        #######################
        ################
        # VOLCANO PLOT #
        ################
        card_wrapper(
            "Association Screening",
            [
                graph_div("volcano", 'volcano-graph'), 
                html.Div(
                    volcano_parameter_dropdown(
                        "volcano-parameter-dropdown",
                    ) + volcano_pvalue_slider(
                        "volcano-pvalue-slider",
                    ) + log_scale_radio_button(
                        "volcano-fdr-radio",
                        label_text="FDR-BH adjustment"
                    ),
                    className="col-sm-4 my-auto",
                ),
                html.Div(id='volcano-selected-cag', style={"display": "none"}),
            ]
        ),
        ##################
        # / VOLCANO PLOT #
        ##################
        ###################
        # CAG DETAIL CARD #
        ###################
        html.Div(
            [
                html.Div(
                    [
                        "CAG Details"
                    ],
                    className="card-header"
                ),
                html.Div(
                    [
                        graph_div("cag-detail-tax", 'cag-detail-tax-graph'),
                        html.Div(
                            basic_slider(
                                "cag-detail-tax-ngenes",
                                "Minimum Number of Genes",
                                included=False,
                                default_value=5,
                            ),
                            className="col-sm-4 my-auto"
                        )
                    ],
                    className="row"
                ),
                # html.Div(
                #     [
                #         # graph_div("cag-detail-heatmap", 'cag-detail-heatmap-graph'),
                #         html.Div(
                #             metadata_field_dropdown(
                #                 "cag-detail-heatmap-color-primary",
                #                 label_text="Color By (primary)",
                #             ) + metadata_field_dropdown(
                #                 "cag-detail-heatmap-color-secondary",
                #                 label_text="Color By (secondary)",
                #             ) + metadata_field_dropdown(
                #                 "cag-detail-heatmap-color-tertiary",
                #                 label_text="Color By (tertiary)",
                #             ),
                #             className="col-sm-4 my-auto",
                #         )
                #     ],
                #     className="row"
                # )

            ],
            className="card"
        ),
        #####################
        # / CAG DETAIL CARD #
        #####################        
        ###################
        # SINGLE CAG PLOT #
        ###################
        card_wrapper(
            "Individual CAG Abundance",
            [
                graph_div("single-cag", 'single-cag-graph'), 
                html.Div(
                    metadata_field_dropdown(
                        "single-cag-xaxis",
                        label_text="X-axis",
                        default_value=metadata_fields[0]
                    ) + plot_type_dropdown(
                        "single-cag-plot-type",
                        options = [
                            {'label': 'Points', 'value': 'scatter'},
                            {'label': 'Line', 'value': 'line'},
                            {'label': 'Boxplot', 'value': 'boxplot'},
                        ]
                    ) + metadata_field_dropdown(
                        "single-cag-color",
                        label_text="Color",
                    ) + metadata_field_dropdown(
                        "single-cag-facet",
                        label_text="Facet",
                    ) + log_scale_radio_button(
                        "single-cag-log"
                    ),
                    className="col-sm-4 my-auto",
                ),
            ]
        ),
        #####################
        # / SINGLE CAG PLOT #
        #####################
    ],
    className="container"
)

# Default figure layout
layout = go.Layout(
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)'
)

# Functions used to render graphs

##################
# RICHNESS GRAPH #
##################
@app.callback(
    Output('richness-graph', 'figure'),
    [
        Input('richness-metric-dropdown', 'value'),
        Input('richness-type-dropdown', 'value'),
    ])
def draw_richness(selected_metric, selected_type):
    assert selected_type in ["hist", "scatter"]
    assert selected_metric in richness_df.columns.values, (selected_metric, richness_df.columns.values)

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
                x=richness_df["n_reads"],
                y=richness_df[selected_metric],
                text=richness_df["specimen"],
                hovertemplate=hovertemplate,
                mode="markers",
            ),
        )
        fig.update_layout(
            xaxis_title="Number of Reads",
            xaxis_range=[0, richness_df["n_reads"].max() * 1.05],
        )

    else:
        assert selected_type == "hist"

        fig = go.Figure(
            data=[
                go.Histogram(
                    y=richness_df[selected_metric],
                    hovertemplate="Range: %{y}<br>Count: %{x}<extra></extra>",
                )
            ],
        )
        fig.update_layout(
            xaxis_title="Number of Specimens"
        )

    fig.update_layout(
        yaxis_range=[0, richness_df[selected_metric].max() * 1.05],
        yaxis_title=metric_names[selected_metric], 
        template="simple_white"
    )

    return fig


#####################
# CAG SUMMARY GRAPH #
#####################
@app.callback(
    Output('cag-summary-graph-hist', 'figure'),
    [
        Input('cag-summary-metric-primary', 'value'),
        Input('cag-summary-size-slider', 'value'),
        Input('cag-summary-entropy-slider', 'value'),
        Input('cag-summary-prevalence-slider', 'value'),
        Input('cag-summary-abundance-slider', 'value'),
        Input('cag-summary-nbinsx-slider', 'value'),
        Input('cag-summary-log', 'value'),
    ])
def draw_cag_summary_graph_hist(
    metric_primary,
    size_range,
    entropy_range,
    prevalence_range,
    abundance_range,
    nbinsx,
    log_scale,
):
    # Apply the filters
    plot_df = cag_summary_df.applymap(
        float
    ).query(
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
            histfunc="sum",
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


@app.callback(
    Output('cag-summary-graph-scatter', 'figure'),
    [
        Input('cag-summary-metric-primary', 'value'),
        Input('cag-summary-metric-secondary', 'value'),
        Input('cag-summary-size-slider', 'value'),
        Input('cag-summary-entropy-slider', 'value'),
        Input('cag-summary-prevalence-slider', 'value'),
        Input('cag-summary-abundance-slider', 'value'),
        Input('global-selected-cag', 'children'),
    ])
def draw_cag_summary_graph_scatter(
    metric_primary,
    metric_secondary,
    size_range,
    entropy_range,
    prevalence_range,
    abundance_range,
    selected_cag_json,
):
    # Apply the filters
    plot_df = cag_summary_df.applymap(
        float
    ).query(
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
            plot_df["CAG"].values == json.loads(selected_cag_json)["id"]
        )
    else:
        selectedpoints = None

    # Draw a scatter plot
    fig = go.Figure(
        go.Scattergl(
            x=plot_df[metric_primary],
            y=plot_df[metric_secondary],
            ids=plot_df["CAG"].values,
            text=plot_df["CAG"].values,
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

@app.callback(
    Output('cag-summary-selected-cag', 'children'),
    [
        Input('cag-summary-graph-scatter', 'clickData'),
    ])
def cag_summary_save_click_data(clickData):
    """Only save the click data when a point in a scatter has been selected"""
    if clickData is not None:
        for point in clickData["points"]:
            if "id" in point:
                # Save the time of the event
                point["time"] = time()
                return json.dumps(point, indent=2)


################
# SELECTED CAG #
################
@app.callback(
    Output('global-selected-cag', 'children'),
    [
        Input('cag-summary-selected-cag', 'children'),
        Input('volcano-selected-cag', 'children'),
    ])
def save_global_selected_cag(cag_summary_click, volcano_click):
    """Save the most recent click"""
    output_time = None
    output_dat = None

    for click_json in [cag_summary_click, volcano_click]:
        if click_json is None:
            continue
        click_dat = json.loads(click_json)
        if output_time is None or click_dat["time"] > output_time:
            output_time = click_dat["time"]
            output_dat = click_dat

    if output_dat is not None:
        return json.dumps(output_dat, indent=2)


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

# LOGIC FOR SHOW/HIDE IN ORDINATION CARD
@app.callback(
    Output('ordination-primary-pc-div', 'style'),
    [
        Input('ordination-algorithm', 'value'),
    ])
def show_hide_ordination_primary_pc(algorithm):
    if algorithm == "pca":
        return {'display': 'block'}
    else:
        return {'display': 'none'}

@app.callback(
    Output('ordination-secondary-pc-div', 'style'),
    [
        Input('ordination-algorithm', 'value'),
    ])
def show_hide_ordination_secondary_pc(algorithm):
    if algorithm == "pca":
        return {'display': 'block'}
    else:
        return {'display': 'none'}

@app.callback(
    Output('ordination-tsne-perplexity-div', 'style'),
    [
        Input('ordination-algorithm', 'value'),
    ])
def show_hide_ordination_perplexity(algorithm):
    if algorithm == "pca":
        return {'display': 'none'}
    else:
        return {'display': 'block'}

@app.callback(
    Output('ordination-graph', 'figure'),
    [
        Input('ordination-algorithm', 'value'),
        Input('ordination-metric', 'value'),
        Input('ordination-primary-pc', 'value'),
        Input('ordination-secondary-pc', 'value'),
        Input('ordination-tsne-perplexity', 'value'),
        Input('ordination-metadata', 'value'),
    ])
def draw_ordination(
    algorithm, 
    metric,
    primary_pc, 
    secondary_pc,
    perplexity,
    metadata,
):
    """Perform ordination and make the display plots."""

    assert metric in distance_df_dict, "Distance metric not found: %s" % metric

    if algorithm == "pca":
        plot_df = run_pca(
            distance_df_dict[metric]
        )

    else:
        assert algorithm == "tsne", "Algorithm not found: %s" % algorithm

        plot_df = run_tsne(
            distance_df_dict[metric],
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
            METADATA = manifest_df[metadata]
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
            sns.color_palette(
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


################
# VOLCANO PLOT #
################
@app.callback(
    Output('volcano-graph', 'figure'),
    [
        Input('volcano-parameter-dropdown', 'value'),
        Input('volcano-pvalue-slider', 'value'),
        Input('volcano-fdr-radio', 'value'),
        Input('global-selected-cag', 'children'),
    ])
def draw_volcano_plot(parameter, neg_log_pvalue_min, fdr_on_off, selected_cag_json):

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
    )

    return fig


@app.callback(
    Output('volcano-selected-cag', 'children'),
    [
        Input('volcano-graph', 'clickData'),
    ])
def volcano_save_click_data(clickData):
    """Only save the click data when a point in a scatter has been selected"""
    if clickData is not None:
        for point in clickData["points"]:
            if "id" in point:
                # Save the time of the event
                point["time"] = time()
                return json.dumps(point, indent=2)


@app.callback(
    Output('volcano-pvalue-slider', 'max'),
    [
        Input('volcano-parameter-dropdown', 'value'),
    ])
def update_volcano_pvalue_slider_max(parameter):
    return max_neg_log_pvalue[parameter]


@app.callback(
    Output('volcano-pvalue-slider', 'marks'),
    [
        Input('volcano-parameter-dropdown', 'value'),
    ])
def update_volcano_pvalue_slider_marks(parameter):
    return {
        str(int(n)): str(int(n))
        for n in np.arange(
            0,
            max_neg_log_pvalue[parameter],
            max(1, int(max_neg_log_pvalue[parameter] / 5))
        )
    }


##################
# CAG DETAIL TAX #
##################
@app.callback(
    [
        Output('cag-detail-tax-graph', 'figure'),
        Output('cag-detail-tax-ngenes', 'max'),
        Output('cag-detail-tax-ngenes', 'marks'),
    ],
    [
        Input('cag-detail-tax-ngenes', 'value'),
        Input('global-selected-cag', 'children'),
    ])
def draw_cag_detail_tax(min_ngenes, selected_cag_json):

    # Parse the selected CAG data
    if selected_cag_json is not None:
        # Set the points which are selected in the scatter plot
        cag_id = json.loads(selected_cag_json)["id"]
    else:
        # Default CAG to plot
        cag_id = 100


    # Read in the taxonomic annotations for this CAG
    cag_df = pd.read_hdf(
        hdf5_fp,
        "/annot/gene/all",
        where="CAG == {}".format(cag_id),
        columns = ["gene", "tax_id"]
    )

    # Format the DataFrame as needed to make a go.Sunburst
    cag_tax_df = make_cag_tax_df(cag_df["tax_id"])
    
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

# ######################
# # CAG DETAIL HEATMAP #
# ######################
# @app.callback(
#     Output('cag-detail-heatmap-graph', 'figure'),
#     [
#         Input('cag-detail-heatmap-color-primary', 'value'),
#         Input('cag-detail-heatmap-color-secondary', 'value'),
#         Input('cag-detail-heatmap-color-tertiary', 'value'),
#     ])
# def draw_cag_detail_heatmap():

#     fig = Go.figure([])
#     return fig

###################
# SINGLE CAG PLOT #
###################
@app.callback(
    Output('single-cag-graph', 'figure'),
    [
        Input('global-selected-cag', 'children'),
        Input('single-cag-xaxis', 'value'),
        Input('single-cag-plot-type', 'value'),
        Input('single-cag-color', 'value'),
        Input('single-cag-facet', 'value'),
        Input('single-cag-log', 'value'),
    ])
def draw_single_cag_plot(selected_cag_json, xaxis, plot_type, color, facet, log_scale):

    # Parse the selected CAG data
    if selected_cag_json is not None:
        # Set the points which are selected in the scatter plot
        cag_id = json.loads(selected_cag_json)["id"]
    else:
        # Default CAG to plot
        cag_id = 100

    plot_df = manifest_df.assign(
        CAG_ABUND = cag_abundance_df.loc[int(cag_id)]
    )

    if plot_type == "scatter":
        plot_func = px.scatter

    elif plot_type == "boxplot":
        plot_func = px.box

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


# @app.callback(
#     Output('cag-membership-table', 'data'),
#     [
#         Input('single-cag-selector', "value"),
#     ])
# def update_cag_membership_table(cag_id):
    
#     if cag_id is None or cag_id == "none":
#         return gene_annot_df.head(
#         ).to_dict(
#             'records'
#         )
#     else:

#         return gene_annot_df.query(
#             "CAG == '{}'".format(cag_id)
#         ).to_dict(
#             'records'
#         )

if __name__ == '__main__':

    app.run_server(
        debug=True
    )

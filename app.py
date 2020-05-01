#!/usr/bin/env python3

import click
import dash
import dash_table
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import plotly.graph_objs as go
from dash.dependencies import Input, Output
import numpy as np
import os
import pandas as pd
import seaborn as sns

# Read in data from the HDF5 defined in `HDF5_fp`
hdf5_fp = os.getenv("HDF5_FP")
assert os.path.exists(hdf5_fp), "Path does not exist: {}".format(hdf5_fp)
with pd.HDFStore(hdf5_fp, "r") as store:

    # Table for RICHNESS GRAPH
    richness_df = pd.read_hdf(store, "/summary/all")

    # Table for CAG SIZE GRAPH
    cag_summary_df = pd.read_hdf(store, "/annot/cag/all")

    # Specimen metadata table (manifest provided by user)
    manifest_df = pd.read_hdf(store, "/manifest")

    # Precomputed PCA (on specimens)
    pca_df = pd.read_hdf(store, "/ordination/pca")

    # Precomputed t-SNE (on specimens)
    tsne_df = pd.read_hdf(store, "/ordination/tsne")

# Precompute some useful metrics

# Precompute the proportion of reads which align
richness_df = richness_df.assign(
    prop_reads_aligned=richness_df["aligned_reads"] / richness_df["n_reads"]
)

# Round to 4 significant figures on the CAG summary metrics
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

# Set index on ordination tables
pca_df.set_index("specimen", inplace=True)
tsne_df.set_index("specimen", inplace=True)


# external CSS stylesheets
external_stylesheets = [
    {
        'href': 'https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/css/bootstrap.min.css',
        'rel': 'stylesheet',
        'integrity': 'sha384-Vkoo8x4CGsO3+Hhxv8T/Q5PaXtkKtu6ug5TOeNV6gBiFeWPGFN9MuhOf23Q9Ifjh',
        'crossorigin': 'anonymous'
    }
]

# Set up the app
app = dash.Dash(
    __name__, 
    external_stylesheets=external_stylesheets
)
app.title = "GLAM Browser"

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
        html.Div(
            [
                ##################
                # RICHNESS GRAPH #
                ##################
                html.Div(
                    [
                        html.Div(
                            [
                                html.A(id="richness"),
                                dcc.Graph(
                                    id='richness-graph'
                                )
                            ],
                            className="col-sm-8",
                        ),
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
                                html.Label('Plot Type'),
                                dcc.Dropdown(
                                    id="richness-type-dropdown",
                                    options=[
                                        {'label': 'Points', 'value': 'scatter'},
                                        {'label': 'Histogram', 'value': 'hist'},
                                    ],
                                    value='scatter'
                                ),
                            ],
                            className="col-sm-4",
                        )
                    ],
                    className="row"
                ),
                ####################
                # / RICHNESS GRAPH #
                ####################
                ##################
                # CAG SIZE GRAPH #
                ##################
                html.Div(
                    [
                        html.Div(
                            [
                                html.A(id="cag-size"),
                                dcc.Graph(
                                    id='cag-size-graph'
                                )
                            ],
                            className="col-sm-8",
                        ),
                        html.Div(
                            [
                                html.Label('CAG Size Window'),
                                dcc.RangeSlider(
                                    id="cag-size-slider",
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
                                html.Br(),
                                html.Label('Number of Bins'),
                                dcc.Slider(
                                    id="cag-nbinsx-slider",
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
                                html.Br(),
                                html.Label('Log Scale'),
                                dcc.RadioItems(
                                    id="cag-size-log",
                                    options=[
                                        {'label': 'On', 'value': 'on'},
                                        {'label': 'Off', 'value': 'off'},
                                    ],
                                    value='off',
                                )
                            ],
                            className="col-sm-4",
                        )
                    ],
                    className="row"
                ),
                ####################
                # / CAG SIZE GRAPH #
                ####################
                ####################
                # ORDINATION GRAPH #
                ####################
                html.Div(
                    [
                        html.Div(
                            [
                                html.A(id="ordination"),
                                dcc.Graph(
                                    id='ordination-graph'
                                )
                            ],
                            className="col-sm-8",
                        ),
                        html.Div(
                            [
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
                                html.Label('Metadata Label'),
                                dcc.Dropdown(
                                    id="ordination-metadata",
                                    options=[
                                        {'label': 'None', 'value': 'none'},
                                    ] + [
                                        {'label': f, "value": f}
                                        for f in metadata_fields
                                    ],
                                    value='none'
                                ),
                            ],
                            className="col-sm-4",
                        )
                    ],
                    className="row"
                ),
                ######################
                # / ORDINATION GRAPH #
                ######################
                ##################
                # CAG DATA TABLE #
                ##################
                html.Div(
                    [
                        html.Div(
                            [
                                html.A(id="cag-table"),
                                dash_table.DataTable(
                                    id='cag-table-datatable',
                                    columns=[
                                        {
                                            "name": {
                                                "mean_abundance": "Mean Abund.",
                                                "prevalence": "Prevalence",
                                                "size": "Size (# of genes)",
                                                "std_abundance": "Std. Abund.",
                                                "CAG": "CAG ID"
                                            }.get(i, i), 
                                            "id": i,
                                            "type": "numeric"
                                        }
                                        for i in cag_summary_df.columns
                                    ],
                                    data=cag_summary_df.query(
                                        "size > 1"
                                    ).to_dict('records'),
                                    page_action='native',
                                    sort_action='native',
                                    filter_action='native',
                                    row_selectable='single',
                                    page_size=10,
                                )
                            ],
                            className="col",
                        )
                    ],
                    className="row"
                ),
                ####################
                # / CAG DATA TABLE #
                ####################
                ######################
                # CAG SUMMARY FIGURE #
                ######################
                html.Div(
                    [
                        html.Div(
                            [
                                html.A(id="cag-summary"),
                                dcc.Graph(
                                    id='cag-summary-graph'
                                )
                            ],
                            className="col-sm-8",
                        ),
                        html.Div(
                            [
                                html.Label('X-axis'),
                                dcc.Dropdown(
                                    id="cag-summary-xaxis-dropdown",
                                    options=[
                                        {'label': 'CAG ID', 'value': 'CAG'},
                                        {'label': 'Size', 'value': 'size'},
                                        {'label': 'Mean Abund.', 'value': 'mean_abundance'},
                                        {'label': 'Prevalence', 'value': 'prevalence'},
                                        {'label': 'Std. Abund', 'value': 'std_abundance'},
                                    ],
                                    value='size'
                                ),
                                html.Br(),
                                html.Label('Y-axis'),
                                dcc.Dropdown(
                                    id="cag-summary-yaxis-dropdown",
                                    options=[
                                        {'label': 'CAG ID', 'value': 'CAG'},
                                        {'label': 'Size', 'value': 'size'},
                                        {'label': 'Mean Abund.', 'value': 'mean_abundance'},
                                        {'label': 'Prevalence', 'value': 'prevalence'},
                                        {'label': 'Std. Abund', 'value': 'std_abundance'},
                                    ],
                                    value='mean_abundance'
                                ),
                                html.Br(),
                                html.Label('Size Filter'),
                                dcc.RangeSlider(
                                    id="cag-summary-size-slider",
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
                                html.Br(),
                                html.Label('Prevalence Filter'),
                                dcc.RangeSlider(
                                    id="cag-summary-prevalence-slider",
                                    min=0,
                                    max=1,
                                    step=0.01,
                                    marks={
                                        "0": "0",
                                        "0.5": "0.5",
                                        "1": "1"
                                    },
                                    value=[0, 1]
                                ),
                                html.Br(),
                                html.Label('Plot Type'),
                                dcc.Dropdown(
                                    id="cag-summary-type-dropdown",
                                    options=[
                                        {'label': 'Points', 'value': 'scatter'},
                                        {'label': 'Histogram', 'value': 'histogram'},
                                    ],
                                    value='scatter'
                                ),
                            ],
                            className="col-sm-4",
                        )
                    ],
                    className="row"
                ),
                ########################
                # / CAG SUMMARY FIGURE #
                ########################
            ],
            className="container"
        )
    ]
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
        
    else:
        assert selected_type == "hist"

        fig = go.Figure(
            data=[
                go.Histogram(
                    y=richness_df[selected_metric],
                )
            ],
        )

    fig.update_layout(
        title={
            'text': "Gene Analysis Summary",
            'y': 0.9,
            'x': 0.5, 
            'xanchor': 'center',
            'yanchor': 'top'
        },
        yaxis_title={
            "prop_reads_aligned": "Prop. Reads Aligned",
            "n_genes_aligned": "Num. Genes Aligned",
            "n_genes_assembled": "Num. Genes Assembled",
        }[selected_metric], 
        template="simple_white"
    )

    return fig


##################
# CAG SIZE GRAPH #
##################
@app.callback(
    Output('cag-size-graph', 'figure'),
    [
        Input('cag-size-slider', 'value'),
        Input('cag-nbinsx-slider', 'value'),
        Input('cag-size-log', 'value'),
    ])
def draw_cag_size(selected_range, nbinsx, log_scale):

    plot_vals = cag_size_log[
        (cag_size_log >= selected_range[0]) & 
        (cag_size_log <= selected_range[1])
    ]

    fig = go.Figure(
        data=[
            go.Histogram(
                x=plot_vals,
                histfunc="sum",
                nbinsx=nbinsx,
            )
        ],
    )

    if log_scale == "on":
        fig.update_layout(yaxis_type="log")

    fig.update_layout(
        title={
            'text': "CAG Size",
            'y': 0.9,
            'x': 0.5, 
            'xanchor': 'center',
            'yanchor': 'top'
        },
        xaxis_title="CAG Size (log10 number of genes per CAG)",
        yaxis_title="Number of Genes (per bin)", 
        template="simple_white"
    )

    return fig


####################
# ORDINATION GRAPH #
####################
@app.callback(
    Output('ordination-graph', 'figure'),
    [
        Input('ordination-algorithm', 'value'),
        Input('ordination-metadata', 'value'),
    ])
def draw_ordination(algorithm, metadata):

    if algorithm == "pca":
        plot_df = pca_df
    else:
        assert algorithm == "tsne"
        plot_df = tsne_df

    if metadata == "none":

        fig = go.Figure(
            data=go.Scatter(
                x=plot_df[plot_df.columns.values[0]],
                y=plot_df[plot_df.columns.values[1]],
                ids=plot_df.index.values,
                text=plot_df.index.values,
                hoverinfo="text",
                mode="markers",
            ),
        )

    else:

        plot_df[metadata] = manifest_df[metadata]

        cmap = dict(zip(
            plot_df[metadata].unique(),
            sns.color_palette(
                "colorblind",
                plot_df[metadata].unique().shape[0],
            ).as_hex()
        ))

        fig = go.Figure()

        for group_name, group_plot_df in plot_df.groupby(metadata):
            fig.add_trace(
                go.Scatter(
                    x=group_plot_df[plot_df.columns.values[0]],
                    y=group_plot_df[plot_df.columns.values[1]],
                    ids=group_plot_df.index.values,
                    text=group_plot_df.index.values,
                    hoverinfo="text",
                    mode="markers",
                    marker={
                        "color": cmap[group_name]
                    },
                    legendgroup=metadata,
                    name=group_name,
                ),
            )

        fig.update_layout(showlegend=True)

    fig.update_layout(
        title={
            'text': "PCA" if algorithm == "pca" else "t-SNE",
            'y': 0.9,
            'x': 0.5, 
            'xanchor': 'center',
            'yanchor': 'top',
        },
        xaxis_title=plot_df.columns.values[0],
        yaxis_title=plot_df.columns.values[1],
        template="simple_white"
    )

    return fig

####################
# CAG SUMMARY GRAPH #
####################
@app.callback(
    Output('cag-summary-graph', 'figure'),
    [
        Input('cag-summary-xaxis-dropdown', 'value'),
        Input('cag-summary-yaxis-dropdown', 'value'),
        Input('cag-summary-size-slider', 'value'),
        Input('cag-summary-prevalence-slider', 'value'),
        Input('cag-summary-type-dropdown', 'value'),
    ])
def draw_ordination(xaxis, yaxis, size_range, prevalence_range, plot_type):

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
    ).apply(
        lambda c: c.apply(np.log10) if c.name == "size" else c
    )

    assert plot_type in ["scatter", "histogram"]
    if plot_type == "scatter":

        fig = go.Figure(
            data=go.Scattergl(
                x=plot_df[xaxis],
                y=plot_df[yaxis],
                ids=plot_df["CAG"].values,
                text=plot_df["CAG"].values,
                hoverinfo="text",
                mode="markers",
            ),
        )

    else:

        fig = go.Figure(
            data=[
                go.Histogram(
                    y=plot_df[xaxis],
                )
            ],
        )

    axis_names = {
        "CAG": "CAG ID",
        "size": "Number of Genes (log10)",
        "mean_abundance": "Mean Abundance",
        "std_abundance": "Std. Abundance",
        "prevalence": "Prevalence",
    }

    fig.update_layout(
        title={
            'text': "CAG Summary Metrics",
            'y': 0.9,
            'x': 0.5, 
            'xanchor': 'center',
            'yanchor': 'top',
        },
        xaxis_title=axis_names.get(xaxis, xaxis),
        yaxis_title=axis_names.get(yaxis, yaxis),
        template="simple_white"
    )

    return fig

if __name__ == '__main__':

    app.run_server(
        debug=True
    )

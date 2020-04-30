#!/usr/bin/env python3

import click
import dash
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

# The limits of CAG sizes
cag_size_min = cag_summary_df["size"].min()
cag_size_max = cag_summary_df["size"].max()
cag_size_range = cag_size_max - cag_size_min
cag_size_step = int(max(1, cag_size_range / 100))

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
                        dbc.DropdownMenuItem("Genes Detected", href="#richness")
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
                                    min=cag_size_min,
                                    max=cag_size_max,
                                    step=cag_size_step,
                                    marks={
                                        str(n): str(n)
                                        for n in np.arange(cag_size_min, cag_size_max, max(1, int(cag_size_range / 5)))
                                    },
                                    value=[
                                        cag_size_min,
                                        cag_size_max
                                    ]
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
        Input('cag-size-log', 'value'),
    ])
def draw_cag_size(selected_range, log_scale):

    plot_vals = cag_summary_df.query(
        "size >= %s" % selected_range[0]
    ).query(
        "size <= %s" % selected_range[1]
    )["size"]

    if log_scale == "on":
        plot_vals = plot_vals.apply(np.log10)

    fig = go.Figure(
        data=[
            go.Histogram(
                x=plot_vals,
                histfunc="sum"
            )
        ],
    )

    fig.update_layout(
        title={
            'text': "CAG Size",
            'y': 0.9,
            'x': 0.5, 
            'xanchor': 'center',
            'yanchor': 'top'
        },
        xaxis_title="CAG Size ({}number of genes per CAG)".format("log10 " if log_scale == "on" else ""), 
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

if __name__ == '__main__':

    app.run_server(
        debug=True
    )

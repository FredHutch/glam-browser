#!/usr/bin/env python3

import click
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import plotly.graph_objs as go
from dash.dependencies import Input, Output
import os
import pandas as pd

# Read in data from the HDF5 defined in `HDF5_fp`
hdf5_fp = os.getenv("HDF5_FP")
assert os.path.exists(hdf5_fp), "Path does not exist: {}".format(hdf5_fp)
with pd.HDFStore(hdf5_fp, "r") as store:
    richness_df = pd.read_hdf(store, "/summary/all")

# Precompute some useful metrics
richness_df = richness_df.assign(
    prop_reads_aligned=richness_df["aligned_reads"] / richness_df["n_reads"]
)

# external CSS stylesheets
external_stylesheets = [
    {
        'href': 'https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/css/bootstrap.min.css',
        'rel': 'stylesheet',
        'integrity': 'sha384-Vkoo8x4CGsO3+Hhxv8T/Q5PaXtkKtu6ug5TOeNV6gBiFeWPGFN9MuhOf23Q9Ifjh',
        'crossorigin': 'anonymous'
    }
]

app = dash.Dash(
    __name__, 
    external_stylesheets=external_stylesheets
)
app.title = "GLAM Browser"

app.layout = html.Div(
    children=[
        html.H1(children='GLAM Browser'),

        html.P(children='''
            Interactive results browser for gene-level metagenomic analysis with geneshot.
        '''),

        html.P(children='''
            Data HDF5: %s
        ''' % os.getenv("HDF5_FP").split("/")[-1]),

        html.Div(
            [
                html.Div([
                    dcc.Graph(
                        id='richness-graph'
                    )
                ]),
                html.Div([
                    html.Label('Display Values'),
                    dcc.Dropdown(
                        id="richness-dropdown",
                        options=[
                            {'label': 'Genes Assembled (#)', 'value': 'n_genes_assembled'},
                            {'label': 'Genes Aligned (#)', 'value': 'n_genes_aligned'},
                            {'label': 'Reads Aligned (%)', 'value': 'prop_reads_aligned'},
                        ],
                        value='prop_reads_aligned'
                    ),
                ])
            ],
            className="row"
        )
    ],
    className="container",
)

# Default figure layout
layout = go.Layout(
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)'
)

# Functions used to render graphs
@app.callback(
    Output('richness-graph', 'figure'),
    [Input('richness-dropdown', 'value')])
def draw_richness(selected_metric):
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
        title={
            'text': "Gene Analysis Summary",
            'y': 0.9,
            'x': 0.5, 
            'xanchor': 'center',
            'yanchor': 'top'
        },
        xaxis_title="Number of Reads",
        yaxis_title={
            "prop_reads_aligned": "Prop. Reads Aligned",
            "n_genes_aligned": "Num. Genes Aligned",
            "n_genes_assembled": "Num. Genes Assembled",
        }[selected_metric], 
        template="simple_white"
    )

    return fig

if __name__ == '__main__':

    app.run_server(
        debug=True
    )

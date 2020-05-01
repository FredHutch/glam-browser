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

    # Corncob results (for volcano plot)
    if "/stats/cag/corncob_wide" in store:
        corncob_df = pd.read_hdf(store, "/stats/cag/corncob_wide")
    elif "/stats/cag/corncob" in store:
        corncob_df = pd.read_hdf(store, "/stats/cag/corncob")
    else:
        corncob_df = None

    # CAG abundance
    cag_abundance_df = pd.read_hdf(store, "/abund/cag/wide")

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

# Set index on CAG abundance table
cag_abundance_df.set_index("CAG", inplace=True)

# Calculate the -log10(p_value)
corncob_df = corncob_df.assign(
    neg_log_pvalue = corncob_df["p_value"].apply(np.log10) * -1
)
max_neg_log_pvalue = corncob_df["neg_log_pvalue"].max()

# external CSS stylesheets
external_stylesheets = [
    {
        'href': 'https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/css/bootstrap.min.css',
        'rel': 'stylesheet',
        'integrity': 'sha384-Vkoo8x4CGsO3+Hhxv8T/Q5PaXtkKtu6ug5TOeNV6gBiFeWPGFN9MuhOf23Q9Ifjh',
        'crossorigin': 'anonymous'
    }
]

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

def graph_div(anchor_id, graph_id):
    """Return a div containing a dcc.Graph and anchor, all within a col-sm-8."""
    return html.Div(
        [
            html.A(id=anchor_id),
            dcc.Graph(
                id=graph_id
            )
        ],
        className="col-sm-8",
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

        return [
            html.Label(label_text),
            dcc.Dropdown(
                id=dropdown_id,
                options=[
                    {'label': l, 'value': l}
                    for l in parameter_list
                ],
                value=parameter_list[0]
            ),
            html.Br(),
        ]


def volcano_pvalue_slider(slider_id, label_text='P-Value Filter'):
    return [
        html.Label(label_text),
        dcc.Slider(
            id=slider_id,
            min=0,
            max=max_neg_log_pvalue,
            step=0.1,
            marks={
                str(int(n)): str(int(n))
                for n in np.arange(
                    0,
                    max_neg_log_pvalue,
                    max(1, int(max_neg_log_pvalue / 5))
                )
            },
            value=1
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
    

def cag_prevalence_slider(slider_id, label_text='CAG Prevalence Filter'):
    return [
        html.Label(label_text),
        dcc.RangeSlider(
            id=slider_id,
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
        html.Br()
    ]
    

def cag_abundance_slider(slider_id, label_text='CAG Abundance Filter'):
    return [
        html.Label(label_text),
        dcc.RangeSlider(
            id=slider_id,
            min=0,
            max=1,
            step=0.001,
            marks={
                "0": "0",
                "1": "1",
            },
            value=[
                0, 
                cag_summary_df["mean_abundance"].apply(float).max()
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
                {'label': 'Mean Abund.', 'value': 'mean_abundance'},
                {'label': 'Prevalence', 'value': 'prevalence'},
                {'label': 'Std. Abund', 'value': 'std_abundance'},
            ],
            value=default_value
        ),
        html.Br()
    ]

def log_scale_radio_button(id_string):
    return [
        html.Br(),
        html.Label('Log Scale'),
        dcc.RadioItems(
            id=id_string,
            options=[
                {'label': 'On', 'value': 'on'},
                {'label': 'Off', 'value': 'off'},
            ],
            value='off',
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
        ##################
        # CAG SIZE GRAPH #
        ##################
        card_wrapper(
            "CAG Size Summary",
            [
                graph_div("cag-size", 'cag-size-graph'),
                html.Div(
                    cag_size_slider(
                        "cag-size-slider"
                    ) + nbins_slider(
                        "cag-nbinsx-slider"
                    ) + log_scale_radio_button(
                        "cag-size-log"
                    ),
                    className="col-sm-4 my-auto",
                )
            ]
        ),
        ####################
        # / CAG SIZE GRAPH #
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
                    ] + metadata_field_dropdown(
                        "ordination-metadata"
                    ) + plot_type_dropdown(
                        "ordination-type-dropdown",
                        label_text="Metadata Overlay",
                        options=[
                            {'label': 'Color', 'value': 'color'},
                            {'label': 'Boxplot', 'value': 'boxplot'},
                            {'label': 'Scatter', 'value': 'scatter'},
                        ],
                        default_value="color"

                    ),
                    className="col-sm-4 my-auto",
                )
            ],
        ),
        ######################
        # / ORDINATION GRAPH #
        ######################
        #########################
        # CAG SUMMARY HISTOGRAM #
        #########################
        card_wrapper(
            "CAG Summary Histogram",
            [
                graph_div("cag-summary-histogram", 'cag-summary-histogram-graph'),
                html.Div(
                    cag_metric_dropdown(
                        "cag-summary-histogram-metric-dropdown",
                        default_value="size"
                    ) + cag_size_slider(
                        "cag-summary-histogram-size-slider"
                    ) + cag_prevalence_slider(
                        "cag-summary-histogram-prevalence-slider"
                    ) + cag_abundance_slider(
                        "cag-summary-histogram-abundance-slider"
                    ) + nbins_slider(
                        "cag-summary-histogram-nbinsx-slider"
                    ) + log_scale_radio_button(
                        "cag-summary-histogram-log"
                    ),
                    className="col-sm-4 my-auto",
                )
            ]
        ),
        ###########################
        # / CAG SUMMARY HISTOGRAM #
        ###########################
        #######################
        # CAG SUMMARY SCATTER #
        #######################
        card_wrapper(
            "CAG Summary Scatterplot",
            [
                graph_div("cag-summary-scatter", 'cag-summary-scatter-graph'), 
                html.Div(
                    cag_metric_dropdown(
                        "cag-summary-scatter-xaxis-dropdown",
                        label_text="X-axis",
                        default_value="size"
                    ) + cag_metric_dropdown(
                        "cag-summary-scatter-yaxis-dropdown",
                        label_text="Y-axis",
                        default_value="mean_abundance"
                    ) + cag_size_slider(
                        "cag-summary-scatter-size-slider"
                    ) + cag_prevalence_slider(
                        "cag-summary-scatter-prevalence-slider"
                    ) + cag_abundance_slider(
                        "cag-summary-scatter-abundance-slider"
                    ),
                    className="col-sm-4 my-auto",
                )
            ]
        ),
        #########################
        # / CAG SUMMARY SCATTER #
        #########################
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
                    ),
                    className="col-sm-4 my-auto",
                )
            ]
        ),
        ##################
        # / VOLCANO PLOT #
        ##################
        ###################
        # SINGLE CAG PLOT #
        ###################
        card_wrapper(
            "Individual CAG Abundance",
            [
                graph_div("single-cag", 'single-cag-graph'), 
                html.Div(
                    [
                        html.Label("CAG"),
                        dcc.Dropdown(
                            id="single-cag-selector",
                            options=[
                                {'label': 'CAG {}'.format(cag_id), 'value': str(cag_id)}
                                for cag_id in cag_summary_df["CAG"].values
                            ],
                            placeholder="Select a CAG",
                        ),
                        html.Br(),
                    ] + metadata_field_dropdown(
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
                    ),
                    className="col-sm-4 my-auto",
                )
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
            xaxis_title="Number of Reads"
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
            xaxis_title="Number of Specimens"
        )

    fig.update_layout(
        title={
            'text': "Gene Analysis Summary",
            'y': 0.9,
            'x': 0.5, 
            'xanchor': 'center',
            'yanchor': 'top'
        },
        yaxis_title=metric_names[selected_metric], 
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
        Input('ordination-type-dropdown', 'value'),
    ])
def draw_ordination(algorithm, metadata, plot_type):

    if algorithm == "pca":
        plot_df = pca_df
    else:
        assert algorithm == "tsne"
        plot_df = tsne_df

    if metadata is None or metadata == "none":

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

        fig.update_layout(
            xaxis_title=plot_df.columns.values[0],
            yaxis_title=plot_df.columns.values[1],
        )

    else:

        # Add the metadata to the table
        plot_df[metadata] = manifest_df[metadata]

        # Overlay the metadata as color
        if plot_type == "color":

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

            fig.update_layout(
                xaxis_title=plot_df.columns.values[0],
                yaxis_title=plot_df.columns.values[1],
                showlegend=True,
            )

        # Make a boxplot by metadata
        elif plot_type == "boxplot":

            fig = px.box(
                plot_df,
                x=metadata, 
                y=plot_df.columns.values[0]
            )

            fig.update_layout(
                xaxis_title=metadata,
                yaxis_title=plot_df.columns.values[0],
            )

        # Make a scatter plot with the metadata
        else:
            assert plot_type == "scatter"

            fig = px.scatter(
                plot_df,
                x=metadata, 
                y=plot_df.columns.values[0]
            )

            fig.update_layout(
                xaxis_title=metadata,
                yaxis_title=plot_df.columns.values[0],
            )

    fig.update_layout(
        title={
            'text': "PCA" if algorithm == "pca" else "t-SNE",
            'y': 0.9,
            'x': 0.5, 
            'xanchor': 'center',
            'yanchor': 'top',
        },
        template="simple_white"
    )

    return fig

#########################
# CAG SUMMARY HISTOGRAM #
#########################
@app.callback(
    Output('cag-summary-histogram-graph', 'figure'),
    [
        Input('cag-summary-histogram-metric-dropdown', 'value'),
        Input('cag-summary-histogram-size-slider', 'value'),
        Input('cag-summary-histogram-prevalence-slider', 'value'),
        Input('cag-summary-histogram-abundance-slider', 'value'),
        Input('cag-summary-histogram-nbinsx-slider', 'value'),
        Input('cag-summary-histogram-log', 'value'),
    ])
def draw_cag_summary_histogram(metric, size_range, prevalence_range, abundance_range, nbinsx, log_scale):

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
    ).apply(
        lambda c: c.apply(np.log10) if c.name == "size" else c
    )

    axis_names = {
        "CAG": "CAG ID",
        "size": "Number of Genes (log10)",
        "mean_abundance": "Mean Abundance",
        "std_abundance": "Std. Abundance",
        "prevalence": "Prevalence",
    }

    fig = go.Figure(
        data=[
            go.Histogram(
                x=plot_df[metric],
                nbinsx=nbinsx,
            )
        ],
    )

    if log_scale == "on":
        fig.update_layout(yaxis_type="log")

    fig.update_layout(
        xaxis_title=axis_names.get(metric, metric),
        yaxis_title="Number of CAGs",
        title={
            'text': "CAG Summary Histogram",
            'y': 0.9,
            'x': 0.5, 
            'xanchor': 'center',
            'yanchor': 'top',
        },
        template="simple_white"
    )

    return fig


#######################
# CAG SUMMARY SCATTER #
#######################
@app.callback(
    Output('cag-summary-scatter-graph', 'figure'),
    [
        Input('cag-summary-scatter-xaxis-dropdown', 'value'),
        Input('cag-summary-scatter-yaxis-dropdown', 'value'),
        Input('cag-summary-scatter-size-slider', 'value'),
        Input('cag-summary-scatter-prevalence-slider', 'value'),
        Input('cag-summary-scatter-abundance-slider', 'value'),
    ])
def draw_cag_summary_scatter(xaxis, yaxis, size_range, prevalence_range, abundance_range):

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
    ).apply(
        lambda c: c.apply(np.log10) if c.name == "size" else c
    )

    axis_names = {
        "CAG": "CAG ID",
        "size": "Number of Genes (log10)",
        "mean_abundance": "Mean Abundance",
        "std_abundance": "Std. Abundance",
        "prevalence": "Prevalence",
    }

    fig = go.Figure(
        data=go.Scattergl(
            x=plot_df[xaxis],
            y=plot_df[yaxis],
            ids=plot_df["CAG"].values,
            text=plot_df["CAG"].values,
            hovertemplate="CAG %{id}<br>X-value: %{x}<br>Y-value: %{y}",
            mode="markers",
        ),
    )

    fig.update_layout(
        xaxis_title=axis_names.get(xaxis, xaxis),
        yaxis_title=axis_names.get(yaxis, yaxis),
        title={
            'text': "CAG Summary Scatterplot",
            'y': 0.9,
            'x': 0.5, 
            'xanchor': 'center',
            'yanchor': 'top',
        },
        template="simple_white"
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
    ])
def draw_volcano_plot(parameter, neg_log_pvalue_min):

    if corncob_df is None:
        return go.Figure()

    # Subset to just this parameter
    plot_df = corncob_df.query(
        "parameter == '{}'".format(parameter)
    ).query(
        "neg_log_pvalue >= {}".format(neg_log_pvalue_min)
    )

    fig = go.Figure(
        data=go.Scattergl(
            x=plot_df["estimate"],
            y=plot_df["neg_log_pvalue"],
            ids=plot_df["CAG"].values,
            text=plot_df["CAG"].values,
            hovertemplate="CAG %{id}<br>Estimate: %{x}<br>p-value (-log10): %{y}<extra></extra>",
            mode="markers",
            opacity=0.5,
        ),
    )

    fig.update_layout(
        xaxis_title="Estimated Coefficient",
        yaxis_title="p-value (-log10)",
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


###################
# SINGLE CAG PLOT #
###################
@app.callback(
    Output('single-cag-graph', 'figure'),
    [
        Input('single-cag-selector', 'value'),
        Input('single-cag-xaxis', 'value'),
        Input('single-cag-plot-type', 'value'),
        Input('single-cag-color', 'value'),
        Input('single-cag-facet', 'value'),
    ])
def draw_single_cag_plot(cag_id, xaxis, plot_type, color, facet):

    if cag_id is None or cag_id == "none":
        fig = go.Figure()
        fig.update_layout(
            template="simple_white",
        )
        return fig

    plot_df = manifest_df.assign(
        CAG_ABUND = cag_abundance_df.loc[int(cag_id)]
    )

    if plot_type == "scatter":
        fig = px.scatter(
            plot_df,
            x = xaxis,
            y = "CAG_ABUND",
            color = None if color == "none" else color,
            facet_col = None if facet == "none" else facet
        )

    elif plot_type == "boxplot":
        fig = px.box(
            plot_df,
            x = xaxis,
            y = "CAG_ABUND",
            color = None if color == "none" else color,
            facet_col = None if facet == "none" else facet
        )

    else:
        assert plot_type == "line"
        fig = px.line(
            plot_df,
            x = xaxis,
            y = "CAG_ABUND",
            color = None if color == "none" else color,
            facet_col = None if facet == "none" else facet
        )

    fig.update_layout(
        template="simple_white",
        yaxis_title="CAG {}".format(cag_id)
    )
    return fig

if __name__ == '__main__':

    app.run_server(
        debug=True
    )

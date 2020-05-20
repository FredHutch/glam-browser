#!/usr/bin/env python3
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_table
import pandas as pd
import numpy as np

# NAVBAR AT THE TOP
def navbar_simple(page_data):
    return dbc.NavbarSimple(
        brand=page_data["page_title"],
        dark=True,
        color="#112345",
        children=[
            dbc.Button(
                'Main Menu',
                id={
                    "type": "open-dataset-button",
                    "index": -1,
                },
                n_clicks=0,
            ),
            html.Div(  # Store the button-press time
                id={
                    "type": "open-dataset-pressed",
                    "index": -1
                },
                style={"display": "none"}
            ),
            html.Div(  # Store which dataset is selected
                id="selected-dataset",
                children=["-1"],
                style={"display": "none"}
            )
        ],
    )


# Summary card for an individual dataset
def dataset_summary_card(ix, dataset):
    return html.Div([
        html.Br(),
        dbc.Card([
            dbc.CardHeader([
                dbc.Row([
                    dbc.Col(
                        html.H4(dataset["name"]),
                        width=6
                    ),
                    dbc.Col(
                        html.Div([
                            dbc.Button(
                                'Open Dataset',
                                id={
                                    "type": "open-dataset-button",
                                    "index": ix,
                                },
                                n_clicks=0,
                                color="primary",
                            ),
                            html.Div( # Hidden div storing whether the dataset has been opened
                                id={
                                    "type": "open-dataset-pressed",
                                    "index": ix
                                },
                                style={"display": "none"}
                            )],
                        style={"text-align": "right"}
                    ),
                        width=6
                    ),
                ])
            ]),
            dbc.CardBody(
                dcc.Markdown(
                    dataset.get("description", "")
                )
            )
        ])
    ])


###########################
# EXPERIMENT SUMMARY CARD #
###########################
def experiment_summary_card():
    return html.Div([
        html.Br(),
        dbc.Card([
            dbc.CardHeader("Experiment"),
            dbc.CardBody(html.Div(id="experiment-summary-card"))
        ])
    ])
def update_experiment_summary_card(dataset_metrics):
    # Make a table with the basic summary of an experiment
    return dbc.Table([
        html.Tbody(
            exp_table_row(  # Wrapper for each row
                "Total Reads",
                "{:,}".format(dataset_metrics.loc["total_reads"]),
                "Aligned Reads",
                "{}%".format(round(
                    100 * \
                    dataset_metrics.loc["aligned_reads"] / \
                    dataset_metrics.loc["total_reads"],
                    1
                ))
            ) + exp_table_row(
                "Genes (#)",
                "{:,}".format(dataset_metrics.loc["num_genes"]),
                "CAGs (#)",
                "{:,}".format(dataset_metrics.loc["num_cags"])
            ) + exp_table_row(
                "Specimens (#)",
                "{:,}".format(dataset_metrics.loc["num_samples"]),
                "Formula",
                "{}".format(dataset_metrics.loc["formula"])
            )
        )
    ], bordered=True, hover=True, responsive=True)
#############################
# / EXPERIMENT SUMMARY CARD #
#############################


#################
# RICHNESS CARD #
#################
def richness_card():
    return card_wrapper(
        "Gene Analysis Summary",
        dbc.Row([
            dbc.Col([
                dbc.Spinner(dcc.Graph(
                    id="richness-graph"
                ))
            ],
                align="center",
                width=8,
            ),
            dbc.Col([
                html.Br(),
                html.Label('Display Values'),
                dcc.Dropdown(
                    id="richness-metric-dropdown",
                    options=[
                        {'label': 'Genes Assembled (#)',
                            'value': 'n_genes_assembled'},
                        {'label': 'Genes Aligned (#)',
                            'value': 'n_genes_aligned'},
                        {'label': 'Reads Aligned (%)',
                            'value': 'prop_reads_aligned'},
                    ],
                    value='prop_reads_aligned'
                ),
                html.Br(),
            ] + plot_type_dropdown(
                "richness-type-dropdown"
            ) + log_scale_radio_button(
                "richness-log-x",
                label_text="Number of Reads - Log Scale"
            ),
                align="center",
                width=4,
            )
        ])
    )
###################
# / RICHNESS CARD #
###################


######################
# SINGLE SAMPLE CARD #
######################
def single_sample_card():
    return card_wrapper(
        "Single Sample Display",
        dbc.Row([
            dbc.Col([
                dbc.Spinner(dcc.Graph(
                            id="single-sample-graph"
                            ))
            ],
                align="center",
                width=8,
            ),
            dbc.Col([
                html.Label('Select Sample'),
                dcc.Dropdown(
                    id="single-sample-primary-dropdown",
                    options=[
                        {'label': 'None', 'value': 'none'}
                    ],
                ),
                html.Br(),
                html.Label('Compare Against'),
                dcc.Dropdown(
                    id="single-sample-secondary-dropdown",
                    options=[
                        {'label': 'CAG Size', 'value': 'cag_size'}
                    ],
                ),
                html.Br(),
                html.Label('Display Metric'),
                dcc.Dropdown(
                    id="single-sample-metric",
                    options=[
                        {'label': 'Relative Abundance', 'value': 'prop'},
                        {'label': 'Centered Log-Ratio', 'value': 'clr'},
                    ],
                    value="prop"
                ),
            ],
                align="center",
                width=4,
            )
        ])
    )
########################
# / SINGLE SAMPLE CARD #
########################


###################
# ORDINATION CARD #
###################
def ordination_card():
    return card_wrapper(
        "Ordination Analysis",
        dbc.Row([
            dbc.Col(
                dbc.Spinner(dcc.Graph(
                            id="ordination-graph"
                            )),
                width=8,
                align="center"
            ),
            dbc.Col([
                html.Label('Distance Metric'),
                dcc.Dropdown(
                    id="ordination-metric",
                    options=[
                        {'label': 'Euclidean', 'value': 'euclidean'},
                        {'label': 'Aitchison', 'value': 'aitchison'},
                        {'label': 'Bray-Curtis', 'value': 'braycurtis'},
                    ],
                    value='braycurtis'
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
            ) + [
                html.Div(
                    id="ordination-anosim-results"
                )
            ],
                width=4,
                align="center"
            )
        ])
    )
#####################
# / ORDINATION CARD #
#####################


####################
# CAG SUMMARY CARD #
####################
def cag_summary_card():
    return card_wrapper(
        "CAG Summary",
        dbc.Row([
            dbc.Col(
                [
                    dbc.Spinner(dcc.Graph(
                                id='cag-summary-graph-hist'
                                )),
                    dbc.Spinner(dcc.Graph(
                                id='cag-summary-graph-scatter'
                                ))
                ],
                width=8,
                align="center"
            ),
            dbc.Col([
                html.Br(),
                html.Label("Histogram Display"),
                dcc.Dropdown(
                    id='cag-summary-histogram-metric',
                    options=[
                        {'label': 'Number of genes', 'value': 'genes'},
                        {'label': 'Number of CAGs', 'value': 'cags'},
                    ],
                    value="genes",
                ),
                html.Br(),
                html.Label("Histogram Log Scale"),
                dcc.Dropdown(
                    id='cag-summary-histogram-log',
                    options=[
                        {'label': 'On', 'value': 'on'},
                        {'label': 'Off', 'value': 'off'},
                    ],
                    value="on",
                ),
                html.Br(),
                html.Div(id='global-selected-cag',
                         style={"display": "none"}),
                html.Div(id='cag-summary-selected-cag',
                         style={"display": "none"}),
            ] + nbins_slider(
                "cag-summary-nbinsx-slider"
            ) + cag_metric_dropdown(
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
            ),
                width=4,
                align="center"
            )
        ])
    )
######################
# / CAG SUMMARY CARD #
######################


####################
# CAG HEATMAP CARD #
####################
def cag_heatmap_card():
    return card_wrapper(
        "CAG Abundance Heatmap",
        dbc.Row([
            dbc.Col(
                dbc.Spinner(dcc.Graph(
                    id='cag-heatmap-graph'
                )),
                width=8,
                align="center"
            ),
            dbc.Col(
                [
                    html.Label("Display CAGs"),
                    dcc.Dropdown(
                        id="cag-heatmap-cag-dropdown",
                        options=[],
                        value=[],
                        multi=True
                    ),
                    html.Br(),
                    html.Label("Display Metadata"),
                    dcc.Dropdown(
                        id="cag-heatmap-metadata-dropdown",
                        options=[],
                        value=[],
                        multi=True
                    ),
                    html.Div(
                        children=[-1],
                        id="cag-heatmap-selected-dataset",
                        style={"display": "none"}
                    ),
                    html.Br(),
                    html.Label("Abundance Metric"),
                    dcc.Dropdown(
                        id='cag-heatmap-abundance-metric',
                        options=[
                            {'label': 'Rel. Abund. (log10)', 'value': 'log10'},
                            {'label': 'Rel. Abund. (log10) (z-score)', 'value': 'zscore'},
                            {'label': 'Rel. Abund.', 'value': 'raw'},
                        ],
                        value="log10",
                    ),
                    html.Br(),
                    html.Label("Group Specimens"),
                    dcc.Dropdown(
                        id='cag-heatmap-cluster',
                        options=[
                            {'label': 'By Metadata', 'value': 'metadata'},
                            {'label': 'By CAG Abundances', 'value': 'cag'},
                        ],
                        value="metadata",
                    ),
                    html.Br(),
                    html.Label("Show Taxonomy"),
                    dcc.Dropdown(
                        id='cag-heatmap-taxa-rank',
                        options=[
                            {'label': 'None', 'value': 'none'},
                            {'label': 'Species', 'value': 'species'},
                            {'label': 'Genus', 'value': 'genus'},
                            {'label': 'Family', 'value': 'family'},
                            {'label': 'Class', 'value': 'class'},
                            {'label': 'Phylum,', 'value': 'Phylum'},
                        ],
                        value="none",
                    ),
                ],
                width=4,
                align="center"
            )
        ])
)
######################
# / CAG HEATMAP CARD #
######################


################
# VOLCANO PLOT #
################
def volcano_card():
    return card_wrapper(
        "Association Screening",
        dbc.Row([
            dbc.Col(
                dbc.Spinner(dcc.Graph(
                            id='volcano-graph'
                            )),
                width=8,
                align="center"
            ),
            dbc.Col(
                corncob_parameter_dropdown(
                    "volcano-parameter-dropdown",
                ) + volcano_pvalue_slider(
                    "volcano-pvalue-slider",
                ) + log_scale_radio_button(
                    "volcano-fdr-radio",
                    label_text="FDR-BH adjustment"
                ) + [
                    html.Br(),
                    html.Label("Compare Against"),
                    dcc.Dropdown(
                        id="corncob-comparison-parameter-dropdown",
                        options=[
                            {'label': 'Estimated Coefficient', 'value': 'coef'},
                        ],
                        value="coef"
                    ),
                    html.Div(id='volcano-selected-cag',
                             style={"display": "none"}),
                ],
                width=4,
                align="center"
            )
        ])
    )
##################
# / VOLCANO PLOT #
##################

#####################
# CAG TAXONOMY CARD #
#####################
def taxonomy_card():
    return card_wrapper(
        "CAG Taxonomy",
        dbc.Row([
            dbc.Col(
                dbc.Spinner(dcc.Graph(
                            id="cag-tax-graph"
                            )),
                width=8,
                align="center",
            ),
            dbc.Col(
                basic_slider(
                    "cag-tax-ngenes",
                    "Minimum Number of Genes",
                    included=False,
                    default_value=5,
                ),
                width=4,
                align="center",
            )
        ])
    )
#######################
# / CAG TAXONOMY CARD #
#######################


###################
# SINGLE CAG CARD #
###################
def single_cag_card():
    return card_wrapper(
        "Individual CAG Abundance",
        dbc.Row([
            dbc.Col(
                dbc.Spinner(dcc.Graph(
                            id="single-cag-graph"
                            )),
                width=8,
                align="center",
            ),
            dbc.Col(
                metadata_field_dropdown(
                    "single-cag-xaxis",
                    label_text="X-axis",
                ) + plot_type_dropdown(
                    "single-cag-plot-type",
                    options=[
                        {'label': 'Points', 'value': 'scatter'},
                        {'label': 'Line', 'value': 'line'},
                        {'label': 'Boxplot', 'value': 'boxplot'},
                        {'label': 'Stripplot', 'value': 'strip'},
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
                width=4,
                align="center",
            )
        ])
    )
#####################
# / SINGLE CAG CARD #
#####################


###############
# GENOME CARD #
###############
def genome_card():
    return card_wrapper(
        "Genome Similarity",
        dbc.Row([
            dbc.Col(
                [
                    dbc.Spinner(dcc.Graph(
                                id="genome-scatter-graph"
                                )),
                    dbc.Spinner(dcc.Graph(
                                id="genome-heatmap-graph"
                                )),
                ],
                width=8,
                align="center",
            ),
            dbc.Col(
                corncob_parameter_dropdown(
                    "genome-parameter-dropdown",
                ) + basic_slider(
                    "genome-scatter-ngenes-slider",
                    "Minimum Size Filter (Num. Genes)",
                    min_value=1,
                    max_value=1000,
                    step_value=10,
                    default_value=100,
                    marks=[1, 500, 1000],
                    included=False
                ),
                width=4,
                align="center",
            )
        ]),
        custom_id="genome-card",
        custom_style={"display": "none"}
    )
#################
# / GENOME CARD #
#################


#################
# MANIFEST CARD #
#################
def manifest_card():
    return html.Div([
        html.Br(),
        dbc.Card([
            dbc.CardHeader("Manifest"),
            dbc.CardBody([
                dbc.Row([
                    dbc.Col([
                        html.Div(  # Hidden div to store the list of rows selected by the user
                            [],
                            id="manifest-rows-selected",
                            style={"display": "none"},
                        ),
                        html.Div(  # Hidden div to store the filtered manifest
                            [],
                            id="manifest-filtered",
                            style={"display": "none"},
                        ),
                    ], width=1),
                    dbc.Col(
                        [
                            dash_table.DataTable(
                                id='manifest-table',
                                columns=[
                                    {"name": "specimen", "id": "specimen"},
                                    {"name": "foo", "id": "foo"},
                                ],
                                data=pd.DataFrame([
                                    {"specimen": "apples", "foo": "bar"},
                                    {"specimen": "carrots", "foo": "baz"},
                                ]).to_dict("records"),
                                row_selectable='multi',
                                style_table={
                                    'minWidth': '100%',
                                },
                                style_header={
                                    "backgroundColor": "rgb(2,21,70)",
                                    "color": "white",
                                    "textAlign": "center",
                                },
                                page_action='native',
                                page_size=20,
                                filter_action='native',
                                sort_action='native',
                                hidden_columns=[],
                                selected_rows=[0, 1],
                                css=[{"selector": ".show-hide",
                                        "rule": "display: none"}],
                            ),
                            html.Br(),
                            html.Label("Show / Hide Columns"),
                            dcc.Dropdown(
                                id="manifest-table-select-columns",
                                options=[{"label": n, "value": n} for n in ["specimen", "foo"]],
                                value=["specimen", "foo"],
                                multi=True,
                            ),
                            html.Br(),
                        ]
                    ),
                    dbc.Col(
                        html.Div(),
                        width=1,
                    ),
                ])
            ])
            ])
        ])
###################
# / MANIFEST CARD #
###################

###############################
# REUSABLE DISPLAY COMPONENTS #
###############################
def card_wrapper(
    card_name,
    card_body,
    custom_id=None,
    custom_style=None
):
    # Make an id for this particular card
    card_id = card_name.lower().replace(" ", "-")
    
    return html.Div([
        html.Br(),
        dbc.Card([
            dbc.CardHeader(
                dbc.Row([
                    dbc.Col(
                        html.Div(
                            card_name,
                            style={"vertical-align": "middle"}
                        ),
                        width=6,
                    ),
                    dbc.Col(
                        html.Div(
                            dbc.Button(
                                "</>",
                                id={
                                    "type": "toggle-collapsable-card",
                                    "parent": card_id
                                },
                                n_clicks=0
                            ),
                            style={"text-align": "right"}
                        ),
                        width=6,
                    )
                ])
            ),
            dbc.CardBody([
                dbc.Collapse(
                    card_body,
                    is_open=True,
                    id={"type": "collapsable-card-body", "parent": card_id}
                )
            ])
        ])
    ],
        id=card_id if custom_id is None else custom_id,
        style=custom_style
    )


def exp_table_row(header1, value1, header2, value2, header_bg="#F4F6F6", value_bg="white", spacer_bg="white"):
    return [     # Table Body
        html.Tr([    # Row
            html.Td(
                header1,
                style={"backgroundColor": header_bg}
            ),
            html.Td(
                value1,
                style={"backgroundColor": value_bg}
            ),
            html.Td(
                "",
                style={"backgroundColor": spacer_bg}
            ),
            html.Td(
                header2,
                style={"backgroundColor": header_bg}
            ),
            html.Td(
                value2,
                style={"backgroundColor": value_bg}
            ),
        ]
        )]


def metadata_field_dropdown(
    dropdown_id,
    label_text='Metadata Label',
    default_value="none",
    multi=False
):
    return [
        html.Label(label_text),
        dcc.Dropdown(
            id={
                "type": "metadata-field-dropdown",
                "name": dropdown_id
            },
            options=[
                {'label': 'None', 'value': 'none'},
            ],
            value=default_value,
            multi=multi
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
        max_value=10,
        default_value=default_value,
        marks=[1, 5, 10],
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


def corncob_parameter_dropdown(
    dropdown_id,
    label_text='Parameter',
):
    return [
        html.Label(label_text),
        dcc.Dropdown(
            id={
                "type": "corncob-parameter-dropdown",
                "name": dropdown_id
            },
            options=[
                {'label': 'None', 'value': 'none'},
            ],
            value="none"
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
            included=False,
        ),
        html.Br()
    ]


def cag_size_slider(slider_id, label_text='CAG Size Filter'):
    return [
        html.Label(label_text),
        dcc.RangeSlider(
            id={
                "type": "cag-size-slider",
                "name": slider_id
            },
            min=0,
            max=3,
            step=0.1,
            marks={
                str(n): str(10**n)
                for n in range(3)
            },
            value=[
                0,
                3
            ]
        ),
        html.Br()
    ]


def cag_metric_slider(slider_id, metric, label_text):
    return [
        html.Label(label_text),
        dcc.RangeSlider(
            id={
                "type": "cag-metric-slider",
                "metric": metric,
                "name": slider_id
            },
            min=0,
            max=1,
            step=0.1,
            marks={
                str(n): str(round(n, 2))
                for n in np.arange(
                    0, 1, 0.2
                )
            },
            value=[0, 1]
        ),
        html.Br()
    ]


def cag_metric_dropdown(slider_id, label_text='Metric', default_value="size"):
    return [
        html.Label(label_text),
        dcc.Dropdown(
            id=slider_id,
            options=[
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
                {'label': 'On   ', 'value': 'on'},
                {'label': 'Off', 'value': 'off'},
            ],
            value=default,
        )
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
            value=20,
            included=False,
        ),
        html.Br()
    ]

#################################
# \ REUSABLE DISPLAY COMPONENTS #
#################################

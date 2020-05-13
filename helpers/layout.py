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
                'Menu',
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
            dbc.CardHeader(dataset["name"]),
            dbc.CardBody(
                dbc.Row([
                    dbc.Col(width=1),
                    dbc.Col(
                        html.Div([
                            html.Br(),
                            dcc.Markdown(
                                dataset.get("description", "")
                            ),
                            html.Br(),
                        ]),
                        width=8,
                    ),
                    dbc.Col(
                        html.Div([
                            html.Br(),
                            dbc.Button(
                                'Open',
                                id={
                                    "type": "open-dataset-button",
                                    "index": ix,
                                },
                                n_clicks=0
                            ),
                            html.Div(
                                id={
                                    "type": "open-dataset-pressed",
                                    "index": ix
                                },
                                style={"display": "none"}
                            ),
                            html.Br(),
                        ]),
                        width=1,
                    ),
                    dbc.Col(width=1),
                ],
                    justify="between"
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
    return html.Div([
        html.Br(),
        dbc.Card([
            dbc.CardHeader("Gene Analysis Summary"),
            dbc.CardBody([
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
                    ),
                        align="center",
                        width=4,
                    )
                ])
            ])
        ])
    ])
###################
# / RICHNESS CARD #
###################


###################
# ORDINATION CARD #
###################
def ordination_card():
    return html.Div([
        html.Br(),
        dbc.Card([
            dbc.CardHeader("Ordination Analysis"),
            dbc.CardBody([
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
                        width=4,
                        align="center"
                    )
                ])
            ]),
        ])
    ])
#####################
# / ORDINATION CARD #
#####################


####################
# CAG SUMMARY CARD #
####################
def cag_summary_card():
    return html.Div([
        html.Br(),
        dbc.Card([
            dbc.CardHeader("CAG Summary"),
            dbc.CardBody([
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
                    dbc.Col(
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
                        ) + [
                            html.Div(id='global-selected-cag',
                                     style={"display": "none"}),
                            html.Div(id='cag-summary-selected-cag',
                                     style={"display": "none"}),

                        ],
                        width=4,
                        align="center"
                    )
                ])
            ]),
        ])
    ])
######################
# / CAG SUMMARY CARD #
######################


################
# VOLCANO PLOT #
################
def volcano_card():
    return html.Div([
        html.Br(),
        dbc.Card([
            dbc.CardHeader("Association Screening"),
            dbc.CardBody([
                dbc.Row([
                    dbc.Col(
                        dbc.Spinner(dcc.Graph(
                            id='volcano-graph'
                        )),
                        width=8,
                        align="center"
                    ),
                    dbc.Col(
                        volcano_parameter_dropdown(
                            "volcano-parameter-dropdown",
                        ) + volcano_pvalue_slider(
                            "volcano-pvalue-slider",
                        ) + log_scale_radio_button(
                            "volcano-fdr-radio",
                            label_text="FDR-BH adjustment"
                        ) + [
                            html.Div(id='volcano-selected-cag', style={"display": "none"}),
                        ],
                        width=4,
                        align="center"
                    )
                ])
            ]),
        ])
    ], id="volcano")
##################
# / VOLCANO PLOT #
##################

#####################
# CAG TAXONOMY CARD #
#####################
def taxonomy_card():
    return html.Div([
        html.Br(),
        dbc.Card([
            dbc.CardHeader("CAG Taxonomy"),
            dbc.CardBody([
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
            ])
        ])
    ])
#######################
# / CAG TAXONOMY CARD #
#######################


###################
# SINGLE CAG CARD #
###################
def single_cag_card():
    return html.Div([
        html.Br(),
        dbc.Card([
            dbc.CardHeader("Individual CAG Abundance"),
            dbc.CardBody([
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
            ])
        ])
    ])
#####################
# / SINGLE CAG CARD #
#####################


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


def volcano_parameter_dropdown(
    dropdown_id,
    label_text='Parameter',
):
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
                {'label': 'On ', 'value': 'on'},
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
            value=20
        ),
        html.Br()
    ]

#################################
# \ REUSABLE DISPLAY COMPONENTS #
#################################

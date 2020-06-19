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
            dbc.CardHeader("Experiment", id="experiment-summary-card-header"),
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
                            'value': 'pct_reads_aligned'},
                    ],
                    value='pct_reads_aligned'
                ),
                html.Br(),
            ] + plot_type_dropdown(
                "richness-type-dropdown"
            ) + metadata_field_dropdown(
                "richness-metadata-dropdown",
                label_text="Color By"
            ) + log_scale_radio_button(
                "richness-log-x",
                label_text="Number of Reads - Log Scale"
            ),
                align="center",
                width=4,
            )
        ]),
        help_text="""
Summary of the detection of genes across samples in this dataset.
For each sample, you may select:

- The proportion of reads from each sample which align uniquely to a single protein-coding gene from the catalog
- The number of genes which were identified by _de novo_ assembly in each sample
- The number of genes which were identified by alignment against the non-redundant gene catalog generated by _de novo_ assembly

Data may be summarized either by a single point per sample or with a frequency histogram.

To mask any sample from this plot, deselect it in the manifest at the bottom of the page.

Note: Click on the camera icon at the top of this plot (or any on this page) to save a PNG to your computer.
        """
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
        ]),
        help_text="""
Summary of the CAGs detected in a single sample. 
You may select any sample to display the relative abundance of every CAG which
was detected, comparing by default against the size (number of genes) in each CAG. 

Optionally, you may choose to compare the abundance of CAGs in a pair of samples.

Abundance metrics:

- Relative Abundance (default): The proportion of gene copies detected in a sample which are assigned to each CAG
- Centered Log-Ratio: The log10 transformed relative abundance of each CAG divided by the geometric mean of relative abundance of all CAGs detected in that sample

Note: Click on the camera icon at the top of this plot (or any on this page) to save a PNG to your computer.
        """
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
        ]),
        help_text="""
Beta-diversity summary of the similarity of community composition across samples.
You may select a distance metric which is used to calculate the similarity of every pair of samples.
Based on that distance matrix, PCA or t-SNE may be used to summarize the groups of samples which
are most similar to each other.

By moving the sliders, you may select different composite indices to display on the plot.

The upper plot shows the frequency histogram of samples across the first composite axis.
The lower scatter plot shows two axes, and metadata from the manifest can be overlaid as a color on each point.

To mask any sample from this plot, deselect it in the manifest at the bottom of the page.

Note: Click on the camera icon at the top of this plot (or any on this page) to save a PNG to your computer.
        """
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
        [
            dbc.Row([
                dbc.Col(
                    cag_metric_dropdown(
                        "cag-summary-metric-primary",
                        default_value="size",
                        label_text="Primary Metric (x-axis)",
                    ) + cag_metric_dropdown(
                        "cag-summary-metric-secondary",
                        default_value="entropy",
                        label_text="Secondary Metric (y-axis)",
                    ),
                    width = 4,
                    align = "center"
                ),
                dbc.Col(
                    cag_size_slider(
                        "cag-summary-size-slider"
                    ) + cag_metric_slider(
                        "cag-summary-entropy-slider",
                        "entropy",
                        "CAG Entropy Filter",
                    ) ,
                    width = 4,
                    align = "center"
                ),
                dbc.Col(
                    cag_metric_slider(
                        "cag-summary-prevalence-slider",
                        "prevalence",
                        "CAG Prevalence Filter"
                    ) + cag_metric_slider(
                        "cag-summary-abundance-slider",
                        "mean_abundance",
                        "CAG Abundance Filter",
                    ),
                    width = 4,
                    align = "center"
                )
            ]),
            dbc.Row([
                dbc.Col(
                    [
                        dbc.Spinner(
                            dcc.Graph(
                                id='cag-summary-graph-scatter'
                            )
                        ),
                    ],
                    width=12,
                    align="center"
                )
            ]),
            dbc.Row([
                dbc.Col(
                    [
                        html.Label("Histogram Display"),
                        dcc.Dropdown(
                            id='cag-summary-histogram-metric',
                            options=[
                                {'label': 'Number of genes', 'value': 'genes'},
                                {'label': 'Number of CAGs', 'value': 'cags'},
                            ],
                            value="genes",
                        )
                    ],
                    width = 4,
                    align = "center"
                ),
                dbc.Col(
                    [
                        html.Label("Histogram Log Scale"),
                        dcc.Dropdown(
                            id='cag-summary-histogram-log',
                            options=[
                                {'label': 'On', 'value': 'on'},
                                {'label': 'Off', 'value': 'off'},
                            ],
                            value="on",
                        ),
                        html.Div(id='global-selected-cag',
                                style={"display": "none"}),
                        html.Div(id='cag-summary-selected-cag',
                                style={"display": "none"}),
                    ],
                    width = 4,
                    align = "center"
                ),
                dbc.Col(
                    nbins_slider(
                        "cag-summary-nbinsx-slider"
                    ),
                    width=4,
                    align="center"
                )
            ]),
            dbc.Row([
                dbc.Col(
                    [
                        dbc.Spinner(
                            dcc.Graph(
                                id='cag-summary-graph-hist'
                            )
                        ),
                    ],
                    width=12,
                    align="center"
                )
            ])
        ],
        help_text="""
Genes were grouped into Co-Abundant Groups (CAGs), and this panel summarizes that set of CAGs on the basis of:

- Size: Number of genes contained in each CAG
- Entropy: The evenness of abundance for a given CAG across all samples
- Mean Abundance across all samples
- Prevalence: The proportion of samples in which a CAG was detected at all
- Standard Deviation of relative abundance values across all samples

*Click on any CAG* to display additional information about that CAG in the panels below.

Note: Masking a sample from the manifest at the bottom of the page does _not_ update the summary CAG metrics displayed here.

You may select to filter which CAGs are displayed using the sliders which are provided.
You may also change the number of bins used to plot the frequency histogram.

By default, the frequency histogram (top) displays the number of _genes_ found in the group of CAGs that fall
within a given range. You may instead choose to display the number of CAGs which fall into each bin.

Note: Click on the camera icon at the top of this plot (or any on this page) to save a PNG to your computer.
        """
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
        [
            dbc.Row([
                dbc.Col(
                    [
                        html.Label("Display CAGs"),
                        dcc.Dropdown(
                            id="cag-heatmap-cag-dropdown",
                            options=[],
                            value=[],
                            multi=True
                        )
                    ],
                    width=4,
                    align="center",
                ),
                dbc.Col(
                    [
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
                        html.Label("Group Specimens"),
                        dcc.Dropdown(
                            id='cag-heatmap-cluster',
                            options=[
                                {'label': 'By Metadata', 'value': 'metadata'},
                                {'label': 'By CAG Abundances', 'value': 'cag'},
                            ],
                            value="metadata",
                        ),
                    ],
                    width=4,
                    align="center",
                ),
                dbc.Col(
                    [
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
                        html.Label("Show Taxonomy"),
                        dcc.Dropdown(
                            id='cag-heatmap-taxa-rank',
                            options=[
                                {'label': 'None', 'value': 'none'},
                                {'label': 'Species', 'value': 'species'},
                                {'label': 'Genus', 'value': 'genus'},
                                {'label': 'Family', 'value': 'family'},
                                {'label': 'Class', 'value': 'class'},
                                {'label': 'Phylum', 'value': 'phylum'},
                            ],
                            value="none",
                        ),
                    ],
                    width=4,
                    align="center"
                )
            ]),
            dbc.Row([
                dbc.Col(
                    dbc.Spinner(dcc.Graph(
                        id='cag-heatmap-graph'
                    )),
                    width=12,
                    align="center"
                )
            ]),
        ],
        help_text="""
The relative abundance of a user-selected group of CAGs is shown in comparison to
specimen metadata as well as the taxonomic annotation of those CAGs.

Type in the box to select CAGs to add to the heatmap. Additionally, clicking on a single CAG in any of
the other displays on this page will add those CAGs to the heatmap.

You may choose to annotate columns by specimen metadata, and you may choose to annotate rows
by the taxonomic annotation of each CAG, when available.
        """
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
        ]),
        help_text="""
The estimated association of each CAG with a specified metadata feature is displayed
as a volcano plot. The values shown in this display must be pre-computed by selecting
the `--formula` flag when running _geneshot_.

Note: Clicking on any CAG in this display will select it for display in the other
panels on this page.
        """
    )
##################
# / VOLCANO PLOT #
##################

###################
# SINGLE CAG CARD #
###################
def single_cag_card():
    return card_wrapper(
        "Individual CAG Abundance",
        [
            dbc.Row([
                dbc.Col(
                    dbc.Spinner(
                        dcc.Graph(id="single-cag-graph")
                    ),
                    width=8,
                    align="center",
                ),
                dbc.Col(
                    [
                        html.Label("Display CAG"),
                        dcc.Dropdown(
                            id="single-cag-dropdown",
                            options=[],
                            value=[],
                        ),
                        html.Br(),
                    ] + metadata_field_dropdown(
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
            ]),
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
            ]),
        ],
        help_text="""
Construct a summary of the abundance of a single CAG in relation to the metadata
assigned to each specimen. By selecting different types of plots, you may flexibly
construct any type of summary display.

The taxonomic annotation of a given CAG is shown as the proportion of
genes which contain a given taxonomic annotation, out of all genes which
were given any taxonomic annotation.

Note: Click on the camera icon at the top of this plot (or any on this page) to save a PNG to your computer.
        """
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
        [
            dbc.Row([
                dbc.Col(
                    corncob_parameter_dropdown(
                        "genome-parameter-dropdown",
                    ),
                    width = 4,
                    align = "center"
                ),
                dbc.Col(
                    basic_slider(
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
            dbc.Row([
                dbc.Col(
                    [
                        dbc.Spinner(
                            dcc.Graph(
                                id="genome-scatter-graph"
                            )
                        )
                    ],
                    width=12,
                    align="center",
                )
            ]),
            dbc.Row([
                dbc.Col(
                    [
                        dbc.Spinner(
                            dcc.Graph(
                                id="genome-heatmap-graph"
                            )
                        ),
                    ],
                    width=12,
                    align="center",
                )
            ])
        ],
        custom_id="genome-card",
        custom_style={"display": "none"},
        help_text="""
Summary of genome alignment against CAGs associated with parameters of interest.

The upper plot summarizes all of the genomes which aligned to the CAGs in this dataset.
Each genome is summarized by (a) the proportion of the genome which aligns to any CAG
which is associated with the indicated parameter (at a fixed p-value threshold), as well as
(b) the average estimated coefficient for those genes which do align against that genome.

*Selecting Genomes*

To select a set of genomes to display in more detail, mouse over the scatter plot and
activate one of the selection tools (box select or lasso select).
Using that selection tool, indicate those genomes which you wish to display more details on.
The lower plot will then be rendered with those set of selected genomes.
Using the controls on the right, you may then filter the CAGs and genomes which are displayed.
        """
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
            ]),
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
    help_text="Missing",
    custom_id=None,
    custom_style=None,
):
    # Must provide help text
    assert help_text is not None, "Must provide help text for every card"

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
                        width=10,
                    ),
                    dbc.Col(
                        html.Div([
                            dbc.Button(
                                "?",
                                id={
                                    "type": "open-help-text",
                                    "parent": card_id
                                },
                                n_clicks=0
                            ),
                            dbc.Modal(
                                [
                                    dbc.ModalHeader(
                                        card_name
                                    ),
                                    dbc.ModalBody(
                                        dcc.Markdown(help_text)
                                    ),
                                    dbc.ModalFooter(
                                        dbc.Button(
                                            "Close",
                                            id={
                                                "type": "close-help-text",
                                                "parent": card_id
                                            },
                                            className = "ml-auto"
                                        )
                                    )
                                ],
                                id = {
                                    "type": "help-text-modal",
                                    "parent": card_id
                                }
                            )
                            ],
                            style={"text-align": "right"}
                        ),
                        width=1,
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
                        width=1,
                    ),
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

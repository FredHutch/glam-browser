#!/usr/bin/env python3
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_table
import pandas as pd
import numpy as np

# NAVBAR AT THE TOP
def navbar_simple():
    return dbc.NavbarSimple(
        brand="GLAM Browser",
        dark=True,
        color="#112345",
        children=[
            dbc.Button(
                'Main Menu',
                id={
                    "type": "open-dataset-button",
                    "index": -1,
                },
                n_clicks=1,
            ),
            html.Div(  # Store the button-press time
                id={
                    "type": "open-dataset-pressed",
                    "index": -1
                },
                style={"display": "none"},
                children=-1
            ),
            html.Div(  # Store which dataset is selected
                id="selected-dataset",
                children=["-1"],
                style={"display": "none"}
            )
        ],
        id="page-title"
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
                                style={"display": "none"},
                                children=-1
                            )],
                        style={"text-align": "right"}
                    ),
                        width=6
                    ),
                ])
            ]),
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
def update_experiment_summary_card(dataset_metrics, description_markdown):
    # Make a table with the basic summary of an experiment
    return html.Div([
        dcc.Markdown(
            description_markdown
        ),
        html.Br(),
        dbc.Table([
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
        ], bordered=True, hover=True, responsive=True),
    ])
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
In order to perform gene-level metagenomic analysis, the first step is to 
estimate the abundance of every microbial gene in every sample.

The analytical steps needed to perform this analysis include:

- Removal of host sequences by subtractive alignment
- _De novo_ assembly of every individual sample
- Deduplication of protein coding sequences across all samples to form a _de novo_ gene catalog
- Alignment of WGS reads from every sample against that gene catalog
- Estimation of the relative abundance of every gene as the proportion of reads which align uniquely (normalized by gene length)

To visualize how genes were quantified in each sample, you may select:

- The proportion of reads from each sample which align uniquely to a single protein-coding gene from the catalog
- The number of genes which were identified by _de novo_ assembly in each sample
- The number of genes which were identified by alignment against the non-redundant gene catalog generated by _de novo_ assembly

Data may be summarized either by a single point per sample or with a frequency histogram.

To mask any sample from this plot, deselect it in the manifest at the bottom of the page.

Note: Click on the camera icon at the top of this plot (or any on this page) to save an image to your computer.
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
**Summary of the CAGs detected in a single sample.**

You may select any sample to display the relative abundance of every CAG which
was detected, comparing by default against the size (number of genes) in each CAG. 

Optionally, you may choose to compare the abundance of CAGs in a pair of samples.

Abundance metrics:

- Relative Abundance (default): The proportion of gene copies detected in a sample which are assigned to each CAG
- Centered Log-Ratio: The log10 transformed relative abundance of each CAG divided by the geometric mean of relative abundance of all CAGs detected in that sample

Note: Click on the camera icon at the top of this plot (or any on this page) to save an image to your computer.
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
**Beta-diversity summary of the similarity of community composition across samples.**

One way to understand the composition of a microbial community is to compare every pair of samples
on the basis of what microbes were detected. This approach is referred to as 'beta-diversity' analysis.

You may select a distance metric which is used to calculate the similarity of every pair of samples.
Based on that distance matrix, PCA or t-SNE may be used to summarize the groups of samples which
are most similar to each other.

By moving the sliders, you may select different composite indices to display on the plot.

The upper plot shows the frequency histogram of samples across the first composite axis.
The lower scatter plot shows two axes, and metadata from the manifest can be overlaid as a color on each point.

To mask any sample from this plot, deselect it in the manifest at the bottom of the page.

Note: Click on the camera icon at the top of this plot (or any on this page) to save an image to your computer.
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
                    dbc.Spinner(
                        dcc.Graph(
                            id='cag-summary-graph-hist'
                        )
                    ),
                    width = 8,
                    align = "center"
                ),
                dbc.Col(
                    cag_metric_dropdown(
                        "cag-summary-metric-primary",
                        default_value="size",
                        label_text="Metric",
                    ) + [
                        html.Label("Histogram Display"),
                        dcc.Dropdown(
                            id='cag-summary-histogram-metric',
                            options=[
                                {'label': 'Number of genes',
                                    'value': 'genes'},
                                {'label': 'Number of CAGs',
                                    'value': 'cags'},
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
                        html.Div(id='global-selected-cag',
                                    style={"display": "none"}),
                        html.Div(id='cag-summary-selected-cag',
                                    style={"display": "none"}),
                        html.Br()
                    ] + nbins_slider(
                        "cag-summary-nbinsx-slider"
                    ),
                    width = 4,
                    align = "center"
                )
            ])
        ],
        help_text="""
A key factor in performing efficient gene-level metagenomic analysis is the grouping of genes by co-abundance.
The term 'co-abundance' is used to describe the degree to which any pair of genes are found at similar relative
abundances in similar samples. The core concept is that any pair of genes which are always found on the same
molecule of DNA are expected to have similar relative abundances as measured by WGS metagenomic sequencing.

In this analysis, genes were grouped into Co-Abundant Groups (CAGs) by average linkage clustering using the
cosine measure of co-abundance across every pair of samples. After constructing CAGs for this dataset, each
of those CAGs can be summarized on the basis of the aggregate abundance of all genes contained in that CAG.
(note: every gene can only belong to a single CAG in this approach).

The biological interpretation of CAGs is that they are expected to correspond to groups of genes which are
consistently found on the same piece of genetic material (chromosome, plasmid, etc.), or that they are found
in organismsm which are highly co-abundant in this dataset.

This panel summarizes that set of CAGs on the basis of:

- Size: Number of genes contained in each CAG
- Mean Abundance across all samples (as the sum of the relative abundance of every gene in that CAG)
- Entropy: The evenness of abundance for a given CAG across all samples
- Prevalence: The proportion of samples in which a CAG was detected at all
- Standard Deviation of relative abundance values across all samples

Note: Masking a sample from the manifest at the bottom of the page does _not_ update the summary CAG metrics displayed here.

By default, this frequency histogram displays the number of _genes_ found in the group of CAGs that fall
within a given range. You may instead choose to display the number of CAGs which fall into each bin.

Note: Click on the camera icon at the top of this plot (or any on this page) to save an image to your computer.
        """
    )
######################
# / CAG SUMMARY CARD #
######################


##############################
# CAG ABUNDANCE HEATMAP CARD #
##############################
def cag_abundance_heatmap_card():
    return card_wrapper(
        "CAG Abundance Heatmap",
        [
            dbc.Row([
                dbc.Col(
                    [
                        html.Label("Display Top CAGs By"),
                        dcc.Dropdown(
                            id={"type": "heatmap-select-cags-by", "parent": "abundance-heatmap"},
                            options=[
                                {"label": "Average Relative Abundance", "value": "abundance"},
                                {"label": "Size (Number of Genes)", "value": "size"},
                            ],
                            value="abundance"
                        ),
                        html.Br()
                    ] + basic_slider(
                        "cag-abundance-heatmap-ncags",
                        "Number of CAGs to Display",
                        min_value=5,
                        max_value=100,
                        default_value=20,
                        marks=[5, 25, 50, 75, 100]
                    ) + cag_size_slider(
                        "cag-abundance-heatmap-size-range"
                    ),
                    width=4,
                    align="center",
                ),
                dbc.Col(
                    [
                        html.Label("Display Metadata"),
                        dcc.Dropdown(
                            id="cag-abundance-heatmap-metadata-dropdown",
                            options=[],
                            value=[],
                            multi=True
                        ),
                        html.Br(),
                        html.Label("Group Specimens"),
                        dcc.Dropdown(
                            id='cag-abundance-heatmap-cluster',
                            options=[
                                {'label': 'By CAG Abundances', 'value': 'cag'},
                                {'label': 'By Metadata', 'value': 'metadata'},
                            ],
                            value="cag",
                        ),
                    ],
                    width=4,
                    align="center",
                ),
                dbc.Col(
                    [
                        html.Label("Abundance Metric"),
                        dcc.Dropdown(
                            id='cag-abundance-heatmap-abundance-metric',
                            options=[
                                {'label': 'Rel. Abund. (log10)', 'value': 'log10'},
                                {'label': 'Rel. Abund. (log10) (z-score)', 'value': 'zscore'},
                                {'label': 'Rel. Abund.', 'value': 'raw'},
                            ],
                            value="log10",
                        ),
                        html.Br(),
                        html.Label("Taxonomic Annotation"),
                        dcc.Dropdown(
                            id='cag-abundance-heatmap-annotate-cags-by',
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
                        id='cag-abundance-heatmap-graph'
                    )),
                    width=12,
                    align="center"
                )
            ]),
        ],
        help_text="""
This display lets you compare the relative abundance of a group of CAGs across all samples.
You may choose to view those CAGs which are most highly abundant, those CAGs containing the
largest number of genes, or those CAGs which are most consistently associated with a parameter
in your formula (if provided).

If you decide to display those CAGs which are most associated with a parameter in the formula,
then you will see the estimated coefficient of association for each CAG against that parameter
displayed to the right of the heatmap.

You may also choose to display the taxonomic annotation of each CAG at a particular taxonomic
level (e.g. species). That will add a color label to each row in the heatmap, and you can see
the name of the organism that represents by moving your mouse over that part of the plot.

The controls at the top of the display help you customize this heatmap. You may choose to include
a greater or smaller number of CAGs; you may choose to filter CAGs based on their size (the
number of genes in each CAG); and you may choose to annotate the samples based on user-defined
metadata from the manifest.

By default, the columns in the heatmap are ordered based on the similarity of CAG abundances
(average linkage clustering), but you may also choose to set the order according to the sorted
metadata for each sample.

Note: Click on the camera icon at the top of this plot (or any on this page) to save an image to your computer.
        """
)
################################
# / CAG ABUNDANCE HEATMAP CARD #
################################


###############################
# CAG ANNOTATION HEATMAP CARD #
###############################
def cag_annotation_heatmap_card():
    return card_wrapper(
        "CAG Annotation Heatmap",
        [
            dbc.Row([
                dbc.Col(
                    [
                        html.Label("Display Top CAGs By"),
                        dcc.Dropdown(
                            id={"type": "heatmap-select-cags-by", "parent": "annotation-heatmap"},
                            options=[
                                {"label": "Average Relative Abundance", "value": "abundance"},
                                {"label": "Size (Number of Genes)", "value": "size"},
                            ],
                            value="abundance"
                        ),
                        html.Br()
                    ],
                    width=4,
                    align="center",
                ),
                dbc.Col(
                    basic_slider(
                        "cag-annotation-heatmap-ncags",
                        "Number of CAGs to Display",
                        min_value=5,
                        max_value=100,
                        default_value=20,
                        marks=[5, 25, 50, 75, 100]
                    ) + cag_size_slider(
                        "cag-annotation-heatmap-size-range"
                    ),
                    width=4,
                    align="center",
                ),
                dbc.Col(
                    [
                        html.Label("Annotation Type"),
                        dcc.Dropdown(
                            id="cag-annotation-heatmap-annotation-type",
                            options=[
                                {'label': 'Functional', 'value': 'eggNOG_desc'},
                                {'label': 'Taxonomic', 'value': 'taxonomic'},
                                {'label': 'Species', 'value': 'species'},
                                {'label': 'Genus', 'value': 'genus'},
                                {'label': 'Family', 'value': 'family'},
                            ],
                            value='taxonomic',
                        ),
                        html.Br(),
                    ] + basic_slider(
                        "cag-annotation-heatmap-nannots",
                        "Max Number of Annotations to Display",
                        min_value=5,
                        max_value=100,
                        default_value=20,
                        marks=[5, 25, 50, 75, 100]
                    ),
                    width=4,
                    align="center",
                ),
            ]),
            dbc.Row([
                dbc.Col(
                    dbc.Spinner(dcc.Graph(
                        id='cag-annotation-heatmap-graph'
                    )),
                    width=12,
                    align="center"
                )
            ]),
        ],
        help_text="""
This display lets you compare the taxonomic or functional annotations of a group of CAGs.

You may choose to view those CAGs which are most highly abundant, those CAGs containing the
largest number of genes, or those CAGs which are most consistently associated with a parameter
in your formula (if provided).

If you decide to display those CAGs which are most associated with a parameter in the formula,
then you will see the estimated coefficient of association for each CAG against that parameter
displayed to the right of the heatmap. In addition, you will see the aggregate association of
each selected annotation against that same parameter from the formula.

The controls at the top of the display help you customize this heatmap. You may choose to include
a greater or smaller number of CAGs; you may choose to filter CAGs based on their size (the
number of genes in each CAG); and you may choose to display either taxonomic or functional annotations.

Note: Click on the camera icon at the top of this plot (or any on this page) to save an image to your computer.
        """
)
#################################
# / CAG ANNOTATION HEATMAP CARD #
#################################


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
                    group="volcano-parameter",
                ) + cag_size_slider(
                    "volcano-cag-size-slider"
                ) + volcano_pvalue_slider(
                    group="volcano-parameter",
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

Note: Click on the camera icon at the top of this plot (or any on this page) to save an image to your computer.
        """
    )
##################
# / VOLCANO PLOT #
##################

###################
# PLOT CAG CARD #
###################
def plot_cag_card():
    return card_wrapper(
        "Plot CAG Abundance",
        [
            dbc.Row([
                dbc.Col(
                    dbc.Spinner(
                        dcc.Graph(id="plot-cag-graph")
                    ),
                    width=8,
                    align="center",
                ),
                dbc.Col(
                    [
                        html.Label("Display CAG(s) By"),
                        dcc.Dropdown(
                            id="plot-cag-selection-type",
                            options=[
                                {"label": "CAG ID", "value": "cag_id"},
                                {"label": "Association & Annotation", "value": "association"},
                            ],
                            value="cag_id",
                        ),
                        html.Br(),
                        html.Div(
                            [
                                html.Label("CAG ID", style={"margin-right": "15px"}),
                                dcc.Input(
                                    id="plot-cag-multiselector",
                                    type="number",
                                    placeholder="<CAG ID>",
                                    debounce=True,
                                    min=0,
                                    max=1, # Will be updated by callbacks
                                    step=1
                                ),
                                html.Br(),
                            ], 
                            id="plot-cag-by-id-div"
                        ),
                        html.Div(
                            corncob_parameter_dropdown(
                                group="plot-cag",
                            ) + [
                                html.Label("Filter by Annotation"),
                                dcc.Dropdown(
                                    id="plot-cag-annotation-multiselector",
                                    options=[],
                                    value=[],
                                    multi=True,
                                    placeholder="None"
                                ),
                                html.Br(),
                                html.Label(
                                    "Number of CAGs", 
                                    style={"margin-right": "15px"}
                                ),
                                dcc.Input(
                                    id="plot-cag-annotation-ncags",
                                    type="number",
                                    placeholder="<NUM CAGs>",
                                    debounce=True,
                                    min=1,
                                    max=1000,
                                    step=1,
                                    value=5,
                                ),
                                html.Br(),
                            ],
                            id="plot-cag-by-association-div"
                        ),
                    ] + metadata_field_dropdown(
                        "plot-cag-xaxis",
                        label_text="X-axis",
                    ) + plot_type_dropdown(
                        "plot-cag-plot-type",
                        options=[
                            {'label': 'Points', 'value': 'scatter'},
                            {'label': 'Line', 'value': 'line'},
                            {'label': 'Boxplot', 'value': 'boxplot'},
                            {'label': 'Stripplot', 'value': 'strip'},
                        ]
                    ) + metadata_field_dropdown(
                        "plot-cag-color",
                        label_text="Color",
                    ) + metadata_field_dropdown(
                        "plot-cag-facet",
                        label_text="Facet",
                    ) + log_scale_radio_button(
                        "plot-cag-log"
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
                        marks=[1, 5, 10],
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

Note: Click on the camera icon at the top of this plot (or any on this page) to save an image to your computer.
        """
    )
#####################
# / PLOT CAG CARD #
#####################


##############################
# ANNOTATION ENRICHMENT CARD #
##############################
def annotation_enrichment_card():
    return card_wrapper(
        "Estimated Coefficients by Annotation",
        [
            dbc.Row([
                dbc.Col(
                    dbc.Spinner(
                        dcc.Graph(id="annotation-enrichment-graph")
                    ),
                    width=8,
                    align="center",
                ),
                dbc.Col(
                    [
                        html.Label("Annotation Group"),
                        dcc.Dropdown(
                            id="annotation-enrichment-type",
                            options=[
                                {'label': 'Functional', 'value': 'eggNOG_desc'},
                                {'label': 'Taxonomic Species', 'value': 'species'},
                                {'label': 'Taxonomic Genus', 'value': 'genus'},
                                {'label': 'Taxonomic Family', 'value': 'family'},
                            ],
                            value='eggNOG_desc',
                        ),
                        html.Br(),
                    ] + corncob_parameter_dropdown(
                        group="annotation-enrichment"
                    ) + basic_slider(
                        "annotation-enrichment-plotn",
                        "Number of Annotations per Plot",
                        min_value=10,
                        max_value=50,
                        default_value=20,
                        marks=[10, 30, 50]
                    ) + [
                        html.Label("Show Positive / Negative"),
                        dcc.Dropdown(
                            id="annotation-enrichment-show-pos-neg",
                            options=[
                                {'label': 'Both', 'value': 'both'},
                                {'label': 'Positive', 'value': 'positive'},
                                {'label': 'Negative', 'value': 'negative'},
                            ],
                            value="both",
                        ),
                        html.Br(),
                        dbc.Button(
                            "Next", 
                            id='annotation-enrichment-button-next', 
                            n_clicks=0,
                            color="light",
                        ),
                        dbc.Button(
                            "Previous", 
                            id='annotation-enrichment-button-previous', 
                            n_clicks=0,
                            color="light",
                        ),
                        dbc.Button(
                            "First", 
                            id='annotation-enrichment-button-first', 
                            n_clicks=0,
                            color="light",
                        ),
                        html.Div(
                            children=[1],
                            id='annotation-enrichment-page-num',
                            style={"display": "none"},
                        ),
                    ],
                    width=4,
                    align="center",
                )
            ])
        ],
        help_text="""
After estimating the coefficient of association for every individual CAG, we are able to
aggregate those estimated coefficients on the basis of the annotations assigned to the genes
found within those CAGs.

The reasoning for this analysis is that a single CAG may be strongly associated with
parameter X, and that CAG may contain a gene which has been taxonomically assigned to
species Y. In that case, we are interested in testing whether _all_ CAGs which contain
any gene which has been taxonomically assigned to the same species Y are estimated to
have an association with parameter X _in aggregate_.

In this analysis we show the estimated coefficient of association for the groups of CAGs
constructed for every unique annotation in the analysis. This may include functional
annotations generated with eggNOG-mapper, as well as taxonomic assignments generated
by alignment against RefSeq.

Note: Click on the camera icon at the top of this plot (or any on this page) to save an image to your computer.
        """
    )
################################
# / ANNOTATION ENRICHMENT CARD #
################################


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
                            dbc.Button("Bulk Select", id="manifest-table-bulk-select-open"),
                            dbc.Modal(
                                [
                                    dbc.ModalHeader("Filter Specimens by Metadata"),
                                    dbc.ModalBody(
                                        [
                                            dcc.Markdown(
"""Enter a formula to filter specimens by the metadata in your table.

To select all samples, simply leave the formula empty.

If the formula cannot be parsed or if no samples pass the filter, then all samples will be selected.
"""
                                            ),
                                            dcc.Input(
                                                id="manifest-table-bulk-select-formula",
                                                placeholder="e.g., day > 1 and participant != 'Jerry'",
                                                type="text",
                                                debounce=True,
                                                size='50',
                                            )
                                        ]
                                    ),
                                    dbc.ModalFooter(
                                        [
                                            dbc.Button("Apply", id="manifest-table-bulk-select-apply"),
                                        ]
                                    ),
                                ],
                                id="manifest-table-bulk-select-modal",
                                centered=True,
                                keyboard=False,
                                backdrop="static"
                            ),
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
                            style={"verticalAlign": "middle"}
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
    label_text='Parameter',
    group="none"
):
    return [
        html.Label(label_text),
        dcc.Dropdown(
            id={
                "type": "corncob-parameter-dropdown",
                "group": group,
            },
            options=[
                {'label': 'None', 'value': 'none'},
            ],
            value="none"
        ),
        html.Br(),
    ]


def volcano_pvalue_slider(label_text='P-Value Filter (-log10)', group="none"):
    """This slider is missing the max and marks, which will be updated by a callback."""
    return [
        html.Label(label_text),
        dcc.Slider(
            id={
                "type": "corncob-pvalue-slider",
                "group": group,
            },
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


def log_scale_radio_button(id_string, default="on", label_text="Log Scale"):
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
            min=10,
            max=100,
            step=1,
            marks={
                str(v): str(v)
                for v in range(10, 101, 10)
            },
            value=50,
            included=False,
        ),
        html.Br()
    ]

#################################
# \ REUSABLE DISPLAY COMPONENTS #
#################################

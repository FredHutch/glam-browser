#!/usr/bin/env python3

from collections import defaultdict
from helpers.io import parse_directory
from helpers.io import hdf5_get_item
from helpers.io import hdf5_taxonomy
from helpers.layout import navbar_simple
from helpers.layout import dataset_summary_card
from helpers.layout import experiment_summary_card
from helpers.layout import update_experiment_summary_card
from helpers.layout import richness_card
from helpers.layout import single_sample_card
from helpers.layout import ordination_card
from helpers.layout import cag_summary_card
from helpers.layout import cag_heatmap_card
from helpers.layout import volcano_card
from helpers.layout import single_cag_card
from helpers.layout import manifest_card
from helpers.plotting import update_richness_graph
from helpers.plotting import run_pca
from helpers.plotting import run_tsne
from helpers.plotting import update_ordination_graph
from helpers.plotting import print_anosim
from helpers.plotting import plot_sample_vs_cag_size
from helpers.plotting import plot_samples_pairwise
from helpers.plotting import draw_cag_summary_graph_hist
from helpers.plotting import draw_cag_summary_graph_scatter
from helpers.plotting import draw_cag_heatmap
from helpers.plotting import draw_volcano_graph
from helpers.plotting import draw_taxonomy_sunburst
from helpers.plotting import draw_single_cag_graph
from helpers.plotting import parse_manifest_json
from helpers.taxonomy import make_cag_tax_df
from flask_caching import Cache
import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State, MATCH, ALL
import json
import numpy as np
import os
import pandas as pd
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from seaborn import color_palette
from time import time

##########
# README #
##########
# The philosophy of the app is as follows. Data will be read in
# from a collection of input files in HDF5 format. By default, all
# of the datasets in a specified folder will be displayed. If a
# manifest.json file is also present, then the 'name' and 'description'
# for each dataset may be read in from that file. The optional
# manifest is a list of dicts, with the 'fp' key used to identify
# which dataset it corresponds to.
# When it comes to organizing the code needed to build the app, all
# of the components can be distinguished as 'layout', 'graphing', and
# 'io'.
# The main app.py will contain the minimum needed to import and
# invoke these elements, but each component will be imported from
# a helper module. Notably, all of the callback decorators must be
# defined in the main app, and the primary layout declaration must
# contain all of the ids described in the callbacks. 


#######################
# INSPECT DATA_FOLDER #
#######################

# Read in data from the folder specified by DATA_FOLDER
data_folder = os.getenv("DATA_FOLDER")
assert os.path.exists(data_folder), "Path does not exist: {}".format(data_folder)

# Parse the contents of the directory
page_data = parse_directory(data_folder)

def parse_fp(selected_dataset):
    """Function to return the file path for whichever dataset is selected."""
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        # No dataset was selected
        return None
    else:
        # Return the path to the indicated HDF5
        return page_data["contents"][selected_dataset[0]]["fp"]

# All of the data in the indicated folder is loaded as a dict
# in the following format:
# {
#   "page_title": "Title for the navbar at the top of the page", 
#   "page_description": "Additional text for landing page",
#   "contents": [
#     {
#       "fp": "path_to_file.hdf5", 
#       "name": "Optional name from manifest.json", 
#       "description": "Optional description (markdown) from manifest.json", 
#     },
#   ] 
# } 

###################
# CSS STYLESHEETS #
###################
external_stylesheets = [
    {
        'href': 'https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/css/bootstrap.min.css',
        'rel': 'stylesheet',
        'integrity': 'sha384-Vkoo8x4CGsO3+Hhxv8T/Q5PaXtkKtu6ug5TOeNV6gBiFeWPGFN9MuhOf23Q9Ifjh',
        'crossorigin': 'anonymous'
    }
]


##################
# SET UP THE APP #
##################

app = dash.Dash(
    __name__, 
    external_stylesheets=external_stylesheets
)
app.title = "GLAM Browser"
app.config.suppress_callback_exceptions = True

####################
# SET UP THE CACHE #
####################
CACHE_CONFIG = {
    # try 'filesystem' if you don't want to setup redis
    'CACHE_TYPE': 'redis',
    'CACHE_REDIS_URL': os.environ.get('REDIS_URL', 'redis://localhost:6379')
}
cache = Cache()
cache.init_app(app.server, config=CACHE_CONFIG)

###################################
# CACHE DATA READ FROM INDEX HDF5 #
###################################
@cache.memoize()
def manifest(fp):
    return hdf5_get_item(
        fp, 
        "/manifest"
    ).set_index(
        "specimen"
    )

@cache.memoize()
def experiment_metrics(fp):
    return hdf5_get_item(
        fp, 
        "/experiment_metrics"
    ).set_index(
        "variable"
    )["value"]

@cache.memoize()
def specimen_metrics(fp):
    return hdf5_get_item(
        fp, 
        "/specimen_metrics"
    ).set_index(
        "specimen"
    )

@cache.memoize()
def analysis_features(fp):
    return hdf5_get_item(
        fp, 
        "/analysis_features"
    )

@cache.memoize()
def cag_annotations(fp):
    return hdf5_get_item(
        fp, 
        "/cag_annotations"
    ).set_index(
        "CAG"
    )

@cache.memoize()
def cag_abundances(fp):
    return hdf5_get_item(
        fp, 
        "/cag_abundances"
    ).set_index(
        "CAG"
    )

@cache.memoize()
def distances(fp, metric):
    # Make sure that this distance metric is included in the analysis
    assert metric in valid_distances(fp)

    return hdf5_get_item(
        fp,
        "/distances/{}".format(metric)
    ).set_index(
        "specimen"
    )

@cache.memoize()
def cag_associations(fp, parameter):
    # Check if any parameter is selected
    if parameter == "none":
        return pd.DataFrame()

    # Make sure that this parameter is valid
    assert parameter in valid_parameters(fp), parameter

    return hdf5_get_item(
        fp,
        "/cag_associations/{}".format(parameter)
    ).set_index(
        "CAG"
    )

@cache.memoize()
def enrichments(fp, parameter, annotation):
    # Make sure that this parameter is valid
    assert parameter in valid_parameters(fp)

    return hdf5_get_item(
        fp,
        "/enrichments/{}/{}".format(parameter, annotation)
    ).set_index(
        "label"
    )

@cache.memoize()
def gene_annotations_by_cag(fp, cag_id):
    return hdf5_get_item(
        fp,
        "/gene_annotations/CAG/{}".format(cag_id)
    )

@cache.memoize()
def gene_annotations_by_parameter(fp, parameter):
    return hdf5_get_item(
        fp,
        "/gene_annotations/parameter/{}".format(parameter)
    )

@cache.memoize()
def taxonomy(fp):
    return hdf5_taxonomy(
        fp
    )

############################################
# CACHE DATA SUMMARIZED FROM LARGER TABLES #
############################################
@cache.memoize()
def cag_summary_describe(fp):
    return cag_annotations(fp).describe()

# Functions to return the pertinent elements of "/analysis_features"
def valid_distances(fp):
    return analysis_features(
        fp
    ).query(
        "group == 'distances'"
    ).query(
        "key == 'metric'"
    )["value"].tolist()

def valid_parameters(fp):
    return analysis_features(
        fp
    ).query(
        "group == 'cag_associations'"
    ).query(
        "key == 'parameters'"
    )["value"].tolist()

def valid_enrichments(fp):
    return analysis_features(
        fp
    ).query(
        "group == 'enrichments'"
    ).query(
        "key == 'annotation'"
    )["value"].tolist()

@cache.memoize()
def cag_taxonomy(cag_id, fp):
    # Skip if there is no taxonomy
    taxonomy_df = taxonomy(fp)
    if taxonomy_df is None:
        return None

    # Read in the taxonomic annotations for this CAG
    cag_df = gene_annotations_by_cag(fp, cag_id)

    # No CAG has been selected (or no annotations are available)
    if cag_df is None:
        return None

    # Format the DataFrame as needed to make a go.Sunburst
    return make_cag_tax_df(cag_df["tax_id"], taxonomy_df)
######################
# PLOTTING FUNCTIONS #
######################
def empty_figure():
    fig = go.Figure()
    fig.update_layout(
        template="simple_white"
    )
    return fig

##########################
# SET UP THE PAGE LAYOUT #
##########################

# The base layout for the app contain:
#   - The Navbar header
#   - A button for the main menu
#   - Hidden div to store the time of button press
#   - Hidden div to store the index of the selected dataset (-1 for none)
#   - Togglable div with a summary of all datasets (summary-display)
#   - Togglable div with the detailed view (detail-display)
app.layout = html.Div(
    children=[
        navbar_simple(page_data),
        html.Div(
            children=[
                dbc.Card(
                    dbc.CardBody([
                        dcc.Markdown(page_data.get("page_description", "")),
                    ])
                )
            ] + [
                dataset_summary_card(
                    ix, dataset
                )
                for ix, dataset in enumerate(page_data["contents"])
            ],
            id="summary-display",
        ),
        html.Div(
            children=[
                # Cards for each of the plotted components
                experiment_summary_card(),
                richness_card(),
                single_sample_card(),
                ordination_card(),
                cag_summary_card(),
                cag_heatmap_card(),
                volcano_card(),
                single_cag_card(),
                manifest_card(),
            ],
            id="detail-display",
            style={"display": "none"}
        )
    ],
    className="container"
)
############################
# / SET UP THE PAGE LAYOUT #
############################


################################
# OPEN DATASET BUTTON CALLBACK #
################################
@app.callback(
    Output(
        {"type": "open-dataset-pressed", "index": MATCH},
        'children'
    ),
    [Input(
        {"type": "open-dataset-button", "index": MATCH},
        'n_clicks'
    )],
)
def open_dataset_button_click(n_clicks):
    return time()
##################################
# / OPEN DATASET BUTTON CALLBACK #
##################################


####################################
# OPEN / CLOSE HELP MODAL CALLBACK #
####################################
@app.callback(
    Output(
        {
            "type": "help-text-modal",
            "parent": MATCH
        },
        "is_open"
    ),
    [
        Input(
            {
                "type": "open-help-text",
                "parent": MATCH
            }, 
            "n_clicks"
        ),
        Input(
            {
                "type": "close-help-text",
                "parent": MATCH
            }, 
            "n_clicks"
        )
    ],
    [State(
        {
            "type": "help-text-modal",
            "parent": MATCH
        },
        "is_open"
    )],
)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open
# / OPEN / CLOSE HELP MODAL CALLBACK #
######################################


##################################
# TOGGLE DISPLAY BUTTON CALLBACK #
##################################
@app.callback(
    Output(
        {"type": "collapsable-card-body", "parent": MATCH},
        'is_open'
    ),
    [Input(
        {"type": "toggle-collapsable-card", "parent": MATCH},
        'n_clicks'
    )],
)
def toggle_card_button_click(n_clicks):
    return n_clicks is None or n_clicks % 2 == 0
####################################
# / TOGGLE DISPLAY BUTTON CALLBACK #
####################################


################################
# OPEN DATASET SWITCH CALLBACK #
################################
@app.callback(
    Output("selected-dataset", 'children'),
    [Input(
        {
            "type": "open-dataset-pressed",
            "index": ALL
        },
        'children'
    )],
    [dash.dependencies.State("selected-dataset", "children")]
)
def open_dataset_switch(button_timestamps, currently_selected):
    # Only switch the dataset if more than 0.1 seconds passed between button press

    # Make a list of the timestamps
    button_timestamps = [
        float(t)
        for t in button_timestamps
    ]

    # Get the rank order
    rank_order = list(range(len(button_timestamps)))
    # Sort by time
    rank_order.sort(key=lambda ix: button_timestamps[ix], reverse=True)

    if button_timestamps[rank_order[0]] - button_timestamps[rank_order[1]] > 0.1:
        return [rank_order[0] - 1]
    else:
        return currently_selected
##################################
# / OPEN DATASET SWITCH CALLBACK #
##################################


####################################
# METADATA FIELD DROPDOWN CALLBACK #
####################################
@app.callback(
    Output({
        "type": "metadata-field-dropdown",
        "name": MATCH
    }, 'options'),
    [
        Input("selected-dataset", "children"),
        Input({
            "type": "metadata-field-dropdown",
            "name": MATCH
        }, "value"),
    ],
)
def metadata_field_dropdown_callback(selected_dataset, dummy_value):
    """Update the metadata field dropdown for the selected dataset."""

    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        # No dataset was selected
        return [{'label': 'None', 'value': 'none'}]
    else:

        return [
            {'label': 'None', 'value': 'none'},
        ] + [
            {
                "label": f,
                "value": f
            }
            for f in manifest(fp).columns.values
        ]
######################################
# / METADATA FIELD DROPDOWN CALLBACK #
######################################


##############################
# CAG METRIC SLIDER CALLBACK #
##############################
@app.callback(
    [
        Output({ "type": "cag-metric-slider", "name": MATCH, "metric": MATCH}, value)
        for value in ["max", "min", "step", "marks", "value"]
    ],
    [
        Input("selected-dataset", "children"),
        Input({
            "type": "cag-metric-slider",
            "name": MATCH,
            "metric": MATCH,
        }, "id"),
    ],
)
def cag_metric_slider_callback_max(selected_dataset, slider_id):
    """Update any CAG metric slider for the selected dataset."""
    metric = slider_id["metric"]

    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        min_val = 0
        max_val = 1
    else:
        # Get the description of the summary metrics
        df = cag_summary_describe(fp)

        min_val = df.loc["min", metric]
        max_val = df.loc["max", metric]

    step = (max_val - min_val) / 20.
    marks = {
        str(n): str(round(n, 2))
        for n in np.arange(
            min_val, max_val, step * 5
        )
    }
    value = [min_val, max_val]
        
    return [max_val, min_val, step, marks, value]
################################
# / CAG METRIC SLIDER CALLBACK #
################################


############################
# CAG SIZE SLIDER CALLBACK #
############################
@app.callback(
    [
        Output({"type": "cag-size-slider","name": MATCH}, value)
        for value in ["max", "min", "step", "marks", "value"]
    ],
    [
        Input("selected-dataset", "children"),
        Input({
            "type": "cag-size-slider",
            "name": MATCH
        }, "id"),
    ],
)
def cag_size_slider_callback(selected_dataset, dummy_value):
    """Update the CAG size slider for the selected dataset."""

    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        max_value = 3
        min_value = 0
    else:
        max_value = cag_summary_describe(fp).loc["max", "size_log10"]
        min_value = cag_summary_describe(fp).loc["min", "size_log10"]

    step_value = (max_value - min_value) / 20.
    marks = {
        str(int(n)): str(10**int(n))
        for n in np.arange(min_value, max_value, 1.)
    }
    value=[np.log10(3), max_value]

    return [max_value, min_value, step_value, marks, value]
##############################
# / CAG SIZE SLIDER CALLBACK #
##############################


###############################
# SHOW / HIDE SUMMARY DISPLAY #
###############################
@app.callback(
    Output("summary-display", 'style'),
    [Input("selected-dataset", "children")],
)
def show_hide_summary_display(selected_dataset):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return {"display": "block"}
    else:
        return {"display": "none"}
#################################
# / SHOW / HIDE SUMMARY DISPLAY #
#################################


##############################
# SHOW / HIDE DETAIL DISPLAY #
##############################
        ##############################        
##############################
@app.callback(
    Output("detail-display", 'style'),
    [Input("selected-dataset", "children")],
)
def show_hide_detail_display(selected_dataset):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return {"display": "none"}
    else:
        return {"display": "block"}
################################
# / SHOW / HIDE DETAIL DISPLAY #
################################


####################################
# EXPERIMENT SUMMARY CARD CALLBACK #
####################################
@app.callback(
    Output("experiment-summary-card", 'children'),
    [Input("selected-dataset", "children")],
)
def experiment_summary_card_callback(selected_dataset):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return
    else:
        return update_experiment_summary_card(
            experiment_metrics(fp)
        )
@app.callback(
    Output("experiment-summary-card-header", 'children'),
    [Input("selected-dataset", "children")],
)
def experiment_summary_card_header_callback(selected_dataset):
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        return "Experiment"
    else:
        # Return the name of the indicated HDF5
        return "Experiment: {}".format(page_data["contents"][selected_dataset[0]]["name"])
######################################
# / EXPERIMENT SUMMARY CARD CALLBACK #
######################################


##########################
# RICHNESS CARD CALLBACK #
##########################
@app.callback(
    Output("richness-graph", 'figure'),
    [
        Input('richness-metric-dropdown', 'value'),
        Input('richness-type-dropdown', 'value'),
        Input({
            "type": "metadata-field-dropdown",
            "name": 'richness-metadata-dropdown'
        }, 'value'),
        Input('richness-log-x', 'value'),
        Input('manifest-filtered', 'children'),
    ],
    [State("selected-dataset", "children")],
)
def richness_graph_callback(
    selected_metric, 
    selected_type, 
    selected_metadata,
    log_x,
    manifest_json,
    selected_dataset, 
):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return empty_figure()
    else:
        return update_richness_graph(
            specimen_metrics(fp),
            selected_metric,
            selected_type,
            selected_metadata,
            log_x,
            manifest_json,
            manifest(fp)
        )
############################
# / RICHNESS CARD CALLBACK #
############################


###############################
# SINGLE SAMPLE CARD CALLBACK #
###############################
@app.callback(
    Output("single-sample-graph", 'figure'),
    [
        Input('single-sample-primary-dropdown', 'value'),
        Input('single-sample-secondary-dropdown', 'value'),
        Input('single-sample-metric', 'value'),
    ],
    [State("selected-dataset", "children")],
)
def single_sample_graph_callback(
    primary_sample, 
    secondary_sample, 
    display_metric,
    selected_dataset, 
):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return empty_figure()
    elif primary_sample == "none" or display_metric is None:
        return empty_figure()
    else:

        if secondary_sample == "cag_size":
            return plot_sample_vs_cag_size(
                cag_abundances(fp)[primary_sample],
                primary_sample,
                cag_annotations(fp),
                display_metric,
            )
        else:
            return plot_samples_pairwise(
                cag_abundances(fp)[primary_sample],
                primary_sample,
                cag_abundances(fp)[secondary_sample],
                secondary_sample,
                display_metric,
                cag_annotations(fp),
            )

@app.callback(
    [
        Output("single-sample-primary-dropdown", value)
        for value in ["options", "value"]
    ],
    [
        Input("selected-dataset", "children"),
    ],
)
def update_single_sample_primary_dropdown(
    selected_dataset, 
):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    options = [
        {'label': 'None', 'value': 'none'}
    ]

    if fp is None:
        return [options, "none"]
    else:

        options = [
            {
                "label": n,
                "value": n
            }
            for n in manifest(fp).index.values
        ]
        return [options, manifest(fp).index.values[0]]

@app.callback(
    [
        Output("single-sample-secondary-dropdown", value)
        for value in ["options", "value"]
    ],
    [
        Input("selected-dataset", "children"),
    ],
)
def update_single_sample_secondary_dropdown(
    selected_dataset, 
):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    options = [
        {'label': 'CAG Size', 'value': 'cag_size'}
    ]
    value = "cag_size"

    if fp is None:
        return [options, value]
    else:

        options.extend([
            {
                "label": n,
                "value": n
            }
            for n in manifest(fp).index.values
        ])
        return [options, value]
#################################
# / SINGLE SAMPLE CARD CALLBACK #
#################################


#############################
# ORDINATION CARD CALLBACKS #
#############################
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
        Input({'type': 'metadata-field-dropdown', 'name': 'ordination-metadata'}, 'value'),
        Input('manifest-filtered', 'children'),
    ],
    [
        State("selected-dataset", "children")
    ])
def ordination_graph_callback(
    algorithm,
    metric,
    primary_pc,
    secondary_pc,
    perplexity,
    metadata,
    manifest_json,
    selected_dataset,
):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return empty_figure()
    else:
        return update_ordination_graph(
            distances(fp, metric),
            algorithm,
            primary_pc,
            secondary_pc,
            perplexity,
            metadata,
            manifest_json,
            manifest(fp),
        )
@app.callback(
    Output('ordination-anosim-results', 'children'),
    [
        Input('ordination-metric', 'value'),
        Input({'type': 'metadata-field-dropdown', 'name': 'ordination-metadata'}, 'value'),
        Input('manifest-filtered', 'children'),
        Input("selected-dataset", "children")
    ])
def ordination_anosim_callback(
    metric,
    metadata,
    manifest_json,
    selected_dataset,
):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return dcc.Markdown("")
    elif metadata == "none":
        return dcc.Markdown("")
    else:

        return print_anosim(
            distances(fp, metric),
            metadata,
            manifest_json,
            manifest(fp),
        )

###############################
# / ORDINATION CARD CALLBACKS #
###############################


##############################
# CAG SUMMARY CARD CALLBACKS #
##############################
@app.callback(
    Output('cag-summary-graph-hist', 'figure'),
    [
        Input("selected-dataset", "children"),
        Input('cag-summary-metric-primary', 'value'),
        Input({"name": 'cag-summary-size-slider', "type": "cag-size-slider"}, 'value'),
        Input({"name": 'cag-summary-entropy-slider', "type": "cag-metric-slider", "metric": "entropy"}, 'value'),
        Input({"name": 'cag-summary-prevalence-slider', "type": "cag-metric-slider", "metric": "prevalence"}, 'value'),
        Input({"name": 'cag-summary-abundance-slider', "type": "cag-metric-slider", "metric": "mean_abundance"}, 'value'),
        Input('cag-summary-nbinsx-slider', 'value'),
        Input('cag-summary-histogram-log', 'value'),
        Input('cag-summary-histogram-metric', 'value'),
    ])
def cag_summary_graph_hist_callback(
    selected_dataset,
    metric_primary,
    size_range,
    entropy_range,
    prevalence_range,
    abundance_range,
    nbinsx,
    log_scale,
    hist_metric,
):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return empty_figure()
    else:
        return draw_cag_summary_graph_hist(
            cag_annotations(fp),
            metric_primary,
            size_range,
            entropy_range,
            prevalence_range,
            abundance_range,
            nbinsx,
            log_scale,
            hist_metric,
        )

@app.callback(
    Output('cag-summary-graph-scatter', 'figure'),
    [
        Input("selected-dataset", "children"),
        Input('cag-summary-metric-primary', 'value'),
        Input('cag-summary-metric-secondary', 'value'),
        Input({"name": 'cag-summary-size-slider', "type": "cag-size-slider"}, 'value'),
        Input({"name": 'cag-summary-entropy-slider', "type": "cag-metric-slider", "metric": "entropy"}, 'value'),
        Input({"name": 'cag-summary-prevalence-slider', "type": "cag-metric-slider", "metric": "prevalence"}, 'value'),
        Input({"name": 'cag-summary-abundance-slider', "type": "cag-metric-slider", "metric": "mean_abundance"}, 'value'),
        Input('global-selected-cag', 'children'),
    ])
def cag_summary_graph_scatter_callback(
    selected_dataset,
    metric_primary,
    metric_secondary,
    size_range,
    entropy_range,
    prevalence_range,
    abundance_range,
    selected_cag_json,
):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return empty_figure()
    else:
        return draw_cag_summary_graph_scatter(
            cag_annotations(fp),
            metric_primary,
            metric_secondary,
            size_range,
            entropy_range,
            prevalence_range,
            abundance_range,
            selected_cag_json,
        )

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
################################
# / CAG SUMMARY CARD CALLBACKS #
################################


#########################
# CAG HEATMAP CALLBACKS #
#########################
@app.callback(
    Output('cag-heatmap-graph', 'figure'),
    [
        Input('cag-heatmap-multiselector', 'value'),
        Input('cag-heatmap-metadata-dropdown', 'value'),
        Input('cag-heatmap-abundance-metric', 'value'),
        Input('cag-heatmap-cluster', 'value'),
        Input('cag-heatmap-taxa-rank', 'value'),
        Input('manifest-filtered', 'children'),
    ],
    [State("selected-dataset", "children")])
def heatmap_graph_callback(
    cags_selected,
    metadata_selected,
    abundance_metric,
    cluster_by,
    taxa_rank,
    manifest_json,
    selected_dataset,
):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return empty_figure()
    
    if len(cags_selected) == 0:
        return empty_figure()

    # Get the abundance of the selected CAGs
    cag_abund_df = cag_abundances(fp).reindex(
        index=cags_selected
    ).T

    # Get the taxonomic hits for the selected CAGs
    if taxa_rank != "none":
        cag_tax_dict = {
            cag_id: cag_taxonomy(cag_id, fp)
            for cag_id in cags_selected
        }
    else:
        cag_tax_dict = {}

    # Draw the figure
    return draw_cag_heatmap(
        cag_abund_df,
        metadata_selected,
        abundance_metric,
        cluster_by,
        taxa_rank,
        manifest_json,
        manifest(fp),
        cag_tax_dict,
    )

@app.callback(
    [
        Output('cag-heatmap-metadata-dropdown', value)
        for value in ['options', 'value']
    ],
    [
        Input("selected-dataset", "children"),
    ])
def update_heatmap_metadata_dropdown(
    selected_dataset,
):
    """When a new dataset is selected, fill in the available metadata."""
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return [], []

    # Get the list of all metadata
    options = [
        {
            "label": f,
            "value": f
        }
        for f in manifest(fp).columns.values
    ]

    return options, []

def get_cag_multiselector_options(
    selected_dataset
):
    """Function to make a list of all of the CAGs for a given dataset."""
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        options = [{
            "label": "None",
            "value": 0
        }]

    else:

        # Get the list of all CAGs
        cag_id_list = cag_annotations(fp).index.values

        # Get the list of CAGs
        options = [
            {
                "label": "CAG {}".format(cag_id),
                "value": cag_id
            }
            for cag_id in cag_id_list
        ]

    return options

@app.callback(
    Output("cag-heatmap-multiselector", 'options'),
    [
        Input("selected-dataset", "children"),
    ])
def update_cag_heatmap_multiselector_options(
    selected_dataset
):
    """When a new dataset is selected, fill in the names of all the CAGs as options."""
    return get_cag_multiselector_options(selected_dataset)

@app.callback(
    [
        Output("cag-heatmap-selected-dataset", "children"),
        Output('cag-heatmap-multiselector', 'value'),
    ],
    [
        Input("selected-dataset", "children"),
        Input('global-selected-cag', 'children'),
    ],
    [
        State("cag-heatmap-selected-dataset", "children"),
        State("cag-heatmap-multiselector", "value"),
    ])
def update_heatmap_cag_dropdown_value(
    selected_dataset,
    clicked_cag_json,
    cag_heatmap_selected_dataset,
    cags_in_dropdown,
):
    # The logic for this callback is a bit involved
    # If a new dataset has been selected, we want to clear the selected values
    # However, if a new CAG has been clicked (global-selected-cag), then add that to the list

    # At the same time, we need to update cag-heatmap-selected-dataset to figure out
    # whether selected-dataset has changed or not

    # If a new dataset is selected, remove all selected values
    if isinstance(selected_dataset, list):
        selected_dataset = selected_dataset[0]
    selected_dataset = int(selected_dataset)

    if isinstance(cag_heatmap_selected_dataset, list):
        cag_heatmap_selected_dataset = cag_heatmap_selected_dataset[0]
    cag_heatmap_selected_dataset = int(cag_heatmap_selected_dataset)

    # A new dataset has been selected
    if cag_heatmap_selected_dataset != selected_dataset:
        
        # Get the path to the selected dataset
        fp = parse_fp([selected_dataset])

        # Pick the top five CAGs to display by mean abundance
        top_five_cags = cag_annotations(fp)["mean_abundance"].sort_values(
            ascending=False
        ).index.values[:5]

        # With a new dataset, select the first five CAGs
        return [selected_dataset], top_five_cags

    else:
        # Reaching this point in the function, a new dataset has _not_ been selected
        # That means that a new CAG has been clicked on somewhere in the display
        if clicked_cag_json is not None:
            # Add the clicked point to the list of selected CAGs in the heatmap
            cags_in_dropdown.append(
                json.loads(clicked_cag_json)["id"]
            )
        return [selected_dataset], cags_in_dropdown


###########################
# / CAG HEATMAP CALLBACKS #
###########################


#####################
# VOLCANO CALLBACKS #
#####################
@app.callback(
    Output('volcano-graph', 'figure'),
    [
        Input("selected-dataset", "children"),
        Input({"type": "corncob-parameter-dropdown", "name": 'volcano-parameter-dropdown'}, 'value'),
        Input("corncob-comparison-parameter-dropdown", "value"),
        Input({"name": 'volcano-cag-size-slider', "type": "cag-size-slider"}, 'value'),
        Input('volcano-pvalue-slider', 'value'),
        Input('volcano-fdr-radio', 'value'),
    ])
def volcano_graph_callback(
    selected_dataset,
    parameter, 
    comparison_parameter,
    cag_size_range,
    neg_log_pvalue_min, 
    fdr_on_off, 
):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    # No dataset has been selected or no parameter has been selected
    if fp is None or parameter == "none":
        return empty_figure()
    else:
        return draw_volcano_graph(
            cag_associations(fp, parameter),
            cag_annotations(fp),
            parameter, 
            comparison_parameter,
            cag_size_range,
            neg_log_pvalue_min, 
            fdr_on_off, 
        )

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
    Output({"type": "corncob-parameter-dropdown","name": MATCH}, "options"),
    [
        Input("selected-dataset", "children"),
    ],
    [State({"type": "corncob-parameter-dropdown","name": MATCH}, "value")])
def update_volcano_parameter_dropdown_options(selected_dataset, dummy):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return [
            {'label': 'None', 'value': 'none'}
        ]
    else:

        # Check to see if the file contains corncob results
        if len(valid_parameters(fp)) == 0:
            return [
                {'label': 'None', 'value': 'none'}
            ]

        else:
            return [
                {'label': l, 'value': l}
                for l in valid_parameters(fp)
            ]

@app.callback(
    [Output("corncob-comparison-parameter-dropdown", value) for value in ["options", "value"]],
    [Input("selected-dataset", "children")]
)
def update_volcano_comparison_dropdown(selected_dataset):
    options = [{'label': 'Estimated Coefficient', 'value': 'coef'}]
    value = "coef"

    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return [options, value]
    else:

        # Get the list of parameters
        parameter_list = valid_parameters(fp)

        # Check to make sure that the file contains corncob results
        if len(parameter_list) > 0:

            options = options + [
                {'label': l, 'value': l}
                for l in parameter_list
            ]

        return [options, value]

@app.callback(
    Output({"type": "corncob-parameter-dropdown","name": MATCH}, "value"),
    [
        Input("selected-dataset", "children"),
    ],
    [State({"type": "corncob-parameter-dropdown", "name": MATCH}, "options")])
def update_volcano_parameter_dropdown_value(selected_dataset, dummy):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return "none"
    else:
        # Get the list of parameters
        parameter_list = valid_parameters(fp)

        if len(parameter_list) == 0:
            return "none"

        if len(parameter_list) > 1:
            return parameter_list[1]
        else:
            return parameter_list[0]

@app.callback(
    [
        Output('volcano-pvalue-slider', value)
        for value in ['max', 'marks']
    ],
    [
        Input("selected-dataset", "children"),
        Input({"type": "corncob-parameter-dropdown", "name": 'volcano-parameter-dropdown'}, 'value'),
    ])
def update_volcano_pvalue_slider(selected_dataset, parameter):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    # No dataset selected, or no parameter selected
    if fp is None or parameter == "none":
        max_value = 1
    else:

        df = cag_associations(
            fp, 
            parameter
        )
        assert "neg_log_pvalue" in df.columns.values, df.columns.values
        max_value = cag_associations(
            fp, 
            parameter
        )["neg_log_pvalue"].max()

    if np.isnan(max_value):
        max_value = 1

    marks = {
        str(int(n)): str(int(n))
        for n in np.arange(
            0,
            max_value,
            max(1, max_value / 5)
        )
    }

    return [max_value, marks]
#######################
# / VOLCANO CALLBACKS #
#######################



#########################
# CAG TAXONOMY CALLBACK #
#########################
@app.callback(
    [
        Output('cag-tax-graph', 'figure'),
        Output('cag-tax-ngenes', 'max'),
        Output('cag-tax-ngenes', 'marks'),
    ],
    [
        Input("selected-dataset", "children"),
        Input('cag-tax-ngenes', 'value'),
        Input('single-cag-multiselector', 'value'),
    ])
def update_taxonomy_graph(selected_dataset, min_ngenes, cag_id):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        marks = {
            n: n
            for n in ["0", "1"]
        }
        return empty_figure(), 1, marks

    # Format the DataFrame as needed to make a go.Sunburst
    cag_tax_df = cag_taxonomy(int(cag_id), fp)

    return draw_taxonomy_sunburst(cag_tax_df, int(cag_id), min_ngenes)
###########################
# / CAG TAXONOMY CALLBACK #
###########################

#######################
# SINGLE CAG CALLBACK #
#######################
@app.callback(
    Output("single-cag-multiselector", 'options'),
    [
        Input("selected-dataset", "children"),
    ])
def update_single_cag_multiselector_options(
    selected_dataset
):
    """When a new dataset is selected, fill in the names of all the CAGs as options."""
    return get_cag_multiselector_options(selected_dataset)
@app.callback(
    Output('single-cag-multiselector', 'value'),
    [
        Input("selected-dataset", "children"),
    ],
    [
        State('single-cag-multiselector', 'value')
    ]
)
def update_single_cag_dropdown_value(
    selected_dataset, selected_cag
):
    # If a new dataset is selected, remove all selected values
    if isinstance(selected_dataset, list):
        selected_dataset = selected_dataset[0]
    selected_dataset = int(selected_dataset)

    # If the Main Menu button has been clicked, just return the existing value
    if selected_dataset == -1:
        return selected_cag

    # Otherwise return the most abundant CAG
    else:
        
        # Get the path to the selected dataset
        fp = parse_fp([selected_dataset])

        # Pick the top CAG to display by mean abundance
        top_cag = cag_annotations(
            fp
        )[
            "mean_abundance"
        ].sort_values(
            ascending=False
        ).index.values[0]

        # With a new dataset, select the first five CAGs
        return top_cag
@app.callback(
    Output('single-cag-graph', 'figure'),
    [
        Input('single-cag-multiselector', 'value'),
        Input({'name': 'single-cag-xaxis',
               "type": "metadata-field-dropdown"}, 'value'),
        Input('single-cag-plot-type', 'value'),
        Input({'name': 'single-cag-color',
               "type": "metadata-field-dropdown"}, 'value'),
        Input({'name': 'single-cag-facet',
               "type": "metadata-field-dropdown"}, 'value'),
        Input('single-cag-log', 'value'),
        Input('manifest-filtered', 'children'),
    ],[
        State("selected-dataset", "children")
    ])
def update_single_cag_graph(
    cag_id,
    xaxis,
    plot_type,
    color,
    facet,
    log_scale,
    manifest_json,
    selected_dataset,
):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return empty_figure()

    # Get the filtered manifest from the browser
    plot_manifest_df = parse_manifest_json(
        manifest_json, 
        manifest(fp)
    )

    plot_df = plot_manifest_df.assign(
        CAG_ABUND = cag_abundances(fp).loc[int(cag_id)]
    )

    return draw_single_cag_graph(
        plot_df,
        int(cag_id),
        xaxis,
        plot_type,
        color,
        facet,
        log_scale
    )
#########################
# / SINGLE CAG CALLBACK #
#########################


#########################
# SELECTED CAG CALLBACK #
#########################
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
###########################
# / SELECTED CAG CALLBACK #
###########################


######################
# MANIFEST CALLBACKS #
######################
@app.callback(
    [
        Output('manifest-table', value)
        for value in ["columns", "data", "selected_rows"]
    ] + [
        Output('manifest-table-select-columns', value)
        for value in ["options", "value"]
    ],
    [
        Input("selected-dataset", "children"),
    ])
def update_manifest_table(selected_dataset):
    """Fill in the values of the manifest with the selected dataset."""

    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        # Deselect and return empty values
        columns = [
            {"name": "specimen", "id": "specimen"},
            {"name": "foo", "id": "foo"},
        ]
        data = pd.DataFrame([
            {"specimen": "apples", "foo": "bar"},
            {"specimen": "carrots", "foo": "baz"},
        ]).to_dict("records")
        selected_rows = [0, 1]
        options = [{"label": n, "value": n} for n in ["specimen", "foo"]]
        value = ["specimen", "foo"]

        return columns, data, selected_rows, options, value

    # Otherwise, fill in with real values for this dataset
    manifest_df = manifest(fp).reset_index()

    columns=[
        {"name": i, "id": i}
        for i in manifest_df.columns.values
    ]

    data = manifest_df.to_dict('records')

    selected_rows = np.arange(0, manifest_df.shape[0])

    options = [
        {
            "label": n,
            "value": n
        }
        for n in manifest_df.columns.values
    ]

    value = [
        i for i in manifest_df.columns.values
    ]

    return columns, data, selected_rows, options, value

    
@app.callback(
    Output('manifest-rows-selected', 'children'),
    [
        Input('manifest-table', 'selected_rows'),
    ])
def manifest_save_rows_selected(selected_rows):
    """Save the list of selected rows to a hidden div."""
    return json.dumps(selected_rows, indent=2)

@app.callback(
    Output('manifest-filtered', 'children'),
    [
        Input("selected-dataset", "children"),
        Input('manifest-rows-selected', 'children'),
    ])
def manifest_save_filtered(selected_dataset, selected_rows_json):
    """Save the filtered manifest to a hidden div.."""
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return None

    # Parse the list of indexes which we should save
    selected_rows = json.loads(selected_rows_json)

    # Filter down the manifest table and save as JSON
    return manifest(
        fp
    ).reset_index(
    ).reindex(
        index=selected_rows
    ).to_json(
        date_format='iso',
        orient='records'
    )

@app.callback(
    Output('manifest-table', 'hidden_columns'),
    [
        Input("selected-dataset", "children"),
        Input('manifest-table-select-columns', 'value'),
    ])
def manifest_update_columns_selected(selected_dataset, selected_columns):
    """Allow the user to hide selected columns in the manifest table."""
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return []

    return [
        n
        for n in manifest(fp).columns.values
        if n not in selected_columns
    ]
########################
# / MANIFEST CALLBACKS #
########################

# Used for gunicorn execution
server = app.server

if __name__ == '__main__':

    app.run_server(
        host='0.0.0.0',
        port=8050,
        debug=True,
    )

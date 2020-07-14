#!/usr/bin/env python3

from collections import defaultdict
from helpers.io import parse_directory
from helpers.io import hdf5_get_item
from helpers.layout import navbar_simple
from helpers.layout import dataset_summary_card
from helpers.layout import experiment_summary_card
from helpers.layout import update_experiment_summary_card
from helpers.layout import richness_card
from helpers.layout import single_sample_card
from helpers.layout import ordination_card
from helpers.layout import cag_summary_card
from helpers.layout import cag_abundance_heatmap_card
from helpers.layout import cag_annotation_heatmap_card
from helpers.layout import volcano_card
from helpers.layout import annotation_enrichment_card
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
from helpers.plotting import draw_cag_abundance_heatmap
from helpers.plotting import draw_cag_annotation_heatmap
from helpers.plotting import draw_volcano_graph
from helpers.plotting import draw_enrichment_graph
from helpers.plotting import draw_taxonomy_sunburst
from helpers.plotting import draw_single_cag_graph
from helpers.plotting import parse_manifest_json
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
    if parameter not in valid_parameters(fp):
        return None

    df = hdf5_get_item(
        fp,
        "/enrichments/{}/{}".format(parameter, annotation)
    )
    if df is not None:
        df.set_index("label", inplace=True)
    return df

@cache.memoize()
def functional_gene_annotations(fp):
    return hdf5_get_item(
        fp,
        "/gene_annotations/functional"
    )

@cache.memoize()
def taxonomic_gene_annotations(fp, rank="all", cag_id=None):
    df = hdf5_get_item(
        fp,
        "/gene_annotations/taxonomic/{}".format(rank)
    )

    if df is None:
        return

    df = df.apply(
        lambda c: c.fillna(0).apply(int) if c.name in ["count", "tax_id", "parent", "total", "CAG"] else c
    )
    
    if cag_id is not None:
        df = df.query(
            "CAG == {}".format(cag_id)
        )

        if df.shape[0] == 0:
            return

        # Sort by number of hits
        df = df.sort_values(
            by="count", 
            ascending=False
        )

    return df
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
        "key == 'parameter'"
    )["value"].tolist()

def valid_enrichments(fp):
    return analysis_features(
        fp
    ).query(
        "group == 'enrichments'"
    ).query(
        "key == 'annotation'"
    )["value"].tolist()

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
                volcano_card(),
                cag_abundance_heatmap_card(),
                cag_annotation_heatmap_card(),
                annotation_enrichment_card(),
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
        Input('cag-summary-nbinsx-slider', 'value'),
        Input('cag-summary-histogram-log', 'value'),
        Input('cag-summary-histogram-metric', 'value'),
    ])
def cag_summary_graph_hist_callback(
    selected_dataset,
    metric_primary,
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
            nbinsx,
            log_scale,
            hist_metric,
        )

################################
# / CAG SUMMARY CARD CALLBACKS #
################################


###################################
# CAG ABUNDANCE HEATMAP CALLBACKS #
###################################
@app.callback(
    [
        Output({"type": "heatmap-select-cags-by", "parent": MATCH}, 'options'),
        Output({"type": "heatmap-select-cags-by", "parent": MATCH}, 'value'),
    ],
    [
        Input("selected-dataset", "children"),
    ],
    [State({"type": "heatmap-select-cags-by", "parent": MATCH}, 'value')])
def abundance_heatmap_graph_select_cags_callback(selected_dataset, _):
    """Add the corncob parameters to the CAG abundance heatmap CAG selection options."""

    # Set the base options
    options = [
        {"label": "Average Relative Abundance", "value": "abundance"},
        {"label": "Size (Number of Genes)", "value": "size"}
    ]

    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        # Return the base options if no dataset is selected
        return options, "abundance"
    else:
        # TODO Add the default parameter from the manifest, if specified
        # Add the parameters
        return options + [
            {"label": parameter, "value": "parameter-{}".format(parameter)}
            for parameter in valid_parameters(fp)
        ], "abundance"


def get_cags_selected_by_criterion(fp, select_cags_by, n_cags, cag_size_range):
    """Function to select the top N CAGs based on a particular criterion."""
    # Either select based on a parameter, or based on some CAG annotation metric

    if select_cags_by.startswith("parameter-"):
        # Parse the parameter name from the value
        parameter = select_cags_by.replace("parameter-", "")
        assert parameter in valid_parameters(fp), \
            "Selecting CAGs by {}, but no parameter found".format(select_cags_by)

        # Get the results for this particular parameter
        ranked_df = cag_associations(
            fp, 
            parameter
        ).assign(
            size_log10 = cag_annotations(fp)["size_log10"]
        ).sort_values(
            by="abs_wald",
            ascending=False,
        )
    else:
        # Selecting either based on abundance or size
        assert select_cags_by in ["abundance", "size"], select_cags_by

        ranked_df = cag_annotations(fp).sort_values(
            by="mean_abundance" if select_cags_by == "abundance" else "size",
            ascending=False
        )

    # Either way, return the set of CAGs in the size range
    return ranked_df.query(
        "size_log10 >= {}".format(cag_size_range[0])
    ).query(
        "size_log10 <= {}".format(cag_size_range[1])
    ).head(
        n_cags
    ).index.values


@app.callback(
    Output('cag-abundance-heatmap-graph', 'figure'),
    [
        Input({"type": "heatmap-select-cags-by", "parent": "abundance-heatmap"}, 'value'),
        Input('cag-abundance-heatmap-ncags', 'value'),
        Input({'name': 'cag-abundance-heatmap-size-range', 'type': 'cag-size-slider'}, 'value'),
        Input('cag-abundance-heatmap-metadata-dropdown', 'value'),
        Input('cag-abundance-heatmap-abundance-metric', 'value'),
        Input('cag-abundance-heatmap-cluster', 'value'),
        Input('cag-abundance-heatmap-annotate-cags-by', 'value'),
        Input('manifest-filtered', 'children'),
    ],
    [State("selected-dataset", "children")])
def abundance_heatmap_graph_callback(
    select_cags_by,
    n_cags,
    cag_size_range,
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

    # Get the top CAGs selected by this criterion
    cags_selected = get_cags_selected_by_criterion(
        fp,
        select_cags_by, 
        n_cags,
        cag_size_range,
    )
    
    # Get the abundance of the selected CAGs
    cag_abund_df = cag_abundances(fp).reindex(
        index=cags_selected
    ).T

    # If we are selecting CAGs by their association with a given parameter,
    # we will then plot the estimated coefficient with that parameter
    if select_cags_by.startswith("parameter-"):
        cag_estimate_dict = cag_associations(
            fp, 
            select_cags_by.replace("parameter-", "")
        ).reindex(
            index=cags_selected
        ).to_dict(
            orient="index" # Index by CAG ID, then column name
        )
    else:
        cag_estimate_dict = None

    # Get the taxonomic annotations for the selected CAGs
    if taxa_rank != "none":
        # Otherwise, get the taxonomic IDs for the CAGs
        cag_taxa_dict = {
            cag_id: taxonomic_gene_annotations(fp, cag_id=cag_id, rank=taxa_rank)
            for cag_id in cags_selected
        }
    else:
        cag_taxa_dict = {}

    # Draw the figure
    return draw_cag_abundance_heatmap(
        cag_abund_df,
        metadata_selected,
        abundance_metric,
        cluster_by,
        taxa_rank,
        manifest_json,
        manifest(fp),
        cag_taxa_dict,
        cag_estimate_dict,
    )

@app.callback(
    [
        Output('cag-abundance-heatmap-metadata-dropdown', value)
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
#####################################
# / CAG ABUNDANCE HEATMAP CALLBACKS #
#####################################


####################################
# CAG ANNOTATION HEATMAP CALLBACKS #
####################################
@app.callback(
    Output('cag-annotation-heatmap-graph', 'figure'),
    [
        Input({"type": "heatmap-select-cags-by", "parent": "annotation-heatmap"}, 'value'),
        Input('cag-annotation-heatmap-ncags', 'value'),
        Input({'name': 'cag-annotation-heatmap-size-range', 'type': 'cag-size-slider'}, 'value'),
        Input('cag-annotation-heatmap-annotation-type', 'value'),
        Input('cag-annotation-heatmap-nannots', 'value'),
    ],
    [State("selected-dataset", "children")])
def annotation_heatmap_graph_callback(
    select_cags_by,
    n_cags,
    cag_size_range,
    annotation_type,
    n_annots,
    selected_dataset,
):
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    if fp is None:
        return empty_figure()

    # Get the top CAGs selected by this criterion
    cags_selected = get_cags_selected_by_criterion(
        fp,
        select_cags_by,
        n_cags,
        cag_size_range,
    )

    # Get the full table of CAG annotations
    if annotation_type == "eggNOG_desc":

        # Read in the functional annotations for all CAGs in the index
        cag_annot_df = functional_gene_annotations(fp)

    elif annotation_type == "taxonomic":

        # Read in the full set of taxonomic annotations
        cag_annot_df = taxonomic_gene_annotations(fp)

    else:
        # Read in the taxonomic annotations at this rank
        cag_annot_df = taxonomic_gene_annotations(fp, rank=annotation_type)

    # Subset to just the CAGs in this list
    cag_annot_df = cag_annot_df.loc[
        cag_annot_df["CAG"].isin(cags_selected)
    ]

    # If the CAGs are selected by parameter, then fetch the annotations by parameter
    if select_cags_by.startswith("parameter-"):

        # Parse the parameter from the `select_cags_by` value
        parameter = select_cags_by.replace("parameter-", "")

        # Read in the enrichments of each annotation according to this parameter

        # To show all taxonomic annotations, read in the enrichments at multiple ranks
        if annotation_type == "taxonomic":
            enrichment_df = pd.concat([
                enrichments(
                    fp, 
                    parameter, 
                    rank
                ).assign(
                    rank=rank
                ).reset_index()
                for rank in ["species", "genus", "family"]
            ]).reset_index(
                drop=True
            )
        else:
            # Otherwise just show the enrichments at this single level
            enrichment_df = enrichments(
                fp, 
                parameter, 
                annotation_type
            )

        # If we are selecting CAGs by their association with a given parameter,
        # we will then plot the estimated coefficient with that parameter
        cag_estimate_dict = cag_associations(
            fp,
            parameter
        ).reindex(
            index=cags_selected
        ).to_dict(
            orient="index"  # Index by CAG ID, then column name
        )
    else:

        parameter = None

        # Do not show any enrichments for annotations
        enrichment_df = None

        # Do not show any estimates for CAGs
        cag_estimate_dict = None

    # Draw the figure
    return draw_cag_annotation_heatmap(
        cag_annot_df,
        annotation_type,
        enrichment_df,
        cag_estimate_dict,
        n_annots,
    )
######################################
# / CAG ANNOTATION HEATMAP CALLBACKS #
######################################


###################################
# ANNOTATION ENRICHMENT CALLBACKS #
###################################
@app.callback(
    Output("annotation-enrichment-graph", "figure"),
    [
        Input("annotation-enrichment-type", "value"),
        Input({
            "type": "corncob-parameter-dropdown",
            "name": "annotation-enrichment-parameter"
        }, "value"),
        Input("annotation-enrichment-plotn", "value"),
        Input("annotation-enrichment-show-pos-neg", "value"),
        Input("annotation-enrichment-page-num", "children"),
    ],
    [State("selected-dataset", "children")]
)
def annotation_enrichment_graph_callback(
    annotation,
    parameter,
    plotn,
    show_pos_neg,
    page_num,
    selected_dataset,
):
    """Render the graph for the annotation enrichment card."""
    
    # Get the path to the selected dataset
    fp = parse_fp(selected_dataset)

    # No dataset has been selected or no parameter has been selected
    if fp is None or parameter == "none":
        return empty_figure()
    else:

        # Read in the enrichments
        enrichment_df = enrichments(
            fp, 
            parameter, 
            annotation
        )

        # If there are no enrichments, skip it
        if enrichment_df is None:
            return empty_figure()

        # Subset to positive / negative estimated coefficients
        if show_pos_neg == "negative":
            enrichment_df = enrichment_df.query(
                "estimate < 0"
            )
        elif show_pos_neg == "positive":
            enrichment_df = enrichment_df.query(
                "estimate > 0"
            )

        # Sort by absolute Wald
        enrichment_df = enrichment_df.sort_values(
            by="abs_wald", 
            ascending=False
        )

        # Remove the Psort annotations
        enrichment_df = enrichment_df.loc[[
            n.startswith("Psort") is False
            for n in enrichment_df.index.values
        ]]

        # Parse the `plotn` and `page_num` parameter
        if isinstance(page_num, list):
            page_num = page_num[0]
        else:
            page_num = 1

        enrichment_df = enrichment_df.head(
            plotn * page_num
        ).tail(
            plotn
        )

        return draw_enrichment_graph(
            enrichment_df,
            annotation,
            parameter, 
        )


@app.callback(
    Output('annotation-enrichment-page-num', 'children'),
    [Input('annotation-enrichment-button-next', 'n_clicks'),
    Input('annotation-enrichment-button-previous', 'n_clicks'),
    Input('annotation-enrichment-button-first', 'n_clicks')],
    [State('annotation-enrichment-page-num', 'children')]
)
def annotation_enrichment_click(btn1, btn2, btn3, page_num):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    if 'annotation-enrichment-button-next' in changed_id:
        return [page_num[0] + 1]
    elif 'annotation-enrichment-button-previous' in changed_id:
        return [max(page_num[0] - 1, 1)]
    else:
        return [1]
#####################################
# / ANNOTATION ENRICHMENT CALLBACKS #
#####################################


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
        assert "neg_log10_pvalue" in df.columns.values, df.columns.values
        max_value = cag_associations(
            fp, 
            parameter
        )["neg_log10_pvalue"].max()

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

    cag_id = int(cag_id)

    # Format the DataFrame as needed to make a go.Sunburst
    cag_tax_df = taxonomic_gene_annotations(fp, cag_id=cag_id)

    if cag_tax_df is None:
        return empty_figure()
    else:
        return draw_taxonomy_sunburst(cag_tax_df, cag_id, min_ngenes)
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
        Input("manifest-table-bulk-select-apply", "n_clicks")
    ],
    [State("manifest-table-bulk-select-formula", "value")]
)
def update_manifest_table(selected_dataset, bulk_select_button, bulk_select_formula):
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

    # By default, select all of the rows
    selected_rows = np.arange(0, manifest_df.shape[0])

    # If the user has used the bulk select feature, apply the provided formula
    bulk_selected_rows = []
    if bulk_select_formula is not None and len(bulk_select_formula) > 0:
        try:
            bulk_selected_rows = manifest_df.query(bulk_select_formula).index.values
        except:
            pass
    # Only use the results of the formula if some samples were selected
    if len(bulk_selected_rows) > 0:
        selected_rows = bulk_selected_rows

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

@app.callback(
    [
        Output('manifest-table-bulk-select-open', 'n_clicks'),
        Output('manifest-table-bulk-select-apply', 'n_clicks'),
    ],
    [
        Input("selected-dataset", "children"),
    ])
def manifest_clear_formula(selected_dataset):
    """Reset the buttons when the dataset changes."""
    return 0, 0

@app.callback(
    Output("manifest-table-bulk-select-modal", "is_open"),
    [
        Input("manifest-table-bulk-select-open", "n_clicks"),
        Input("manifest-table-bulk-select-apply", "n_clicks"),
    ]
)
def manifest_bulk_select_toggle_modal(open_n, apply_n):
    if open_n is None:
        return False

    dismiss_clicks = 0
    if apply_n is not None:
        dismiss_clicks += apply_n
    
    return open_n > dismiss_clicks
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

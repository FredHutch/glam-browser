#!/usr/bin/env python3

from collections import defaultdict
from helpers.io import parse_directory
from helpers.io import hdf5_manifest
from helpers.io import hdf5_richness
from helpers.io import hdf5_cag_summary
from helpers.io import hdf5_metrics
from helpers.io import hdf5_distances
from helpers.io import hdf5_corncob
from helpers.io import hdf5_taxonomy
from helpers.io import hdf5_get_item
from helpers.layout import navbar_simple
from helpers.layout import dataset_summary_card
from helpers.layout import experiment_summary_card
from helpers.layout import update_experiment_summary_card
from helpers.layout import richness_card
from helpers.plotting import update_richness_graph
from helpers.plotting import run_pca
from helpers.plotting import run_tsne
from helpers.plotting import update_ordination_graph
from helpers.plotting import draw_cag_summary_graph_hist
from helpers.plotting import draw_cag_summary_graph_scatter
from helpers.plotting import draw_volcano_graph
from helpers.plotting import draw_taxonomy_sunburst
from helpers.plotting import draw_single_cag_graph
from helpers.plotting import parse_manifest_json
from helpers.layout import ordination_card
from helpers.layout import cag_summary_card
from helpers.layout import volcano_card
from helpers.layout import taxonomy_card
from helpers.taxonomy import make_cag_tax_df
from helpers.layout import single_cag_card
from helpers.layout import manifest_card
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
# CACHE FUNCTIONS TO READ IN DATA #
###################################
@cache.memoize()
def manifest(fp):
    return hdf5_manifest(fp)

@cache.memoize()
def richness(fp):
    return hdf5_richness(fp)

@cache.memoize()
def cag_summary(fp):
    return hdf5_cag_summary(fp)

@cache.memoize()
def distances(fp, metric):
    return hdf5_distances(fp, metric)

@cache.memoize()
def metrics(fp):
    return hdf5_metrics(fp)

@cache.memoize()
def corncob(fp):
    return hdf5_corncob(fp)

@cache.memoize()
def taxonomy(fp):
    return hdf5_taxonomy(fp)

@cache.memoize()
def read_cag_annotations(fp, cag_id):
    return hdf5_get_item(
        fp,
        "/annot/gene/all",
        where="CAG == {}".format(cag_id),
        columns = ["gene", "tax_id"]
    )

@cache.memoize()
def read_cag_abundance(fp, cag_id):
    return hdf5_get_item(
        fp,
        "/abund/cag/wide",
        where="CAG == {}".format(cag_id),
    ).loc[cag_id]


#####################################
# CACHE FUNCTIONS TO SUMMARIZE DATA #
#####################################
@cache.memoize()
def cag_summary_describe(fp):
    return cag_summary(fp).describe()

@cache.memoize()
def metadata_fields(fp):
    return [
        f
        for f in manifest(fp).columns.values
        if f not in ["R1", "R2", "I1", "I2"]
    ]

@cache.memoize()
def max_neg_log_pvalue(fp):
    df = corncob(fp)
    if df is None:
        return defaultdict(1)
    else:
        return df.groupby(
            "parameter"
        )["neg_log_pvalue"].max()

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
                        html.Div(id="load-data-output", style={"display": "none"}),
                        html.Div(id="load-data-input", style={"display": "none"}),
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
                ordination_card(),
                cag_summary_card(),
                volcano_card(),
                taxonomy_card(),
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
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        return [{'label': 'None', 'value': 'none'}]
    else:
        # Get the path to the indicated HDF5
        fp = page_data["contents"][selected_dataset[0]]["fp"]
        options = [
            {'label': 'None', 'value': 'none'},
        ] + [
            {
                "label": f,
                "value": f
            }
            for f in metadata_fields(fp)
        ]
        return options
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

    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        min_val = 0
        max_val = 1
    else:
        # Get the path to the indicated HDF5
        fp = page_data["contents"][selected_dataset[0]]["fp"]
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
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        max_value = 3
        min_value = 0
    else:
        # Get the path to the indicated HDF5
        fp = page_data["contents"][selected_dataset[0]]["fp"]
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
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
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
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
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
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        return
    else:
        # Get the path to the indicated HDF5
        fp = page_data["contents"][selected_dataset[0]]["fp"]
        return update_experiment_summary_card(
            metrics(fp)
        )
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
        Input('manifest-filtered', 'children'),
    ],
    [State("selected-dataset", "children")],
)
def richness_graph_callback(
    selected_metric, 
    selected_type, 
    manifest_json,
    selected_dataset, 
):
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        return go.Figure()
    else:
        fp = page_data["contents"][selected_dataset[0]]["fp"]
        return update_richness_graph(
            richness(fp),
            selected_metric,
            selected_type,
            manifest_json,
            manifest(fp)
        )
############################
# / RICHNESS CARD CALLBACK #
############################


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
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        return go.Figure()
    else:
        fp = page_data["contents"][selected_dataset[0]]["fp"]
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
        Input('cag-summary-log', 'value'),
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
):
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        return go.Figure()
    else:
        # Get the path to the indicated HDF5
        fp = page_data["contents"][selected_dataset[0]]["fp"]
        return draw_cag_summary_graph_hist(
            cag_summary(fp),
            metric_primary,
            size_range,
            entropy_range,
            prevalence_range,
            abundance_range,
            nbinsx,
            log_scale,
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
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        return go.Figure()
    else:
        # Get the path to the indicated HDF5
        fp = page_data["contents"][selected_dataset[0]]["fp"]
        return draw_cag_summary_graph_scatter(
            cag_summary(fp),
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


#####################
# VOLCANO CALLBACKS #
#####################
@app.callback(
    Output('volcano-graph', 'figure'),
    [
        Input("selected-dataset", "children"),
        Input('volcano-parameter-dropdown', 'value'),
        Input('volcano-pvalue-slider', 'value'),
        Input('volcano-fdr-radio', 'value'),
        Input('global-selected-cag', 'children'),
    ])
def volcano_graph_callback(
    selected_dataset,
    parameter, 
    neg_log_pvalue_min, 
    fdr_on_off, 
    selected_cag_json
):
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        return go.Figure()
    else:
        # Get the path to the indicated HDF5
        fp = page_data["contents"][selected_dataset[0]]["fp"]
        return draw_volcano_graph(
            corncob(fp),
            parameter, 
            neg_log_pvalue_min, 
            fdr_on_off, 
            selected_cag_json
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
    [
        Output('volcano-parameter-dropdown', value)
        for value in ["options", "value"]
    ],
    [
        Input("selected-dataset", "children"),
    ])
def update_volcano_parameter_dropdown(selected_dataset):
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        options =[
            {'label': 'None', 'value': 'none'},
        ]
        default_value = "none"
    else:
        # Get the path to the indicated HDF5
        fp = page_data["contents"][selected_dataset[0]]["fp"]

        # Get the list of parameters
        parameter_list = corncob(fp)[
            "parameter"
        ].drop_duplicates(
        ).sort_values(
        ).tolist(
        )

        options = [
            {'label': l, 'value': l}
            for l in parameter_list
        ]

        if len(parameter_list) > 1:
            default_value = parameter_list[1]
        else:
            default_value = parameter_list[0]

    return options, default_value

@app.callback(
    [
        Output('volcano-pvalue-slider', value)
        for value in ['max', 'marks']
    ],
    [
        Input("selected-dataset", "children"),
        Input('volcano-parameter-dropdown', 'value'),
    ])
def update_volcano_pvalue_slider(selected_dataset, parameter):
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        max_value = 1
    else:
        # Get the path to the indicated HDF5
        fp = page_data["contents"][selected_dataset[0]]["fp"]
        max_value = corncob(fp).query(
            "parameter == '{}'".format(parameter)
        )["neg_log_pvalue"].max()

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
        Input('global-selected-cag', 'children'),
    ])
def update_taxonomy_graph(selected_dataset, min_ngenes, selected_cag_json):
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        marks = {
            n: n
            for n in ["0", "1"]
        }
        return go.Figure(), 1, marks

    # Get the path to the indicated HDF5
    fp = page_data["contents"][selected_dataset[0]]["fp"]

    # Parse the selected CAG data
    if selected_cag_json is not None:
        # Set the points which are selected in the scatter plot
        cag_id = json.loads(selected_cag_json)["id"]
    else:
        # Default CAG to plot
        cag_id = 1000

    # Read in the taxonomic annotations for this CAG
    cag_df = read_cag_annotations(fp, cag_id)

    # Format the DataFrame as needed to make a go.Sunburst
    cag_tax_df = make_cag_tax_df(cag_df["tax_id"], taxonomy(fp))

    return draw_taxonomy_sunburst(cag_tax_df, cag_id, min_ngenes)
###########################
# / CAG TAXONOMY CALLBACK #
###########################

#######################
# SINGLE CAG CALLBACK #
#######################
@app.callback(
    Output('single-cag-graph', 'figure'),
    [
        Input('global-selected-cag', 'children'),
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
    selected_cag_json,
    xaxis,
    plot_type,
    color,
    facet,
    log_scale,
    manifest_json,
    selected_dataset,
):
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        return go.Figure()

    # Get the path to the indicated HDF5
    fp = page_data["contents"][selected_dataset[0]]["fp"]

    # Get the filtered manifest from the browser
    plot_manifest_df = parse_manifest_json(manifest_json, manifest(fp))

    # Parse the selected CAG data
    if selected_cag_json is not None:
        # Set the points which are selected in the scatter plot
        cag_id = json.loads(selected_cag_json)["id"]
    else:
        # Default CAG to plot
        cag_id = 1000

    plot_df = plot_manifest_df.assign(
        CAG_ABUND = read_cag_abundance(fp, int(cag_id))
    )

    return draw_single_cag_graph(
        plot_df,
        cag_id,
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

    if selected_dataset == [-1] or selected_dataset == ["-1"]:
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

    # Get the path to the indicated HDF5
    fp = page_data["contents"][selected_dataset[0]]["fp"]

    # Get the manifest for this dataset
    manifest_df = manifest(fp).reset_index()

    columns=[
        {"name": i, "id": i}
        for i in manifest_df.columns
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
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        return None

    # Get the path to the indicated HDF5
    fp = page_data["contents"][selected_dataset[0]]["fp"]

    # Parse the list of indexes which we should save
    selected_rows = json.loads(selected_rows_json)

    # Filter down the manifest table and save as JSON
    return manifest(fp).reset_index(
    ).reindex(
        index=selected_rows
    ).to_json(date_format='iso', orient='records')

@app.callback(
    Output('manifest-table', 'hidden_columns'),
    [
        Input("selected-dataset", "children"),
        Input('manifest-table-select-columns', 'value'),
    ])
def manifest_update_columns_selected(selected_dataset, selected_columns):
    """Allow the user to hide selected columns in the manifest table."""
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        return []

    # Get the path to the indicated HDF5
    fp = page_data["contents"][selected_dataset[0]]["fp"]

    return [
        n
        for n in manifest(fp).columns.values
        if n not in selected_columns
    ]
########################
# / MANIFEST CALLBACKS #
########################

######################
# LOAD DATA CALLBACK #
######################
@app.callback(
    Output('load-data-output', 'children'),
    [Input('load-data-input', 'children')])
def load_data_callback(trigger):
    for dataset in page_data["contents"]:
        _ = manifest(dataset["fp"])
        _ = richness(dataset["fp"])
        _ = cag_summary(dataset["fp"])
        _ = distances(dataset["fp"], "euclidean")
        _ = distances(dataset["fp"], "aitchison")
        _ = distances(dataset["fp"], "braycurtis")
        _ = metrics(dataset["fp"])
        _ = corncob(dataset["fp"])
        _ = taxonomy(dataset["fp"])
        _ = read_cag_annotations(dataset["fp"], 1000)
    return 1
########################
# / LOAD DATA CALLBACK #
########################



if __name__ == '__main__':

    app.run_server(
        host='0.0.0.0',
        port=8050,
    )

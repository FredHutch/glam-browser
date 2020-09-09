#!/usr/bin/env python3

from collections import defaultdict
from helpers.io import Manifest
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
from helpers.layout import plot_cag_card
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
from filelock import Timeout, FileLock
import zarr


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

# Get the Google Tag Manager container ID, if specified
gtm_container = os.getenv("GTM_CONTAINER")

# Parse the contents of the directory
# See the documentation for options on how to structure the data folder
page_data = Manifest(data_folder)


###################
# CSS STYLESHEETS #
###################
external_stylesheets = [
    {
        'href': 'https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/css/bootstrap.min.css',
        'rel': 'stylesheet',
        'integrity': 'sha384-Vkoo8x4CGsO3+Hhxv8T/Q5PaXtkKtu6ug5TOeNV6gBiFeWPGFN9MuhOf23Q9Ifjh',
        'crossorigin': 'anonymous'
    },
    "assets/page_click.js"
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
    ).apply(
        lambda c: c.apply(str) if c.name == "specimen" else c
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
    df = hdf5_get_item(
        fp, 
        "/specimen_metrics"
    )
    # Fix the mismatch in index variable type, if any
    df = pd.DataFrame({
        col_name: {
            str(k): v
            for k, v in zip(
                df["specimen"].values,
                df[col_name].values
            )
            if pd.isnull(v) is False
        }
        for col_name in df.columns.values
        if col_name != "specimen"
    })
    return df

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
def cag_abundances(fp, cag_id=None, specimen=None, retry=True, timeout=60):
    assert cag_id is not None or specimen is not None, "Must specify cag_id or specimen"
    assert cag_id is None or specimen is None, "Cannot specify cag_id with specimen"

    # Read in the list of specimens and index values
    specimen_list = cache.get("{}/specimen_index".format(fp))
    cag_id_list = cache.get("{}/cag_id_index".format(fp))

    if specimen_list is None or cag_id_list is None:
        if retry:
            cache_cag_abundances(fp)
            return cag_abundances(fp, cag_id=cag_id, specimen=specimen, retry=False)
        else:
            return

    # Make sure that the CAG ID or specimen is in the index
    if cag_id is not None:
        cag_id_ix = cag_id_list.index(int(cag_id))
        assert cag_id_ix is not None, "Key {} not found in index".format(cag_id)
    elif specimen is not None:
        specimen_ix = specimen_list.index(str(specimen))
        assert specimen_ix is not None, "Key {} not found in index".format(specimen)

    # Set up a file lock to prevent multiple concurrent access attempts
    lock = FileLock("{}.cag_abundances.zarr.lock".format(fp), timeout=timeout)

    # Read in the abundances from the zarr index
    with lock:
        z = zarr.open(
            "{}.cag_abundances.zarr".format(fp),
            mode="r"
        )
        if cag_id is not None:
            values = z[cag_id_ix, :]
            return pd.Series(
                values,
                index=specimen_list
            )

        elif specimen is not None:
            values = z[:, specimen_ix]
            return pd.Series(
                values,
                index=cag_id_list
            )

def cache_cag_abundances(fp, timeout=60):

    # If the zarr cache has already been made, don't write it again
    if cag_abundances(fp, cag_id=0, retry=False) is not None:
        return

    # Read in the complete set of CAG abundances
    df = hdf5_get_item(
        fp, 
        "/cag_abundances"
    ).set_index(
        "CAG"
    )

    # Set up a file lock to prevent multiple concurrent access attempts
    lock = FileLock("{}.cag_abundances.zarr.lock".format(fp), timeout=timeout)

    # Write out in zarr format
    with lock:
        zarr.save(
            "{}.cag_abundances.zarr".format(fp),
            df
        )

    # Save the list of specimens and CAG IDs to redis
    cache.set(
        "{}/specimen_index".format(fp),
        [str(v) for v in df.columns.values]
    )
    cache.set(
        "{}/cag_id_index".format(fp),
        [int(v) for v in df.index.values]
    )


# Populate the cache for all files in this folder
for fp in page_data.all_filepaths:
    cache_cag_abundances(fp)


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

    df = hdf5_get_item(
        fp,
        "/cag_associations/{}".format(parameter)
    )
    
    if df is None:

        return

    else:

        # Set the index
        df.set_index("CAG", inplace=True)

        # Calculate the Wald
        df = df.assign(
            wald = df["estimate"] / df["std_error"]
        )

        # Calculate the absolute value of the Wald
        df = df.assign(
            abs_wald = df["wald"].abs()
        )

        # Sort by abs_wald (descending)
        return df.sort_values(
            by="abs_wald",
            ascending=False,
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
        if 'wald' not in df.columns.values:
            df = df.assign(
                wald = df["estimate"] / df["std_error"]
            )
        if 'abs_wald' not in df.columns.values:
            df = df.assign(
                abs_wald = df["wald"].abs()
            )

        # Protect against duplicate labels by taking the most associated one
        df = df.sort_values(
            by="abs_wald",
            ascending=False,
        ).groupby(
            "label"
        ).head(
            1
        ).set_index(
            "label"
        )
    return df

def functional_gene_annotations(fp, cag_id):
    df = functional_gene_annotations_shard(fp, cag_id % 1000)
    if df is None:
        return
    # Filter down to the annotations for this CAG
    df = df.query(
        "CAG == {}".format(cag_id)
    )
    if df.shape[0] == 0:
        return
    # Explicitly mask the "Psort" annotations
    df = df.loc[
        df["label"].apply(lambda n: n.startswith("Psort") is False)
    ]
    # Don't return a DataFrame with zero rows
    if df.shape[0] == 0:
        return None
    else:
        return df

@cache.memoize()
def functional_gene_annotations_shard(fp, group_ix):
    return hdf5_get_item(
        fp,
        "/gene_annotations/functional/{}".format(group_ix)
    )


def taxonomic_gene_annotations(fp, cag_id, rank="all"):
    df = taxonomic_gene_annotations_shard(fp, cag_id % 1000, rank=rank)

    if df is None:
        return

    # Filter down to the annotations for this CAG
    df = df.query(
        "CAG == {}".format(cag_id)
    )

    if df.shape[0] == 0:
        return

    # Format required  columns as integers
    df = df.apply(
        lambda c: c.fillna(0).apply(int) if c.name in ["count", "tax_id", "parent", "total", "CAG"] else c
    )
    
    # Sort by number of hits
    df = df.sort_values(
        by="count", 
        ascending=False
    )

    return df

@cache.memoize()
def taxonomic_gene_annotations_shard(fp, group_ix, rank="all"):
    return hdf5_get_item(
        fp,
        "/gene_annotations/taxonomic/{}/{}".format(rank, group_ix)
    )

@cache.memoize()
def genome_manifest(fp):
    return hdf5_get_item(
        fp,
        "/genome_manifest"
    )
def genomic_alignment_annotations(fp, cag_id):
    df = genomic_alignment_annotations_shard(fp, cag_id % 1000)

    if df is None:
        return
    
    else:

        # Filter down to the annotations for this CAG
        df = df.query(
            "CAG == {}".format(cag_id)
        )

        if df.shape[0] == 0:
        
            return
        
        else:

            return df
@cache.memoize()
def genomic_alignment_annotations_shard(fp, group_ix):
    return hdf5_get_item(
        fp, 
        "/genome_containment/{}".format(group_ix)
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
        dcc.Location(id='url', refresh=False),
        navbar_simple(),
        html.Div(
            children=[],
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
                plot_cag_card(),
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


############################
# LOGIN / LOGOUT CALLBACKS #
############################
# Open and close the login modal based on button clicks
@app.callback(
    Output("login-modal", "is_open"),
    [
        Input("login-button", "n_clicks"),
        Input("login-modal-apply-button", "n_clicks"),
    ]
)
def open_close_login_modal(open_nclicks, close_nclicks):
    if open_nclicks > close_nclicks:
        return True
    else:
        return False
# Logging out will clear the username and password
@app.callback(
    [
        Output("login-username", "value"),
        Output("login-key", "value"),
        Output({"type": "open-dataset-button", "index": -1}, "n_clicks")
    ],
    [
        Input("logout-button", "n_clicks"),
    ],
    [
        State({"type": "open-dataset-button", "index": -1}, "n_clicks")
    ]
)
def clear_username_password(logout_nclicks, main_menu_nclicks):
    return None, None, main_menu_nclicks + 1
# Logging in will change the login / logout button appearance
@app.callback(
    [
        Output("login-button", "style"),
        Output("logout-button", "style"),
        Output("logout-button", "children"),
    ],
    [
        Input("login-modal", "is_open"),
        Input("logout-button", "n_clicks"),
    ],
    [
        State("login-username", "value"),
    ]
)
def login_logout_callback(login_is_open, logout_nclicks, username):
    # Style for buttons being shown or hidden
    show_style = {"margin": "10px"}
    hide_style = {"margin": "10px", "display": "none"}

    # Set up the context for the callback
    ctx = dash.callback_context

    # Figure out which element changed
    if not ctx.triggered:

        # Nothing has changed yet
        return show_style, hide_style, "Logout"

    elif ctx.triggered[0]["prop_id"] == "logout-button.n_clicks":

        # The logout button was clicked
        return show_style, hide_style, "Logout"

    elif username is None:
        # If there is no username, show the user the login button
        return show_style, hide_style, "Logout"

    # implicit else
    # The login model must have just opened or closed
    # If there is a username, use that to update the logout button text
    logout_text = "Logout ({})".format(username)
    return hide_style, show_style, logout_text
##############################
# / LOGIN / LOGOUT CALLBACKS #
##############################


#############################
# SUMMARY DISPLAY CALLBACKS #
#############################
@app.callback(
    Output("page-title", "brand"),
    [
        Input("selected-dataset", 'children'),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ],
)
def page_title_callback(selected_dataset, page, key):
    """Return the title of the page."""
    page_title = page_data.page_title(page=page, key=key)

    if page_title is None:
        return "GLAM Browser"
    else:
        return page_title

@app.callback(
    Output("summary-display", "children"),
    [
        Input({
            "type": "open-dataset-pressed",
            "index": -1
        }, "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ],
    [
        State("selected-dataset", 'children'),
    ],
)
def summary_display_callback(_, page, key, selected_dataset):
    """Fill out the summary display div."""
    page_description = page_data.page_description(page=page, key=key)
    if page_description is None:
        page_description = """
This page may have been reached in error, please check your username/password and try again.

For assistance, please consult the GLAM Browser documentation or contact its maintainers.
"""

    # The summary page will start with a description panel
    children = [
        dbc.Card(
            dbc.CardBody([
                dcc.Markdown(page_description),
            ])
        )
    ]

    # Add the summary cards for each dataset
    for ix, dataset in enumerate(page_data.dataset_list(page=page, key=key)):
        children.append(
            dataset_summary_card(
                ix, dataset
            )
        )

    return children

###############################
# / SUMMARY DISPLAY CALLBACKS #
###############################



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
    if n_clicks > 0:
        return time()
    else:
        return -1
@app.callback(
    Output({"type": "open-dataset-button", "index": -1}, "style"),
    [Input("selected-dataset", "children")]
)
def show_hide_main_menu_button(selected_dataset):
    # Only show the Main Menu button if a dataset has been selected
    if selected_dataset[0] == -1:
        return {"margin": "10px", "display": "none"}
    else:
        return {"margin": "10px"}
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
######################################
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
    )]
)
def open_dataset_switch(button_timestamps):
    """Control whether a dataset is open, or the summary page."""
    
    # Make a list of the timestamps
    button_timestamps = [
        float(t)
        for t in button_timestamps
    ]

    # Get the most recent timestamp
    most_recent_timestamp = max(button_timestamps)

    # Switch to whichever button was pressed most recently
    for button_ix, button_timestamp in enumerate(button_timestamps):

        # If this button was pressed most recently
        if button_timestamp == most_recent_timestamp:

            # Return the index position (adjusted -1 to include the main menu)
            return [button_ix - 1]
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
        Input("login-username", "value"),
        Input("login-key", "value"),
        Input({
            "type": "metadata-field-dropdown",
            "name": MATCH
        }, "value"),
    ],
)
def metadata_field_dropdown_callback(selected_dataset, page, key, dummy_value):
    """Update the metadata field dropdown for the selected dataset."""

    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
        Input("login-username", "value"),
        Input("login-key", "value"),
        Input({
            "type": "cag-metric-slider",
            "name": MATCH,
            "metric": MATCH,
        }, "id"),
    ],
)
def cag_metric_slider_callback_max(selected_dataset, page, key, slider_id):
    """Update any CAG metric slider for the selected dataset."""
    metric = slider_id["metric"]

    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
        Input("login-username", "value"),
        Input("login-key", "value"),
        Input({
            "type": "cag-size-slider",
            "name": MATCH
        }, "id"),
    ],
)
def cag_size_slider_callback(selected_dataset, page, key, dummy_value):
    """Update the CAG size slider for the selected dataset."""

    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
    value=[1., max_value]

    return [max_value, min_value, step_value, marks, value]
##############################
# / CAG SIZE SLIDER CALLBACK #
##############################


###############################
# SHOW / HIDE SUMMARY DISPLAY #
###############################
@app.callback(
    Output("summary-display", 'style'),
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ],
)
def show_hide_summary_display(selected_dataset, page, key):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ],
)
def show_hide_detail_display(selected_dataset, page, key):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ],
)
def experiment_summary_card_callback(selected_dataset, page, key):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    if fp is None:
        return
    else:
        # Get the description of this dataset
        dataset_description = {
            d["fp"]: d
            for d in page_data.dataset_list(
                page=page, 
                key=key
            )
        }[
            fp.split("/")[-1]
        ].get(
            "description", ""
        )

        return update_experiment_summary_card(
            experiment_metrics(fp),
            dataset_description
        )
@app.callback(
    Output("experiment-summary-card-header", 'children'),
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ],
)
def experiment_summary_card_header_callback(selected_dataset, page, key):
    if selected_dataset == [-1] or selected_dataset == ["-1"]:
        return "Experiment"
    else:
        # Get the list of datasets on this page
        dataset_list = page_data.dataset_list(
            page=page,
            key=key
        )

        # If we cannot access the data, return a safe value
        if dataset_list is None or len(dataset_list) == 0:
            return "Experiment"

        # Return the name of the indicated HDF5
        return "Experiment: {}".format(dataset_list[
            selected_dataset[0]
        ][
            "name"
        ])
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
    [
        State("selected-dataset", "children"),
        State("login-username", "value"),
        State("login-key", "value"),
    ],
)
def richness_graph_callback(
    selected_metric, 
    selected_type, 
    selected_metadata,
    log_x,
    manifest_json,
    selected_dataset, 
    page,
    key,
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
    [
        State("selected-dataset", "children"),
        State("login-username", "value"),
        State("login-key", "value"),
    ],
)
def single_sample_graph_callback(
    primary_sample, 
    secondary_sample, 
    display_metric,
    selected_dataset, 
    page,
    key,
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    if fp is None:
        return empty_figure()
    elif primary_sample == "none" or display_metric is None:
        return empty_figure()
    else:

        if secondary_sample == "cag_size":
            return plot_sample_vs_cag_size(
                cag_abundances(fp, specimen=primary_sample),
                primary_sample,
                cag_annotations(fp),
                display_metric,
            )
        else:
            return plot_samples_pairwise(
                cag_abundances(fp, specimen=primary_sample),
                primary_sample,
                cag_abundances(fp, specimen=secondary_sample),
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
        Input("login-username", "value"),
        Input("login-key", "value"),
    ],
)
def update_single_sample_primary_dropdown(
    selected_dataset, page, key
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
        Input("login-username", "value"),
        Input("login-key", "value"),
    ],
)
def update_single_sample_secondary_dropdown(
    selected_dataset, page, key
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
        State("selected-dataset", "children"),
        State("login-username", "value"),
        State("login-key", "value"),
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
    page,
    key,
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
    Output({'type': 'metadata-field-dropdown', 'name': 'ordination-metadata'}, 'value'),
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ]
)
def ordination_graph_metadata(
    selected_dataset,
    page,
    key,
):
    # Get the default value, if specified
    default_value = page_data.default(selected_dataset, "overlay_label", page=page, key=key)

    if default_value is None:
        return "none"
    else:
        return default_value

@app.callback(
    Output('ordination-anosim-results', 'children'),
    [
        Input('ordination-metric', 'value'),
        Input({'type': 'metadata-field-dropdown', 'name': 'ordination-metadata'}, 'value'),
        Input('manifest-filtered', 'children'),
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ])
def ordination_anosim_callback(
    metric,
    metadata,
    manifest_json,
    selected_dataset,
    page,
    key,
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
        Input("login-username", "value"),
        Input("login-key", "value"),
        Input('cag-summary-metric-primary', 'value'),
        Input('cag-summary-nbinsx-slider', 'value'),
        Input('cag-summary-histogram-log', 'value'),
        Input('cag-summary-histogram-metric', 'value'),
    ])
def cag_summary_graph_hist_callback(
    selected_dataset,
    page,
    key,
    metric_primary,
    nbinsx,
    log_scale,
    hist_metric,
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
        Input("login-username", "value"),
        Input("login-key", "value"),
    ],
    [State({"type": "heatmap-select-cags-by", "parent": MATCH}, 'value')])
def abundance_heatmap_graph_select_cags_callback(selected_dataset, page, key, _):
    """Add the corncob parameters to the CAG abundance heatmap CAG selection options."""

    # Set the base options
    options = [
        {"label": "Average Relative Abundance", "value": "abundance"},
        {"label": "Size (Number of Genes)", "value": "size"}
    ]

    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    if fp is None:
        # Return the base options if no dataset is selected
        return options, "abundance"
    else:
        # Get the default parameter, if defined in the manifest
        default_parameter = page_data.default(
            selected_dataset,
            "parameter",
            key=key,
            page=page
        )
        
        # Get the list of available parameters
        parameter_list = valid_parameters(fp)
        
        # Add the parameters to the list of options
        options = options + [
            {"label": parameter, "value": "parameter-{}".format(parameter)}
            for parameter in parameter_list
        ]

        # If the default in the manifest is valid, use that
        if default_parameter in parameter_list:
            return options, "parameter-{}".format(default_parameter)
        else:
            return options, "abundance"


def get_cags_selected_by_criterion(fp, select_cags_by, cag_size_range):
    """Function to select the top CAGs sorted on a particular criterion."""
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
    ).index.values


@app.callback(
    Output('abundance-heatmap-selected-cags', 'value'),
    [
        Input({"type": "heatmap-select-cags-by", "parent": "abundance-heatmap"}, 'value'),
        Input('cag-abundance-heatmap-ncags', 'value'),
        Input({'name': 'cag-abundance-heatmap-size-range', 'type': 'cag-size-slider'}, 'value'),
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ])
def abundance_heatmap_select_cags_callback(
    select_cags_by,
    n_cags,
    cag_size_range,
    selected_dataset,
    page,
    key,
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    if fp is None:
        return ""

    # Return the top CAGs selected by this criterion
    return ",".join([
        str(cag_id)
        for cag_id in get_cags_selected_by_criterion(
            fp,
            select_cags_by, 
            cag_size_range,
        )[
            :n_cags
        ]
    ])
@app.callback(
    Output('cag-abundance-heatmap-graph', 'figure'),
    [
        Input({"type": "heatmap-select-cags-by", "parent": "abundance-heatmap"}, 'value'),
        Input('abundance-heatmap-selected-cags', 'value'),
        Input('cag-abundance-heatmap-metadata-dropdown', 'value'),
        Input('cag-abundance-heatmap-abundance-metric', 'value'),
        Input('cag-abundance-heatmap-cluster', 'value'),
        Input('cag-abundance-heatmap-annotate-cags-by', 'value'),
        Input('manifest-filtered', 'children'),
    ],
    [
        State("selected-dataset", "children"),
        State("login-username", "value"),
        State("login-key", "value"),
    ])
def abundance_heatmap_graph_callback(
    select_cags_by,
    cags_selected,
    metadata_selected,
    abundance_metric,
    cluster_by,
    taxa_rank,
    manifest_json,
    selected_dataset,
    page,
    key,
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    if fp is None:
        return empty_figure()

    # Format the list of CAGs selected
    cags_selected = [
        int(cag_id.strip(" "))
        for cag_id in cags_selected.split(",")
    ]
    
    # Get the abundance of the selected CAGs
    cag_abund_df = pd.DataFrame({
        cag_id: cag_abundances(fp, cag_id=cag_id)
        for cag_id in cags_selected
    })

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
            cag_id: taxonomic_gene_annotations(fp, cag_id, rank=taxa_rank)
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
        Input("login-username", "value"),
        Input("login-key", "value"),
    ])
def update_heatmap_metadata_dropdown(
    selected_dataset, page, key
):
    """When a new dataset is selected, fill in the available metadata."""
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
    Output('annotation-heatmap-selected-cags', 'value'),
    [
        Input({"type": "heatmap-select-cags-by", "parent": "annotation-heatmap"}, 'value'),
        Input('cag-annotation-heatmap-ncags', 'value'),
        Input({'name': 'cag-annotation-heatmap-size-range', 'type': 'cag-size-slider'}, 'value'),
    ],
    [
        State('cag-annotation-heatmap-annotation-type', 'value'),
        State("selected-dataset", "children"),
        State("login-username", "value"),
        State("login-key", "value"),
    ])
def annotation_heatmap_select_cags(
    select_cags_by,
    n_cags,
    cag_size_range,
    annotation_type,
    selected_dataset,
    page,
    key,
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    if fp is None:
        return ""

    # Get the top CAGs selected by this criterion, ordered by priority
    cags_selected = get_cags_selected_by_criterion(
        fp,
        select_cags_by,
        cag_size_range,
    )

    # Record the final set of selected CAGs as a set
    final_cags_selected = set([])

    # Get the sizes (number of genes) of these CAGs
    cag_sizes = cag_annotations(fp)["size"].reindex(index=cags_selected)

    # Get the full table of CAG annotations
    if annotation_type == "eggNOG_desc":

        # Check CAGs in the order of priority
        for cag_id in cags_selected:

            # Check to see if it has any valid annotations
            if functional_gene_annotations(fp, cag_id) is not None:
                final_cags_selected.add(cag_id)

                if len(final_cags_selected) >= n_cags:
                    break

    elif annotation_type == "taxonomic":

        # Check CAGs in the order of priority
        for cag_id in cags_selected:

            # Check to see if it has any valid annotations
            if taxonomic_gene_annotations(fp, cag_id) is not None:
                final_cags_selected.add(cag_id)

                if len(final_cags_selected) >= n_cags:
                    break

    elif annotation_type == "genomes":

        # Check CAGs in the order of priority
        for cag_id in cags_selected:

            # Check to see if it has any valid annotations
            df = genomic_alignment_annotations(fp, cag_id)
            # Only keep this CAG if it has >50% containment to at least one genome
            if df is not None and df["cag_prop"].max() >= 0.5:

                final_cags_selected.add(cag_id)

                if len(final_cags_selected) >= n_cags:
                    break

    else:
        # Check CAGs in the order of priority
        for cag_id in cags_selected:

            # Check to see if it has any valid annotations
            if taxonomic_gene_annotations(fp, cag_id, rank=annotation_type) is not None:

                final_cags_selected.add(cag_id)

                if len(final_cags_selected) >= n_cags:
                    break

    return ",".join([
        str(cag_id)
        for cag_id in final_cags_selected
    ])

@app.callback(
    Output('cag-annotation-heatmap-graph', 'figure'),
    [
        Input({"type": "heatmap-select-cags-by", "parent": "annotation-heatmap"}, 'value'),
        Input('annotation-heatmap-selected-cags', 'value'),
        Input('cag-annotation-heatmap-annotation-type', 'value'),
        Input('cag-annotation-heatmap-nannots', 'value'),
    ],
    [
        State("selected-dataset", "children"),
        State("login-username", "value"),
        State("login-key", "value"),
    ])
def annotation_heatmap_graph_callback(
    select_cags_by,
    selected_cags,
    annotation_type,
    n_annots,
    selected_dataset,
    page,
    key,
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    if fp is None:
        return empty_figure()

    # Format the list of CAGs selected
    selected_cags = [
        int(cag_id.strip(" "))
        for cag_id in selected_cags.split(",")
    ]

    # Get the sizes (number of genes) of these CAGs
    cag_sizes = cag_annotations(fp)["size"].reindex(index=selected_cags)

    # Get the full table of CAG annotations
    if annotation_type == "eggNOG_desc":

        # Read in the functional annotations for the set of selected CAGs
        cag_annot_dict = {
            cag_id: functional_gene_annotations(fp, cag_id)
            for cag_id in selected_cags
            if functional_gene_annotations(fp, cag_id) is not None
        }

    elif annotation_type == "taxonomic":

        # Read in the taxonomic annotations for the set of selected CAGs
        cag_annot_dict = {
            cag_id: taxonomic_gene_annotations(fp, cag_id)
            for cag_id in selected_cags
            if taxonomic_gene_annotations(fp, cag_id) is not None
        }

    elif annotation_type == "genomes":

        # Read in the genomic alignment annotations for the set of selected CAGs
        cag_annot_dict = {
            cag_id: genomic_alignment_annotations(fp, cag_id).rename(
                columns = {
                    "genome": "name",
                    "n_genes": "count",
                }
            )
            for cag_id in selected_cags
            if genomic_alignment_annotations(fp, cag_id) is not None
        }

    else:
        # Read in the selected taxonomic annotations for the set of selected CAGs
        cag_annot_dict = {
            cag_id: taxonomic_gene_annotations(fp, cag_id, rank=annotation_type)
            for cag_id in selected_cags
            if taxonomic_gene_annotations(fp, cag_id, rank=annotation_type) is not None
        }

    # Check to see if the selected annotation information is available for these CAGs
    if len(cag_annot_dict) > 0:

        cag_annot_df = pd.concat([
            df.assign(
                CAG = cag_id
            )
            for cag_id, df in cag_annot_dict.items()
        ])
    
    else:
        return empty_figure()

    # If the CAGs are selected by parameter, then fetch the annotations by parameter
    if select_cags_by.startswith("parameter-"):

        # Parse the parameter from the `select_cags_by` value
        parameter = select_cags_by.replace("parameter-", "")

        # Read in the enrichments of each annotation according to this parameter

        # To show all taxonomic annotations, read in the enrichments at multiple ranks
        if annotation_type == "taxonomic":
            # In the combined DataFrame, make sure to include the rank of interest
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
        
        # The annotation type is genome alignment
        elif annotation_type == "genomes":

            # There are no annotations computed genome-by-genome
            enrichment_df = None

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
            index=selected_cags
        ).to_dict(
            orient="index"  # Index by CAG ID, then column name
        )
    else:

        parameter = None

        # Do not show any enrichments for annotations
        enrichment_df = None

        # Do not show any estimates for CAGs
        cag_estimate_dict = None

    # If the annotation type is "genomes", read in the table of genome names
    if annotation_type == "genomes":

        # Read in the table linking the genome 'id' to its name
        annotation_names = genome_manifest(fp)

        # Make sure that the data exists
        if annotation_names is not None:

            # Format as a Series, indexed by 'id'
            annotation_names = annotation_names.set_index("id")["name"]
        else:

            annotation_names = None
    else:
        annotation_names = None

    # Draw the figure
    return draw_cag_annotation_heatmap(
        cag_annot_df,
        annotation_type,
        enrichment_df,
        cag_estimate_dict,
        n_annots,
        annotation_names,
        cag_sizes,
    )

@app.callback(
    Output('cag-annotation-heatmap-annotation-type', 'options'),
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ])
def annotation_heatmap_type_options(
    selected_dataset,
    page,
    key,
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    # Set the default options
    options = [
        {'label': 'Functional', 'value': 'eggNOG_desc'},
        {'label': 'Taxonomic', 'value': 'taxonomic'},
        {'label': 'Species', 'value': 'species'},
        {'label': 'Genus', 'value': 'genus'},
        {'label': 'Family', 'value': 'family'},
    ]

    # Return the default options if no dataset is selected
    if fp is None:
        return options

    # Request genome alignment information for this dataset
    df = genome_manifest(fp)

    # Add the option to display genome information, if present
    if df is not None:
        options.append(
            {'label': 'Genomes', 'value': 'genomes'},
        )

    return options
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
            "group": "annotation-enrichment",
        }, "value"),
        Input("annotation-enrichment-plotn", "value"),
        Input("annotation-enrichment-show-pos-neg", "value"),
        Input("annotation-enrichment-page-num", "children"),
    ],
    [
        State("selected-dataset", "children"),
        State("login-username", "value"),
        State("login-key", "value"),
    ]
)
def annotation_enrichment_graph_callback(
    annotation,
    parameter,
    plotn,
    show_pos_neg,
    page_num,
    selected_dataset,
    page,
    key,
):
    """Render the graph for the annotation enrichment card."""
    
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
        Input("login-username", "value"),
        Input("login-key", "value"),
        Input({
            "type": "corncob-parameter-dropdown", 
            "group": "volcano-parameter",
        }, 'value'),
        Input("corncob-comparison-parameter-dropdown", "value"),
        Input({"name": 'volcano-cag-size-slider', "type": "cag-size-slider"}, 'value'),
        Input({
            'group': 'volcano-parameter',
            'type': 'corncob-pvalue-slider'
        }, 'value'),
        Input('volcano-fdr-radio', 'value'),
    ])
def volcano_graph_callback(
    selected_dataset,
    page,
    key,
    parameter, 
    comparison_parameter,
    cag_size_range,
    neg_log_pvalue_min, 
    fdr_on_off, 
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    # No dataset has been selected or no parameter has been selected
    if fp is None or parameter == "none":
        return empty_figure()
    else:

        # Get the associations for this parameter
        plot_df = cag_associations(fp, parameter)

        # If the comparison is to be made against a second parameter, add it
        if comparison_parameter != "coef":
            plot_df = pd.concat([
                plot_df.assign(parameter = parameter),
                cag_associations(
                    fp, 
                    comparison_parameter
                ).assign(
                    parameter=comparison_parameter
                )
            ])

        return draw_volcano_graph(
            plot_df,
            cag_annotations(fp),
            parameter, 
            comparison_parameter,
            cag_size_range,
            neg_log_pvalue_min, 
            fdr_on_off, 
        )

@app.callback(
    Output({"type": "corncob-parameter-dropdown", "group": MATCH}, "options"),
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ],
    [State({"type": "corncob-parameter-dropdown", "group": MATCH}, "value")])
def update_volcano_parameter_dropdown_options(selected_dataset, page, key, dummy):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ]
)
def update_volcano_comparison_dropdown(selected_dataset, page, key):
    options = [{'label': 'Estimated Coefficient', 'value': 'coef'}]
    value = "coef"

    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
    Output({"type": "corncob-parameter-dropdown", "group": MATCH}, "value"),
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ],
    [State({"type": "corncob-parameter-dropdown", "group": MATCH}, "options")])
def update_volcano_parameter_dropdown_value(selected_dataset, page, key, dummy):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    if fp is None:
        return "none"
    else:
        # Get the list of parameters
        parameter_list = valid_parameters(fp)

        if len(parameter_list) == 0:
            return "none"

        # Get the default parameter, if defined in the manifest
        default_parameter = page_data.default(
            selected_dataset, 
            "parameter", 
            key=key, 
            page=page
        )
        if default_parameter is not None and default_parameter in parameter_list:
            return default_parameter

        elif len(parameter_list) > 1:
            return parameter_list[1]
        else:
            return parameter_list[0]

@app.callback(
    Output({
        'type': 'corncob-pvalue-slider',
        'group': MATCH,
    }, "max"),
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
        Input({
            "type": "corncob-parameter-dropdown", 
            "group": MATCH}, 'value'),
    ])
def update_volcano_pvalue_slider_max(selected_dataset, page, key, parameter):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    # No dataset selected, or no parameter selected
    if fp is None or parameter == "none":
        return 1
    else:

        df = cag_associations(
            fp, 
            parameter
        )
        assert "neg_log10_pvalue" in df.columns.values, df.columns.values
        return cag_associations(
            fp, 
            parameter
        )["neg_log10_pvalue"].max()

@app.callback(
    Output({
        'type': 'corncob-pvalue-slider',
        'group': MATCH,
    }, "marks"),
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
        Input({
            "type": "corncob-parameter-dropdown", 
            "group": MATCH}, 'value'),
    ])
def update_volcano_pvalue_slider_marks(selected_dataset, page, key, parameter):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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

    return marks
#######################
# / VOLCANO CALLBACKS #
#######################



#########################
# CAG TAXONOMY CALLBACK #
#########################
@app.callback(
    Output('cag-tax-graph', 'figure'),
    [
        Input("plot-cag-selection-type", "value"),
        Input({"type": "corncob-parameter-dropdown", "group": "plot-cag"}, "value"),
        Input("plot-cag-annotation-ncags", "value"),
        Input("plot-cag-annotation-multiselector", "value"),
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
        Input('plot-cag-multiselector', 'value'),
    ])
def update_taxonomy_graph(
    selection_type,
    parameter,
    annotation_ncags,
    annotations,
    selected_dataset, 
    page,
    key,
    cag_id
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    if fp is None or cag_id is None:
        return empty_figure()

    # If a single CAG has been selected, add that to the plot
    if selection_type == "cag_id":
        cag_id = int(cag_id)

        # Format the DataFrame as needed to make a go.Sunburst
        cag_tax_df = taxonomic_gene_annotations(fp, cag_id)

        plot_title = "CAG {}".format(cag_id)

    else:

        # Pick the CAGs to plot based on their association with
        # the selected parameter, as well as their annotation with
        # the selected taxa and/or functions 
        cags_selected = select_cags_by_association_and_annotation(
            fp,
            parameter,
            annotations,
            annotation_ncags
        )
        
        if len(cags_selected) == 0:
            return empty_figure()

        # Format the DataFrame as needed to make a go.Sunburst
        cag_tax_df = pd.concat([
            taxonomic_gene_annotations(
                fp, cag_id
            )
            for cag_id in cags_selected
        ])
        cag_tax_df = cag_tax_df.groupby(
            ["name", "rank", "tax_id", "parent"]
        ).sum(
        ).reset_index(
        )

        if len(cags_selected) > 1:
            plot_title = "{:,} CAGs with annotations for {}".format(
                len(cags_selected),
                " / ".join([
                    a.split("-", 1)[1]
                    for a in annotations
                ])
            )
        else:
            plot_title = "CAG {} - with annotations for {}".format(
                int(cags_selected[0]),
                " / ".join([
                    a.split("-", 1)[1]
                    for a in annotations
                ])
            )

    if cag_tax_df is None:
        return empty_figure()
    else:
        return draw_taxonomy_sunburst(cag_tax_df, plot_title)
###########################
# / CAG TAXONOMY CALLBACK #
###########################

#####################
# PLOT CAG CALLBACK #
#####################
@app.callback(
    [
        Output("plot-cag-by-id-div", "style"),
        Output("plot-cag-by-association-div", "style"),
    ],
    [Input("plot-cag-selection-type", "value")]
)
def plot_cag_show_hide_selection_controls(
    selection_type
):
    """When the user selects 'Association & Annotation', show the appropriate controls."""
    if selection_type == "cag_id":
        return {"display": "block"}, {"display": "none"}
    else:
        return {"display": "none"}, {"display": "block"}

@app.callback(
    [
        Output("plot-cag-annotation-multiselector", "options"),
        Output("plot-cag-annotation-multiselector", "placeholder")
    ],
    [
        Input("plot-cag-selection-type", "value"),
        Input({"type": "corncob-parameter-dropdown", "group": "plot-cag"}, "value"),
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ]
)
def plot_cag_annotation_options(
    selection_type,
    parameter,
    selected_dataset,
    page,
    key,
):
    """When the user selects 'Association & Annotation', show the appropriate annotations."""
    options = [{"label": "None", "value": "None"}]
    if selection_type == "cag_id":
        return options, "Loading..."

    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    # Make sure that this parameter is valid
    if parameter not in valid_parameters(fp):
        return options, "Loading..."

    # Get the table of CAG associations
    corncob_df = cag_associations(fp, parameter)
    if corncob_df is None:
        return options, "Loading..."

    # Get the annotations for no more than 500 CAGs
    cags_checked = set([])
    cag_id_list = corncob_df.index.values

    # Record the number of times that each annotation is seen
    all_annot = defaultdict(int)

    # Iterate over each CAG
    for cag_id in cag_id_list:

        # Get the annotations for no more than 500 CAGs
        if len(cags_checked) >= 500:
            break

        # Include all of the taxa annotated for these CAGs
        for rank in ["family", "genus", "species"]:
            # Read in the taxonomic annotations for this CAG
            df = taxonomic_gene_annotations(
                fp, cag_id, rank=rank
            )
            # Add the labels
            if df is not None:
                cags_checked.add(cag_id)
                for label in df["name"].values:
                    # Record the label along with the category
                    all_annot[
                        "{}-{}".format(
                            rank,
                            label
                        )
                    ] += 1

        # Include all of the functional annotations for these CAGs
        df = functional_gene_annotations(fp, cag_id)
        if df is not None:
            cags_checked.add(cag_id)
            for label in df["label"].values:
                all_annot[
                    "functional-{}".format(
                        label
                    )
                ] += 1
    
    if len(all_annot) == 0:
        return options, "Loading..."

    # Sort the annotations by the number of CAGs found
    # Don't display the category in the 'label', but store it in the 'value'
    return [
            {
                "label": label.split("-", 1)[1], 
                "value": label
            }
        for label in pd.Series(
            all_annot
        ).sort_values(
            ascending=False
        ).index.values
    ], "Select..."


@app.callback(
    Output("plot-cag-multiselector", 'max'),
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ])
def update_single_cag_multiselector_options(
    selected_dataset, page, key
):
    """When a new dataset is selected, fill in maximum allowed CAG ID value."""

    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    if fp is None:
        return 1

    else:

        # Get the list of all CAGs
        cag_id_list = cag_annotations(fp).index.values

        return max(cag_id_list)
    
@app.callback(
    Output('plot-cag-multiselector', 'value'),
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ]
)
def update_single_cag_dropdown_value(
    selected_dataset, page, key
):
    # If a new dataset is selected, update the default CAG selected
    if isinstance(selected_dataset, list):
        selected_dataset = selected_dataset[0]
    selected_dataset = int(selected_dataset)

    # If the Main Menu button has been clicked, just return 0
    if selected_dataset == -1:
        return 0

    # Otherwise return the most abundant CAG
    else:
        
        # Get the path to the selected dataset
        fp = page_data.parse_fp(
            [selected_dataset], 
            page=page, 
            key=key
        )

        # See if there is a default parameter
        default_parameter = page_data.default(
            selected_dataset,
            "parameter",
            page=page,
            key=key
        )

        # Get the table of CAG annotations
        annot_df = cag_annotations(fp)

        if default_parameter is not None:

            # Get the table of CAG associations with this parameter
            assoc_df = cag_associations(fp, default_parameter)

            if assoc_df is not None:

                # Pick the top CAG based on the wald and size
                cag_scoring = assoc_df["abs_wald"] * annot_df["size_log10"]
                cag_scoring = cag_scoring.dropna().sort_values()
                if cag_scoring.shape[0] > 1:
                    return cag_scoring.index.values[-1]

        # Otherwise, pick the top CAG to display by mean abundance
        top_cag = annot_df[
            "mean_abundance"
        ].sort_values(
            ascending=False
        ).index.values[0]

        # With a new dataset, select the first five CAGs
        return top_cag

@app.callback(
    Output({'name': 'plot-cag-xaxis',
            "type": "metadata-field-dropdown"}, 'value'),
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ])
def update_single_cag_default_x(
    selected_dataset,
    page,
    key,
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    if fp is None:
        return "none"

    # See if there is a default value specified in the manifest
    default_value = page_data.default(
        selected_dataset, 
        "x",
        page=page,
        key=key,
    )
    if default_value is not None and default_value in manifest(fp).columns.values:
        return default_value
    else:
        return "none"

@app.callback(
    Output('plot-cag-plot-type', 'value'),
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ])
def update_single_cag_default_plot_type(
    selected_dataset,
    page,
    key,
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    if fp is None:
        return "scatter"

    # See if there is a default value specified in the manifest
    default_value = page_data.default(
        selected_dataset, 
        "plot_type",
        page=page,
        key=key,
    )

    if default_value is not None and default_value in ["scatter", "line", "boxplot", "strip"]:
        return default_value
    else:
        return "scatter"

@app.callback(
    Output({'name': 'plot-cag-color',
            "type": "metadata-field-dropdown"}, 'value'),
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ])
def update_single_cag_default_color(
    selected_dataset,
    page,
    key,
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    if fp is None:
        return "none"

    # See if there is a default value specified in the manifest
    default_value = page_data.default(
        selected_dataset, 
        "color",
        page=page,
        key=key,
    )
    if default_value is not None and default_value in manifest(fp).columns.values:
        return default_value
    else:
        return "none"

@app.callback(
    Output({'name': 'plot-cag-facet',
            "type": "metadata-field-dropdown"}, 'value'),
    [
        Input("selected-dataset", "children"),
        Input("login-username", "value"),
        Input("login-key", "value"),
    ])
def update_single_cag_default_facet(
    selected_dataset,
    page,
    key,
):
        # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    if fp is None:
        return "none"

    # See if there is a default value specified in the manifest
    default_value = page_data.default(
        selected_dataset, 
        "facet",
        page=page,
        key=key,
    )
    if default_value is not None and default_value in manifest(fp).columns.values:
        return default_value
    else:
        return "none"

@app.callback(
    Output('plot-cag-graph', 'figure'),
    [
        Input("plot-cag-selection-type", "value"),
        Input({"type": "corncob-parameter-dropdown", "group": "plot-cag"}, "value"),
        Input("plot-cag-annotation-ncags", "value"),
        Input("plot-cag-annotation-multiselector", "value"),
        Input('plot-cag-multiselector', 'value'),
        Input({'name': 'plot-cag-xaxis',
               "type": "metadata-field-dropdown"}, 'value'),
        Input('plot-cag-plot-type', 'value'),
        Input({'name': 'plot-cag-color',
               "type": "metadata-field-dropdown"}, 'value'),
        Input({'name': 'plot-cag-facet',
               "type": "metadata-field-dropdown"}, 'value'),
        Input('plot-cag-log', 'value'),
        Input('manifest-filtered', 'children'),
    ],[
        State("selected-dataset", "children"),
        State("login-username", "value"),
        State("login-key", "value"),
    ])
def update_single_cag_graph(
    selection_type,
    parameter,
    annotation_ncags,
    annotations,
    cag_id,
    xaxis,
    plot_type,
    color,
    facet,
    log_scale,
    manifest_json,
    selected_dataset,
    page,
    key,
):
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

    if fp is None:
        return empty_figure()

    # Get the filtered manifest from the browser
    plot_manifest_df = parse_manifest_json(
        manifest_json, 
        manifest(fp)
    )

    # If a single CAG has been selected, add that to the plot
    if selection_type == "cag_id":
        plot_df = plot_manifest_df.assign(
            CAG_ABUND = cag_abundances(fp, cag_id=cag_id)
        )
        plot_title = "CAG {}".format(cag_id)

    else:

        # Pick the CAGs to plot based on their association with
        # the selected parameter, as well as their annotation with
        # the selected taxa and/or functions 
        cags_selected = select_cags_by_association_and_annotation(
            fp,
            parameter,
            annotations,
            annotation_ncags
        )

        if len(cags_selected) == 0:
            return empty_figure()

        # Now that we have a list of CAGs, add the abundance and format the title
        plot_df = plot_manifest_df.assign(
            CAG_ABUND = pd.DataFrame({
                cag_id: cag_abundances(fp, cag_id=cag_id)
                for cag_id in cags_selected
            }).sum(axis=1)
        )

        if len(cags_selected) > 1:
            plot_title = "{:,} CAGs with annotations for {}".format(
                len(cags_selected),
                " / ".join([
                    a.split("-", 1)[1]
                    for a in annotations
                ])
            )
            axis_label = "Relative Abundance"
        else:
            plot_title = "CAG {} - with annotations for {}".format(
                int(cags_selected[0]),
                " / ".join([
                    a.split("-", 1)[1]
                    for a in annotations
                ])
            )

    axis_label = "Relative Abundance"

    return draw_single_cag_graph(
        plot_df,
        plot_title,
        axis_label,
        xaxis,
        plot_type,
        color,
        facet,
        log_scale
    )

def select_cags_by_association_and_annotation(
    fp,
    parameter,
    annotations,
    annotation_ncags
):
    # Return none if no annotations are selected
    if len(annotations) == 0:
        return []
    # Read the association metrics for each CAG against this parameter
    corncob_df = cag_associations(fp, parameter)
    if corncob_df is None:
        return []

    # Get the list of the top CAGs, most consistently associated with this paramter
    cag_id_list = corncob_df.index.values

    # The `annotations` is a list which looks something like:
    # ["genus-Escherichia", "functional-MatE"]
    # To parse this, we will select CAGs which match _any_ of the
    # taxonomic labels (if any), and _all_ of the functional labels (if any)

    taxon_labels = {
        rank: set([
            a.split("-", 1)[1]
            for a in annotations
            if a.startswith(rank + "-")
        ])
        for rank in ["species", "genus", "family"]
    }
    taxon_labels = {
        k: v
        for k, v in taxon_labels.items()
        if len(v) > 0
    }

    functional_labels = set([
        a.split("-", 1)[1]
        for a in annotations
        if a.startswith("functional-")
    ])
    # Make sure that all labels are present and accounted for
    assert len(annotations) == sum([len(v) for v in taxon_labels.values()]) + len(functional_labels)

    # Get the list of CAGs based on these criteria
    selected_cags = set([])

    # Iterate over each CAG
    for cag_id in cag_id_list:

        # Stop if we have reached the maximum number of CAGs
        if len(selected_cags) >= annotation_ncags:
            return list(selected_cags)

        # Iterate over each taxonomic level
        for rank, rank_annotations in taxon_labels.items():

            # If we have already added this CAG, skip it
            if cag_id in selected_cags:
                continue

            # Only read over the CAGs if there are annotations to look for
            if len(rank_annotations) > 0:

                # Read in the assignments at this level for this CAG
                rank_df = taxonomic_gene_annotations(
                    fp, cag_id, rank=rank
                )

                # If this CAG has no annotations, skip it
                if rank_df is None:
                    continue

                # Check if this CAG has the annotation of interest
                elif rank_df["name"].isin(rank_annotations).any():

                    # If a functional label was provided, see if that is present
                    if len(functional_labels) > 0:

                        df = functional_gene_annotations(
                            fp,
                            cag_id
                        )
                        if df is not None:
                            if df["label"].isin(functional_labels).any():
                                selected_cags.add(cag_id)

                    else:
                        selected_cags.add(cag_id)

        # If no taxonomic labels are provided, just check the functional labels
        if len(taxon_labels) == 0:

            # Get the functional labels
            df = functional_gene_annotations(
                fp,
                cag_id
            )
            if df is not None:
                if df["label"].isin(functional_labels).any():
                    selected_cags.add(cag_id)

    return list(selected_cags)
#######################
# / PLOT CAG CALLBACK #
#######################


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
        Input("login-username", "value"),
        Input("login-key", "value"),
        Input("manifest-table-bulk-select-apply", "n_clicks")
    ],
    [State("manifest-table-bulk-select-formula", "value")]
)
def update_manifest_table(
    selected_dataset, 
    page,
    key,
    bulk_select_button, 
    bulk_select_formula
):
    """Fill in the values of the manifest with the selected dataset."""

    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
        Input("login-username", "value"),
        Input("login-key", "value"),
        Input('manifest-rows-selected', 'children'),
    ])
def manifest_save_filtered(selected_dataset, page, key, selected_rows_json):
    """Save the filtered manifest to a hidden div.."""
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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
        Input("login-username", "value"),
        Input("login-key", "value"),
        Input('manifest-table-select-columns', 'value'),
    ])
def manifest_update_columns_selected(selected_dataset, page, key, selected_columns):
    """Allow the user to hide selected columns in the manifest table."""
    # Get the path to the selected dataset
    fp = page_data.parse_fp(selected_dataset, page=page, key=key)

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

# ADD THE GTM CONTAINER, IF PROVIDED
if gtm_container is not None and isinstance(gtm_container, str) and gtm_container.startswith("GTM-"):
    app.index_string = '''<!DOCTYPE html>
<html>
    <head>
    {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <script>(function(w,d,s,l,i){w[l]=w[l]||[];w[l].push({'gtm.start':new Date().getTime(),event:'gtm.js'});var f=d.getElementsByTagName(s)[0],j=d.createElement(s),dl=l!='dataLayer'?'&l='+l:'';j.async=true;j.src='https://www.googletagmanager.com/gtm.js?id='+i+dl;f.parentNode.insertBefore(j,f);})(window,document,'script','dataLayer','GTM-N4XKFR5');</script>
    </head>
    <body>
<noscript><iframe src="https://www.googletagmanager.com/ns.html?id=GTM-N4XKFR5"
height="0" width="0" style="display:none;visibility:hidden"></iframe></noscript>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
'''

if __name__ == '__main__':

    app.run_server(
        host='0.0.0.0',
        port=8050,
    )

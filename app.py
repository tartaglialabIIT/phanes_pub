"""PHANES – a Dash web app for querying and visualising a protein dataset."""

import dash
import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
from dash import Input, Output, State, callback, ctx, dash_table, dcc, html, no_update, callback_context, ALL
from dash.exceptions import PreventUpdate
import io
import base64
import datetime
import time
import re
import urllib.parse
import smtplib
from email.message import EmailMessage
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from pathlib import Path
from api_v6 import *

import base64
import os
import uuid
import subprocess
import threading
import json
import shutil

from dotenv import load_dotenv
from itsdangerous import URLSafeTimedSerializer, SignatureExpired, BadSignature
from flask import abort, send_file


# Local dev: .env in WorkingDirectory. Production: systemd EnvironmentFile
# (e.g. /etc/profet/profet-db-dash.env) — see profetdb.service.
load_dotenv()

SECRET_KEY = os.environ["SECRET_KEY"]
PROFETDB_BASE_URL = os.environ["PROFETDB_BASE_URL"].rstrip("/")
DOWNLOAD_TTL_DAYS = int(os.environ.get("DOWNLOAD_TTL_DAYS", "5"))
SMTP_PASSWORD = os.environ.get("SMTP_PASSWORD", "")
MAX_PDB = 10

# --- feature calculation -------------------------------------------------
# Conda env holding the physico-chemical / AlphaFold / graph calculator.
# calculate_all_features.sh calls this env "catg"; keep the two in sync.
FEATURE_ENV = os.environ.get("PROFET_FEATURE_ENV", "catgranule2")
#FEATURE_ENV = os.environ.get("PROFET_FEATURE_ENV", "catg")
# Conda env holding Thermometer (R + bio3d/igraph/reshape2/svMisc, biopython).
THERMO_ENV = os.environ.get("PROFET_THERMO_ENV", "thermo")
#THERMO_ENV = os.environ.get("PROFET_THERMO_ENV", "catg")
# Root of the Thermometer tool (scripts/, files/, bin/). Ships with the repo;
# override THERMO_HOME if it is installed elsewhere on the host.
THERMO_HOME = os.environ.get(
    "THERMO_HOME", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "thermometer"))
# Thermometer runs a NAMD minimisation per structure and is slow: two ~100-residue
# structures took ~9 min on 8 cores. Bound it, and bound how many jobs run at once
# (each job forks cpu//NAMD_PROCS workers, each of which forks NAMD).
PCAF_TIMEOUT = int(os.environ.get("PROFET_PCAF_TIMEOUT", str(2 * 3600)))
THERMO_TIMEOUT = int(os.environ.get("PROFET_THERMO_TIMEOUT", str(12 * 3600)))
MAX_CONCURRENT_JOBS = int(os.environ.get("PROFET_MAX_CONCURRENT_JOBS", "1"))
# The file users actually receive: scalars as CSV + per-residue profiles as a
# compressed pickle, zipped. Defined once so the in-app download, the emailed
# link and the "job finished" checks can never drift apart. cmb_df.pkl stays in
# the job folder as the intermediate.
RESULT_ARCHIVE = "phanes_results.zip"


def result_archive_path(folder):
    return os.path.join(folder, RESULT_ARCHIVE)
_job_slots = threading.Semaphore(MAX_CONCURRENT_JOBS)

from data import *
from descriptions import build_feature_descriptions, build_analysis_descriptions

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

serializer = URLSafeTimedSerializer(SECRET_KEY, salt='pdb-download')
numeric_df = generate_scalar_df()
numeric_df_scaled = generate_preproccesed_scalar_df()
numeric_df_scaled_noScaling = generate_preproccesed_scalar_noScaling_df()
profiles_df = generate_profiles_df()
mapping_df = generate_mapping_df()
gene_to_uniprot = (
    mapping_df.groupby("gene_names")["uniprot_id"]
    .apply(list)
    .to_dict()
)
uniprot_to_gene = (
    mapping_df.groupby("uniprot_id")["gene_names"]
    .apply(list)
    .to_dict()
)
valid_uniprot_ids = set(mapping_df["uniprot_id"])


UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
#"STRUCTURAL": ['n_contacts', 'RG_protein', 'rmsd_prot', 'Length', 'fullCharge', 'extCharge', 'Coil_FullSeq', 'AlphaHelix_FullSeq', 'BetaBridge_FullSeq', 'Strand_FullSeq', '3_10Helix_FullSeq', 'PiHelix_FullSeq', 'Turn_FullSeq', 'Bend_FullSeq', 'Coil_ExtSeq', 'AlphaHelix_ExtSeq', 'BetaBridge_ExtSeq', 'Strand_ExtSeq', '3_10Helix_ExtSeq', 'PiHelix_ExtSeq', 'Turn_ExtSeq', 'Bend_ExtSeq', 'asa_mean', 'asa_std', 'RBD_Full_min_0', 'RBD_Full_min_1', 'RBD_Full_min_2', 'RBD_int_min_0', 'RBD_int_min_1', 'RBD_int_min_2', 'RBD_ext_min_0', 'RBD_ext_min_1', 'RBD_ext_min_2', 'average_plddt', 'stddev_plddt', 'n_contacts_norm', 'RG_protein_norm'],
FEATURE_GROUPS = {
	"PHYSICAL/CHEMICAL" : ['BetaSheet_1', 'Aggregation_1', 'Disorder_1', 'Aggregation_2', 'Membrane_1', 'Hydrophobicity_1', 'BetaSheet_2', 'BetaSheet_3', 'AlphaHelix_1', 'Membrane_2', 'BetaSheet_4', 'Burial_1', 'Burial_2', 'Hydrophobicity_2', 'NucleicAcidBinding_1', 'AlphaHelix_2', 'BetaSheet_5', 'charge', 'LowFolding_1', 'Membrane_3', 'Turn_1', 'fg', 'BetaSheet_6', 'NucleicAcidBinding_2', 'BetaSheet_7', 'Disorder_2', 'Membrane_4', 'Burial_3', 'Hydrophobicity_3', 'BetaSheet_8', 'Burial_4', 'Disorder_3', 'Aggregation_3', 'Membrane_5', 'Hydrophobicity_4', 'Hydrophobicity_5', 'Hydrophobicity_6', 'Aggregation_4', 'Aggregation_5', 'Disorder_4', 'Aggregation_6', 'Hydrophobicity_7', 'Membrane_6', 'AlphaHelix_3', 'Disorder_5', 'Disorder_6', 'AlphaHelix_4', 'Burial_5', 'NucleicAcidBinding_3', 'AlphaHelix_5', 'NucleicAcidBinding_4', 'BetaSheet_9', 'Aggregation_7', 'Hydrophobicity_8', 'rg', 'Disorder_7', 'NucleicAcidBinding_5', 'Disorder_8', 'Disorder_9', 'Aggregation_8', 'Burial_6', 'AlphaHelix_6', 'NucleicAcidBinding_6', 'AlphaHelix_7', 'Aggregation_9', 'BetaSheet_10', 'Burial_7', 'Disorder_10', 'Hydrophobicity_9', 'AlphaHelix_8', 'Burial_8', 'AlphaHelix_9', 'NucleicAcidBinding_7', 'NucleicAcidBinding_8', 'NucleicAcidBinding_9', 'AlphaHelix_10', 'NucleicAcidBinding_10', 'Hydrophobicity_10', 'Membrane_7', 'Membrane_8', 'Membrane_9', 'Burial_9'],
	"STRUCTURAL": ['n_contacts', 'RG_protein', 'rmsd_prot', 'Length', 'fullCharge', 'extCharge', 'Coil_FullSeq', 'AlphaHelix_FullSeq', 'BetaBridge_FullSeq', 'Strand_FullSeq', '3_10Helix_FullSeq', 'PiHelix_FullSeq', 'Turn_FullSeq', 'Bend_FullSeq', 'Coil_ExtSeq', 'AlphaHelix_ExtSeq', 'BetaBridge_ExtSeq', 'Strand_ExtSeq', '3_10Helix_ExtSeq', 'PiHelix_ExtSeq', 'Turn_ExtSeq', 'Bend_ExtSeq', 'asa_mean', 'asa_std', 'RBD_Full_min_0', 'RBD_Full_min_1', 'RBD_Full_min_2', 'RBD_int_min_0', 'RBD_int_min_1', 'RBD_int_min_2', 'RBD_ext_min_0', 'RBD_ext_min_1', 'RBD_ext_min_2', 'average_plddt', 'stddev_plddt', 'n_contacts_norm', 'RG_protein_norm'],
	"GRAPH": ['number_of_nodes', 'number_of_edges', 'Density', 'size_of_largest_connected_component', 'avgerage_degree', 'degree_assortativity', 'betweeness_centrality', 'eigenvector_centrality', 'closeness_centrality', 'katz_centrality', 'harmonic_centrality', 'pagerank', 'hub_scores', 'authority_scores', 'local_clustering_coefficient', 'average_clustering_coefficient', 'global_transitivity', 'degree_centrality', 'average_shortest_path_length', 'diameter', 'entropy'],
	"THERMOSTABILITY": ['Ts_tot_PN', 'Ts_tot_N', 'Ts_tot_P', 'Ts_C_PN', 'Ts_C_N', 'Ts_C_P', 'Ts_lj_PN', 'Ts_lj_N', 'Ts_lj_P'],
	"EXPRESSION": ['mean_expression','tau']
}
FEATURE_GROUPS["ALL"] = (
    FEATURE_GROUPS["PHYSICAL/CHEMICAL"] +
    FEATURE_GROUPS["STRUCTURAL"] +
    FEATURE_GROUPS["GRAPH"] +
    FEATURE_GROUPS["THERMOSTABILITY"] +
    FEATURE_GROUPS["EXPRESSION"]
)
BASE_DIR = Path(__file__).resolve().parent
CONTENT = {
    "lnk-about": ("About",               BASE_DIR / "content" / "about.md"),
    "lnk-docs":  ("Documentation",       BASE_DIR / "content" / "Documentation.md"),
    "lnk-tut":   ("Tutorial",            BASE_DIR / "content" / "tutorial.md"),
    "lnk-feat":  ("Feature Description", BASE_DIR / "content" / "feature_description.md"),
}

length_list = ["Length", "number_of_nodes"]
FEATURE_GROUPS["ALL - Length"] = (
    list(set(FEATURE_GROUPS["ALL"]) - set(length_list))
)

PER_RESIDUE_FEATURES = ['Membrane_8', 'Aggregation_4', 'Aggregation_3', 'LowFolding_1', 'Aggregation_5', 'Aggregation_2', 'Aggregation_7', 'Aggregation_6', 'Aggregation_8', 'Aggregation_9', 'Aggregation_1', 'AlphaHelix_8', 'Turn_1', 'Disorder_8', 'BetaSheet_2', 'AlphaHelix_9', 'BetaSheet_10', 'Burial_1', 'Disorder_5', 'BetaSheet_8', 'charge', 'Disorder_9', 'Membrane_2', 'Disorder_10', 'Burial_4', 'fg', 'Hydrophobicity_8', 'Hydrophobicity_3', 'Hydrophobicity_10', 'Hydrophobicity_5', 'Hydrophobicity_6', 'Hydrophobicity_7', 'Hydrophobicity_4', 'Hydrophobicity_1', 'Hydrophobicity_9', 'Hydrophobicity_2', 'AlphaHelix_5', 'BetaSheet_6', 'Membrane_9', 'Membrane_4', 'Burial_2', 'Burial_8', 'Membrane_1', 'Membrane_5', 'AlphaHelix_3', 'AlphaHelix_7', 'AlphaHelix_6', 'rg', 'AlphaHelix_10', 'AlphaHelix_4', 'BetaSheet_3', 'BetaSheet_4', 'BetaSheet_7', 'BetaSheet_5', 'BetaSheet_1', 'Disorder_1', 'Disorder_2', 'AlphaHelix_1', 'Disorder_7', 'NucleicAcidBinding_2', 'NucleicAcidBinding_9', 'NucleicAcidBinding_8', 'NucleicAcidBinding_4', 'NucleicAcidBinding_1', 'NucleicAcidBinding_5', 'NucleicAcidBinding_7', 'NucleicAcidBinding_6', 'NucleicAcidBinding_10', 'NucleicAcidBinding_3', 'Burial_5', 'Burial_7', 'Burial_9', 'Burial_6', 'Burial_3', 'AlphaHelix_2', 'BetaSheet_9', 'Disorder_6', 'Disorder_3', 'Membrane_7', 'Membrane_3', 'Membrane_6', 'Disorder_4', 'local_clustering_coefficient', 'degrees', 'degree_centrality', 'betweeness_centrality', 'eigenvector_centrality', 'closeness_centrality', 'katz_centrality', 'harmonic_centrality', 'pagerank', 'hub_scores', 'authority_scores', 'eccentricity', 'triangles', 'eigenvalues']
all_options = [{"label":v, "value":v} for v in numeric_df.index] + [{"label":v, "value":v} for v in gene_to_uniprot.keys()]

FEATURE_DESCRIPTIONS = build_feature_descriptions(list(numeric_df.columns) + PER_RESIDUE_FEATURES)
ANALYSIS_DESCRIPTIONS = build_analysis_descriptions()

# ---------------------------------------------------------------------------
# App
# ---------------------------------------------------------------------------

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.FLATLY],
    title="PHANES",
    suppress_callback_exceptions=True,
    requests_pathname_prefix="/phanes/",
    routes_pathname_prefix="/phanes/",
)
server = app.server


@app.server.route("/phanes/download/<token>")
def download_token(token):
    try:
        data = serializer.loads(token, max_age=DOWNLOAD_TTL_DAYS * 24 * 3600)
    except SignatureExpired:
        abort(410, "This download link has expired.")
    except BadSignature:
        abort(403, "Invalid download link.")
        
    folder = folder_for_job(data["job_id"])
    file_path = result_archive_path(folder)

    if not os.path.exists(file_path):
        abort(404, "These results are no longer available")

    # Flask's send_file streams the file; no base64, no whole-file buffering.
    return send_file(file_path, as_attachment=True, download_name=RESULT_ARCHIVE)

# ---------------------------------------------------------------------------
# Layout helpers
# ---------------------------------------------------------------------------

def _make_input_card() -> dbc.Card:
    return dbc.Card(
        dbc.CardBody([
            html.H5("Selection Panel"),
            
            html.Label("Protein Selection"),
            dcc.Input(
                id="input-1a-text",
                type="text",
                placeholder="Paste UniprotIDs or Gene names (comma separated)",
                debounce=True,
                style={"width": "100%", "marginBottom": "10px"}
            ),

            dcc.Dropdown(
                id="input-1a",
                #options=[{"label":v, "value":v} for v in numeric_df.index] + [{"label":v, "value":v} for v in gene_to_uniprot.keys()],
                options=[],
                multi=True,
                placeholder="Or select proteins/genes...",
                #style={
                #    "maxHeight": "90px",
                #    "overflowY": "auto",
                #    "padding-right": "35px",
                #}
            ),
                      
            html.Label("Secondary Input. Check against"),
            dcc.Input(
                id="input-1b-text",
                type="text",
                placeholder="Paste UniprotIDs or Gene names (comma separated - Optional)",
                debounce=True,   
                style={"width": "100%", "marginBottom": "10px"}
            ),
            
            dcc.Dropdown(
                id="input-1b",
                #options=[{"label":v, "value":v} for v in numeric_df.index] + [{"label":v, "value":v} for v in gene_to_uniprot.keys()],
                options=[],
                multi=True,
                placeholder="Or select proteins/genes to check against... (Optional)",
                #style={
                #    "maxHeight": "90px",
                #    "overflowY": "auto",
                #}
            ),
            
            html.Label("Feature Selection"),
            dcc.Input(
                id="input-columns-text",
                type="text",
                placeholder="Aggregation_8, n_contacts, Membrane_2 (comma separated)",
                debounce=True,   
                style={"width": "100%", "marginBottom": "10px"}
            ),
            dcc.Dropdown(
                id="input-columns",
                options=[{"label": f"🔹 {k}", "value": f"__group__{k}"} for k in FEATURE_GROUPS.keys()]
                    +
                    [{"label":c, "value":c} for c in numeric_df.columns]
                ,
                multi=True,
                placeholder="Or select features... (Optional)"
            ),
        ])
    )
    
    
def _make_plot_buttons() -> dbc.Card():
    return dbc.Card(
        dbc.CardBody([
            html.H5("Select Analysis Plot"),
            
            dbc.Button("PCA", id="btn-pca", className="mb-2 w-100"),
            dbc.Button("UMAP", id="btn-umap", className="mb-2 w-100"),
            dbc.Button("SpiderPlot", id="btn-spider", className="mb-2 w-100"),
            dbc.InputGroup([
                        dbc.InputGroupText("Number of Features"),
                        dbc.Input(
                            id="input-spider-topn",
                            type="number",
                            min=1, max=10, step=1, value=5,
                        ),
                    ], className="mb-2"),
            
            dbc.Button("Heatmap", id="btn-heatmap", className="mb-2 w-100"),
            dbc.Button("Protein Similarity Network", id="btn-pin", className="mb-2 w-100"),
            dbc.InputGroup([
                        dbc.InputGroupText("Similarity Threshold"),
                        dbc.Input(
                            id="input-psn-thr",
                            type="number",
                            min=0.1, max=1.0, step=0.05, value=0.5,
                        ),
                    ], className="mb-2"),
            dbc.Button("ViolinPlot", id="btn-violin", className="mb-2 w-100"),
            dbc.Button("Protein Group Comparison", id="btn-pgc", className="mb-2 w-100"),
            dbc.InputGroup([
                        dbc.InputGroupText("Number of Features"),
                        dbc.Input(
                            id="input-pgc-topn",
                            type="number",
                            min=1, max=len(numeric_df.columns), step=1, value=10,
                        ),
                    ], className="mb-2"),
            dbc.Button("Profile Plot", id="btn-profile", className="mb-2 w-100"),
            dbc.Button("Feature Correlation", id="btn-corr", className="mb-2 w-100"),
            dbc.Tooltip("Right-click for info", target="btn-pca",     placement="right"),
            dbc.Tooltip("Right-click for info", target="btn-umap",    placement="right"),
            dbc.Tooltip("Right-click for info", target="btn-spider",  placement="right"),
            dbc.Tooltip("Right-click for info", target="btn-heatmap", placement="right"),
            dbc.Tooltip("Right-click for info", target="btn-pin",     placement="right"),
            dbc.Tooltip("Right-click for info", target="btn-violin",  placement="right"),
            dbc.Tooltip("Right-click for info", target="btn-pgc",     placement="right"),
            dbc.Tooltip("Right-click for info", target="btn-profile", placement="right"),
            dbc.Tooltip("Right-click for info", target="btn-corr",    placement="right"),
        ])
    )

def _make_download_buttons() -> dbc.Card():
    return dbc.Card(
        dbc.CardBody([
            html.H5("Download Data"),
            
            dbc.Button("Download Scalar Data", "btn-scl", className="mb-2 w-100"),
            dbc.Button("Download per-residue Data", "btn-prf", className="mb-2 w-100"),
            dcc.Download(id="dwn-all"),
        ])
    )
    
def _make_upload_pdb_interface() -> dbc.Card():
    return dbc.Card(
        dbc.CardBody([
            html.H5("Calculate Data"),
            dbc.Input(
                id="input-email",
                type="email",
                placeholder="Enter email (optional)",
                style={"marginBottom": "10px"}
            ),
            dcc.Upload(
                id="upload-pdb",
                children=html.Div(["Drag & Drop or ", html.A("Select PDB Files")]),
                multiple=True,
                accept=".pdb",
                style={
                    "width": "100%",
                    "height": "60px",
                    "lineHeight": "60px",
                    "borderWidth": "1px",
                    "borderStyle": "dashed",
                    "borderRadius": "5px",
                    "textAlign": "center",
                    "marginBottom": "10px",
                },
            ),
            dcc.Store(id="pdb-files-store", data=[]), # accumulate pdb files
            html.Div(id="upload-pdb-feedback", style={"marginBottom": "10px", "fontSize": "0.9rem"}),

            dbc.Button("Run Analysis", id="btn-run-pdb", color="primary", className="w-100"),
            dbc.Progress(
                id="pdb-progress",
                value=0,
                max=100,
                striped=True,
                animated=True,
                className="mt-2"
            ),

            html.Div(id="pdb-status", className="mt-2"),

            # Shown when a finished job is restored on page load. The automatic
            # download only fires once per job, so this is how you get the file
            # again -- and being a real click it also survives the browser
            # policies that block downloads with no user gesture behind them.
            dbc.Button(
                "Download results",
                id="btn-download-again",
                color="secondary",
                outline=True,
                size="sm",
                className="mt-2",
                style={"display": "none"},
            ),

            dcc.Download(id="download-pdb-result"),

            dcc.Store(id="pdb-job-id", storage_type="local"),   # track job
            # job id whose file has already been handed to this browser
            dcc.Store(id="pdb-delivered", storage_type="local"),
            dcc.Interval(id="resume-trigger", interval=500, max_intervals=1),
            dcc.Interval(id="pdb-poll", interval=2000, disabled=True),  # polling
        
        ])
        
    )
    
def _render_file_list(files, message=""):
    rows = []
    for f in files:
        rows.append(
            html.Div([
                html.Span(f["filename"], style={"flex": "1", "overflow": "hidden",
                                                    "textOverflow": "ellipsis"}),
                html.Button("X",
                            id={"type": "remove-pdb", "index": f["filename"]},
                            n_clicks=0,
                            style={"border": "none", "background": "none", "color": "#c0392b",
                                   "cursor": "pointer", "fontWeight": "bold"},
                ),
            ],
            style={"display": "flex", "alignItems": "center", "gap": "8px",
                   "borderBottom": "1px solid #eee", "padding": "2px 0"}
            )
        )
    children = []
    if message:
        children.append(html.Div(message, style={"color": "#c0392b", "marginBottom": "4px"}))
    if files:
        children.append(html.Small(f"{len(files)}/{MAX_PDB} file(s) selected:"))
        children.append(html.Div(rows, style={"marginTop": "4px"}))
    return children


# ---------------------------------------------------------------------------
# App layout
# ---------------------------------------------------------------------------
app.layout = dbc.Container(
    [
        dcc.Store(
            id="feature-description",
            data={
                "features": FEATURE_DESCRIPTIONS,
                "analyses": ANALYSIS_DESCRIPTIONS
            }
        ),
        html.Div(id="feature-info-dummy", style={"display": "none"}),
        dbc.Navbar(
            dbc.Container([
                dbc.NavbarBrand(
                    [html.I(className="bi bi-search me-2"), "PHANES"],
                    className="fs-4 fw-bold",
                ),
                dbc.Nav([
                    dbc.NavLink("About", id="lnk-about", href="#", className="small text-white-50 me-2"),
                    dbc.NavLink("Documentation", id="lnk-docs", href="#", className="small text-white-50 me-2"),
                    dbc.NavLink("Tutorial", id="lnk-tut", href="#", className="small text-white-50 me-2"),
                    dbc.NavLink("Feature Description", id="lnk-feat", href="#", className="small text-white-50 me-2"),
                    ],
                    className="ms-3", navbar=True,
                    )],
                fluid=True, className="d-flex align-items-center"
            ),
            color="primary",
            dark=True,
            className="mb-4",
        ),
        #dbc.Offcanvas(
        #    dcc.Markdown(id="doc-oc-body"),
        #    id="doc-oc",
        #    title="",
        #    placement="end",
        #    is_open=False,
        #    scrollable=True,
        #    backdrop=True,
        #    style={"width": "680px"},
        #),
        dbc.Modal([
            dbc.ModalHeader(id="doc-modal-title"),
            dbc.ModalBody(html.Div(id="doc-modal-body"))
            ],
            id="doc-modal", size="xl", scrollable=True, is_open=False),
        dbc.Row(
            [
                # ---- sidebar ----
                dbc.Col(
                    [_make_input_card(), _make_plot_buttons(), _make_download_buttons(), _make_upload_pdb_interface()],
                    md=3,
                    className="mb-4",
                ),
                # ---- main area ----
                dbc.Col(
                    [
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    dcc.Loading(
                                        html.Div(
                                                dcc.Graph(id="main-plot", 
                                                    style={"height": "1000px", "width": "100%", "minHeight": "1000px"},
                                                    #config={"scrollZoom": True}
                                                ),     
                                                id="plt-scroll-container",
                                                style={
                                                    "overflowY": "auto",
                                                    "overflowX": "auto",
                                                    "width": "100%",
                                                    "height": "1000px",
                                                }
                                        )
                                        
                                    ),
                                    html.Div([
                                        dbc.Button("Download Figure", id="btn-dwn-fig", color="primary"),
                                        dbc.Button("Fit to Screen", id="btn-autoscale", color="secondary", className="ms-2"),
                                        dbc.Button("Add selected to input", id="btn-add-selected", color="secondary", className="ms-2"),
                                        dcc.Download(id="dwn-fig")
                                    ]),
                                    html.Hr(),
                                    dash_table.DataTable(
                                        id="protein-table",
                                        columns=[{"name": c.replace("_", " ").title(), "id": c} for c in numeric_df.columns],
                                        page_size=15,
                                        sort_action="native",
                                        sort_mode="multi",
                                        filter_action="native",
                                        style_table={"overflowX": "auto"},
                                        style_header={
                                            "backgroundColor": "#2c3e50",
                                            "color": "white",
                                            "fontWeight": "bold",
                                        },
                                        style_cell={
                                            "padding": "8px 12px",
                                            "fontSize": "13px",
                                            "textAlign": "left",
                                        },
                                        style_data_conditional=[
                                            {
                                                "if": {"row_index": "odd"},
                                                "backgroundColor": "#f8f9fa",
                                            }
                                        ],
                                        tooltip_header={c: "Right click for info" for c in numeric_df.columns},
                                        tooltip_delay=200,
                                        tooltip_duration=None,
                                    ),
                                    #dbc.Button("Render Table", id="btn-rnd-tbl", color="primary"),
                                    #dbc.Button("Download Table", id="btn-dwn-tbl", color="primary"),
                                    #dcc.Download(id="dwn-tbl"),
                                    html.Div([
                                        dbc.Button("Render Table", id="btn-rnd-tbl", className="mb-2", color="primary")
                                    ]),
                                    html.Div([
                                        dbc.Button("Download Table", id="btn-dwn-tbl", className="mb-2", color="primary"),
                                        dcc.Download(id="dwn-tbl")
                                    ]),
                                ]
                            ),
                            className="shadow-sm mb-4",
                        ),
                    ],
                    md=9,
                ),
            ]
        ),
    ],
    fluid=True,
)

app.clientside_callback(
    "function(data){ window.PROFET_DESCRIPTIONS = data || {features:{},analyses:{}}; return ''; }",
    Output("feature-info-dummy", "children"),
    Input("feature-description", "data"),
)

# Auto-LaTeX equation image-links -> $...$  (LaTeX taken from the ?input= URL param,
# tolerating empty alt text and **bold** wrappers).
_EQ_IMG = re.compile(
    r'\[\s*\**\s*!\[.*?\]\[.*?\]\s*\**\s*\]'
    r'\((https?://[^\s]*?mathjax-viewer[^\s]*?)#\d+\)'
)
_MD_UNESCAPE = re.compile(r'\\([\\`*_{}\[\]()#+.!\-])')

def _to_math(m):
    mm = re.search(r'[?&]input=([^#&]*)', m.group(1))
    if not mm:
        return m.group(0)
    latex = _MD_UNESCAPE.sub(r'\1', urllib.parse.unquote(mm.group(1))).strip()
    return f"${latex}$" if latex else ""

# In-text Paperpile citations [(N)]/[(N,M)]/[(N-M)](paperpile.com/c/...) -> link each number to #ref-N
_CITELINK = re.compile(r'\[([^\]]*?)\]\(https?://paperpile\.com/c/[^)]*\)')
_NUM = re.compile(r'\d+')

def _fix_cite(m):
    label = m.group(1)
    if not _NUM.search(label):
        return label
    return _NUM.sub(lambda d: f'[{d.group(0)}](#ref-{d.group(0)})', label)

# Bibliography entries "N. [text ... Available from: DOI](paperpile.com/b/...)"
# -> add an #ref-N anchor, drop the broken wrapper, make the DOI clickable.
_REF = re.compile(r'^(\d+)\\?\.\s*\[(.*?)\]\(https?://paperpile\.com/b/[^)]*\)', re.M)
_AVAIL = re.compile(r'(Available from:\s*)(https?://\S+)')

def _fix_ref_md(m):
    n, txt = m.group(1), m.group(2)
    txt = _AVAIL.sub(lambda d: f'{d.group(1)}[{d.group(2)}]({d.group(2)})', txt)
    return n, f"{n}. {txt}"

def render_doc(path):
    """Body + references as Dash components. References are html.Div(id='ref-N') so the
    in-text [(N)](#ref-N) links have real DOM scroll targets — dcc.Markdown strips
    raw-HTML ids, so anchors must be real Dash components."""
    text = Path(path).read_text(encoding="utf-8")
    text = _EQ_IMG.sub(_to_math, text)      # equation images -> $...$
    text = _CITELINK.sub(_fix_cite, text)   # in-text citations -> [(N)](#ref-N)
    first = _REF.search(text)               # start of the bibliography (consecutive N. entries)
    if not first:
        return [dcc.Markdown(text, mathjax=True)]
    comps = [dcc.Markdown(text[:first.start()], mathjax=True)]   # body + "References" heading
    for m in _REF.finditer(text[first.start():]):
        n, entry = _fix_ref_md(m)
        comps.append(html.Div(dcc.Markdown(entry), id=f"ref-{n}",
                              style={"scrollMarginTop": "10px"}))
    return comps

def folder_for_job(job_id):
    return os.path.join(os.getcwd(), UPLOAD_FOLDER, job_id)
    
def make_download_token(job_id):
    return serializer.dumps({"job_id": job_id})

def sns_to_base64(g):
    buf = io.BytesIO()
    g.savefig(buf, format="png", bbox_inches="tight")
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode("utf-8")
    return f"data:image/png;base64,{encoded}"

def is_valid_email(email):
    if not email:
        return False
    pattern = r'^[\w\.-]+@[\w\.-]+\.\w+$'
    return re.match(pattern, email) is not None

def resolve_proteins(input_list):
    resolved = []
    if not input_list:
        return resolved
    for item in input_list:
        if not item:
            continue
        
        item = item.strip().upper()
        
        if item in numeric_df.index:
            resolved.append(item)
        if item in gene_to_uniprot.keys():
            resolved.extend(gene_to_uniprot[item])
            
    #return list(set(resolved))
    return list(dict.fromkeys(resolved))

def validate_inputs(input_1a, input_1b, columns):
    input_1a = input_1a or []
    input_1b = input_1b or []
    columns= columns or []
    
    valid_dimensionality = (len(columns) == 0 or len(columns) > 2)
    per_residue_cols = [c for c in columns if c in PER_RESIDUE_FEATURES]
    
    return {
        "pca": len(input_1a) >= 1 and valid_dimensionality,
        "umap": len(input_1a) >= 1 and valid_dimensionality,
        "spider": len(input_1a) >= 1,
        "heatmap": len(input_1a) >= 1 and len(columns) >= 1 and 'Length' not in columns,
        "pin": len(input_1a) >= 2 and 'Length' not in columns and (len(columns) == 0 or len(columns) >= MIN_SIMILARITY_FEATURES),
        #"pin": len(input_1a) >= 2,
        "violin": len(input_1a) >= 1 and len(columns) >= 1,
        #"pgc": len(input_1a) >= 2 and len(input_1b) >= 2,
        "pgc": len(input_1a) >= 2,
        "profile": (
            (len(input_1a) == 1 and 1 <= len(per_residue_cols) <= 12)
            or (len(per_residue_cols) == 1 and 1<= len(input_1a) <= 8)
        ),
        "render": len(input_1a) >= 1,
        "corr": len(columns) == 0 or len(columns) >= 2,
    }
    

def send_email_results(email, job_id):
    
    folder = folder_for_job(job_id)
    file_path = result_archive_path(folder)

    if not os.path.exists(file_path):
        return
     
    sender = "webservice@tartaglialab.com"
    recipient = email
    #link = f"download/{make_download_token(job_id)}"
    link = f"{PROFETDB_BASE_URL}/download/{make_download_token(job_id)}"
    body = ("Your PHANES submission has successfully finished.\n\n"
            f"Download your results within {DOWNLOAD_TTL_DAYS} days here \n"
            f"{link}\n\n"
            )
    msg = MIMEText(body, "plain")
    msg["Subject"] = "PHANES Feature Calculation Results"
    msg["From"] = sender
    msg["To"] = recipient
    
    # do this if we opt not to use link
    #with open(file_path, "rb") as f:
    #    msg.add_attachment(
    #        f.read(),
    #        maintype="application",
    #        subtype="octet-stream",
    #        filename="results.pkl"
    #    )
    
    with smtplib.SMTP("smtp.gmail.com", 587) as smtp:
        smtp.starttls()
        if SMTP_PASSWORD:
            smtp.login(sender, SMTP_PASSWORD)
        smtp.sendmail(sender, recipient, msg.as_string())
        
    

def _conda_bin():
    conda_exe = os.environ.get("CONDA_EXE", "")
    if os.path.isfile(conda_exe):
        return conda_exe
    return "/home/aarmaos-iit.local/miniconda3/bin/conda"


def _run_step(label, env_name, argv, cwd, folder, timeout, extra_env=None):
    """Run one calculation step in a conda env, appending its output to the
    job log. Returns (completed_process, combined_output)."""
    cmd = [_conda_bin(), "run", "-n", env_name] + argv
    run_env = dict(os.environ, **(extra_env or {}))
    try:
        proc = subprocess.run(
            cmd, cwd=cwd, env=run_env, capture_output=True, text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"{label} timed out after {timeout}s") from None

    output = (proc.stdout or "") + (proc.stderr or "")
    entry = ["", f"===== {label} (exit {proc.returncode}) =====",
             "$ " + " ".join(cmd), output]
    with open(os.path.join(folder, "feature_calculation.log"), "a") as log:
        log.write("\n".join(entry) + "\n")
    return proc, output


def run_feature_calculation(job_id, folder, file_paths, email=None):

    def job():
        try:
            save_progress(folder, 1, "Queued...")
            with _job_slots:
                save_progress(folder, 2, "Launching calculation...")
                code_dir = os.path.join(os.path.dirname(__file__), "code")

                # --- 1. physico-chemical + AlphaFold + graph  (0-50%) ------
                proc, output = _run_step(
                    "FeatureCalculator", FEATURE_ENV,
                    ["python", "FeatureCalculatorPcAfGraphs.py",
                     "-D", folder, "-O", folder,
                     "--progress-start", "0", "--progress-end", "50"],
                    code_dir, folder, PCAF_TIMEOUT,
                )
                if proc.returncode != 0:
                    raise RuntimeError(
                        f"FeatureCalculator failed (exit {proc.returncode}):"
                        f"\n{output}")

                # --- 2. thermostability  (50-92%) --------------------------
                # The driver owns its own work dir and removes it on success, so
                # nothing is left behind in the Thermometer tool directory.
                thermo_out = os.path.join(folder, "thermo")
                proc, output = _run_step(
                    "Thermometer", THERMO_ENV,
                    ["python", os.path.join(THERMO_HOME, "thermometer_new_parallel.py"),
                     "-D", folder, "-O", thermo_out,
                     "--progress-file", os.path.join(folder, "progress.json"),
                     "--progress-start", "50", "--progress-end", "92"],
                    THERMO_HOME, folder, THERMO_TIMEOUT,
                    extra_env={"THERMO_HOME": THERMO_HOME},
                )
                # Exit 1 means "some structures failed", not "the job failed":
                # those rows get NaN Ts_* values in step 3. Only a missing
                # RESULTS dir is fatal.
                if not os.path.isdir(os.path.join(thermo_out, "RESULTS")):
                    raise RuntimeError(
                        f"Thermometer produced no RESULTS (exit {proc.returncode}):"
                        f"\n{output}")

                # --- 3. merge thermostability into cmb_df.pkl  (95%) -------
                save_progress(folder, 95, "Combining results...")
                proc, output = _run_step(
                    "CombineThermo", FEATURE_ENV,
                    ["python", os.path.join(os.path.dirname(__file__),
                                            "combine_thermo_features.py"),
                     "-D", folder, "--thermo-dir", thermo_out],
                    os.path.dirname(__file__), folder, 600,
                )
                if proc.returncode != 0:
                    raise RuntimeError(
                        f"Merging thermostability features failed "
                        f"(exit {proc.returncode}):\n{output}")

                # --- 4. package into the file users download ---------------
                save_progress(folder, 98, "Packaging results...")
                proc, output = _run_step(
                    "PackageResults", FEATURE_ENV,
                    ["python", os.path.join(os.path.dirname(__file__),
                                            "package_results.py"),
                     "-D", folder],
                    os.path.dirname(__file__), folder, 900,
                )
                if proc.returncode != 0:
                    raise RuntimeError(
                        "Packaging results failed "
                        f"(exit {proc.returncode}): {output}")

            # mark completion
            save_progress(folder, 100, "Done")

            with open(os.path.join(folder, "done.txt"), "w") as f:
                f.write("done")

            if email:
                send_email_results(email, job_id)

        except Exception as e:
            save_progress(folder, 100, "Error")

            with open(os.path.join(folder, "error.txt"), "w") as f:
                f.write(str(e))

    threading.Thread(target=job).start()
    
def save_progress(folder, progress, status):
    with open(os.path.join(folder, "progress.json"), "w") as f:
        json.dump({"progress": progress, "status": status}, f)
# ---------------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------------
@callback(
    Output("pdb-files-store", "data"),
    Output("upload-pdb-feedback", "children", allow_duplicate=True),
    Input("upload-pdb", "contents"),
    Input({"type": "remove-pdb", "index": ALL}, "n_clicks"),
    State("upload-pdb", "filename"),
    State("pdb-files-store", "data"),
    prevent_initial_call=True,
)
def manage_pdb_files(contents, remove_clicks, filenames, store):
    store = store or []
    trigger = ctx.triggered_id
    message = ""
    
    if isinstance(trigger, dict) and trigger.get("type") == "remove-pdb":
        if not any(remove_clicks or []):
            raise PreventUpdate
        store = [f for f in store if f["filename"] != trigger["index"]]
        return store, _render_file_list(store)
        
    if trigger == "upload-pdb" and contents:
        existing = {f["filename"] for f in store}
        for content, name in zip(contents, filenames):
            if not name.lower().endswith(".pdb"):
                continue
            if name in existing:
                continue
            store.append({"filename": name, "content": content})
            existing.add(name)
            
        if len(store) > MAX_PDB:
            store = store[:MAX_PDB]
            message = f"Exceeds {MAX_PDB} - PDB limit"
            
        return store, _render_file_list(store, message)
        
    raise PreventUpdate

#@callback(
#    Output("doc-oc", "is_open"),
#    Output("doc-oc", "title"),
#    Output("doc-oc-body", "children"),
#    
#    Input("lnk-about", "n_clicks"),
#    Input("lnk-docs", "n_clicks"),
#    Input("lnk-tut", "n_clicks"),
#    Input("lnk-feat", "n_clicks"),
#    prevent_initial_call=True
#)
#def show_doc(*_):
#    key = ctx.triggered_id
#    if not key:
#        raise PreventUpdate
#    title, path = CONTENT[key]
#    text = Path(path).read_text(encoding="utf-8")
#    return True, title, text

@callback(
    Output("doc-modal", "is_open"),
    Output("doc-modal-title", "children"),
    Output("doc-modal-body", "children"),
    
    Input("lnk-about", "n_clicks"),
    Input("lnk-docs", "n_clicks"),
    Input("lnk-tut", "n_clicks"),
    Input("lnk-feat", "n_clicks"),
    prevent_initial_call=True
)
def show_doc(*_):
    key = ctx.triggered_id
    if not key:
        raise PreventUpdate
    title, path = CONTENT[key]
    #text = Path(path).read_text(encoding="utf-8")
    return True, title, render_doc(path)


@callback(
    Output("main-plot", "figure", allow_duplicate=True),
    Input("main-plot", "relayoutData"),
    State("main-plot", "figure"),
    prevent_initial_call=True
)
def sync_dendro_on_zoom(relayout_data, figure):
    
    if not relayout_data or not figure:
        raise PreventUpdate

    # Heatmap is xaxis3/yaxis3 (row=2, col=2 in a 2x2 grid → axes 3 and 4)
    # You may need to verify the axis keys by printing relayout_data when you zoom
    x_start = relayout_data.get("xaxis4.range[0]")
    x_end   = relayout_data.get("xaxis4.range[1]")
    y_start = relayout_data.get("yaxis4.range[0]")
    y_end   = relayout_data.get("yaxis4.range[1]")

    layout = figure["layout"]
    
    if relayout_data.get("xaxis4.autorange") or relayout_data.get("autosize"):
        heatmap_trace = next(t for t in figure["data"] if t.get("type") == "heatmap")
        n_cols = len(heatmap_trace["x"])
        n_rows = len(heatmap_trace["y"])
        layout["xaxis2"]["range"] = [-0.5, n_cols - 0.5]
        layout["xaxis2"]["autorange"] = False
        layout["yaxis3"]["range"] = [-0.5, n_rows - 0.5]
        layout["yaxis3"]["autorange"] = False
        figure["layout"] = layout
        return figure

    # sync col dendrogram x-axis (xaxis2) with heatmap x zoom
    if x_start is not None and x_end is not None:
        layout["xaxis2"]["range"] = [x_start, x_end]
        layout["xaxis2"]["autorange"] = False

    # sync row dendrogram y-axis (yaxis4) with heatmap y zoom
    if y_start is not None and y_end is not None:
        layout["yaxis3"]["range"] = [y_start, y_end]
        layout["yaxis3"]["autorange"] = False

    # handle reset (double-click)
    if "xaxis3.autorange" in relayout_data:
        n_cols = len(figure["data"][-1]["x"])
        n_rows = len(figure["data"][-1]["y"])
        layout["xaxis2"]["range"] = [-0.5, n_cols - 0.5]
        layout["yaxis3"]["range"] = [-0.5, n_rows - 0.5]

    figure["layout"] = layout
    return figure

# ---------------------------------------------------------------------------
# Autoscale plot
# ---------------------------------------------------------------------------
@callback(
    Output("main-plot", "figure", allow_duplicate=True),
    Input("btn-autoscale", "n_clicks"),
    State("main-plot", "figure"),
    prevent_initial_call=True
)
def autoscale_plot(n, figure):
    if not figure:
        raise PreventUpdate
        
    layout = figure["layout"]
        
    layout["height"] = 1000
    layout["width"] = None   # reset to container width
    layout["autosize"] = True
    
    heatmap_trace = next((t for t in figure["data"] if t.get("type") == "heatmap"), None)
    if figure.get("layout", {}).get("meta", {}).get("plot_info", {}) == "corr":
        return figure
    if heatmap_trace:
        n_cols = len(heatmap_trace["x"])
        n_rows = len(heatmap_trace["y"])
        
        layout["xaxis4"]["autorange"] = True
        layout["yaxis4"]["autorange"] = True
        
        layout["xaxis2"]["range"] = [-0.5, n_cols - 0.5]
        layout["xaxis2"]["autorange"] = False
        
        layout["yaxis3"]["range"] = [-0.5, n_rows - 0.5]
        layout["yaxis3"]["autorange"] = False
    
    figure["layout"] = layout
    return figure

# ---------------------------------------------------------------------------
# Download Table Contents
# ---------------------------------------------------------------------------
@callback(
    Output("dwn-tbl", "data"),
    Input("btn-dwn-tbl", "n_clicks"),
    State("protein-table", "data"),
    prevent_initial_call=True
)
def download_table(n, table_data):
    
    if not table_data:
        return dcc.send_string("No data available", "empty.txt")
        
    df = pd.DataFrame(table_data)
    
    return dcc.send_data_frame(
        df.to_csv,
        f"profet_data_table{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
        index="Uniprot ID"
    )

# ---------------------------------------------------------------------------
# Download Main Figure
# ---------------------------------------------------------------------------
@callback(
    Output("dwn-fig", "data"),
    
    Input("btn-dwn-fig", "n_clicks"),
    
    State("main-plot", "figure"),
    prevent_initial_call=True
)
def download_current_figure(n, fig_dic):
    if not fig_dic:
        return dash.no_update

    #fig = go.Figure(fig_dic)
    filename = "profet_figure.html"
    html_str = pio.to_html(fig_dic, full_html=True, include_plotlyjs="cdn")
    
    return dcc.send_string(html_str, filename)

# ---------------------------------------------------------------------------
# Download Data
# ---------------------------------------------------------------------------
@callback(
    Output("dwn-all", "data"),
    
    Input("btn-scl", "n_clicks"),
    Input("btn-prf", "n_clicks"),
    
    prevent_initial_call=True
)
def download_pkl(n_scalar, n_profiles):
    trigger_id = ctx.triggered_id
    
    
    if trigger_id == "btn-scl":
        buffer = io.StringIO()
        numeric_df.to_csv(buffer, index=True)
        filename = "profet_numeric_v6.csv"
        return dcc.send_string(buffer.getvalue(), filename)
        
    elif trigger_id == "btn-prf":
        buffer = io.BytesIO()
        profiles_df.to_pickle(buffer)
        filename = "profet_per_residue_v6.pkl"
        buffer.seek(0)
        return dcc.send_bytes(buffer.read(), filename)     
    
    return dash.no_update
        

# ---------------------------------------------------------------------------
# INPUT 1A
# ---------------------------------------------------------------------------

@callback(
    Output("input-1a", "options"),
    Input("input-1a", "search_value"),
    State("input-1a", "value"),
    prevent_initial_call=True
)
def search_1a(search, selected):
    selected = selected or []
    selected_opts = [{"label": v, "value": v} for v in selected]
    if not search or len(search) < 2:
        return selected_opts
    matches = [o for o in all_options if search.lower() in o["label"].lower()]
    seen = {o["value"] for o in matches}
    return matches + [o for o in selected_opts if o["value"] not in seen]


@callback(
    Output("input-1a", "value"),
    Output("input-1a", "options", allow_duplicate=True),
    Output("input-1a-text", "value"),
    Input("input-1a-text", "value"),
    Input("input-1a", "value"),
    State("input-1a", "options"),
    prevent_initial_call=True
)
def update_ids_1a(text_value, dropdown_value, current_options):
    trigger = ctx.triggered_id
    current_options = current_options or []

    if trigger == "input-1a":
        dropdown_value = dropdown_value or []
        resolved = resolve_proteins(dropdown_value)

        # stop the re-trigger loop
        if set(resolved) == set(dropdown_value):
            return no_update, no_update, no_update

        # inject ALL resolved ids into options before setting value
        # this is critical — Dash drops any value not in options
        opts = {o["value"] for o in current_options}
        new_opts = list(current_options) + [
            {"label": v, "value": v} for v in resolved if v not in opts
        ]
        return resolved, new_opts, no_update

    if trigger == "input-1a-text":
        if not text_value:
            return no_update, no_update, no_update
        pasted_ids = [s.strip().upper() for s in text_value.split(",") if s.strip()]
        resolved = resolve_proteins(pasted_ids)
        updated = list(dict.fromkeys((dropdown_value or []) + resolved))
        opts = {o["value"] for o in current_options}
        new_opts = list(current_options) + [
            {"label": v, "value": v} for v in updated if v not in opts
        ]
        return updated, new_opts, ""

    return no_update, no_update, no_update
    
@callback(
    Output("input-1b", "options"),
    Input("input-1b", "search_value"),
    State("input-1b", "value"),
    prevent_initial_call=True
)
def search_1b(search, selected):
    selected = selected or []
    selected_opts = [{"label": v, "value": v} for v in selected]
    if not search or len(search) < 2:
        return selected_opts
        #return all_options
        #return []
    matches = [o for o in all_options if search.lower() in o["label"].lower()]
    seen = {o["value"] for o in matches}
    return matches + [o for o in selected_opts if o["value"] not in seen]

@callback(
    Output("input-1b", "value"),
    Output("input-1b", "options", allow_duplicate=True),
    Output("input-1b-text", "value"),
    Input("input-1b-text", "value"),
    Input("input-1b", "value"),
    State("input-1b", "options"),
    prevent_initial_call=True
)
def update_ids_1b(text_value, dropdown_value, current_options):

    trigger = ctx.triggered_id

    # --- Dropdown selection ---
    if trigger == "input-1b":
        dropdown_value = dropdown_value or []
        resolved = resolve_proteins(dropdown_value)
        
        if set(resolved) == set(dropdown_value):
            return no_update, no_update, no_update
            
        opts = {o["value"] for o in current_options}
        new_opts = list(current_options) + [
            {"label": v, "value": v} for v in resolved if v not in opts
        ]
        return resolved, new_opts, no_update

    # --- Text paste ---
    if trigger == "input-1b-text":
        if not text_value:
            return no_update, no_update, no_update

        pasted_ids = [s.strip().upper() for s in text_value.split(",") if s.strip()]
        resolved = resolve_proteins(pasted_ids)
        updated = list(dict.fromkeys((dropdown_value or []) + resolved))
        opts = {o["value"] for o in current_options}
        new_opts = list(current_options) + [
            {"label": v, "value": v} for v in updated if v not in opts
        ]
       
        return updated, new_opts, ""

    return no_update, no_update, no_update
    
    
@callback(
   Output("input-columns", "value"),
   Output("input-columns-text", "value"),
   Input("input-columns-text", "value"),
   Input("input-columns", "value"),
   prevent_initial_call=True
)
def update_ids(text_value, dropdown_value):

    trigger = ctx.triggered_id
    
    all_feats = set(FEATURE_GROUPS["ALL"])
    all_nonlength = set(FEATURE_GROUPS["ALL - Length"])

    # --- Dropdown selection ---
    if trigger == "input-columns":
        dropdown_value = dropdown_value or []
        
        group_keys = [v for v in dropdown_value if isinstance(v,str) and v.startswith('__group__')]
        individual_feats = [v for v in dropdown_value if isinstance(v,str) and not v.startswith('__group__')]
        if not group_keys:
            return individual_feats, no_update
            
        current_set = set(individual_feats)
        last = group_keys[-1].replace("__group__", "")
        new_group_feats = FEATURE_GROUPS.get(last, [])
        
        if current_set == all_feats or current_set == all_nonlength:
            return list(dict.fromkeys(new_group_feats)), no_update
        else:
            return list(dict.fromkeys(individual_feats + new_group_feats)), no_update

    # --- Text paste ---
    if trigger == "input-columns-text":
        if not text_value:
            return no_update, no_update
            
        
        search_list = [v.strip() for v in text_value.split(",") if v.strip()]
        cols_lower = {c.lower(): c for c in numeric_df.columns}
        valid = [cols_lower[v.lower()] for v in search_list if v.lower() in cols_lower]
        current = [v for v in (dropdown_value or []) if not (isinstance(v, str) and v.startswith('__group__'))]
        updated = list(dict.fromkeys(current + valid))

        return updated, ""

    return no_update, no_update


@callback(
    Output("input-columns", "value", allow_duplicate=True),
    Input("main-plot", "clickData"),
    State("input-columns", "value"),
    State("main-plot", "figure"),
    prevent_initial_call=True
)
def add_feature_from_plot(click_data, current_cols, fig):
    if not click_data or not fig:
        return no_update
    
    info = fig.get("layout", {}).get("meta", {}).get("plot_info", {})
    
    if info != "pgc":
        return no_update
    
    point = click_data["points"][0]
    feat = point.get("text")
    
    if not feat:
        return no_update
        
    current = current_cols or []
    if feat in current:
        return no_update
        
    return list(dict.fromkeys(current + [feat]))
    
@callback(
    Output("input-1a", "value", allow_duplicate=True),
    Output("input-1a", "options", allow_duplicate=True),
    Input("main-plot", "clickData"),
    State("input-1a", "value"),
    State("main-plot", "figure"),
    prevent_initial_call=True
)
def add_id_from_plot(click_data, current_ids, fig):
    if not click_data or not fig:
        return no_update, no_update
        
    info = fig.get("layout", {}).get("meta", {}).get("plot_info", {})
    selected_trace = fig["data"][-1]
    print(selected_trace.get("hovertext", []))
    if info != "pca" and info != "umap":
        return no_update, no_update
    
    point = click_data["points"][0]
    id = point.get("hovertext")
    
    if not id:
        return no_update, no_update
        
    current = current_ids or []
    
    if id in current:
        return no_update, no_update
    opts = [{"label": v, "value": v} for v in current + [id]]    
    return list(dict.fromkeys(current + [id])), opts
    
@callback(
    Output("input-1a", "value", allow_duplicate=True),
    Output("input-1a", "options", allow_duplicate=True),
    Output("input-columns", "value", allow_duplicate=True),
    Input("btn-add-selected", "n_clicks"),
    State("main-plot", "relayoutData"),
    State("main-plot", "figure"),
    State("input-1a", "value"),
    prevent_initial_call=True
)
def add_plot_elements_to_input_1a(n, relayout_data, fig, current_ids):
    trigger = ctx.triggered_id
    
    #if trigger != "btn-add-selected":
    #    return no_update, no_update
     
    print(relayout_data)
    info = fig.get("layout", {}).get("meta", {}).get("plot_info", "")
    x_min = relayout_data.get("xaxis.range[0]", None)
    x_max = relayout_data.get("xaxis.range[1]", None)
    y_min = relayout_data.get("yaxis.range[0]", None)
    y_max = relayout_data.get("yaxis.range[1]", None)
    zoomed_in = all(v is not None for v in [x_min, x_max, y_min, y_max])
    
    selected_points_trace = fig["data"][-1]
    print(selected_points_trace)
    xs = selected_points_trace.get("x", [])
    shape = selected_points_trace.get("x", []).get("_inputArray", {}).get("shape", 0)
    xs = []
    ys = []
    
    for i in range(int(shape)):
        xs.append(float(selected_points_trace.get("x", []).get("_inputArray", {}).get(f"{i}", 0)))
        ys.append(float(selected_points_trace.get("y", []).get("_inputArray", {}).get(f"{i}", 0)))
    print(f'zoomed in : {zoomed_in}')
    print(f'info : {info}')
    ids = selected_points_trace.get("text", []) if info == "pgc" else selected_points_trace.get("hovertext", [])
    print(f'ids before zoom: {ids}')
    
    # ids that are in the plot area
    if zoomed_in:
        ids = [
            id for x, y, id in zip(xs, ys, ids)
            if x_min <= x <= x_max and y_min <= y <= y_max
        ]

    print(ids)
    if not ids:
        return no_update, no_update, no_update
    
    if info == "pca" or info == "umap":
        opts = [{"label": v, "value": v} for v in ids]
        return list(dict.fromkeys(ids)), opts, no_update
        
    if info == "pgc":
        to_return = [{"label": v, "value":v} for v in ids]
        return no_update, no_update, list(dict.fromkeys(ids)) 


@callback(
    Output("main-plot", "figure"),
    Output("main-plot", "style"),
    #Output("heatmap-img", "src"),
    #Output("heatmap-img", "style"),
    Output("protein-table", "data"),
    Output("protein-table", "columns"),
    Input("btn-pca", "n_clicks"),
    Input("btn-umap", "n_clicks"),
    Input("btn-spider", "n_clicks"),
    Input("btn-heatmap", "n_clicks"),
    Input("btn-pin", "n_clicks"),
    Input("btn-violin", "n_clicks"),
    Input("btn-pgc", "n_clicks"),
    Input("btn-profile", "n_clicks"), 
    Input("btn-rnd-tbl", "n_clicks"),
    Input("btn-corr", "n_clicks"),
    State("input-1a", "value"),
    State("input-1b", "value"),
    State("input-columns", "value"),
    State("input-spider-topn", "value"),
    State("input-psn-thr", "value"),
    State("input-pgc-topn", "value"),
    prevent_initial_call=True,
)
def render_main_plot(pca, umap, spider, heatmap, pin, violin, pgc, profile, render, corr,
    input_1a, input_1b, columns, spider_topn, psn_threshold, pgc_topn):

    trigger = ctx.triggered_id
    
    ids1 = input_1a
    ids2 = input_1b
    feats = columns
    
    if trigger == "btn-corr":
        return plot_correlation(numeric_df, feats), {"display": "block"}, [], []
    
    selected_proteins = ids1 + ids2 if ids2 else ids1  # combine lists if you want both
    selected_proteins = [p for p in selected_proteins if p in numeric_df.index]

    table_df = numeric_df.loc[selected_proteins, feats] if feats else numeric_df.loc[selected_proteins]
    table_df = table_df.reset_index()
    table_df.rename(columns={'index': "Uniprot ID"}, inplace=True)
    cols = [{"name": "Uniprot ID", "id": "Uniprot ID"}] + [{"name": c.replace("_", " ").title(), "id": c} for c in table_df.columns if "Uniprot ID" not in c]
    
    valid = validate_inputs(input_1a, input_1b, columns)
        
    if trigger == "btn-rnd-tbl" and valid["render"]:
        return {"data": [], "layout": {"title": "Invalid Selection"}}, {"display": "block"}, table_df.to_dict("records"), cols
    
    if trigger == "btn-pca" and valid["pca"]:
        return plot_pca(numeric_df_scaled, ids1, ids2, feats), {"display": "block"}, table_df.to_dict("records"), cols
    
    elif trigger == "btn-umap" and valid["umap"]:
        return plot_umap(numeric_df_scaled, ids1, ids2, feats), {"display": "block"}, table_df.to_dict("records"), cols
        
    elif trigger == "btn-spider" and valid["spider"]:
        top_n = int(spider_topn) if spider_topn is not None else 5
        if feats:
            return plot_spider_plot_features(numeric_df_scaled_noScaling, ids1, feats[:top_n], ids2), {"display": "block"}, table_df.to_dict("records"), cols
        else:
            return plot_spider_plot(numeric_df_scaled_noScaling, ids1, ids2, top_n=top_n), {"display": "block"}, table_df.to_dict("records"), cols
        
    elif trigger == "btn-heatmap" and valid["heatmap"]:
        return plot_heatmap_new(numeric_df_scaled, ids1,feats), {"display": "block"}, table_df.to_dict("records"), cols
        #img_src = plot_heatmap_with_dendrogram(numeric_df_scaled, ids1,feats)
        ##img_src = sns_to_base64(g)
        #return {}, {"display": "none"}, img_src, {"display": "block"}, table_df.to_dict("records"), cols
        
    elif trigger == "btn-pin" and valid["pin"]:
        thr = float(psn_threshold) if psn_threshold is not None else 0.5
        thr = max(0.1, min(1.0, thr))
        # Colour-code only when a second protein set (input 2) is supplied; otherwise single colour.
        return plot_protein_similarity_network(numeric_df_scaled, ids1, feats, thr, ids2=ids2), {"display": "block"}, table_df.to_dict("records"), cols
        
    elif trigger == "btn-violin" and valid["violin"]:
        return plot_violinplot_w_selection(numeric_df, feats[0], ids1), {"display": "block"}, table_df.to_dict("records"), cols
        
    elif trigger == "btn-pgc" and valid["pgc"]:
        print(pgc_topn)
        pgc_top_n = int(pgc_topn) if pgc_topn is not None else 10
        temp_df = transform_features(numeric_df)
        return plot_feature_difference(numeric_df, ids1, ids2, top_n=pgc_top_n), {"display": "block"}, table_df.to_dict("records"), cols
        
    elif trigger == "btn-profile" and valid["profile"]:
        profile_feats = [f + '_profile' for f in feats if f + '_profile' in profiles_df.columns]
        if len(ids1) == 1:
            return plot_profile_one_seq_multiple_features(profiles_df, ids1[0], profile_feats), {"display": "block"}, table_df.to_dict("records"), cols
        else:
            return plot_profile_multiple_seq_one_feature(profiles_df, ids1, profile_feats[0]), {"display": "block"}, table_df.to_dict("records"), cols
    
    return {"data": [], "layout": {"title": "Invalid Selection"}}, {"display": "block"}, [], []
    
    
       
    
@callback(
    Output("btn-pca", "disabled"),
    Output("btn-umap", "disabled"),
    Output("btn-spider", "disabled"),
    Output("btn-heatmap", "disabled"),
    Output("btn-pin", "disabled"),
    Output("btn-violin", "disabled"),
    Output("btn-pgc", "disabled"),
    Output("btn-profile", "disabled"),
    Output("btn-rnd-tbl", "disabled"),
    Output("btn-corr", "disabled"),

    Input("input-1a", "value"),
    Input("input-1b", "value"),
    Input("input-columns", "value"),
)
def toggle_buttons(input_1a, input_1b, columns):
    v = validate_inputs(input_1a, input_1b, columns)

    return tuple(not v[k] for k in ["pca","umap","spider","heatmap","pin","violin","pgc","profile", "render", "corr"])

  
@callback(
    Output("pdb-status", "children"),
    Output("pdb-progress", "value"),
    Output("pdb-progress", "label"),
    Output("pdb-progress", "style"),
    Output("pdb-poll", "disabled"),   
    Output("pdb-job-id", "data"),

    Input("btn-run-pdb", "n_clicks"),
    Input("pdb-poll", "n_intervals"),

    #State("upload-pdb", "contents"), contents
    #State("upload-pdb", "filename"), filenames
    State("pdb-files-store", "data"),
    State("pdb-job-id", "data"),
    State("input-email", "value"),

    prevent_initial_call=True
)
def handle_pdb(n_clicks, n_intervals, store, job_id, email):

    trigger = ctx.triggered_id
    SHOW = {"display": "block"}
    HIDE = {"display": "none"}

    # -------------------------
    # 1. USER CLICKED BUTTON
    # -------------------------
    if trigger == "btn-run-pdb":
        
        if not store:
            return "No files uploaded", 0, "", HIDE, True, None

        #if not contents:
        #    return "No files uploaded", 0, "", True, None

        #if len(contents) > 10:
        #    return "Max 10 files allowed", 0, "", True, None
            
        if email and not is_valid_email(email):
            return "Invalid email address", 0, "", HIDE, True, None 

        job_id = str(uuid.uuid4())
        folder = os.path.join(os.getcwd(), UPLOAD_FOLDER, job_id)
        os.makedirs(folder, exist_ok=True)

        file_paths = []
        for f in store:
            data = base64.b64decode(f["content"].split(",")[1])
            path = os.path.join(folder, f["filename"])
            
            with open(path, "wb") as fh:
                fh.write(data)
                
            file_paths.append(path)
        #for content, name in zip(contents, filenames):
        #    data = base64.b64decode(content.split(",")[1])
        #    path = os.path.join(folder, name)
        #
        #    with open(path, "wb") as f:
        #        f.write(data)
        #
        #    file_paths.append(path)

        save_progress(folder, 0, "Starting...")
        run_feature_calculation(job_id, folder, file_paths, email)
        
        if email:
            with open(os.path.join(folder, "email.txt"), "w") as f:
                f.write(email)
            return "Analysis submitted. When it's finished you will be notified via email.", 0, "", HIDE, True, job_id

        return "Processing started...", 0, "0%", SHOW, False, job_id

    # -------------------------
    # 2. POLLING
    # -------------------------
    if trigger == "pdb-poll":

        if not job_id:
            return "No job", 0, "", HIDE, True, job_id

        folder = os.path.join(os.getcwd(), UPLOAD_FOLDER, job_id)

        done_file = os.path.join(folder, "done.txt")
        error_file = os.path.join(folder, "error.txt")
        progress_file = os.path.join(folder, "progress.json")

        if os.path.exists(error_file):
            return "❌ Error occurred", 100, "Error", SHOW, True, job_id

        if os.path.exists(done_file):
            return "✅ Done!", 0, "0%", SHOW, True, job_id
            #return "✅ Done!", 100, "100%", True, job_id

        if os.path.exists(progress_file):
            with open(progress_file) as f:
                data = json.load(f)

            progress = data.get("progress", 0)
            status = data.get("status", "")

            return status, progress, f"{progress}%", SHOW, False, job_id

        return "Starting...", 0, "0%", SHOW, False, job_id

    return (dash.no_update,) * 6
    
@callback(
    Output("pdb-status", "children", allow_duplicate=True),
    Output("pdb-progress", "value", allow_duplicate=True),
    Output("pdb-progress", "label", allow_duplicate=True),
    Output("pdb-progress", "style", allow_duplicate=True),
    Output("pdb-poll", "disabled", allow_duplicate=True),
    Input("resume-trigger", "n_intervals"),
    State("pdb-job-id", "data"),
    prevent_initial_call=True,
)
def resume_job_on_load(_, job_id):
    HIDE = {"display": "none"}
    SHOW = {"display": "block"}
    # On page (re)load, the local-storage Store rehydrates job_id. If a job is
    # still running, re-enable polling; if it finished while away, restore "Done".
    if not job_id:
        raise PreventUpdate

    folder = os.path.join(os.getcwd(), UPLOAD_FOLDER, job_id)
    if not os.path.isdir(folder):          # cleaned up or never existed
        raise PreventUpdate

    email_marker = os.path.join(folder, "email.txt")
    done_file = os.path.join(folder, "done.txt")
    error_file = os.path.join(folder, "error.txt")
    progress_file = os.path.join(folder, "progress.json")
    
    if os.path.exists(email_marker):
        if os.path.exists(error_file):
            return "❌ Error occurred", 0, "", HIDE, True
        return ("Analysis submitted. When it's finished you will be notified via email.",
                0, "", HIDE, True)

    if os.path.exists(error_file):
        return "❌ Error occurred", 100, "Error", SHOW, True

    if os.path.exists(done_file):
        # This restores the "Done" state; it no longer re-triggers the download.
        # download_result auto-sends each job once and then only on an explicit
        # click of btn-download-again, which this status also reveals.
        return "✅ Done!", 0, "0%", SHOW, True

    if os.path.exists(progress_file):      # still running -> resume polling
        with open(progress_file) as f:
            data = json.load(f)
        progress = data.get("progress", 0)
        status = data.get("status", "")
        return status, progress, f"{progress}%", SHOW, False

    return "Resuming...", 0, "0%", SHOW, False
  
@callback(
    Output("download-pdb-result", "data"),
    Output("pdb-delivered", "data"),
    Input("pdb-status", "children"),
    Input("btn-download-again", "n_clicks"),
    State("pdb-job-id", "data"),
    State("pdb-delivered", "data"),
    prevent_initial_call=True
)
def download_result(status, n_clicks, job_id, delivered):

    clicked = ctx.triggered_id == "btn-download-again"

    if not clicked and "Done" not in (status or ""):
        return dash.no_update, dash.no_update

    if not job_id:
        return dash.no_update, dash.no_update

    # pdb-job-id lives in localStorage, so a finished job is restored on every
    # page load and used to re-emit "Done". Without this guard that replayed a
    # download of the same file forever. Auto-send each job exactly once; after
    # that the user asks for it explicitly via btn-download-again.
    if not clicked and delivered == job_id:
        return dash.no_update, dash.no_update

    folder = folder_for_job(job_id)
    file_path = result_archive_path(folder)

    # Gate on done.txt, not on the status text. The status the poller shows comes
    # straight from progress.json, and an intermediate stage writing "Done" there
    # used to fire this callback early and hand the user an archive with no
    # thermostability columns.
    if not os.path.exists(os.path.join(folder, "done.txt")):
        return dash.no_update, dash.no_update

    if not os.path.exists(file_path):
        return dash.no_update, dash.no_update

    return dcc.send_file(file_path), job_id


@callback(
    Output("btn-download-again", "style"),
    Input("pdb-status", "children"),
    State("pdb-job-id", "data"),
    prevent_initial_call=True,
)
def toggle_download_button(status, job_id):
    """Offer the file whenever a finished job is on disk -- in particular after
    a page reload, where the automatic download deliberately no longer fires."""
    if not job_id or "Done" not in (status or ""):
        return {"display": "none"}

    folder = folder_for_job(job_id)
    if not os.path.exists(os.path.join(folder, "done.txt")):
        return {"display": "none"}
    if not os.path.exists(result_archive_path(folder)):
        return {"display": "none"}

    return {"display": "block"}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    app.run(debug=True, port=8080, use_reloader=True)

import pandas as pd
import numpy as np
import joblib
import re
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
from math import pi
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, MinMaxScaler, PowerTransformer, RobustScaler
from sklearn.impute import SimpleImputer, KNNImputer
from sklearn.linear_model import LinearRegression
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram as scipy_dendrogram
from scipy.stats import mannwhitneyu, ttest_ind
from scipy.spatial.distance import pdist, squareform
import umap
import networkx as nx
import sys
import os
from glob import glob
from api_v6 import *
from scipy.stats import skew, kurtosis, describe
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
import io
import base64


pc = ['BetaSheet_1', 'Aggregation_1', 'Disorder_1', 'Aggregation_2',       'Membrane_1', 'Hydrophobicity_1', 'BetaSheet_2', 'BetaSheet_3',       'AlphaHelix_1', 'Membrane_2', 'BetaSheet_4', 'Burial_1',       'Burial_2', 'Hydrophobicity_2', 'NucleicAcidBinding_1',       'AlphaHelix_2', 'BetaSheet_5', 'charge', 'LowFolding_1',       'Membrane_3', 'Turn_1', 'fg', 'BetaSheet_6',       'NucleicAcidBinding_2', 'BetaSheet_7', 'Disorder_2', 'Membrane_4',       'Burial_3', 'Hydrophobicity_3', 'BetaSheet_8', 'Burial_4',       'Disorder_3', 'Aggregation_3', 'Membrane_5', 'Hydrophobicity_4',       'Hydrophobicity_5', 'Hydrophobicity_6', 'Aggregation_4',       'Aggregation_5', 'Disorder_4', 'Aggregation_6', 'Hydrophobicity_7',       'Membrane_6', 'AlphaHelix_3', 'Disorder_5', 'Disorder_6',       'AlphaHelix_4', 'Burial_5', 'NucleicAcidBinding_3', 'AlphaHelix_5',       'NucleicAcidBinding_4', 'BetaSheet_9', 'Aggregation_7',       'Hydrophobicity_8', 'rg', 'Disorder_7', 'NucleicAcidBinding_5',       'Disorder_8', 'Disorder_9', 'Aggregation_8', 'Burial_6',       'AlphaHelix_6', 'NucleicAcidBinding_6', 'AlphaHelix_7',       'Aggregation_9', 'BetaSheet_10', 'Burial_7', 'Disorder_10',       'Hydrophobicity_9', 'AlphaHelix_8', 'Burial_8', 'AlphaHelix_9',       'NucleicAcidBinding_7', 'NucleicAcidBinding_8',       'NucleicAcidBinding_9', 'AlphaHelix_10', 'NucleicAcidBinding_10',       'Hydrophobicity_10', 'Membrane_7', 'Membrane_8', 'Membrane_9',       'Burial_9']
af = ['n_contacts', 'RG_protein', 'rmsd_prot', 'Length',       'fullCharge', 'extCharge', 'Coil_FullSeq', 'AlphaHelix_FullSeq',       'BetaBridge_FullSeq', 'Strand_FullSeq', '3_10Helix_FullSeq',       'PiHelix_FullSeq', 'Turn_FullSeq', 'Bend_FullSeq', 'Coil_ExtSeq',       'AlphaHelix_ExtSeq', 'BetaBridge_ExtSeq', 'Strand_ExtSeq',       '3_10Helix_ExtSeq', 'PiHelix_ExtSeq', 'Turn_ExtSeq', 'Bend_ExtSeq',       'asa_mean', 'asa_std', 'RBD_Full_min_0', 'RBD_Full_min_1',       'RBD_Full_min_2', 'RBD_int_min_0', 'RBD_int_min_1',       'RBD_int_min_2', 'RBD_ext_min_0', 'RBD_ext_min_1', 'RBD_ext_min_2',       'average_plddt', 'stddev_plddt', 'n_contacts_norm',       'RG_protein_norm']
graph = ['number_of_nodes', 'number_of_edges', 'Density',       'size_of_largest_connected_component', 'avgerage_degree',       'degree_assortativity', 'betweeness_centrality',       'eigenvector_centrality', 'closeness_centrality',       'katz_centrality', 'harmonic_centrality', 'pagerank', 'hub_scores',       'authority_scores', 'local_clustering_coefficient',       'average_clustering_coefficient', 'global_transitivity',       'degree_centrality', 'average_shortest_path_length', 'diameter',       'entropy']
th = ['Ts_tot_PN', 'Ts_tot_N', 'Ts_tot_P', 'Ts_C_PN', 'Ts_C_N', 'Ts_C_P', 'Ts_lj_PN', 'Ts_lj_N', 'Ts_lj_P']
expr = ['mean_expression','tau']

# ------PREPROCESSING------
def allometric_residual(feature_series, length_series):
    log_f = np.log1p(feature_series.values)
    log_L = np.log1p(length_series.values).reshape(-1,1)
    
    reg = LinearRegression().fit(log_L, log_f)
    b = reg.coef_[0]
    return log_f - reg.predict(log_L)

def divide_by_std(df):
    '''Assumes 1 dof'''
    return df / df.std()

def transform_features(
    df,
    selected_features=None,
    kurtosis_threshold=10,
    skewness_threshold=2
):
    
    data = df[selected_features].copy() if selected_features else df.copy()
    
    heavy_features = []
    for c in data.columns:
        k = kurtosis(data[c].loc[~data[c].isna()], fisher=False)
        s = skew(data[c].loc[~data[c].isna()])
        if k > kurtosis_threshold or s > skewness_threshold:
            heavy_features.append(c)

    if len(heavy_features) > 0:
        pt = PowerTransformer(method='yeo-johnson', standardize=False)

        data[heavy_features] = pt.fit_transform(data[heavy_features])
        
    return data
    

def preprocess_noscaling_features(
    df,
    selected_features=None,
    knn_neighbors=15,
    kurtosis_threshold=10,
    skewness_threshold=2
):
    '''A preprocessing function that uses standard scaler and powertransformer for heavy tailed features'''
    data = df[selected_features].copy() if selected_features else df.copy()
    length_col = 'Length'
    extensive_features = ['n_contacts','number_of_nodes','number_of_edges','size_of_largest_connected_component']
    residual_features = ['harmonic_centrality','diameter','average_shortest_path_length','RG_protein','rmsd_prot','entropy']
    features_to_drop = ['authority_scores','number_of_nodes']
    
    info = {}
    # detect extensive features
    info['extensive_features'] = extensive_features

    # Normalize extensive by length
    for f in extensive_features:
        if f in data.columns:
            data[f] = df[f] / df[length_col]
    
    # normalize with allometric residuals log(Y/Ypred) weakly extensive features      
    for f in residual_features:
        if f in data.columns:
            data[f] = allometric_residual(data[f], data[length_col])


    # detect heavy tailed features
    heavy_features = []
    for c in data.columns:
        k = kurtosis(data[c].loc[~data[c].isna()], fisher=False)
        s = skew(data[c].loc[~data[c].isna()])
        if k > kurtosis_threshold or s > skewness_threshold:
            heavy_features.append(c)

    info['yeojohnson_features'] = heavy_features

    if len(heavy_features) > 0:
        pt = PowerTransformer(method='yeo-johnson', standardize=False)

        data[heavy_features] = pt.fit_transform(data[heavy_features])
        
    knn_scaler = StandardScaler()
    X_scaled = knn_scaler.fit_transform(data)

    imputer = KNNImputer(n_neighbors=knn_neighbors, weights='distance')
    X_imputed = imputer.fit_transform(X_scaled)

    X_final = knn_scaler.inverse_transform(X_imputed)

    data_final = pd.DataFrame(X_final, index=data.index, columns=data.columns)

    return data_final, info

def preprocess_features(
    df,
    selected_features=None,
    knn_neighbors=15,
    kurtosis_threshold=10,
    skewness_threshold=2,
    pca_componets=30,
    for_umap=False
):
    '''A preprocessing function that uses standard scaler and powertransformer for heavy tailed features'''
    data = df[selected_features].copy() if selected_features else df.copy()
    length_col = 'Length'
    extensive_features = ['n_contacts','number_of_nodes','number_of_edges','size_of_largest_connected_component']
    residual_features = ['harmonic_centrality','diameter','average_shortest_path_length','RG_protein','rmsd_prot','entropy']
    features_to_drop = ['authority_scores','number_of_nodes']
    
    info = {}
    # detect extensive features
    info['extensive_features'] = extensive_features

    # Normalize extensive by length
    for f in extensive_features:
        if f in data.columns:
            data[f] = df[f] / df[length_col]
    
    # normalize with allometric residuals log(Y/Ypred) weakly extensive features      
    for f in residual_features:
        if f in data.columns:
            data[f] = allometric_residual(data[f], data[length_col])

    # detect heavy tailed features
    heavy_features = []
    for c in data.columns:
        k = kurtosis(data[c].loc[~data[c].isna()], fisher=False)
        s = skew(data[c].loc[~data[c].isna()])
        if k > kurtosis_threshold or s > skewness_threshold:
            heavy_features.append(c)

    info['yeojohnson_features'] = heavy_features

    if len(heavy_features) > 0:
        pt = PowerTransformer(method='yeo-johnson', standardize=False)

        data[heavy_features] = pt.fit_transform(data[heavy_features])

    # KNN imputation
    knn_scaler = StandardScaler()
    X_scaled = knn_scaler.fit_transform(data)

    imputer = KNNImputer(n_neighbors=knn_neighbors, weights='distance')
    X_imputed = imputer.fit_transform(X_scaled)
    
    scaler = StandardScaler()
    X_final = scaler.fit_transform(X_imputed)

    data_final = pd.DataFrame(X_final, index=data.index, columns=data.columns)

    # PCA reduction for UMAP
    if for_umap:
        n_comp = min(pca_componets, X_final.shape[1])
        pca = PCA(n_components=n_comp)
        X_final = pca.fit_transform(X_final)

        info['pca_variance_explained'] = np.sum(pca.explained_variance_ratio_)
        data_final = pd.DataFrame(X_final, index=data_imputed.index)

    return data_final, info


def normalize_profile_feature(df, feature):
    """df: the profiles_df"""
    all_values = np.concatenate(df[feature].values).reshape(-1,1)
    scaler = MinMaxScaler()
    scaler.fit(all_values)

    return df[feature].apply(lambda x: scaler.transform(np.array(x).reshape(-1,1)).ravel())
    
    
    
 # ------CLUSTERING------
def perform_umap(df_scaled, selected_features=None, n_components=2, n_neighbors=15, min_dist=0.1, metric='euclidean'):
    """
    Perform UMAP on descriptor data and return UMAP results for visualization.

    Args:
    - df_scaled (pd.DataFrame): Preprocessed DataFrame with protein descriptors (rows = proteins, cols = descriptors).
    - n_components (int): Number of UMAP components to compute (default=2 for visualization).
    - n_neighbors (int): Determines the balance between local and global structure.
    - min_dist (float): Minimum distance between embedded points (smaller = more clustered).
    - metric (str): Distance metric for UMAP (default='euclidean').

    Returns:
    - umap_df (pd.DataFrame): DataFrame of UMAP components for plotting.
    """
    #df_scaled, info = preprocess_features(df, selected_features=selected_features) if selected_features else preprocess_features(df)
    original_df_index = df_scaled.index
    if selected_features:
        df_scaled = df_scaled[selected_features]
    # Perform UMAP
    reducer = umap.UMAP(n_components=n_components, n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, random_state=42)
    umap_result = reducer.fit_transform(df_scaled)

    # Create a DataFrame with UMAP results
    umap_df = pd.DataFrame(umap_result, columns=[f"UMAP{i+1}" for i in range(n_components)], index=original_df_index)
    return umap_df


def plot_umap(df, ids1, ids2=None, feats=None, save_file=None, labels=None):
    """
    New UMAP function that uses plotly for an interranctive scatterplot of selected group vs another group
    """
    # checks both for not None and empty list (both False)
    if ids2:
        df = df.loc[df.index.to_series().isin(list(ids1) + list(ids2))].copy()
    umap_df = perform_umap(df, selected_features=feats) if feats else perform_umap(df)

    label_selected = 'selected'
    label_unselected = 'Rest'
    if labels:
        label_selected = labels[0]
        label_unselected = labels[1]
        
    umap_df['is_selected'] = umap_df.index.to_series().isin(ids1)
    umap_df['group'] = umap_df['is_selected'].apply(
        lambda x: label_selected if x else label_unselected
    )
    umap_df['label'] = umap_df.index
        
    color_rest = "rgba(99, 110, 250, 0.25)"   
    color_sel = "rgba(239, 85, 59, 0.8)"     

    fig = px.scatter(
        umap_df,
        x="UMAP1",
        y="UMAP2",
        color="group",
        hover_name="label",
        title="UMAP of Protein Group",
        color_discrete_map={
            label_unselected: color_rest,
            label_selected: color_sel
        },
        category_orders={"group": [label_unselected, label_selected]}
    )
    
    fig.update_traces(
        selector={"name": [label_selected]}
    )
    fig.data[-1].marker.size = 8
    fig.data[-1].marker.line = dict(width=1, color="black")
    
    fig.update_layout(
        meta=dict(
            plot_info="umap",
            selected=ids1
        )
    )
    
    if save_file:
        fig.write_html(save_file)
    return fig
    
    
def perform_pca(df_scaled, selected_features=None, n_components=2):
    """
    Perform PCA on descriptor data and return PCA results for visualization.

    Args:
    - df_scaled (pd.DataFrame): Preprocessed DataFrame with protein descriptors (rows = proteins, cols = descriptors).
    - n_components (int): Number of PCA components to compute (default=2 for visualization).

    Returns:
    - pca_df (pd.DataFrame): DataFrame of UMAP components for plotting.
    """
    #df_scaled, info = preprocess_features(df, selected_features=selected_features) if selected_features else preprocess_features(df)
    original_df_index = df_scaled.index
    if selected_features:
        df_scaled = df_scaled[selected_features]
    # Perform PCA
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(df_scaled)

    # Create a DataFrame with PCA results
    pca_df = pd.DataFrame(pca_result, columns=[f"PCA{i+1}" for i in range(n_components)], index=original_df_index)
    return pca_df



def plot_pca(df, ids1, ids2=None, feats=None, save_file=None, labels=None):
    if ids2:
        df = df.loc[df.index.to_series().isin(list(ids1) + list(ids2))].copy()

    pca_df = perform_pca(df, selected_features=feats) if feats else perform_pca(df)

    label_selected = 'selected'
    label_unselected = 'Rest'
    if labels:
        label_selected = labels[0]
        label_unselected = labels[1]

    pca_df['is_selected'] = pca_df.index.to_series().isin(ids1)
    pca_df['group'] = pca_df['is_selected'].apply(
        lambda x: label_selected if x else label_unselected
    )
    pca_df['label'] = pca_df.index
    
    color_rest = "rgba(99, 110, 250, 0.25)"   
    color_sel = "rgba(239, 85, 59, 0.8)"
    
    fig = px.scatter(
        pca_df,
        x="PCA1",
        y="PCA2",
        color="group",
        hover_name="label",
        title="Principal Component Analysis",
        color_discrete_map={
            label_unselected: color_rest,
            label_selected: color_sel
        },
        category_orders={"group": [label_unselected, label_selected]}
    )
    
    fig.update_traces(
        selector={"name": [label_selected]},
        marker=dict(size=11, line=dict(width=1, color="black"))
    )
    fig.data[-1].marker.size = 8
    fig.data[-1].marker.line = dict(width=1, color="black")
    
    fig.update_layout(
        meta=dict(plot_info="pca")
    )
    
    return fig
    
def plot_pca_v1(df, ids1, ids2=None, feats=None, save_file=None, labels=None):

    if ids2 is not None:
        df = df.loc[df.index.to_series().isin(list(ids1) + list(ids2))].copy()

    pca_df = perform_pca(df, selected_features=feats) if feats else perform_pca(df)

    label_selected = 'selected'
    label_unselected = 'Rest'
    if labels:
        label_selected = labels[0]
        label_unselected = labels[1]

    pca_df['is_selected'] = pca_df.index.to_series().isin(ids1)
    pca_df['group'] = pca_df['is_selected'].apply(
        lambda x: label_selected if x else label_unselected
    )
    pca_df['label'] = pca_df.index

    # --- split data ---
    df_rest = pca_df[pca_df["group"] == label_unselected]
    df_sel = pca_df[pca_df["group"] == label_selected]

    color_rest = "rgba(99, 110, 250, 0.25)"   
    color_sel = "rgba(239, 85, 59, 1.0)"

    fig = px.scatter(
        df_rest,
        x="PCA1",
        y="PCA2",
        hover_name="label",
        title="PCA of Protein Group"
    )

    # add unselected first (bottom layer)
    fig.data[0].marker.color = color_rest
    fig.data[0].name = label_unselected

    # add selected
    fig.add_trace(
        px.scatter(
            df_sel,
            x="PCA1",
            y="PCA2",
            hover_name="label"
        ).data[0]
    )
    
    fig.data[1].marker.color = color_sel
    fig.data[1].name = label_selected
    fig.data[1].marker.size = 11
    fig.data[1].marker.line = dict(width=1, color="black")

    #sel_trace = px.scatter(
    #    df_sel,
    #    x="PCA1",
    #    y="PCA2",
    #    hover_name="label"
    #).data[0]
    #sel_trace.marker.color = color_sel
    #sel_trace.name = label_selected
    #sel_trace.marker.size = 11
    #sel_trace.marker.line = dict(width=1, color="black")
    #fig.add_trace(sel_trace)
    
    #fig.data = fig.data[::-1]

    if save_file:
        fig.write_html(save_file)

    return fig


# ------SPIDERPLOT------
def is_correlated(feat_1, feat_2, corr_matrix, thr=.9):
    if corr_matrix.loc[feat_1, feat_2] > thr:
        return True
    return False

    
def compute_statistics(df):
    means = df.mean(axis=0, numeric_only=True)
    std_devs = df.std(axis=0, numeric_only=True)
    return means, std_devs


def compute_z_scores(df, ids, means, std_devs):
    # Guard against constant / near-constant features. A feature whose spread is
    # negligible produces astronomically large z-scores for the few proteins that
    # deviate at all (e.g. number_of_nodes/Length ~= 1 for nearly every protein,
    # so its std is tiny). Because get_top_features ranks by |z|, such a feature
    # would dominate the spider plot and blow it off-scale. We set the std of
    # degenerate features to NaN so their z becomes NaN and they drop to the
    # bottom of the ranking instead.
    feats = means.index
    pop = df[feats]
    degenerate = (std_devs == 0)
    if len(pop) >= 50:
        # robust, unit-free flatness test: inter-quartile range tiny relative to
        # the feature's typical magnitude => effectively constant.
        iqr = (pop.quantile(0.75) - pop.quantile(0.25)).abs()
        scale = pop.abs().median() + 1e-12
        degenerate = degenerate | ((iqr / scale) < 1e-3)
    std_devs = std_devs.mask(degenerate, np.nan)

    z_scores = {}
    for id in ids:
        prot = df.loc[df.index.to_series().isin([id])]
        if prot.empty:
            print(f'{id} is not contained in PROFET')
            continue
        z = (prot[feats] - means) / std_devs
        z_scores[id] = z.iloc[0]
    return z_scores


def get_top_features(z_scores, corr_matrix,  top_n=5):
    first = next(iter(z_scores))
    #top_features = z_scores[first].abs().sort_values(ascending=False).head(top_n)
    top_features = z_scores[first].abs().sort_values(ascending=False)
    feats = []
    vals = []
    for i, f in top_features.items():
        if len(feats) == top_n:
            break
        if i == top_features.index[0]:
            feats.append(i)
            vals.append(f)
            continue
        if i in pc:
            name = re.findall(r'^(.*)_\d+$', i)
            if name:
                name=name[0]
            else:
                feats.append(i)
                vals.append(f)
                continue
            is_in_ls = [k for k in feats if name in k]
            if len(is_in_ls) > 0:
                continue
            else:
                feats.append(i)
                vals.append(f)
        is_correlated_series = pd.Series(data=[is_correlated(i, feat, corr_matrix) for feat in feats])
        if is_correlated_series.any():
            continue
        else:
            feats.append(i)
            vals.append(f)
    top_features = pd.Series(data=vals, index=feats)
    return top_features



def plot_spider_helper(df, ids, means, top_features, save_path=None):
    """
    Plot a radar (spider) chart for selected proteins and top features using Plotly.
    Each protein has its own color.

    Parameters:
    - df: DataFrame with protein data
    - ids: list of protein IDs to plot (first one is anchor)
    - means: Series/DataFrame with mean values for features
    - top_features: Series with top features to include (index = feature names)
    - save_path: optional, file path to save static image (requires kaleido)

    Returns:
    - fig: Plotly figure
    """
    
    labels = list(top_features.index)
    anchor_protein = ids[0]
    prots = df.loc[df.index.to_series().isin(ids)]
    
    # Prepare mean values
    mean_values = means[labels].values
    mean_values = np.append(mean_values, mean_values[0])
    
    # Prepare feature labels for closing the loop
    labels_loop = labels + [labels[0]]
    
    # Set the lengend annotation
    mean_name = "Mean - Human Proteome" if df.shape[0] == 20335 else "Group Mean"
    
    fig = go.Figure()
    
    # Add mean line (dashed red) with fill
    fig.add_trace(go.Scatterpolar(
        r=mean_values,
        theta=labels_loop,
        mode='lines',
        line=dict(color='red', dash='dash', width=2),
        fill='toself',
        name=mean_name,
        opacity=0.2
    ))
    
    # Assign a color for each protein using Plotly qualitative palette
    protein_colors = px.colors.qualitative.Plotly
    num_colors = len(protein_colors)
    
    # Add each protein line with its own color
    for idx_counter, (idx, r) in enumerate(prots.iterrows()):
        prot_values = r[labels].values
        prot_values = np.append(prot_values, prot_values[0])
        color = protein_colors[idx_counter % num_colors]  # now idx_counter is integer
        fig.add_trace(go.Scatterpolar(
            r=prot_values,
            theta=labels_loop,
            mode='lines',
            line=dict(color=color, dash='solid', width=2),
            fill='toself',
            name=str(idx),  # use actual protein ID for legend
            opacity=0.2
    ))
    
    # Layout
    fig.update_layout(
        margin=dict(l=100, r=100, t=100, b=100),  # left, right, top, bottom
        polar=dict(
            radialaxis=dict(visible=True),
        ),
        title=dict(
            text=f'{anchor_protein} - Top {len(labels)} Features with Significant Deviation',
            y=0.95,
            x=0.5,
            xanchor='center',
            yanchor='top'
        ),
        showlegend=True,
        width=1000,
        height=800,
    )
    
    # Optional: save static image
    if save_path:
        fig.write_image(save_path, engine="kaleido", scale=2)
    
    return fig
    
def plot_spider_helper_z(z_scores, ids, top_features, save_path=None):

    labels = list(top_features.index)
    anchor_protein = ids[0]

    labels_loop = labels + [labels[0]]

    fig = go.Figure()

    # Mean line = zero
    mean_values = np.zeros(len(labels))
    mean_values = np.append(mean_values, mean_values[0])

    fig.add_trace(go.Scatterpolar(
        r=mean_values,
        theta=labels_loop,
        mode='lines',
        line=dict(color='red', dash='dash', width=2),
        name="Mean (Z=0)",
        opacity=0.2
    ))

    protein_colors = px.colors.qualitative.Plotly

    for idx_counter, pid in enumerate(ids):

        z = z_scores[pid][labels].values
        z = np.append(z, z[0])

        fig.add_trace(go.Scatterpolar(
            r=z,
            theta=labels_loop,
            mode='lines',
            line=dict(color=protein_colors[idx_counter % len(protein_colors)], width=2),
            fill='toself',
            name=pid,
            opacity=0.3
        ))

    # --- Dynamic radial range ---------------------------------------------
    # Baseline is [-3.5, 3.5] (roughly the range of |z| for a normal feature),
    # but grow the axis symmetrically when a protein deviates further so genuine
    # outliers stay on-scale instead of being clipped to the rim.
    max_abs = 3.5
    for pid in ids:
        vals = z_scores[pid][labels].values.astype(float)
        vals = np.abs(vals[~np.isnan(vals)])
        if vals.size:
            max_abs = max(max_abs, float(vals.max()))
    axis_limit = 3.5 if max_abs <= 3.5 else float(np.ceil(max_abs))

    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[-axis_limit, axis_limit]
            )
        ),
        title=f"{anchor_protein} - Top Deviating Features (Z-score)",
        showlegend=True
    )

    return fig
    

def plot_spider_helper_old(df, ids, means, top_features, save_path):
    labels = list(top_features.index)

    anchor_protein = ids[0]
    prots = df.loc[df.index.to_series().isin(ids)]

    mean_values = means[labels].values
    
    num_vars = len(labels)
    angles = [n / float(num_vars) * 2 * pi for n in range(num_vars)]
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(16, 9), subplot_kw=dict(polar=True))

    mean_values = np.concatenate((mean_values, [mean_values[0]]))
    ax.plot(angles, mean_values, linewidth=2, linestyle='dashed', label='Mean')
    ax.fill(angles, mean_values, 'r', alpha=0.1)

    for i, r  in prots.iterrows():
        name = i
        prot_values = r[labels].values
        prot_values = np.concatenate((prot_values, [prot_values[0]]))
        ax.plot(angles, prot_values, linewidth=2, linestyle='solid', label=f'{name}')
        ax.fill(angles, prot_values, 'b', alpha=0.1)
           
    
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels)

    plt.title(f'{anchor_protein} - Top 5 Features with Significant Deviation', y=1.08)
    plt.legend(loc='upper right', bbox_to_anchor=(1.5,1))
    if save_path:
        plt.savefig(save_path)

    #plt.show()
    return fig


def plot_spider_plot(df_scaled, ids, ids_to_check_against=None, top_n=5, save_path=None):
    '''A function to check a group of proteins against the human proteome or another group of proteins'''

    if ids_to_check_against:
        df_scaled = df_scaled.loc[df_scaled.index.to_series().isin(list(ids) + list(ids_to_check_against))].copy()
        #mask = df.index.to_series().str.contains('|'.join(list(ids_to_check_against) + list(ids)))
        #df = df[mask]
    
    means, std_devs = compute_statistics(df_scaled)
    
    z_scores = compute_z_scores(df_scaled, ids, means, std_devs)
    
    # Get top 5 features with significant deviation
    top_features = get_top_features(z_scores, df_scaled.corr(method='spearman').abs(), top_n)
    
    # Plot the spider plot
    #fig = plot_spider_helper(df_scaled, ids, means, top_features, save_path)
    fig = plot_spider_helper_z(z_scores, ids, top_features)
    return fig
  
  
def plot_spider_plot_features(df_scaled, ids, feats, ids_to_check_against=None, save_path=None):
    '''A function to check a group of proteins against the human proteome or another group of proteins'''

    if ids_to_check_against:
        df_scaled = df_scaled.loc[df_scaled.index.to_series().isin(list(ids) + list(ids_to_check_against))].copy()
        #mask = df.index.to_series().str.contains('|'.join(list(ids_to_check_against) + list(ids)))
        #df = df[mask]
    
    df_scaled = df_scaled[feats]
    
    means, std_devs = compute_statistics(df_scaled)
    
    z_scores = compute_z_scores(df_scaled, ids, means, std_devs)
    print(z_scores)
    print(type(z_scores[next(iter(z_scores))]))
    print(z_scores[next(iter(z_scores))])
    
    feats = pd.Series(data=z_scores[ids[0]].values, index=z_scores[ids[0]].index)

    fig = plot_spider_helper_z(z_scores, ids, feats)
    return fig
    
    
# ------PROFILEPLOT------
def plot_profile_one_seq_multiple_features(
    df, target_id, feats, save_path=None
):
    """
    Plot one protein sequence with multiple per-residue features.
    
    - Each feature uses its own native y-axis
    - Only the first y-axis is visible
    - Other y-axes are invisible but fully functional
    - Sliding window + AA annotations preserved
    """

    fig = go.Figure()

    window_size = 100

    if isinstance(target_id, str):
        list_target_id = [target_id]
    
    if df.loc[df.index.to_series().isin(list_target_id)].empty:
        print(f'{target_id} not in PROFET')
        return -1
    sequence = df.loc[df.index.to_series().isin(list_target_id)]['seq'].values[0]
    
    df = df[feats]
    for f in feats:
        df.loc[:, f] = normalize_profile_feature(df, f)
        
    df = df.loc[df.index.to_series().isin(list_target_id)]
    
    max_len = len(sequence)
    initial_range = [1, min(window_size, max_len)]

    x = list(range(1, max_len + 1))
    aa = list(sequence)

    aa_row_y = 0.005
    yaxis_layout = {}

    # ---------- add traces ----------
    for i, f in enumerate(feats):
        name = target_id
        values = df.loc[target_id, f]

        axis_name = "y" if i == 0 else f"y{i+1}"

        fig.add_trace(go.Scatter(
            x=x,
            y=values,
            mode='lines+markers',
            name=f"{name}_{f}",
            yaxis=axis_name,
            text=aa,
            hovertemplate=(
                f"<b>{f}</b><br>"
                "AA: %{text}<br>"
                "Value: %{y}<br>"
                "Pos: %{x}"
            )
        ))

    for i in range(1, len(feats) + 1):
        axis_name = f'yaxis{i}'
        side = 'left' if i % 2 == 1 else 'right'
        position = 0.05 * (i - 1) if side == 'left' else 1 - 0.05 * (i - 1)

        if i == 1:
            yaxis_layout[axis_name] = dict(
                title=feats[i-1],
                showticklabels=True,
                showgrid=True,
                side='left',
                domain=[0.09, 1]
            )
        else:
            yaxis_layout[axis_name] = dict(
                title='',
                showticklabels=False,
                showgrid=False,
                overlaying="y" if i > 1 else None,
                side=side,
                position=position
            )
        
    annotations = []
    for j, aa_char in enumerate(sequence):
        if initial_range[0] <= j + 1 <= initial_range[1]:
            annotations.append(dict(
                x=j + 1,
                y=aa_row_y,
                text=aa_char,
                showarrow=False,
                xref="x",
                yref="paper",
                font=dict(size=10)
            ))
            
    # ---------- slider steps ----------
    steps = []
    for start in range(1, max_len - window_size + 1, 10):
        end = min(start + window_size - 1, max_len)

        step_annotations = []
        for j in range(start - 1, end):
            step_annotations.append(dict(
                x=j + 1,
                y=aa_row_y,
                text=sequence[j],
                showarrow=False,
                xref="x",
                yref="paper",
                font=dict(size=10)
            ))

        steps.append(dict(
            label=f"{start}-{end}",
            method="update",
            args=[
                {
                    "x": [x[start-1:end]] * len(feats),
                    "y": [
                        df.loc[target_id, f][start-1:end]
                        for f in feats
                    ],
                },
                {
                    "annotations": step_annotations,
                    "xaxis.range": [start, end]
                    #"xaxis": {"range": [start, end],
                    #    "title": "Amino Acid Position"}
                }
            ]
        ))

    # ---------- layout ----------
    fig.update_layout(
        title=f"{name} per-residue features (native scales)",
        title_y=0.97,
        xaxis=dict(
            title="Amino Acid Position",
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            range=initial_range
        ),
        annotations=annotations,
        sliders=[dict(
            active=0,
            currentvalue={"prefix": "Window: "},
            steps=steps,
            x=0.03,
            y=1.2,
            xanchor="left",
            len=0.9
        )],
        margin=dict(t=120, b=150),
        height=800,
        showlegend=True,
        updatemenus=[
            dict(
                type="buttons",
                direction="left",
                x=.5,
                y=1.2,
                xanchor="center",
                buttons=[
                    dict(
                        label="Full sequence",
                        method="update",
                        args=[
                            {
                                "x": [x] * len(feats),
                                "y": [df.loc[target_id, f] for f in feats],
                            },
                            {
                                "xaxis.range": [1, max_len],
                                "annotations": []
                            }
                        ]
                    )
                ]
            )
        ],
        **yaxis_layout
    )

    if save_path:
        fig.write_html(save_path)

    #fig.show()
    return fig


def plot_profile_multiple_seq_one_feature(df, target_ids, feat):
    """
    Plot multiple sequences in the same figure with synced x-axis, 
    and show a separate amino acid label row (tick-style) for each sequence.
    """
    fig = go.Figure()
    window_size = 100
    df = df.loc[df.index.to_series().isin(target_ids)]
    if df.empty:
        print(f'None of the ids were found in the database')
        return -1
        
    sequences = df['seq'].values
    max_len = max(len(seq) for seq in sequences)
    initial_range = [1, min(window_size, max_len)]
    
    row_spacing = 0.08
    base_offset = 0.15

    all_values = []
    # ----------- traces -----------
    for target_id in target_ids:
        if target_id not in df.index:
            print(f'{name} not contained in db')
            continue
    #for target_id, sequence in zip(target_ids, sequences):
        sequence = df.loc[target_id, 'seq']
        name = target_id
        
        profile_df = df.loc[df.index.to_series().isin([target_id])]
        #if profile_df.empty:
        #    print(f'{name} not contained in db')
        #    continue
            
        values = profile_df[feat].iloc[0]
        all_values.extend(values)

        x = list(range(1, len(sequence) + 1))
        y = values
        aa = list(sequence)

        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode='lines+markers',
            name=target_id,
            text=aa,
            hovertemplate="AA: %{text}<br>Value: %{y:.3f}<br>Pos: %{x}"
        ))

    ymin = np.nanmin(all_values)
    ymax = np.nanmax(all_values)
    pad = 0.05
    y_range = [ymin - pad, ymax + pad]
    # ---------- first frame annotations ----------
    annotations = []
    for i, sequence in enumerate(sequences):
        #for j, aa in enumerate(sequence):
        for j in range(initial_range[0]-1, min(initial_range[1], len(sequence))):
            annotations.append(
                dict(
                    x=j + 1,
                    #y=-(base_offset + i * row_spacing),  
                    y=0.15 - i * row_spacing,
                    text=sequence[j],
                    showarrow=False,
                    xref="x",
                    yref="paper",
                    font=dict(size=10)
                )
            )

    # Slider steps
    steps = []
    for start in range(1, max_len - window_size + 1, 10):
        end = min(start + window_size - 1, max_len)

        step_annotations = []
        for i, sequence in enumerate(sequences):
            for j in range(start - 1, min(end, len(sequence))):
                step_annotations.append(
                    dict(
                        x=j + 1,
                        #y=-(base_offset + i * row_spacing),
                        y=0.15 - i * row_spacing,
                        text=sequence[j],
                        showarrow=False,
                        xref="x",
                        yref="paper",
                        font=dict(size=10)
                    )
                )

        steps.append(dict(
            label=f"{start}-{end}",
            method="relayout",
            args=[
                {
                    "xaxis.range": [start, end],
                    "annotations": step_annotations
                }
            ]
        ))
        
    full_annotations = []
    
    updatemenus = [
        dict(
            type="buttons",
            buttons=[
                dict(
                    label="Full Sequence",
                    method="relayout",
                    args=[{
                        "xaxis.range":[1, max_len],
                        "annotations":full_annotations
                    }]
                )
            ],
            direction="left",
            x=.5,
            y=1.2,
            xanchor="center"
        )
    ]

    fig.update_layout(
        title=f"{feat} profile plot",
        xaxis=dict(
            title="Amino Acid Position",
            showticklabels=False,  
            showgrid=False,
            zeroline=False,
            range=initial_range
        ),
        yaxis=dict(
            title=f"{feat}",
            range=y_range,
            domain=[0.25,1]
        ),
        annotations=annotations,
        sliders=[dict(
            active=0,
            currentvalue={"prefix": "Window: "},
            steps=steps,
            x=0.1, 
            y=1.1,  
            xanchor="left",
            len=0.9
        )],
        updatemenus=updatemenus,
        margin=dict(t=100, b=150 + int(len(sequences)) * 20), 
        autosize=False,
        height=1000,
        #height=600 + int(len(sequences)) * 30,
        showlegend=True,
        uirevision="constant"
    )

    #fig.show()
    return fig
    
# ------VIOLINPLOT------
def plot_violinplot(df, col, save_file=None):
    plt.figure(figsize=(12, 8))
    sns.violinplot(y=df[col])
    plt.title(f'Violin Plot for {col}')
    plt.ylabel("Values")
    plt.xlabel(f"{col}")
    plt.tight_layout()
    if save_file:
        plt.savefig(save_file, dpi=300)
    #plt.show()
    return plt



def plot_violinplot_w_selection(df, col, ids: list, save_file=None):
    """
    Create a Plotly violin plot for a column with selected IDs highlighted.
    
    Parameters:
    - df: DataFrame with numeric data
    - col: column to plot
    - ids: list of row indices (or protein IDs) to highlight
    - save_file: optional, path to save static image (requires kaleido)
    
    Returns:
    - Plotly figure
    """
    
    # Create base violin plot
    fig = px.violin(
        df,
        y=col,
        box=True,  # Show box inside violin
        points="all",  # Show all data points
        hover_data=df.columns
    )
    
    # Add horizontal lines for selected IDs
    colors = px.colors.qualitative.Pastel
    j=0
    for i, idx in enumerate(ids):
        if idx not in df.index:
            print(f"{idx} not in DataFrame")
            continue
        y_val = df.loc[idx, col]
        if pd.isna(y_val):
            fig.add_annotation(
                text=f"{idx}: No Value",
                x=0.05,  
                y=0.55 + (j * 0.05),  
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(color=colors[i % len(colors)], size=12)
                )
            j += 1
        else:
            fig.add_hline(
                y=y_val,
                line_dash="dash",
                line_color=colors[i % len(colors)],
                annotation_text=str(idx),
                annotation_position="top right",
                annotation_font_size=12
            )
    
    # Update layout
    fig.update_layout(
        title=f"Violin Plot for {col}",
        yaxis_title="Values",
        xaxis_title=col,
        template="plotly_white",
        width=800,
        height=600,
        showlegend=False
    )
    
    # Optional: save static image
    if save_file:
        fig.write_image(save_file, engine="kaleido", scale=2)
    
    return fig



def plot_violinplot_w_selection_old(df, col, ids: list, save_file=None):
    '''
    ids is a list
    '''
    plt.figure(figsize=(12, 8))
    sns.violinplot(y=df[col])
    plt.title(f'Violin Plot for {col}')
    sub_df = df.loc[df.index.to_series().isin(ids)]
    colors = sns.color_palette("pastel", len(ids))
    for i, id in enumerate(ids):
        index = id
        selected_protein_value = df.loc[df.index.to_series().isin([index])][col]
        if selected_protein_value.empty:
            print(f'{id} not in PROFET')
        else:
            plt.axhline(selected_protein_value.iloc[0], color=colors[i], linestyle="--", linewidth=1.5, label=f"{id}")
    plt.ylabel("Values")
    plt.xlabel(f"{col}")
    plt.legend()
    plt.tight_layout()
    if save_file:
        plt.savefig(save_file, dpi=300)
    #plt.show()
    return plt


# ------HEATMAP------
def plot_heatmap_new(df_scaled, ids, feats, save_path=None):
    mask = df_scaled.index.to_series().isin(ids)
    index_series = df_scaled.index.to_series()
    df_scaled = df_scaled.loc[mask]
    df_scaled = df_scaled[feats]

    not_included_list = list(set(ids) - set(index_series))
    if not_included_list:
        print(f'{", ".join(not_included_list)} not in PROFET')

    # --- hierarchical clustering ---
    row_linkage = linkage(df_scaled.values, method="average", metric="euclidean")
    col_linkage = linkage(df_scaled.values.T, method="average", metric="euclidean")
    row_order = leaves_list(row_linkage)
    col_order = leaves_list(col_linkage)
    df_clustered = df_scaled.iloc[row_order, col_order]
    
    n_rows = len(df_clustered.index)
    n_cols = len(df_clustered.columns)
    
    def make_dendro_traces(lm, orientation='top'):
        dend = scipy_dendrogram(lm, no_plot=True)
        traces = []
        for xs, ys in zip(dend['icoord'], dend['dcoord']):
            xs_scaled = [(x - 5) / 10 for x in xs]
            if orientation == 'top':
                t = go.Scatter(
                    x=xs_scaled, y=ys,
                    mode='lines', line=dict(color='#444', width=1),
                    showlegend=False, hoverinfo='none'
                )
            else:
                t = go.Scatter(
                    x=ys, y=xs_scaled,
                    mode='lines', line=dict(color='#444', width=1),
                    showlegend=False, hoverinfo='none'
                )
            traces.append(t)
        return traces
    
    col_dencro_traces = make_dendro_traces(col_linkage, orientation='top')
    row_dencro_traces = make_dendro_traces(row_linkage, orientation='right')

    # --- build subplot grid: top-left empty, top-right col dendro, bottom-left row dendro, bottom-right heatmap ---
    fig = make_subplots(
        rows=2, cols=2,
        column_widths=[0.15, 0.85],
        row_heights=[0.15, 0.85],
        horizontal_spacing=0.06,
        vertical_spacing=0.01,
    )

    # row dendrogram (left)
    for trace in row_dencro_traces:
        trace.showlegend = False
        fig.add_trace(trace, row=2, col=1)

    # col dendrogram (top)
    for trace in col_dencro_traces:
        trace.showlegend = False
        fig.add_trace(trace, row=1, col=2)

    # heatmap
    fig.add_trace(
        go.Heatmap(
            z=df_clustered.values,
            x=list(df_clustered.columns),
            y=list(df_clustered.index),
            colorscale="RdBu_r",
            zmin=-3.5, zmax=3.5,
            colorbar=dict(title="Std Value", x=1.02),
            hovertemplate="Protein: %{y}<br>Feature: %{x}<br>Value: %{z:.3f}<br>Right-click for info<extra></extra>"
        ),
        row=2, col=2
    )
 
    # --- layout ---
    plot_height = max(600, 20 * len(df_clustered))
    plot_width  = max(1000, 60 * len(df_clustered.columns))

    fig.update_layout(
        title="Hierarchical Clustering Heatmap",
        height=plot_height,
        width=plot_width,
        showlegend=False,
        paper_bgcolor="white",
        plot_bgcolor="white",
        margin=dict(l=20, r=120, t=80, b=120),
    )

    # hide dendrogram axes ticks
    fig.update_xaxes(
        range=[-0.5, n_cols - 0.5],
        showticklabels=False, showgrid=False, zeroline=False,
        fixedrange=True,
        row=1, col=2
    )
    fig.update_yaxes(
        showticklabels=False, showgrid=False, zeroline=False,
        fixedrange=True,
        row=1, col=2
    )
    
    fig.update_xaxes(
        showticklabels=False, showgrid=False, zeroline=False,
        autorange='reversed',   # distances grow left from heatmap
        fixedrange=True,
        row=2, col=1
    )
    fig.update_yaxes(
        range=[-0.5, n_rows - 0.5],
        showticklabels=False, showgrid=False, zeroline=False,
        fixedrange=True,
        row=2, col=1
    )

    # --- heatmap axes (row=2, col=2) ---
    fig.update_xaxes(tickangle=-45, showgrid=False, row=2, col=2)
    fig.update_yaxes(showgrid=False, row=2, col=2)

    # --- empty top-left cell ---
    fig.update_xaxes(showticklabels=False, showgrid=False, zeroline=False, row=1, col=1)
    fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False, row=1, col=1)

    
    if save_path:
        fig.write_html(save_path)

    return fig

def plot_heatmap_with_dendrogram(df_scaled, ids, feats, save_path=None):
    """
    Clustered heatmap with dendrograms using seaborn.

    Args:
        df_scaled (pd.DataFrame): Preprocessed dataframe.
        ids (list[str]): Protein IDs to include.
        feats (list[str]): Features to include.
        save_path (str, optional): Path to save figure (PNG/PDF).
    """

    # --- filter data ---
    mask = df_scaled.index.to_series().isin(ids)
    index_series = df_scaled.index.to_series()

    df_scaled = df_scaled.loc[mask, feats]

    not_included_list = list(set(ids) - set(index_series))
    if len(not_included_list) > 0:
        print(f'{", ".join(not_included_list)} not in PROFET')

    if df_scaled.empty:
        raise ValueError("No matching IDs found in dataframe.")

    # --- clustering ---
    row_linkage = linkage(df_scaled, method="average", metric="euclidean")
    col_linkage = linkage(df_scaled.T, method="average", metric="euclidean")

    # --- seaborn clustermap ---
    g = sns.clustermap(
        df_scaled,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        cmap="RdBu_r",
        figsize=(12, 10),
        linewidths=0.5,
        #annot=True,           
        fmt=".2f",            
        vmin=-3, vmax=3,
        cbar_kws={"label": "Standardized Value"}
    )

    # --- labels ---
    g.ax_heatmap.set_xlabel("Features")
    g.ax_heatmap.set_ylabel("Proteins")

    #g.fig.canvas.draw()
    #
    ## Now rotate labels
    g.ax_heatmap.set_yticklabels(
        g.ax_heatmap.get_yticklabels(),
        rotation=0,
        va="center"
    )
    #row_order = g.dendrogram_row.reordered_ind
    ## apply labels manually
    #g.ax_heatmap.set_yticks(range(len(row_order)))
    #g.ax_heatmap.set_yticklabels(
    #    [df_scaled.index[i] for i in row_order],
    #    rotation=90
    #)
    g.ax_heatmap.set_xticklabels(
        g.ax_heatmap.get_xticklabels(),
        rotation=45
    )
    #plt.setp(g.ax_heatmap.get_yticklabels(), rotation=90)
    plt.title("Clustered Heatmap", pad=80)
    
    buf = io.BytesIO()
    g.savefig(buf, format="png", bbox_inches="tight")
    plt.close()
    
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode("utf-8")

    # --- save ---
    if save_path:
        plt.savefig(save_path, bbox_inches="tight")

    return f"data:image/png;base64,{encoded}"

def plot_heatmap(df_scaled, ids, feats, save_path=None):
    """
    Plotly based heatmap

    Args:
        df_scaled (pd.DataFrame): Preprocessed Input dataframe.
        ids (list[str]): List of protein IDs (or substrings to match index).
        feats (list[str]): List of columns to include.
        save_path (str, optional): Path to save interactive HTML file.
    """

    mask = df_scaled.index.to_series().isin(ids)
    index_series = df_scaled.index.to_series()
    #df_scaled, info = preprocess_features(df, selected_features=feats)
    df_scaled = df_scaled.loc[mask]
    df_scaled = df_scaled[feats]
    not_included_list = list(set(ids) - set(index_series))
    if len(not_included_list) > 0:
        print(f'{", ".join(not_included_list)} not in PROFET')
    
    # --- hierarchical clustering to reorder rows & columns ---
    row_linkage = linkage(df_scaled, method="average", metric="euclidean")
    col_linkage = linkage(df_scaled.T, method="average", metric="euclidean")

    row_order = leaves_list(row_linkage)
    col_order = leaves_list(col_linkage)

    df_clustered = df_scaled.iloc[row_order, col_order]

    # --- build interactive plotly heatmap ---
    fig = px.imshow(
        df_clustered,
        color_continuous_scale="RdBu_r",
        labels=dict(x="Features", y="Proteins", color="Standardized Value"),
        zmin=0, zmax=1,
        aspect="auto"
    )

    # --- hover info ---
    fig.update_traces(
        hovertemplate="Protein: %{y}<br>Feature: %{x}<br>Value: %{z:.3f}<br>Right-click for info<extra></extra>"
    )

    # --- layout ---
    fig.update_layout(
        title="Clustered Heatmap",
        xaxis_title="Features",
        yaxis_title="Proteins",
        width=max(1000, 60 * len(df_scaled.columns)),
        height=max(600, 20 * len(df_scaled)),
        margin=dict(l=120, r=120, t=80, b=80)
    )

    # --- save optional HTML ---
    if save_path:
        fig.write_html(save_path)

    #fig.show()
    return fig

# ------PROTEIN SIMILARITY NETWORKS------
# Minimum number of variance-carrying features required for a meaningful cosine network.
# With 1-2 features cosine similarity is degenerate (it can only take values near +/-1).
MIN_SIMILARITY_FEATURES = 3


def _empty_message_figure(message):
    '''Return an empty figure that just displays a centered message. Used by the safety lock so
    the Dash callback always receives a valid figure instead of a broken / empty / NaN network.'''
    fig = go.Figure()
    fig.add_annotation(text=message, xref="paper", yref="paper", x=0.5, y=0.5,
                       showarrow=False, font=dict(size=14), align="center")
    fig.update_layout(
        title="Protein Similarity Network",
        title_x=0.5,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        margin=dict(l=20, r=20, t=60, b=20)
    )
    return fig


def plot_protein_similarity_network(df_scaled, ids, feats=None, similarity_threshold=0.5, save_file=None, ids2=None):
    '''Expects an already-standardized dataframe (features centered to mean 0 and scaled to unit
    variance across all proteins, e.g. numeric_df_scaled). Standardization removes the shared
    feature offset so cosine similarity reflects directional (profile) differences instead of
    collapsing toward 1; on standardized data cosine is equivalent to the Pearson correlation.

    If ids2 is provided, nodes in `ids` are coloured as one group and nodes in `ids2` as another,
    with a legend; otherwise all nodes share a single colour (backward-compatible).'''
    if feats:
        df_scaled = df_scaled[feats]

    # --- SAFETY LOCK ---------------------------------------------------------------------------
    # (1) Drop zero / near-zero-variance features (global std over all proteins). They carry no
    #     information and are numerically unstable for cosine. Variance is assessed before subsetting.
    stds = df_scaled.std()
    df_scaled = df_scaled.loc[:, stds > 1e-8]
    if df_scaled.shape[1] == 0:
        return _empty_message_figure(
            "Selected features have no variance, so a similarity network cannot be computed.<br>"
            f"Please select at least {MIN_SIMILARITY_FEATURES} informative features."
        )

    # Subset to the requested proteins (standardization is already done over the full proteome).
    group2 = list(ids2) if ids2 else []
    all_ids = list(ids) + group2
    df_scaled = df_scaled.loc[df_scaled.index.to_series().isin(all_ids)]

    # (2) Drop features that are (near-)constant *within the selected proteins*; the global filter
    #     above cannot see this. A within-subset constant column only injects a shared offset.
    sub_stds = df_scaled.std()
    df_scaled = df_scaled.loc[:, sub_stds > 1e-8]

    # (3) Require a minimum number of features: cosine similarity is degenerate (only ~+/-1) for
    #     1-2 dimensions and would yield a meaningless, over-connected network.
    if df_scaled.shape[1] < MIN_SIMILARITY_FEATURES:
        return _empty_message_figure(
            f"Only {df_scaled.shape[1]} informative feature(s) for the selected proteins.<br>"
            f"Cosine similarity needs at least {MIN_SIMILARITY_FEATURES}."
        )

    # (4) Drop proteins whose feature vector has ~zero norm: cosine is undefined (0/0 -> NaN) for a
    #     zero vector, which would silently corrupt the network.
    norms = np.linalg.norm(df_scaled.values, axis=1)
    df_scaled = df_scaled.loc[norms > 1e-12]
    if df_scaled.shape[0] < 2:
        return _empty_message_figure(
            "Fewer than two proteins have a non-degenerate feature vector; cannot build a network."
        )

    similarity_matrix = 1 - squareform(pdist(df_scaled, metric='cosine'))
    # (5) NaN backstop: never let a NaN similarity reach the edge-threshold logic below.
    similarity_matrix = np.nan_to_num(similarity_matrix, nan=0.0)
    node_labels = df_scaled.index.to_series()

    G = nx.Graph()
    for protein in node_labels:
        G.add_node(protein)

    for i in range(len(node_labels)):
        for j in range(i + 1, len(node_labels)):
            if similarity_matrix[i, j] >= similarity_threshold:
                G.add_edge(node_labels.iloc[i], node_labels.iloc[j], weight=similarity_matrix[i, j])

    # GRAPH
    pos = nx.spring_layout(G, seed=42)

    edge_x, edge_y = [], []
    for u, v in G.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        mode='lines',
        hoverinfo='none',
        showlegend=False,
        line=dict(width=0.5, color='#888')
    )

    if group2:
        # Two coloured groups with a legend.
        group2_set = set(group2)
        node_traces = []
        for name, color, members in [
            ("Group 0", "#1f77b4", [n for n in G.nodes() if n not in group2_set]),
            ("Group 1", "#d62728", [n for n in G.nodes() if n in group2_set]),
        ]:
            node_traces.append(go.Scatter(
                x=[pos[n][0] for n in members],
                y=[pos[n][1] for n in members],
                mode='markers+text',
                text=[str(n) for n in members],
                name=name,
                hoverinfo='text',
                textposition="top center",
                marker=dict(size=12, color=color, line_width=2)
            ))
        fig = go.Figure(data=[edge_trace] + node_traces)
        show_legend = True
    else:
        node_x, node_y, node_text = [], [], []
        for node in G.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            node_text.append(f"{node}")

        node_trace = go.Scatter(
            x=node_x,
            y=node_y,
            mode='markers+text',
            text=node_text,
            hoverinfo='text',
            textposition="top center",
            marker=dict(size=12, color='skyblue', line_width=2)
        )
        fig = go.Figure(data=[edge_trace, node_trace])
        show_legend = False

    fig.update_layout(
        title="Protein Similarity Network",
        title_x=0.5,
        showlegend=show_legend,
        hovermode='closest',
        margin=dict(l=20, r=20, t=60, b=20)
    )

    if save_file:
        fig.write_html(save_file)

    #fig.show()
    return fig


# ------PROTEIN GROUP SCATTERPLOT------
def compare_feature_distributions(df1, df2, name1="Group1", name2="Group2"):
    """
    Compare feature distributions between two DataFrames with same numeric columns.

    Returns a summary DataFrame with:
      - means, stds for each group
      - absolute mean difference
      - relative difference
      - Cohen's d (using weighted pooled std)
      - rank by relative difference
      - rank by effect size
    """
    # --- compute stats ---
    n1, n2 = df1.shape[0], df2.shape[0]
    mean1, mean2 = df1.mean(axis=0), df2.mean(axis=0)
    std1, std2 = df1.std(axis=0), df2.std(axis=0)
    n = n1 + n2

    # absolute and relative differences
    abs_diff = (mean1 - mean2).abs()
    relative_diff = abs_diff / ((mean1.abs() + mean2.abs()) / 2)

    # pooled std (weighted)
    pooled_std = np.sqrt(
        ((n1 - 1) * std1**2 + (n2 - 1) * std2**2) / (n1 + n2 - 2)
    )

    correction = 1 - (3 / (4 * (n1 + n2) - 9))
    correction_durlak = ((n-3)/(n-2.25)) * np.sqrt((n-2)/n)
    pooled_std = pooled_std.apply(lambda x: np.nan if x==0 else x)

    # Cohen's d (effect size)
    cohen_d = (mean1 - mean2) / pooled_std
    hedge_g = cohen_d * correction

    # summary table
    summary = pd.DataFrame({
        f"{name1}_mean": mean1,
        f"{name2}_mean": mean2,
        f"{name1}_std": std1,
        f"{name2}_std": std2,
        "abs_mean_diff": abs_diff,
        "relative_diff": relative_diff,
        "cohen_d": cohen_d,
        "hedge_g": hedge_g
    })

    metric = "cohen_d" if n >= 50 else "hedge_g"
    # ranks (smaller rank = larger difference)
    summary["rank_by_relative_diff"] = summary["relative_diff"].abs().rank(ascending=False)
    summary["rank_by_effect_size"] = summary["cohen_d"].abs().rank(ascending=False) if n >= 50 else summary["hedge_g"].abs().rank(ascending=False)

    # order by effect size (most different features first)
    summary.sort_values("cohen_d", key=lambda x: x.abs(), ascending=False, inplace=True) if n >= 50 else summary.sort_values("hedge_g", key=lambda x: x.abs(), ascending=False, inplace=True)

    return summary, metric



def feature_difference_helper(
    summary,
    metric,
    name1="Group1",
    name2="Group2",
    top_n=10,
    save_path=None
):
    """
    Interactive scatter plot of relative difference vs Cohen's d.
    Top-N features are highlighted and labeled.
    """

    # Sort by absolute effect size
    summary_sorted = summary.copy()
    metric_title = "Cohen's d" if metric == 'cohen_d' else "Hedge's g"
    if metric == "cohen_d":
        summary_sorted["abs_d"] = summary_sorted["cohen_d"].abs()
    else:
        summary_sorted["abs_d"] = summary_sorted["hedge_g"].abs()
    summary_sorted = summary_sorted.sort_values("abs_d", ascending=False)

    # Split data
    top_df = summary_sorted.head(top_n)
    rest_df = summary_sorted.iloc[top_n:]

    fig = go.Figure()

    # --- Background points (hover only) ---
    fig.add_trace(go.Scatter(
        x=rest_df["relative_diff"],
        y=rest_df[metric],
        mode="markers",
        name="Other features",
        marker=dict(
            size=7,
            #color="lightgray",
            color="rgba(99, 110, 250, 0.25)",
            opacity=0.6
        ),
        text=rest_df.index,
        hovertemplate=(
            f"<b>%{{text}}</b><br>"
            f"Relative diff: %{{x:.3f}}<br>"
            f"{metric_title}: %{{y:.3f}}<br>Right-click for info<extra></extra>"
        )
    ))

    # highlighted
    fig.add_trace(go.Scatter(
        x=top_df["relative_diff"],
        y=top_df[metric],
        mode="markers+text",
        name=f"Top {top_n} features",
        marker=dict(
            size=10,
            color="crimson",
            opacity=0.9
        ),
        text=top_df.index,
        textposition="top center",
        hovertemplate=(
            f"<b>%{{text}}</b><br>"
            f"Relative diff: %{{x:.3f}}<br>"
            f"{metric_title}: %{{y:.3f}}<br>Right-click for info<extra></extra>"
        )
    ))

    # Reference lines
    fig.add_hline(y=0, line_dash="dash", line_color="gray")
    fig.add_vline(x=0, line_dash="dash", line_color="gray")

    fig.update_layout(
        title=f"Feature Differences: {name1} vs {name2}",
        xaxis_title="Relative Difference",
        yaxis_title=f"{metric_title} (Effect Size)",
        template="plotly_white",
        width=850,
        height=650,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
        meta=dict(
            plot_info="pgc",
            selected=top_df.index.to_list()
        )
    )

    if save_path:
        fig.write_html(save_path)

    #fig.show()
    return fig


def plot_feature_difference(df, ids1, ids2=None, top_n=10, name1="Group1", name2="Group2", save_path=None):
    data = df.copy()
            
    g1_df = data.loc[data.index.to_series().isin(ids1)].copy()
    
    if ids2 is not None:
        g2_df = data.loc[data.index.to_series().isin(ids2)].copy()
    else:
        g2_df = data.loc[~data.index.to_series().isin(ids1)].copy()
        name2 = "Human Proteome"

    summary, metric = compare_feature_distributions(g1_df, g2_df, name1, name2)

    fig = feature_difference_helper(summary, metric, name1, name2, top_n=top_n, save_path=save_path)
    return fig
    
# CORRELATION

def plot_correlation(df, features, thr=0.9):
    if features:
        df = df[features].copy()
    c_mat = df.corr(method='spearman').abs()
    
    mask = np.triu(np.ones_like(c_mat, dtype=bool), k=1)
    c_mat_masked = c_mat.where(~mask)
    
    fig = go.Figure(
        go.Heatmap(
            z=c_mat_masked.values,
            x=list(c_mat.columns),
            y=list(c_mat.index),
            zmin=0, zmax=1,
            zmid=0.5,
            colorbar=dict(
                title=dict(text="Spearman ρ", side="right"),
                thickness=15,
                len=0.6,
                ),
            hovertemplate="Feature X: %{x}<br>Feature Y: %{y}<br>ρ: %{z:.3f}<br>Right-click for info<extra></extra>",
        )
    )
    fig.update_layout(
        meta=dict()
    )
    fig.update_layout(
        title=dict(
            text=f"Spearman Rank Correlation Matrix<br>Redundancy threshold: |ρ| > {thr}",
            font=dict(size=14),
            x=0.5,
        ),
        width=1200,
        height=1100,
        xaxis=dict(
            tickangle=45,
            showgrid=False,
            tickfont=dict(size=8),
        ),
        yaxis=dict(
            showgrid=False,
            tickfont=dict(size=8),
            autorange="reversed",
        ),
        paper_bgcolor="white",
        plot_bgcolor="white",
        meta={"plot_info": "corr"}
    )
    
    return fig
    

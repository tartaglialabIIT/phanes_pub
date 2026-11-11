"""Loads the dataframes"""

import joblib
import pandas as pd


def generate_dfs() -> tuple[pd.DataFrame, pd.DataFrame]:
    numeric_pkl_path = './profet_numeric_v6_clean.pkl'
    profiles_pkl_path = './profet_numeric_v6_profiles.pkl'
    
    numeric_df = joblib.load(numeric_pkl_path)
    profiles_df = joblib.load(profiles_pkl_path)
    
    return numeric_df, profiles_df


def generate_scalar_df() -> pd.DataFrame:
    numeric_pkl_path = './profet_numeric_v6_clean.pkl'
    
    numeric_df = joblib.load(numeric_pkl_path)
    
    return numeric_df
    
    
def generate_preproccesed_scalar_df() -> pd.DataFrame:
    pkl_path = './profet_numeric_v6_preprocessed.pkl'
    
    scaled_df = joblib.load(pkl_path)
    
    return scaled_df
    
    
def generate_preproccesed_scalar_noScaling_df() -> pd.DataFrame:
    pkl_path = './profet_numeric_v6_preprocessed_noScaling.pkl'
    
    scaled_df = joblib.load(pkl_path)
    
    return scaled_df
    
    
def generate_profiles_df() -> pd.DataFrame:
    pkl_path = './profet_numeric_v6_profiles.pkl'
    
    profiles_df = joblib.load(pkl_path)
    
    return profiles_df
    
def generate_mapping_df() -> pd.DataFrame:
    pkl_path = './uniprot_id_gene_name.pkl'
    
    mapping_df = joblib.load(pkl_path)
    
    return mapping_df
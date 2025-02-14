#!/usr/bin/env python3

"""
Name = Jori de Leuw
Description = this script learns and test machine learning models
"""
import pandas as pd
import joblib
from scipy.stats import zscore

def feature_table_merger(df1, df2):
    """This function combines two pandas dataframes based on index

    Args:
        df1 (pd df): a pandas dataframe with a CLASS column
        df2 (pd df): a pandas dataframe with a CLASS column

    Returns:
        pd df: the merged df
    """
    #merge dfs based on index
    merged_df = pd.merge(df1, df2, left_index=True, right_index=True)
    return merged_df

def orphan_classifier(ml_model, data):
    """This function loads the machine learning model and makes predictions

    Args:
        ml_model (str): path to the machine learning model
        data (pd df): pandas dataframe containing the features

    Returns:
        pred: predictions made by the model
    """
    loaded_model = joblib.load(ml_model)
    pred = loaded_model.predict(data)
    return pred

def main(diamond_file, sequence_file, model):
    """This function reads the feature tables, normalizes it and makes predictions

    Args:
        diamond_file (str): path to diamond feature table
        sequence_file (str): path to sequence feature table
        model (str): path to machine learning model

    Returns:
        predictions: predictions made by the model
    """
    #load features
    diamond_features = pd.read_table(diamond_file, index_col=0)
    sequence_features = pd.read_table(sequence_file, index_col=0)
    #normalize features (this has to be done manually, to prevent error with nan)
    diamond_features = diamond_features.apply(zscore, axis = 1)
    sequence_features = sequence_features.apply(zscore, axis = 1)
    #merge features
    combined_features = pd.merge(diamond_features, sequence_features, left_index=True, right_index=True)
    predictions = orphan_classifier(model, combined_features)
    return predictions
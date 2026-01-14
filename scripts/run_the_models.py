#!/usr/bin/env python3
import pandas as pd
import joblib
from scipy.stats import zscore
import sys
import numpy as np

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

def main():
    """This function reads the feature tables, normalizes it and makes predictions

    Args:
        diamond_file (str): path to diamond feature table
        model (str): path to machine learning model

    Returns:
        predictions: predictions made by the model
    """
    diamond_file = sys.argv[1]
    model = sys.argv[2]
    model_name = sys.argv[2].split('/')[-1].split('.')[0].split('_',1)[1]
    #load features
    diamond_features = pd.read_table(diamond_file, index_col=0)
    #sequence_features = pd.read_table(sequence_file, index_col=0)
    #normalize features (this has to be done manually, to prevent error with nan)
    diamond_features = diamond_features.apply(zscore, axis = 0)
    #sequence_features = sequence_features.apply(zscore, axis = 1)
    #merge features
    #combined_features = pd.merge(diamond_features, sequence_features, left_index=True, right_index=True)
    predictions = orphan_classifier(model, diamond_features)
    print(predictions)
    np.savetxt(model_name + '_out.txt', predictions)
    #with open(model_name + '_output.txt', 'w') as f:
     #   f.write(str(predictions))
if __name__ == "__main__":
   main()
	

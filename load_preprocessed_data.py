"""
loader for preprocessed data (computed via analysis_main.py)
just point it at your data folder and get everything back as a dict (see implementation in figures.ipynb)
"""

data_dir='preprocessed/' # put dir where you wrote the intermediate data to 

import os
import json
import pickle
import pandas as pd
import xarray as xr
import geopandas as gpd
import numpy as np


def load_data(data_dir=data_dir):
    """
    load all the preprocessed data into one big dictionary
    handles all the common file types automatically
    """
    
    if not os.path.exists(data_dir):
        raise FileNotFoundError(f"can't find {data_dir} - did you run main_analysis.py yet?")
    
    data = {}
    files = [f for f in os.listdir(data_dir) if not f.startswith('.')]
    
    for file in files:
        name = os.path.splitext(file)[0]
        path = os.path.join(data_dir, file)
        
        try:
            # nc files 
            if file.endswith('.nc'):
                ds = xr.open_dataset(path)
                # if it's just one variable, extract the array
                if len(ds.data_vars) == 1:
                    data[name] = ds[list(ds.data_vars)[0]]
                else:
                    data[name] = ds
            
            # csv files (malaria data, correlations)
            elif file.endswith('.csv'):
                data[name] = pd.read_csv(path)
            
            # parquet files (gdfs w/ geometry)
            elif file.endswith('.parquet'):
                data[name] = gpd.read_parquet(path)
            
            # numpy arrays
            elif file.endswith('.npy'):
                data[name] = np.load(path)
            
            # json files (some metadata, labels)
            elif file.endswith('.json'):
                with open(path) as f:
                    data[name] = json.load(f)
            
            # pickle files (random items like year strings)
            elif file.endswith('.pkl'):
                with open(path, 'rb') as f:
                    data[name] = pickle.load(f)
                    
        except Exception as e:
            print(f"couldn't load {file}: {e}\npls check file names")
    
    print(f"loaded {len(data)} datasets")
    return data
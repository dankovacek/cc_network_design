# calculate_distances.py

import os
import time

import pandas as pd
import numpy as np
import xarray as xr

import scipy.stats as st

from itertools import combinations

import shapely.geometry as geom
import geopandas as gpd

import multiprocessing as mp
from multiprocessing import Pool


from notification_alert import client


# hysets_folder = '/media/danbot/T7 Touch/hysets_series/'
# if production:
#     hysets_folder = 'data/hysets_series/'

hysets_df = pd.read_csv('data/HYSETS_watershed_properties.txt', sep=';', dtype={'Official_ID': str})

# create the centroid geometry
hysets_df['centroid_geom'] = hysets_df.apply(lambda xy: geom.Point((xy['Centroid_Lon_deg_E'], xy['Centroid_Lat_deg_N'])), axis=1)


# create a dictionary of identifying information to facilitate
# selection of specific watersheds
basin_metadata = ['Watershed_ID', 'Official_ID', 'Name']

basin_centroid_geom = ['centroid_geom']

basin_characteristics_cols = ['Centroid_Lat_deg_N', 'Centroid_Lon_deg_E',
              'Drainage_Area_GSIM_km2', 'Drainage_Area_km2',
              'Elevation_m', 'Gravelius', 'Perimeter', 'Slope_deg', 'Aspect_deg', 
              'Land_Use_Crops_frac', 'Land_Use_Forest_frac', 'Land_Use_Grass_frac',
              'Land_Use_Shrubs_frac', 'Land_Use_Snow_Ice_frac', 'Land_Use_Urban_frac',
              'Land_Use_Water_frac', 'Land_Use_Wetland_frac', 
              'Permeability_logk_m2', 'Porosity_frac']

hysets_dict = hysets_df[basin_metadata + basin_centroid_geom + basin_characteristics_cols].set_index('Official_ID').to_dict(orient='index')


results_fnames = os.listdir('results/')
results_fnames = [f'results/{e}' for e in results_fnames if "COMBINED_pairs" in e]

def calculate_pair_centroid_distance(pair):
    pair_df = hysets_df[hysets_df['Official_ID'].isin(pair)]
    hdf = gpd.GeoDataFrame(pair_df, geometry=foo['centroid_geom'], crs='EPSG:4326')
    hdf = hdf.to_crs(3005)
    hdf.reset_index(inplace=True)
    return hdf.loc[0, 'geometry'].distance(hdf.loc[1, 'geometry']) / 1000


def create_line(row):
    return geom.LineString([hysets_dict[row['b1']]['centroid_geom'], hysets_dict[row['b2']]['centroid_geom']])


all_dfs = []

for r in results_fnames:

    print(r)
    df = pd.read_pickle(r)
        
    geometry = gpd.GeoDataFrame({'geometry': df.apply(lambda row: create_line(row), axis=1)}, crs='EPSG:4326')

    geometry = geometry.to_crs(3005)
    geometry['centroid_distance_km'] = geometry.length / 1000  # convert to km

    df['distance_btwn_centroids_km'] = geometry['centroid_distance_km']
    df = df[df['distance_btwn_centroids_km'] < 1000]
    
    all_dfs.append(df)

combined_df = pd.concat(all_dfs)
print(combined_df.head())
print(len(combined_df))
combined_df.to_pickle('results/FINAL_SET_HYSETS_PAIRS.csv')

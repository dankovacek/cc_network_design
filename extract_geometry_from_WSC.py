import os

import pandas as pd
import numpy as np
import fiona

from multiprocessing import Pool

import sqlite3

from shapely.geometry import Point
import geopandas as gpd


# import basin characteristics
WSC_db_folder = '/media/danbot/T7 Touch/hydat_db/'
metadata_fn = 'WSC_Stations_Master.csv'

df = pd.read_csv(WSC_db_folder + metadata_fn)
df.head()

df['geometry'] = df.apply(lambda row: Point(row['Longitude'], row['Latitude']), axis=1)
df['num_years_record'] = df['Year To'] - df['Year From']
# filter for stations in BC and Alberta
df = df[df['Province'].isin(['BC', 'AB'])]

stn_pairs_list = np.random.choice(df['Station Number'].to_numpy(), size=(int(1E4), 2), replace=True)
stn_pairs_list = [list(sorted(e)) for e in stn_pairs_list]

def check_time_periods(p):
    stn1, stn2 = p[0], p[1]
    s1 = df[df['Station Number']==stn1]
    s2 = df[df['Station Number']==stn2]
    start1, end1 = s1['Year From'].to_numpy()[0], s1['Year To'].to_numpy()[0]
    start2, end2 = s2['Year From'].to_numpy()[0], s2['Year To'].to_numpy()[0]
    if end1 < end2:
        overlap_duration = end1 - start2
    else:
        overlap_duration = end2 - start1
    if overlap_duration > 50:
        return p
    else:
        return None

# filter for pairs that have minimum 50 years of concurrent data

pool = Pool()
overlapping_records = pool.map(check_time_periods, stn_pairs_list)
pool.close()
pool.join()



points_df = gpd.GeoDataFrame(df, crs='EPSG:4326')

# convert to BC albers for geometric distance operation
points_df = points_df.to_crs(3005)

Hysets_series_folder = '/media/danbot/T7 Touch/hysets_series/'
hysets_files = sorted(os.listdir(Hysets_series_folder))

filtered_pairs = [e for e in overlapping_records if e is not None]
filtered_pairs = [pair for pair in filtered_pairs if (pair[0] + '.csv' in hysets_files) & (pair[1] + '.csv' in hysets_files)]
print('n filtered pairs:')
print(len(filtered_pairs))
print('')

WSC_basin_polygons_path = '/media/danbot/T7 Touch/hydat_db/'

# get all wsc_catchment data into its own dataframe
gdb_path = os.path.join(WSC_basin_polygons_path, 'WSC_Basins.gdb.zip')
all_layers = fiona.listlayers(gdb_path)
all_layer_names = [e.split('_')[1].split('_')[0] for e in all_layers]

geometry_dict = {}
with fiona.open(gdb_path) as gdb_file:
    for l in all_layers[:1]:
        foo = next(iter(gdb_file))
        stn_id = foo['properties']['Station']
        geom = foo['geometry']['coordinates']
        print(geom)
        gdf = gpd.GeoDataFrame()
        geometry_dict[stn_id] = {'centroid': geom.centroid} 
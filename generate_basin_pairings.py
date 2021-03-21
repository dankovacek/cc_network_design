import os
import time

import pandas as pd
import numpy as np
import xarray as xr

import scipy.stats as st

from itertools import combinations

from shapely.geometry import Point
import geopandas as gpd

from multiprocessing import Pool

from notification_alert import client

# import basin characteristics
WSC_db_folder = '/media/danbot/T7 Touch/hydat_db/'
# metadata_fn = 'WSC_Stations_Master.csv'
hysets_folder = '/media/danbot/T7 Touch/hysets_series/'
camels_folder = '/media/danbot/T7 Touch/camels_db/usgs_streamflow/'

# create a list of all pairs of basins in the camels database
# note that all filenames are 8 digits long
camels_files, camels_files_stations = [], []
for f in os.listdir(camels_folder):
    camels_files += [camels_folder + f + '/' + e for e in os.listdir(camels_folder + f)]
    camels_files_stations += [e.split('_')[0] for e in os.listdir(camels_folder + f)]

hysets_df = pd.read_csv('data/HYSETS_watershed_properties.txt', sep=';', dtype={'Official_ID': str})

# check drainage area flag information and filter for basins < threshold area
basin_area_threshold = 2000
print(f' {len(hysets_df)} hysets basins')
hysets_df = hysets_df[hysets_df['Drainage_Area_km2'] <= basin_area_threshold]
print(f' {len(hysets_df)} hysets basins < {basin_area_threshold} km^2')

# filter out basins with large area discrepancies between area calc methods
def check_area(row):
    """If there is no area for either method, return false
    If the relative area difference between methods is > 0.15, return False,
    Otherwise return True."""
    basin = row['Official_ID']
    da1_check = np.isnan(row['Drainage_Area_km2'])
    da2_check = np.isnan(row['Drainage_Area_GSIM_km2'])
    if da1_check & da2_check:
        return False
    if ~da2_check:
        da1 = row['Drainage_Area_km2']
        da2 = row['Drainage_Area_GSIM_km2']
        pct_diff = np.abs(da1 - da2) / ((da1 + da2)/2)
        
        if pct_diff > 0.15:
            print(f'{basin}: {pct_diff} area difference')
            
            return False

    return True

hysets_df['area_check'] = hysets_df.apply(lambda row: check_area(row), axis=1)
hysets_df = hysets_df[hysets_df['area_check']]
hysets_df.drop(labels=['area_check'], axis=1, inplace=True)

# strip leading zeros from usgs official IDs
def format_usgs_id(row):
    basin_id = row['Official_ID']
    if row['Source'] != 'USGS':
        return basin_id
    if len(basin_id) == 8:
        return basin_id

    basin_id = row['Official_ID'].lstrip()
    if len(basin_id) == 8:
        return basin_id
    return False

# for all USGS stations, only keep the ones that are in the camels set
# which correspond to all those with 8 digit Official ID.
hysets_df['Official_ID'] = hysets_df.apply(lambda row: format_usgs_id(row), axis=1)
hysets_df = hysets_df[hysets_df['Official_ID'] != False]

# usgs_ids = hysets_df.loc[hysets_df['Source'] == 'USGS', ['Official_ID']].to_numpy().flatten()
# print(usgs_ids[:10])
# print(list(set([len(e) for e in usgs_ids])))
# print('')
# print('_____________-')
# print(asdfsd)

# extract the HYDAT subset
hydat_df = hysets_df[hysets_df['Source'] == 'HYDAT'].copy()
# extract the USGS subset
usgs_df = hysets_df[hysets_df['Source'] == 'USGS'].copy()
# extract the Mexico subset
# mex_df = hysets_df[hysets_df['Source'] == 'Mexico'].copy()

# create a centroid shapely Point
hydat_df['centroid_geom'] = hydat_df.apply(lambda xy: Point((xy['Centroid_Lon_deg_E'], xy['Centroid_Lat_deg_N'])), axis=1)
usgs_df['centroid_geom'] = usgs_df.apply(lambda xy: Point((xy['Centroid_Lon_deg_E'], xy['Centroid_Lat_deg_N'])), axis=1)
# mex_df['centroid_geom'] = mex_df.apply(lambda xy: Point((xy['Centroid_Lon_deg_E'], xy['Centroid_Lat_deg_N'])), axis=1)

# create a dictionary of identifying information to facilitate
# selection of specific watersheds
basin_metadata = ['Watershed_ID', 'Official_ID', 'Name']

basin_centroid_geom = ['centroid_geom']

basin_characteristics_cols = ['Centroid_Lat_deg_N', 'Centroid_Lon_deg_E',
              'Drainage_Area_km2',
              'Elevation_m', 'Gravelius', 'Perimeter', 'Slope_deg', 'Aspect_deg', 
              'Land_Use_Crops_frac', 'Land_Use_Forest_frac', 'Land_Use_Grass_frac',
              'Land_Use_Shrubs_frac', 'Land_Use_Snow_Ice_frac', 'Land_Use_Urban_frac',
              'Land_Use_Water_frac', 'Land_Use_Wetland_frac', 
              'Permeability_logk_m2', 'Porosity_frac']

area_flag_col = 'Drainage_Area_GSIM_km2'


hydat_dict = hydat_df[basin_metadata + basin_centroid_geom + basin_characteristics_cols].set_index('Official_ID').to_dict(orient='index')

usgs_dict = usgs_df[basin_metadata + basin_centroid_geom + basin_characteristics_cols].set_index('Official_ID').to_dict(orient='index')
# mex_dict = mex_df[basin_metadata + basin_centroid_geom + basin_characteristics_cols].set_index('Official_ID').to_dict(orient='index')


hydat_stns = list(hydat_dict.keys()) 
camels_stns = np.intersect1d(list(usgs_dict.keys()), camels_files_stations)

# mex_stns = list(mex_dict.keys())
# filter the usgs list for just the stations in CAMELS
print(f'There are {len(hydat_stns)} HYDAT station records in the HYSETS database.')
print(f'There are {len(camels_stns)} CAMELS station records in the HYSETS database.')
# print(f'There are {len(mex_stns)} Mexico station records in the HYSETS database.')



all_properties_dict = {'hydat': hydat_dict,
                      'camels': usgs_dict,
                    #    'mex': mex_dict,
                      }


def check_pair_properties(row, basin_dict):
    """
    
    """
    bid = str(row['Official_ID'])

    all_char_status = [np.isnan(all_properties_dict[basin_dict][bid][c]) for c in basin_characteristics_cols]
    
    missing_chars = [c for c in basin_characteristics_cols if np.isnan(all_properties_dict[basin_dict][bid][c])]
    
    # if any characteristics are True, the corresponding value is nan,
    # so return False, otherwise True because the set is complete
    if sum(all_char_status) != 0:
        # print(f'{sum(all_char_status)} characteristics are missing for {bid}')
        # print(f'    {missing_chars}')
        return False
    else:
        return True


# Get just the CAMELS subset from the USGS set
usgs_df['camels_check'] = [True if e[-8:] in camels_stns else False for e in usgs_df['Official_ID']]
camels_df = usgs_df[usgs_df['camels_check']].copy()
camels_df.drop(labels=['camels_check'], axis=1, inplace=True)

## Filter out basins with incomplete characteristics
hydat_df['char_check'] = hydat_df.apply(lambda row: check_pair_properties(row, 'hydat'), axis=1)
hydat_df = hydat_df[hydat_df['char_check']]
hydat_df.drop(labels=['char_check'], axis=1, inplace=True)

camels_df['char_check'] = camels_df.apply(lambda row: check_pair_properties(row, 'camels'), axis=1)
camels_df = camels_df[camels_df['char_check']]
camels_df.drop(labels=['char_check'], axis=1, inplace=True)

# hydat_stn_pairs_list = list(combinations(hydat_stns, 2))
# camels_stn_pairs_list = list(combinations(camels_stns, 2))
# mex_stn_pairs_list = list(combinations(mex_stns, 2))
combined_stations = list(camels_df['Official_ID'].to_numpy()) + list(hydat_df['Official_ID'].to_numpy())
print(f'After filtering, there are {len(combined_stations)} basins. {len(camels_df)} + {len(hydat_df)}')

# generate all pairs of stations
combined_pairs_list = list(combinations(combined_stations, 2))
print(f'{len(combined_pairs_list)} basin pairs will be generated for comparison.')

def extract_streamflow_series(stn):
    df = pd.read_csv(f'{hysets_folder}{stn}.csv', index_col=['time'])
    df.dropna(inplace=True)
    return df

def get_concurrence_length_and_COD(pair):
    df1 = extract_streamflow_series(pair[0])
    df1.rename(mapper={'discharge': f'{pair[0]}'}, inplace=True, axis=1)
    
    df2 = extract_streamflow_series(pair[1])
    df2.rename(mapper={'discharge': f'{pair[1]}'}, inplace=True, axis=1)
    concurrent_df = pd.concat([df1, df2], join='inner', axis=1)
    
    if len(concurrent_df) > 364:
        try:
            out = st.linregress(concurrent_df.to_numpy()) 
            return len(concurrent_df), out[2]**2   
        except Exception as ex:
            print('regression attempt failed')
            return len(concurrent_df), np.nan
    else:
        return len(concurrent_df), np.nan


print('')
print('#######################')
print('')

combined_pair_df = pd.DataFrame(combined_pairs_list, columns=['b1', 'b2'])

try:
    t0 = time.time()
    pool = Pool()
    result = pool.map(get_concurrence_length_and_COD, 
                combined_pairs_list[:100])
    pool.close()
    pool.join()
    t1 = time.time()

    concurrent_length_array = [e[0] for e in result]
    cod_array = [e[1] for e in result]

    print(f'Time to calculate concurrent period lengths: {t1 - t0:.1f}')

    combined_pair_df['concurrent_length_days'] = concurrent_length_array
    combined_pair_df['similarity'] = cod_array

    combined_pair_df = combined_pair_df[combined_pair_df['similarity'].isna()]

    # combined_pair_df = combined_pair_df[combined_pair_df['concurrent_days'] > 364]
    print(f'{len(combined_pair_df)} basin pairs meet the concurrence length, basin area, and characteristic information criteria.')

    # write the list of unique pairs to disk so you 
    # don't have to go through that process again
    combined_pair_df.to_pickle('results/COMBINED_pairs_min365d_output.csv')

    t_hours = (t1 - t0) / 3600


    message = client.messages \
                    .create(
                        body=f"Code Run completed in {t_hours:.1f}. Huzzah!",
                        from_='+16048006923',
                        to='+16048420619'
                    )


except Exception as ex:
    msg = str(ex)[:25]
    
    message = client.messages \
                    .create(
                        body=f"Code run failed. {msg}",
                        from_='+16048006923',
                        to='+16048420619'
                    )


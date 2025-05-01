#%% 
'''
Code Purpose: Cross-reference the Sydney Radio Stars Catalogue with BL opendata observations.
Author: Owen A. Johnson 
Date: 2025-04-30 
'''

import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import scienceplots; plt.style.use(['science', 'ieee'])
import astropy.units as u
from astropy.coordinates import SkyCoord
from tqdm import tqdm

# --- User Input ---
verbose = False 
plot = False

# --- Functions --- 
def beam_FWHM(freq, d):
    """
    Calculate beam FWHM.

    Parameters:
    freq : float
        Observing frequency in MHz
    d : float
        Dish diameter in meters

    Returns:
    float
        Beam FWHM in degrees (using factor ~1.2 for typical tapering)
    """
    wavelength = 300 / freq  # Wavelength in meters
    fwhm_rad = 1.2 * wavelength / d  # FWHM in radians
    return fwhm_rad * (180 / np.pi)  # Convert to degrees

# --- Sydney Radio Stars Catalogue ---
SRSC_df = pd.read_csv('SRSC.csv')
print('Number of sources in the SRSC:', len(SRSC_df))

# --- BL opendata observations ---
BL_df = pd.read_csv('BL_opendata_metadata.csv')
print('Number of entries in the BL opendata set:', len(BL_df))

BL_unique_df = BL_df.drop_duplicates(subset='ra_deg')
BL_unique_df = BL_unique_df[(BL_unique_df['dec_deg'] >= -90) & (BL_unique_df['dec_deg'] <= 90)]
BL_unique_df = BL_unique_df[(BL_unique_df['ra_deg'] >= 0) & (BL_unique_df['ra_deg'] < 360)]
BL_unique_df.to_csv('BL_unique_metadata.csv', index=False)

print('Number of unique sources in the BL opendata set:', len(BL_unique_df))

SRSC_coords = SkyCoord(ra=SRSC_df['Archival_ra'].values * u.deg, dec=SRSC_df['Archival_dec'].values * u.deg)
BL_coords = SkyCoord(ra=BL_unique_df['ra_deg'].values * u.deg, dec=BL_unique_df['dec_deg'].values * u.deg)

if verbose:
    print(BL_unique_df['ra_deg'].describe())
    print(BL_unique_df['dec_deg'].describe())

if plot: 
    # --- Aitoff Projection of SRSC and BL opendata sources ---
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection='aitoff')
    ax.grid(True)

    # Plot SRSC sources
    ax.scatter(SRSC_coords.ra.wrap_at(180 * u.deg).radian, SRSC_coords.dec.radian, s=1, color='blue', label='SRSC Sources', zorder = 2)

    # Plot BL opendata sources
    ax.scatter(BL_coords.ra.wrap_at(180 * u.deg).radian, BL_coords.dec.radian, s=1, color='red', label='BL Opendata Sources', zorder = 1, alpha=0.2)

    plt.xlabel('RA (deg)')
    plt.ylabel('Dec (deg)')
    plt.legend()
    
# --- Find SRSC sources that are in the FWHM of BL opendata sources ---
output_rows = []

for bl_source in tqdm(range(len(BL_unique_df))):
    freq = BL_unique_df['center_freq'].iloc[bl_source]
    telescope = str(BL_unique_df['telescope'].iloc[bl_source]).replace(' ', '')

    if telescope == 'Parkes':
        diameter = 64
    elif telescope == 'GBT':
        diameter = 100
    else: 
        continue  

    fwhm = beam_FWHM(freq, diameter)
    separation = SRSC_coords.separation(BL_coords[bl_source])
    
    # --- Check if SRSC source is within FWHM of BL source ---
    mask = SRSC_coords.separation(BL_coords[bl_source]).deg < fwhm
    matched_SRSC = SRSC_df[mask].copy()
    matched_SRSC['separation_deg'] = separation[mask]
    
    for idx, row in matched_SRSC.iterrows():
        output_rows.append({
            'bl_source_id': BL_unique_df['source_name'].iloc[bl_source],
            'bl_ra_deg': BL_unique_df['ra_deg'].iloc[bl_source],
            'bl_dec_deg': BL_unique_df['dec_deg'].iloc[bl_source],
            'telescope': telescope,
            'center_freq': freq,
            'FWHM_deg': fwhm,
            'separation_deg': row['separation_deg'],
            'data_url': BL_unique_df['data_url'].iloc[bl_source],
            'srsc_source_id': row['star_name'],
            'srsc_ra_deg': row['Archival_ra'], 
            'srsc_dec_deg': row['Archival_dec'],
            'detection_telescope': row['Radio_telescope'],
            'DOI': row['Reference'],
            'Radio_flux_peak': row['Radio_I_flux_peak'],
            'Survey_ID': row['Radio_component_id']
        })

CROSSREF_df = pd.DataFrame(output_rows)
CROSSREF_df = CROSSREF_df.sort_values(by='Radio_flux_peak', ascending=False) # sort by flux decending
CROSSREF_df.to_csv('BL_SRSC_crossref.csv', index=False)
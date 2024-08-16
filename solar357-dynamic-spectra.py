'''
Code Purpose: Plot out 357 dynamic spectra from .fil 
Author: Owen A. Johnson 
'''
import argparse
import numpy as np
import your
import matplotlib.pyplot as plt
import scienceplots
from sigpyproc.readers import FilReader
from astropy.time import Time
import os
from scipy.interpolate import interp1d
import matplotlib.dates as mdates
from solar_rfi import *

plt.style.use(['science','ieee', 'no-latex'])

def parse_args():
    """
    Purpose: Parse command line arguments
    """
    parser = argparse.ArgumentParser(description='Plot out 357 dynamic spectra from .fil')
    parser.add_argument('fil', type=str, help='Path to .fil file')
    parser.add_argument('start_time', type=str, help='Start time in HH:MM format')
    parser.add_argument('end_time', type=str, help='End time in HH:MM format')
    args = parser.parse_args()
    return args

def get_time(fil):
    sppFil = FilReader(fil)
    nsamples = sppFil.header.nsamples
    tsamp = sppFil.header.tsamp
    startmjd = sppFil.header.tstart
    time = startmjd + np.linspace(0, nsamples * tsamp, nsamples) / (24 * 60 * 60)
    return time

def mjd_to_date(mjd):
    t = Time(mjd, format='mjd')
    return t.to_datetime()

def get_data(fil, start_idx, end_idx):
    """
    Purpose: Retrieve the data from the filterbank file for the specified range of samples.
    """
    your_object = your.Your(fil)
    data = your_object.get_data(nstart=start_idx, nsamp=(end_idx - start_idx + 1))
    return data

def assign_frequencies(data, mode):
    """
    Assign the frequency range based on the RCU mode.
    """
    if mode == 3:
        freq_start, freq_end = 10, 90  # Mode 3: 30-90 MHz
    elif mode == 5:
        freq_start, freq_end = 110, 190  # Mode 5: 110-190 MHz
    elif mode == 7:
        freq_start, freq_end = 210, 270  # Mode 7: 210-270 MHz
    elif mode == 'gap1':
        freq_start, freq_end = 90, 110
    elif mode == 'gap2':
        freq_start, freq_end = 190, 210
    else:
        raise ValueError("Unsupported RCU mode")

    nchans = data.shape[1]
    freqs = np.linspace(freq_start, freq_end, nchans)
    return freqs

def resample_data(data, target_nchans, original_nchans):
    """
    Resample the 2D data array to match the target number of channels.
    Applies interpolation along the frequency axis for each time sample.
    """
    original_freqs = np.linspace(0, original_nchans, original_nchans)
    target_freqs = np.linspace(0, original_nchans, target_nchans)

    # Interpolate each row (time sample) separately
    resampled_data = np.array([np.interp(target_freqs, original_freqs, row) for row in data])

    return resampled_data

def bkg_subtraction(data): 
    """
    Take the first 10 seconds of data and average it to get the background and divide data by mean 
    """
    bkg_data = data[:10]
    return data / np.mean(bkg_data, axis=0)


def main():
    args = parse_args()
    fil_path = args.fil
    time_mjd = get_time(fil_path) 
    time_utc = mjd_to_date(time_mjd)
    
    header = your.Your(fil_path).your_header 
   
    start_time = Time(time_utc[0].strftime('%Y-%m-%d') + ' ' + args.start_time).mjd
    end_time = Time(time_utc[0].strftime('%Y-%m-%d') + ' ' + args.end_time).mjd

    # Find the index of the start and end time
    start_idx = np.argmin(np.abs(time_mjd - start_time))
    end_idx = np.argmin(np.abs(time_mjd - end_time))

    print('======  Time ======')
    print('Obs Start Time :', time_utc[0])
    print('Obs End Time   :', time_utc[-1])
    print('Plot Start Time:', time_utc[start_idx])
    print('Plot End Time  :', time_utc[end_idx])

    data = get_data(fil_path, start_idx, end_idx)
    print('Data Shape:', data.shape)

    masked_data = rfi_masking(data, your.Your(fil_path), fil_path, rfi_plot_cond=False)

    # Split data based on the modes
    mode3_data = data[:, :200]
    mode5_data = data[:, 200:400]
    mode7_data = data[:, 400:]

    # Background subtraction
    mode3_data = bkg_subtraction(mode3_data)
    mode5_data = bkg_subtraction(mode5_data)
    mode7_data = bkg_subtraction(mode7_data)

    mode7_data = resample_data(mode7_data, 200, 88)  # resample 

    gap_data = np.full((data.shape[0], 50), np.nan) # 0.4 MHz frequency res

    print('=== Data Shapes ===')
    print('Mode 3 Shape:', mode3_data.shape)
    print('Mode 5 Shape:', mode5_data.shape)
    print('Mode 7 Shape:', mode7_data.shape)

    mode3_freqs = assign_frequencies(mode3_data, 3)
    mode5_freqs = assign_frequencies(mode5_data, 5)
    mode7_freqs = assign_frequencies(mode7_data, 7)
    
    gap1_freqs = assign_frequencies(gap_data, 'gap1')
    gap2_freqs = assign_frequencies(gap_data, 'gap2')

    # Print min and max freq for each mode
    print('Mode 3: Min Freq:', mode3_freqs[0], 'Max Freq:', mode3_freqs[-1])
    print('Mode 5: Min Freq:', mode5_freqs[0], 'Max Freq:', mode5_freqs[-1])
    print('Mode 7: Min Freq:', mode7_freqs[0], 'Max Freq:', mode7_freqs[-1])

    # --- Combine Data --- 
    combined_data = np.concatenate((mode3_data, gap_data, mode5_data, gap_data), axis=1)
    combined_freqs = np.concatenate((mode3_freqs, gap1_freqs, mode5_freqs, gap2_freqs))

    # --- Plotting ---
    plt.figure(figsize=(10, 5))
    vmax, vmin = np.nanpercentile(combined_data, (90, 30))
    plt.imshow(combined_data.T, aspect='auto', vmax=vmax, vmin=vmin, extent=[mdates.date2num(time_utc[start_idx]), mdates.date2num(time_utc[end_idx]), combined_freqs[-1], combined_freqs[0]])

    # Format x-axis for UTC time
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    plt.gca().xaxis.set_major_locator(mdates.AutoDateLocator())

    plt.xlabel('Time (UTC)')
    plt.ylabel('Frequency (MHz)')
    plt.yticks([10, 90, 110, 190])
    plt.savefig('dynamic_spectra.png', dpi=300)

if __name__ == '__main__':
    main()

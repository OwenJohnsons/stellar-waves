"""
Script for Plotting Solar Dynamic Spectra from I-LOFAR 357 .fil Files

This script processes and visualizes dynamic spectra from .fil files, specifically targeting
solar observations. The script handles multiple RCU modes (frequency bands), applies background 
subtraction, and masks Radio Frequency Interference (RFI) using the IQRM method. It then generates 
plots for the dynamic spectra across different frequency bands, including handling gaps in the 
spectrum caused by the FM band and TV broadcasting.

Dependencies:
- argparse
- numpy
- your
- matplotlib
- scienceplots
- sigpyproc
- astropy
- scipy
- iqrm

Usage:
    python dynamic_spectra.py <fil_file_path> <start_time_HH:MM> <end_time_HH:MM>

Example:
    python dynamic_spectra.py /path/to/357_2024-08-21T09:00:00.fil 09:00 09:30

Author:
    Owen A. Johnson
    Date: 2024-08-21
"""

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
# from solar_rfi import *
from iqrm import iqrm_mask

# Set the plotting style to use scientific formatting with IEEE standards, without LaTeX
plt.style.use(['science', 'ieee', 'no-latex'])

def parse_args():
    """
    Parses command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments containing the file path to the .fil file,
                            start time, and end time.
    """
    parser = argparse.ArgumentParser(description='Plot out 357 dynamic spectra from .fil')
    parser.add_argument('fil', type=str, help='Path to .fil file')
    parser.add_argument('start_time', type=str, help='Start time in HH:MM format')
    parser.add_argument('end_time', type=str, help='End time in HH:MM format')
    args = parser.parse_args()
    return args

def get_time(fil):
    """
    Extracts the observation time from the .fil file.

    Args:
        fil (str): Path to the .fil file.

    Returns:
        np.ndarray: Array of time values in Modified Julian Date (MJD) format.
    """
    sppFil = FilReader(fil)
    nsamples = sppFil.header.nsamples
    tsamp = sppFil.header.tsamp
    startmjd = sppFil.header.tstart
    time = startmjd + np.linspace(0, nsamples * tsamp, nsamples) / (24 * 60 * 60)
    return time

def mjd_to_date(mjd):
    """
    Converts MJD time to a datetime object.

    Args:
        mjd (float): Modified Julian Date.

    Returns:
        datetime.datetime: Converted datetime object.
    """
    t = Time(mjd, format='mjd')
    return t.to_datetime()

def get_data(fil, start_idx, end_idx):
    """
    Retrieves the data from the filterbank file for the specified range of samples.

    Args:
        fil (str): Path to the .fil file.
        start_idx (int): Starting index for data extraction.
        end_idx (int): Ending index for data extraction.

    Returns:
        np.ndarray: Extracted data from the .fil file.
    """
    your_object = your.Your(fil)
    data = your_object.get_data(nstart=start_idx, nsamp=(end_idx - start_idx + 1))
    return data

def assign_frequencies(data, mode):
    """
    Assigns the frequency range based on the RCU mode.

    Args:
        data (np.ndarray): Data array with shape (time, frequency).
        mode (int or str): RCU mode identifier.

    Returns:
        np.ndarray: Frequency values corresponding to the given mode.
    
    Raises:
        ValueError: If an unsupported RCU mode is provided.
    """
    if mode == 3:
        freq_start, freq_end = 30, 90  # Mode 3: 30-90 MHz
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
    Resamples the 2D data array to match the target number of channels.

    Args:
        data (np.ndarray): Original data array with shape (time, original_nchans).
        target_nchans (int): Target number of frequency channels after resampling.
        original_nchans (int): Original number of frequency channels.

    Returns:
        np.ndarray: Resampled data array with shape (time, target_nchans).
    """
    original_freqs = np.linspace(0, original_nchans, original_nchans)
    target_freqs = np.linspace(0, original_nchans, target_nchans)

    # Interpolate each row (time sample) separately
    resampled_data = np.array([np.interp(target_freqs, original_freqs, row) for row in data])

    return resampled_data

def bkg_subtraction(data): 
    """
    Performs background subtraction on the data.

    Args:
        data (np.ndarray): Data array with shape (time, frequency).

    Returns:
        np.ndarray: Background-subtracted data.
    """
    bkg_data = data[:10]
    return data / np.mean(bkg_data, axis=0)

def rfi_masking(data): 
    """ 
    Masks RFI (Radio Frequency Interference) using the IQRM method.

    Args:
        data (np.ndarray): Data array with shape (time, frequency).

    Returns:
        np.ndarray: RFI-masked data with RFI affected regions replaced by NaNs.
    """
    spectral_std = np.std(data, axis=0)
    mask, votes = iqrm_mask(spectral_std, radius=2)
    return np.where(mask, np.nan, data)

def main():
    """
    Main function to process the .fil file and generate dynamic spectra plots.

    Steps:
        - Parse command-line arguments.
        - Extract observation time and data.
        - Perform background subtraction and RFI masking.
        - Assign frequencies based on the RCU mode.
        - Plot the dynamic spectra with appropriate frequency bands.
    """
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

    # Split data based on the modes
    mode7_data = data[:, :88]
    mode5_data = data[:, 88:288]
    mode3_data = data[:, 288:]

    # Background subtraction
    mode3_data = bkg_subtraction(mode3_data)
    mode5_data = bkg_subtraction(mode5_data)
    mode7_data = bkg_subtraction(mode7_data)

    # RFI masking
    mode3_data = rfi_masking(mode3_data)
    mode5_data = rfi_masking(mode5_data)
    mode7_data = rfi_masking(mode7_data)

    # mode7_data = resample_data(mode7_data, 200, 88)  # resample 

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

    time_utc_num = mdates.date2num(time_utc)

    upper_per = mode3_data.mean() + 1 * mode3_data.std()
    lower_per = mode3_data.mean() - 1 * mode3_data.std()

    # Plotting 
    fig = plt.figure(figsize=(16, 4 * 5))  # Increase figure height for additional subplots
    gs = fig.add_gridspec(5, 1, hspace=0, height_ratios=[1, 0.1, 1, 0.1, 0.44])

    axes = [fig.add_subplot(gs[i]) for i in range(5)]

    # Plot Mode 3
    mode3_vim, mode3_vmax = np.nanpercentile(mode3_data, [5, upper_per])
    im = axes[0].imshow(mode3_data.T, aspect='auto', cmap='Oranges', vmin=mode3_vim, vmax=mode3_vmax, extent=[time_utc_num[start_idx], time_utc_num[end_idx], mode3_freqs[0], mode3_freqs[-1]])
    axes[0].get_xaxis().set_ticks([])
    axes[0].set_yticks([30, 50, 70, 90])
    axes[0].set_title('Solar Dynamic Spectra: %s to %s' % (time_utc[start_idx].strftime('%Y-%m-%d %H:%M:%S'), time_utc[end_idx].strftime('%Y-%m-%d %H:%M:%S')))

    # Plot Gap1 (NaN values between mode 3 and mode 5)
    gap1_vim, gap1_vmax = np.nanpercentile(gap_data, [lower_per, upper_per])
    axes[1].imshow(gap_data.T, aspect='auto', cmap='Greys', vmin=gap1_vim, vmax=gap1_vmax, extent=[time_utc_num[start_idx], time_utc_num[end_idx], gap1_freqs[0], gap1_freqs[-1]])
    axes[1].get_xaxis().set_ticks([])
    axes[1].set_yticks([90, 100, 110])
    axes[1].text(0.5, 0.5, 'FM Band', horizontalalignment='center', verticalalignment='center', transform=axes[1].transAxes, fontsize=12)

    # Plot Mode 5
    # mode5_vim, mode5_vmax = np.nanpercentile(mode5_data, [lower_per, upper_per])
    mode5_vim, mode5_vmax = lower_per, upper_per
    im = axes[2].imshow(mode5_data.T, aspect='auto', cmap='Oranges', vmin=mode5_vim, vmax=mode5_vmax, extent=[time_utc_num[start_idx], time_utc_num[end_idx], mode5_freqs[0], mode5_freqs[-1]])
    axes[2].get_xaxis().set_ticks([])
    axes[2].set_yticks([120, 140, 160, 180])
    axes[2].set_ylabel('Frequency (MHz)')

    # Plot Gap2 (NaN values between mode 5 and mode 7)
    gap2_vim, gap2_vmax = np.nanpercentile(gap_data, [lower_per, upper_per])
    axes[3].imshow(gap_data.T, aspect='auto', cmap='Greys', vmin=gap2_vim, vmax=gap2_vmax, extent=[time_utc_num[start_idx], time_utc_num[end_idx], gap2_freqs[0], gap2_freqs[-1]])
    axes[3].get_xaxis().set_ticks([])
    axes[3].set_yticks([190, 200, 210])
    axes[3].text(0.5, 0.5, 'TV Broadcasting', horizontalalignment='center', verticalalignment='center', transform=axes[3].transAxes, fontsize=12)

    # Plot Mode 7
    mode7_vim, mode7_vmax = np.nanpercentile(mode7_data, [lower_per, upper_per])
    im = axes[4].imshow(mode7_data.T, aspect='auto', cmap='Oranges', vmin=mode7_vim, vmax=mode7_vmax, extent=[time_utc_num[start_idx], time_utc_num[end_idx], mode7_freqs[0], mode7_freqs[-1]])
    axes[4].set_yticks([210, 230, 250, 270])
    axes[4].set_xlabel('Time (UTC)')

    # Format the x-axis as date
    axes[4].xaxis_date()
    axes[4].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    # flip the y-axis
    axes[0].invert_yaxis()
    axes[1].invert_yaxis()
    axes[2].invert_yaxis()
    axes[3].invert_yaxis()
    axes[4].invert_yaxis()

    plt.savefig('dynamic_spectra/Sun357_dynamic_spectra_%s_%s_%s.png' % (time_utc[0].strftime('%Y-%m-%d'), args.start_time, args.end_time), dpi=300)

if __name__ == '__main__':
    main()

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
    return t.to_value('iso', 'date_hms')

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
        freq_start, freq_end = 10, 90  # Mode 3: 10-90 MHz
    elif mode == 5:
        freq_start, freq_end = 110, 190  # Mode 5: 110-190 MHz
    elif mode == 7:
        freq_start, freq_end = 210, 270  # Mode 7: 210-270 MHz
    else:
        raise ValueError("Unsupported RCU mode")

    nchans = data.shape[1]
    freqs = np.linspace(freq_start, freq_end, nchans)
    return freqs

def main():
    args = parse_args()
    fil_path = args.fil
    time = get_time(fil_path) 
    date = mjd_to_date(time[0])[0:10]
    header = your.Your(fil_path).your_header 
   

    start_time = date + ' ' + args.start_time
    end_time = date + ' ' + args.end_time 

    # Convert to MJD
    start_time_mjd = Time(start_time).mjd
    end_time_mjd = Time(end_time).mjd

    # Find the index of the start and end time
    start_idx = np.argmin(np.abs(time - start_time_mjd))
    end_idx = np.argmin(np.abs(time - end_time_mjd))

    print('======  Time ======')
    print('Obs Start Time :', mjd_to_date(time[0]))
    print('Obs End Time   :', mjd_to_date(time[-1]))
    print('Plot Start Time:', mjd_to_date(time[start_idx]))
    print('Plot End Time  :', mjd_to_date(time[end_idx]))

    data = get_data(fil_path, start_idx, end_idx)
    print('Data Shape:', data.shape)

    # Split data based on the modes
    mode3_data = data[:, :200]
    mode5_data = data[:, 200:400]
    mode7_data = data[:, 400:]

    mode3_freqs = assign_frequencies(mode3_data, 3)
    mode5_freqs = assign_frequencies(mode5_data, 5)
    mode7_freqs = assign_frequencies(mode7_data, 7)

    # NaN arrays with the gaps between 90 and 110 MHz and 190 and 210 MHz
    gap1 = np.full((data.shape[0], 20), np.nan)
    gap2 = np.full((data.shape[0], 20), np.nan)

    # Concatenate the data
    final_data = np.concatenate((mode3_data, gap1, mode5_data, gap2, mode7_data), axis=1)
    final_freqs = np.concatenate((mode3_freqs, np.full(200, np.nan), mode5_freqs, np.full(200, np.nan), mode7_freqs))

    # Plot the dynamic spectra
    plt.figure(figsize=(12, 6))
    plt.imshow(final_data.T, aspect='auto', extent=[time[start_idx], time[end_idx], final_freqs[0], final_freqs[-1]], cmap='viridis')
    plt.xlabel('Time')
    plt.ylabel('Frequency (MHz)')
    plt.savefig('test.png')


if __name__ == '__main__':
    main()

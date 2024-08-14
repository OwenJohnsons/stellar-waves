'''
Code Purpose: Plot out 357 dynamic spectra from .fil 
Author: Owen A. Johnson 
'''
import argparse
import numpy as np
import your
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science','ieee', 'no-latex'])
from sigpyproc.readers import FilReader
from astropy.time import Time
import os

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

def main():
    args = parse_args()
    fil_path = args.fil
    time = get_time(fil_path) 
    date = mjd_to_date(time[0])[0:10]

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

    # --- Plotting Section --- #
    header = your.Your(fil_path).your_header
    top_freq = 270
    bottom_freq = 10
    time_res = (time[end_idx] - time[start_idx]) * 24 * 3600 / (end_idx - start_idx)

    print('====== Frequency Info ======')
    print('Top Frequency   :', top_freq)
    print('Bottom Frequency:', bottom_freq)
    print('Time Resolution (s) :', time_res)

    subdir = '%s-%s-%sbit' % (header.source_name, header.tstart_utc, header.nbits)

    if not os.path.exists('dynamic_spectra'):
        os.makedirs('dynamic_spectra')

    if not os.path.exists('dynamic_spectra/%s' % (subdir)):
        os.makedirs('dynamic_spectra/%s' % (subdir))

    # Prepare data for plotting (entire frequency range)
    freq_range_start = bottom_freq
    freq_range_end = top_freq
    start_idx_freq = 0
    end_idx_freq = header.nchans

    # Slicing data for plotting
    sliced_data = data[:, start_idx_freq:end_idx_freq]

    y = np.linspace(freq_range_start, freq_range_end, sliced_data.shape[1])  # Frequency range
    x = (np.arange(sliced_data.shape[0]) * time_res / (60**2))  # Time in hours

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)

    x_extent = [0, x.max()]
    y_extent = [freq_range_start, freq_range_end]

    vmax, vmin = np.nanpercentile(sliced_data, (90, 30))
    ax.imshow(sliced_data.T, extent=x_extent + y_extent, aspect='auto', vmax=vmax, vmin=vmin, interpolation='nearest', cmap='viridis')
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Frequency (MHz)')
    ax.set_title('Dynamic Spectra for %s starting at %s : Time Resolution %s (ms)' % (header.source_name, start_time, np.round(time_res*1000, 1)))
    ax.invert_yaxis()

    plt.savefig('dynamic_spectra/%s/%s-%s-%sbit-FullBand.png' % (subdir, header.source_name, header.tstart_utc, header.nbits))
    plt.close()

if __name__ == '__main__':
    main()

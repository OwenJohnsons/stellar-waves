import os
import your
from your.candidate import Candidate
import glob as glob 
import matplotlib.pyplot as plt 
from matplotlib import gridspec
from your.utils.plotter import save_bandpass
import numpy as np
import scienceplots; import scienceplots; plt.style.use(['science','ieee', 'no-latex'])
from tqdm import tqdm 
from sigpyproc.readers import FilReader
import argparse

def factor(number, target):
    factors = []

    for i in range(1, int(number**0.5) + 1):
        if number % i == 0:
            factors.append(i)
            factors.append(number // i)

    closest_factor = min(factors, key=lambda x: abs(x - target))    
    return closest_factor

def rfi_masking(your_data, your_object, fil, rfi_plot_cond=True):
    rfi_mask = your.utils.rfi.savgol_filter(your_data.mean(0), your_object.foff, frequency_window=90, sigma=3)
    removed_data = your_data[:, ~rfi_mask]
    masked_data = np.where(rfi_mask, np.round(np.median(removed_data)), your_data)

    print('=== RFI Statistics ===')
    print('Number of masked channels for %s: %s' % (your_object.your_header.filename, np.sum(rfi_mask)))
    print('Percentage of masked channeles: %s' % float((np.sum(rfi_mask)/your_object.your_header.nchans)*100))
    print('Data: %s ± %s' % (np.mean(your_data), np.std(your_data)))
    print('Data with removed channels: %s ± %s' % (np.mean(removed_data), np.std(removed_data)))
    print('Data with masked channels: %s ± %s' % (np.mean(masked_data), np.std(masked_data)))

    if rfi_plot_cond == True: 
        plot_rfi(your_data, fil, mask=rfi_mask)

    return masked_data

def downsampling(array, nsamp, channels, original_sampling_rate, target_sampling_rate): 

    tobs = nsamp * original_sampling_rate

    print('=== Downsampling Information ===')
    print('Observation Duration: %s (hr)' % (tobs/60**2))
    print('Original Sampling Rate: %s (s)' % original_sampling_rate)
    print('Target Sampling Rate: %s (s)' % target_sampling_rate)

    downsampled_nsamp = int(nsamp/(target_sampling_rate/original_sampling_rate))
    downsampled_nsamp = factor(nsamp, downsampled_nsamp)
    print('Downsampling Factor: %s' % downsampled_nsamp)
    print('New time resolution: %s (s)' % (tobs/downsampled_nsamp))
    print('Ingested Array Shape: %s, %s' % (array.shape[0], array.shape[1]))

    if target_sampling_rate < original_sampling_rate:
        raise ValueError('Target sampling rate must be greater than original sampling rate.')
    
    reshaped_array = array.reshape(downsampled_nsamp, -1, channels)
    # downsampled_array = np.mean(reshaped_array, axis=1)
    downsampled_array = np.median(reshaped_array, axis=1)

    print('Downsampled Array Shape: %s, %s' % (downsampled_array.shape[0], downsampled_array.shape[1]))
    
    return downsampled_array, tobs/downsampled_nsamp

def data_scrunch(fil_fnames, time_res): 

    data_dictionary = {}; labels=['I', 'Q', 'U', 'V']; i = 0 

    if len(fil_fnames) != 4: 
        raise ValueError('Please provide 4 filterbanks for 4 polarisations')
    
    for fil in fil_fnames: 
        your_object = your.Your(fil)
        header = your_object.your_header
        tsamp = header.tsamp; top_freq = header.fch1

        sppFil = FilReader(fil)
        sppshape = (sppFil.header.nsamples, sppFil.header.nchans, sppFil.header.nifs)
        nsamples = sppshape[0]

        data = your_object.get_data(nstart=0, nsamp=nsamples)
        # masked_data = rfi_masking(data, your_object, fil, rfi_plot_cond=False)
        ds_data, new_tsamp = downsampling(data, nsamples, header.nchans, tsamp, time_res)
        data_dictionary[labels[i]] = ds_data; i += 1

    return data_dictionary, new_tsamp

def data2plot(data, start_index, end_index, freq_range_start, freq_range_end, time_res):

    # Slice the data based on frequency range
    # print('Indexes for slicing: %s, %s' % (start_index, end_index))
    sliced_data = data[:, start_index:end_index]

    y = np.linspace(freq_range_start, freq_range_end, sliced_data.shape[1])  # Frequency range
    x = (np.arange(sliced_data.shape[0])*time_res/(60**2)) 

    xo, yo = np.meshgrid(x, y)

    return xo, yo, sliced_data

def plot_rfi(data, fil, mask=None):
    print('Plotting RFI masks for %s' % fil)
    channels = np.arange(data.shape[1])
    c = Candidate(fp=str(fil))
    c.data = data
    c.dedisperse(target='GPU')
    data = c.dedispersed
    bandpass = data.mean(0)
    ts = data.mean(1)
    
    if mask is not None:
        data[:, mask] = 0

    fig = plt.figure(figsize=(8, 6))
    gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 3])

    ax0 = plt.subplot(gs[1, 0])
    ax0.imshow(data.T, aspect="auto", interpolation=None)
    ax0.set_xlabel("Time Samples")
    ax0.set_ylabel("Frequency Channels")

    ax1 = plt.subplot(gs[1, 1], sharey=ax0)
    plt.setp(ax1.get_yticklabels(), visible=False)
    ax1.plot(bandpass, channels, "k")
    ax1.set_xlabel("Flux (Arb. Units)")
    
    ax2 = plt.subplot(gs[0, 0], sharex=ax0)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.plot(ts, "k")
    ax2.set_ylabel("Flux (Arb. Units)")

    if mask is not None:
        for channel, val in zip(channels[mask], bandpass[mask]):
            ax0.axhline(channel, color="r", xmin=0, xmax=0.03, lw=0.1)
            ax1.scatter(val, channel, color="r", marker="o")

            plt.tight_layout()
            plt.savefig("rfi_masks/%s.png" % fil[0:-4])

# --- Main --- 
parser = argparse.ArgumentParser(description='Dynamic Spectra Plotter for 4 polarisations')
parser.add_argument('-b', '--bandwidth', type=int, help='Bandwidth in MHz', required=False, default=10)
parser.add_argument('-t', '--tres', type=int, help='Time resolution in seconds', required=False, default=4)
args = parser.parse_args()

# --- User Inputs --- #
bandwidth = args.bandwidth # MHz 
tres_user = args.tres # seconds 

filterbank_list = sorted(glob.glob('/mnt/ucc4_data2/data/Owen/raw-radio-stars/CR_Draconis_2023-12-12T14:09:51_S*'))
print('Filterbanks in directory....')
for file in filterbank_list:
    print(file)

if not os.path.exists('dynamic_spectra'):
    os.makedirs('dynamic_spectra')

if not os.path.exists('rfi_masks'):
    os.makedirs('rfi_masks')

dict, new_tsamp = data_scrunch(filterbank_list, time_res=tres_user)
print('Data Scrunching Complete....')

# --- Reading Header Information --- #
header = your.Your(filterbank_list[0]).your_header # Using one header as reference

tsamp = header.tsamp; top_freq = header.fch1
channel_bw = header.foff; channel_num = header.nchans 
nsamples = header.nspectra
bottom_freq = top_freq - abs(channel_bw*channel_num)

for i in tqdm(np.arange(bottom_freq, (top_freq - bandwidth), bandwidth)):
    # Define frequency range for slicing
    freq_range_start = i  # MHz
    freq_range_end = i + bandwidth   

    # Find indices corresponding to the frequency range
    start_idx = int((freq_range_start - bottom_freq) / abs(channel_bw))
    end_idx = int((freq_range_end - bottom_freq) / abs(channel_bw))

    # - Stokes Data Conversion for Plotting - 
    xI, yI, sliced_dataI = data2plot(dict['I'], start_index=start_idx, end_index=end_idx, freq_range_start=freq_range_start, freq_range_end=freq_range_end, time_res=new_tsamp)
    xQ, yQ, sliced_dataQ = data2plot(dict['Q'], start_index=start_idx, end_index=end_idx, freq_range_start=freq_range_start, freq_range_end=freq_range_end, time_res=new_tsamp)
    xU, yU, sliced_dataU = data2plot(dict['U'], start_index=start_idx, end_index=end_idx, freq_range_start=freq_range_start, freq_range_end=freq_range_end, time_res=new_tsamp)
    xV, yV, sliced_dataV = data2plot(dict['V'], start_index=start_idx, end_index=end_idx, freq_range_start=freq_range_start, freq_range_end=freq_range_end, time_res=new_tsamp)
  
    fig = plt.figure(figsize=(8, 4 * 2))
    gs = fig.add_gridspec(4, 1, hspace=0)

    axes = [fig.add_subplot(gs[i, 0]) for i in range(4)]
    vmax_multiple = 1

    x_extent = [0, xI.max()]  
    y_extent = [freq_range_start, freq_range_end] 

    # Plot Stokes I
    axes[0].imshow(sliced_dataI, extent=x_extent + y_extent, aspect='auto', vmax=np.mean(sliced_dataI.T) + vmax_multiple*np.std(sliced_dataI.T))
    axes[0].set_ylabel('Frequency (MHz)')
    axes[0].set_title('Dynamic Spectra for %s starting at %s : Time Resolution %s (s)' % (header.source_name, header.tstart_utc, np.round(new_tsamp,1)))
    axes[0].annotate('Stokes I', xy=(1, 1), xytext=(-5, -5), xycoords='axes fraction', textcoords='offset points', horizontalalignment='right', verticalalignment='top', fontsize=12, color='black')

    # Plot Stokes Q
    axes[1].imshow(sliced_dataQ, extent=x_extent + y_extent, aspect='auto', vmax=np.mean(sliced_dataQ.T) + vmax_multiple*np.std(sliced_dataQ.T))
    axes[1].set_ylabel('Frequency (MHz)')
    axes[1].annotate('Stokes Q', xy=(1, 1), xytext=(-5, -5), xycoords='axes fraction', textcoords='offset points', horizontalalignment='right', verticalalignment='top', fontsize=12, color='black')

    # Plot Stokes U
    axes[2].imshow(sliced_dataU, extent=x_extent + y_extent, aspect='auto', vmax=np.mean(sliced_dataU.T) + vmax_multiple*np.std(sliced_dataU.T))
    axes[2].set_ylabel('Frequency (MHz)')
    axes[2].annotate('Stokes U', xy=(1, 1), xytext=(-5, -5), xycoords='axes fraction', textcoords='offset points', horizontalalignment='right', verticalalignment='top', fontsize=12, color='black')

    # Plot Stokes V
    axes[3].imshow(sliced_dataV, extent=x_extent + y_extent, aspect='auto', vmax=np.mean(sliced_dataV.T) + vmax_multiple*np.std(sliced_dataV.T))
    axes[3].set_ylabel('Frequency (MHz)')
    axes[3].annotate('Stokes V', xy=(1, 1), xytext=(-5, -5), xycoords='axes fraction', textcoords='offset points', horizontalalignment='right', verticalalignment='top', fontsize=12, color='black')

    axes[-1].set_xlabel('Time (Hours)')

    plt.savefig('dynamic_spectra/%s-%s-%sbit-4pol-Freq[%s,%s].png' % (header.source_name, header.tstart_utc, header.nbits, i, i + bandwidth))
    plt.close()

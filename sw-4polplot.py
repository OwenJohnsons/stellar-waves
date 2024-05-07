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

def factor(number, target):
    factors = []

    for i in range(1, int(number**0.5) + 1):
        if number % i == 0:
            factors.append(i)
            factors.append(number // i)

    closest_factor = min(factors, key=lambda x: abs(x - target))    
    return closest_factor

def rfi_masking(your_data, your_object, fil, rfi_plot_cond=True):
    # composite_mask = your.utils.rfi.sk_sg_filter(your_data, your_object, your_object.your_header.nchans, sigma=3)
    rfi_mask = your.utils.rfi.savgol_filter(your_data.mean(0), your_object.foff, frequency_window=10, sigma=3)
    print('=== RFI Statistics ===')
    print('Number of masked channels for %s: %s' % (your_object.your_header.filename, np.sum(rfi_mask)))
    print('Percentage of masked channeles: %s' % float((np.sum(rfi_mask)/your_object.your_header.nchans)*100))
    removed_data = your_data[:, ~rfi_mask]
    print('Data: %s ± %s' % (np.mean(your_data), np.std(your_data)))
    print('Data with removed channels: %s ± %s' % (np.mean(removed_data), np.std(removed_data)))
    masked_data = np.where(rfi_mask, np.round(np.mean(removed_data)), your_data)
    print('Data with masked channels: %s ± %s' % (np.mean(masked_data), np.std(masked_data)))
    print('======================')

    if rfi_plot_cond == True: 
        plot_rfi(your_data, fil, mask=rfi_mask)

    return masked_data

def data_scrunch(fil_fnames, time_res): 

    data_dictionary = {}; labels=['I', 'Q', 'U', 'V']; i = 0 

    if len(fil_fnames) != 4: 
        raise ValueError('Please provide 4 filterbanks for 4 polarisations')
    
    for fil in fil_fnames: 
        your_object = your.Your(fil)
        header = your_object.your_header
        tsamp = header.tsamp; top_freq = header.fch1
        channel_bw = header.foff; channel_num = header.nchans 
        nsamples = header.nspectra
        bottom_freq = top_freq - abs(channel_bw*channel_num)

        new_tobs = 50
        tobs = nsamples*tsamp
        new_tobs_factor = int(new_tobs/time_res)
        scrunch = int(tobs/time_res)

        data = your_object.get_data(nstart=0, nsamp=new_tobs_factor*scrunch)

        masked_data = rfi_masking(data, your_object, fil, rfi_plot_cond=False)
        reshaped_data = data.reshape((scrunch, -1, masked_data.shape[1]))
        averaged_data = np.mean(reshaped_data, axis=1)

        data_dictionary[labels[i]] = averaged_data; i += 1


    return data_dictionary

def data2plot(data, start_index, end_index, freq_range_start, freq_range_end, time_res):

    # Slice the data based on frequency range
    # print('Indexes for slicing: %s, %s' % (start_index, end_index))
    sliced_data = data[:, start_index:end_index]

    y = np.linspace(freq_range_start, freq_range_end, sliced_data.shape[1])  # Frequency range
    x = (np.arange(sliced_data.shape[0]) * time_res) / (60 * 60)

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


# --- User Inputs --- #
bandwidth = 10 # MHz
tres_user = 8 # Time resolution in seconds

filterbank_list = sorted(glob.glob('/mnt/ucc4_data2/data/Owen/raw-radio-stars/CR_Draconis_2023-12-12T14:09:51_S*'))
print('Filterbanks in directory....')
for file in filterbank_list:
    print(file)

if not os.path.exists('dynamic_spectra'):
    os.makedirs('dynamic_spectra')

if not os.path.exists('rfi_masks'):
    os.makedirs('rfi_masks')

dict = data_scrunch(filterbank_list, time_res=10)
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
    xI, yI, sliced_dataI = data2plot(dict['I'], start_index=start_idx, end_index=end_idx, freq_range_start=freq_range_start, freq_range_end=freq_range_end, time_res=tres_user)
    xQ, yQ, sliced_dataQ = data2plot(dict['Q'], start_index=start_idx, end_index=end_idx, freq_range_start=freq_range_start, freq_range_end=freq_range_end, time_res=tres_user)
    xU, yU, sliced_dataU = data2plot(dict['U'], start_index=start_idx, end_index=end_idx, freq_range_start=freq_range_start, freq_range_end=freq_range_end, time_res=tres_user)
    xV, yV, sliced_dataV = data2plot(dict['V'], start_index=start_idx, end_index=end_idx, freq_range_start=freq_range_start, freq_range_end=freq_range_end, time_res=tres_user)
  
    fig = plt.figure(figsize=(8, 4 * 2))
    gs = fig.add_gridspec(4, 1, hspace=0)

    axes = [fig.add_subplot(gs[i, 0]) for i in range(4)]

    # Plot stokes I
    axes[0].contourf(xI, yI, sliced_dataI.T, cmap='magma', vmax=np.mean(sliced_dataI.T) + 2*np.std(sliced_dataI.T)) 
    axes[0].set_ylabel('Frequency (MHz)')
    axes[0].set_title('Dynamic Spectra for %s starting at %s' % (header.source_name, header.tstart_utc))
    axes[0].annotate('Stokes I',xy=(1, 1),  xytext=(-5, -5),  xycoords='axes fraction',  textcoords='offset points',horizontalalignment='right',  verticalalignment='top',  fontsize=12,  color='black')
    
    # Plot stokes Q
    axes[1].contourf(xQ, yQ, sliced_dataQ.T, cmap='magma', vmax=np.mean(sliced_dataQ.T) + 2*np.std(sliced_dataQ.T))
    axes[1].set_ylabel('Frequency (MHz)')
    axes[1].annotate('Stokes Q',xy=(1, 1),  xytext=(-5, -5),  xycoords='axes fraction',  textcoords='offset points',horizontalalignment='right',  verticalalignment='top',  fontsize=12,  color='black')

    # Plot stokes U
    axes[2].contourf(xU, yU, sliced_dataU.T, cmap='magma', vmax=np.mean(sliced_dataU.T) + 2*np.std(sliced_dataU.T))
    axes[2].set_ylabel('Frequency (MHz)')
    axes[2].annotate('Stokes U',xy=(1, 1),  xytext=(-5, -5),  xycoords='axes fraction',  textcoords='offset points',horizontalalignment='right',  verticalalignment='top',  fontsize=12,  color='black')

    # Plot stokes V
    axes[3].contourf(xV, yV, sliced_dataV.T, cmap='magma', vmax=np.mean(sliced_dataV.T) + 2*np.std(sliced_dataV.T))
    axes[3].set_ylabel('Frequency (MHz)')
    axes[3].annotate('Stokes V',xy=(1, 1),  xytext=(-5, -5),  xycoords='axes fraction',  textcoords='offset points',horizontalalignment='right',  verticalalignment='top',  fontsize=12,  color='black')
    axes[-1].set_xlabel('Time (Hours)')

    plt.savefig('dynamic_spectra/%s-%s-%sbit-4pol-Freq[%s,%s].png' % (header.source_name, header.tstart_utc, header.nbits, i, i + bandwidth))
    plt.close()

    break 

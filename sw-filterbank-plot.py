"""
Code Purpose: Plot dynamic spectra from filterbank or HDF5 files.
Author: Owen A. Johnson 
Date: 2025-05-01
Usage: python sw-filterbank-plot.py <filterbank_file> -t <time_resolution> -b <bandwidth>
Description: Here you just need to provide a path to the fil or h5, a time resolution in seconds, and a bandwidth slice width in MHz, this will downsample and then plot for the give slice width. 
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scienceplots; plt.style.use(['science', 'ieee', 'no-latex'])
from your import Your
from sigpyproc.readers import FilReader
from tqdm import tqdm
from blimpy import Waterfall


# ---- Functions ----
def factor(number, target):
    factors = [i for i in range(1, int(number**0.5) + 1) if number % i == 0]
    factors += [number // i for i in factors]
    return min(factors, key=lambda x: abs(x - target))

def downsample(data, nsamp, channels, orig_tsamp, target_tsamp):
    if target_tsamp < orig_tsamp:
        raise ValueError("Target tsamp must be >= original tsamp")

    downsample_factor = factor(nsamp, int(nsamp * orig_tsamp / target_tsamp))
    reshaped = data.reshape(downsample_factor, -1, channels)
    downsampled = np.median(reshaped, axis=1)
    new_tsamp = (nsamp * orig_tsamp) / downsample_factor
    
    print(f"\n=== Downsampling Information ===")
    print(f"Original Sampling Rate: {orig_tsamp} (s)")
    print(f"Target Sampling Rate: {target_tsamp} (s)")
    print(f"Downsampling Factor: {downsample_factor}")
    print(f"New Time Resolution: {new_tsamp} (s)")
    print(f"Ingested Array Shape: {data.shape[0]}, {data.shape[1]}")
    print(f"Downsampled Array Shape: {downsampled.shape[0]}, {downsampled.shape[1]}")
    print(f"Observation Duration: {nsamp * orig_tsamp / 3600} (hr)\n")
    
    return downsampled, new_tsamp

def data2plot(data, start_index, end_index, freq_start, freq_end, tsamp):
    sliced = data[:, start_index:end_index]
    y = np.linspace(freq_start, freq_end, sliced.shape[1])
    x = np.arange(sliced.shape[0]) * tsamp / 3600  # hours
    return np.meshgrid(x, y), sliced

# --- Class for Loading --- 
class DynamicSpectrumLoader:
    def __init__(self, filename):
        self.filename = filename
        self.ext = os.path.splitext(filename)[-1].lower()

        if self.ext == ".fil":
            self.loader = self._load_fil()
        elif self.ext == ".h5":
            self.loader = self._load_h5()
        else:
            raise ValueError("Unsupported file type. Only .fil and .h5 are supported.")

    def _load_fil(self):
        print(f"Reading .fil file: {self.filename}")
        self.your_obj = Your(self.filename)
        self.header = self.your_obj.your_header
        self.source = self.header.source_name
        self.tsamp = self.header.tsamp
        self.nchans = self.header.nchans
        self.foff = self.header.foff
        self.fch1 = self.header.fch1
        self.tstart_utc = self.header.tstart_utc

        fil = FilReader(self.filename)
        self.nsamp = fil.header.nsamples
        self.data = self.your_obj.get_data(nstart=0, nsamp=self.nsamp)

        return "fil"

    def _load_h5(self):
        print(f"Reading .h5 file: {self.filename}")
        self.wf = Waterfall(self.filename)
        self.header = self.wf.header
        self.source = self.header.get("source_name", "unknown_source")
        self.tsamp = self.header["tsamp"]
        self.nchans = self.header["nchans"]
        self.foff = self.header["foff"]
        self.fch1 = self.header["fch1"]
        self.tstart_utc = str(self.header.get("tstart", "unknown_start"))
        self.nsamp = self.wf.data.shape[0]
        self.data = np.squeeze(self.wf.data)

        return "h5"

    def downsample_data(self, target_tsamp):
        self.data, self.new_tsamp = downsample(
            self.data, self.nsamp, self.nchans, self.tsamp, target_tsamp
        )

    def get_freq_bounds(self):
        top_freq = self.fch1
        bottom_freq = top_freq - abs(self.foff) * self.nchans
        return bottom_freq, top_freq

# ---- Argument Parser ----
parser = argparse.ArgumentParser(description="Dynamic Spectrum Plotter for a Single Filterbank or HDF5")
parser.add_argument("filterbank", type=str, help="Path to the filterbank or .h5 file")
parser.add_argument("-t", "--tres", type=float, default=4, help="Desired time resolution in seconds")
parser.add_argument("-b", "--bandwidth", type=float, default=10, help="Bandwidth slice in MHz")
args = parser.parse_args()

# ---- Load data ----
ds = DynamicSpectrumLoader(args.filterbank)
ds.downsample_data(args.tres)
bottom_freq, top_freq = ds.get_freq_bounds()

# ---- Output Directory ----
outdir = f"dynamic_spectra/{ds.source}-{ds.tstart_utc}"
os.makedirs(outdir, exist_ok=True)

# ---- Plotting ----
for f0 in tqdm(np.arange(bottom_freq, top_freq - args.bandwidth, args.bandwidth)):
    f1 = f0 + args.bandwidth
    idx0 = int((f0 - bottom_freq) / abs(ds.foff))
    idx1 = int((f1 - bottom_freq) / abs(ds.foff))

    (X, Y), Z = data2plot(ds.data, idx0, idx1, f0, f1, ds.new_tsamp)

    # --- Spectrum ---
    spectrum = np.mean(Z, axis=0)              # shape (n_freq,)
    freqs = np.linspace(f0, f1, Z.shape[1])    # shape (n_freq,)

    # Setup figure with two panels, sharing y-axis
    fig = plt.figure(figsize=(12, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1], wspace=0.05)

    # --- Dynamic spectrum panel ---
    ax0 = fig.add_subplot(gs[0])
    vmax, vmin = np.nanpercentile(Z, (90, 30))
    ax0.pcolormesh(X, Y, Z.T, shading="auto", cmap="viridis", vmax=vmax, vmin=vmin)
    ax0.set_ylabel("Frequency (MHz)")
    ax0.set_xlabel("Time (hr)")
    ax0.set_title(f"{ds.source} Dynamic Spectrum\n{f0:.2f}–{f1:.2f} MHz, {ds.new_tsamp:.2f}s resolution")

    # --- Spectrum panel (right) ---
    ax1 = fig.add_subplot(gs[1], sharey=ax0)
    ax1.plot(spectrum, freqs, color='black')
    ax1.set_xlabel("Arb. Power")
    ax1.tick_params(labelleft=False, labelright=True)
    ax1.grid(True, alpha=0.3)

    # --- Save ---
    outname = f"{outdir}/{ds.source}_{ds.tstart_utc}_Freq[{f0:.2f},{f1:.2f}].png"
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()


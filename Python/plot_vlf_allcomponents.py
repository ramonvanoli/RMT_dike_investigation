import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import detrend
import xml.etree.ElementTree as ET
# This scripts plots the amplitude spectra of RMT data for all components.

plt.style.use('ggplot')

def read_ats_direct(ats_file, xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    sample_rate = float(root.find('.//sample_freq').text)
    start_date = root.find('.//start_date').text
    start_time = root.find('.//start_time').text
    datetime_str = f"{start_date}T{start_time}.0"

    with open(ats_file, 'rb') as f:
        f.read(1024)
        raw = np.frombuffer(f.read(), dtype=np.int32).astype(np.float64)

    chan = {
        'sample_rate': sample_rate,
        'datetime': datetime_str,
    }
    return chan, raw

# ── base path & component definitions ────────────────────────────────────────
base_dir = "C:/Users/Baar/Masterthesis_RMT/Proc_Potsdam/Data_Dikeproject/Dike_1_Beerenpolderse/Line_3/ADU_GFZ/meas_2026-03-10_15-25-52/"

components = {
    'Ex': ('C00', 'TEx'),
    'Ey': ('C01', 'TEy'),
    'Hx': ('C02', 'THx'),
    'Hy': ('C03', 'THy'),
    'Hz': ('C04', 'THz'),
}

xml_file = base_dir + "036_2026-03-10_15-25-52_2026-03-10_15-26-12_R000_524288H.xml"

window_size = 524288

# ── process all components ────────────────────────────────────────────────────
results = {}

for comp_name, (chan_code, type_code) in components.items():
    ats_file = base_dir + f"036_V01_{chan_code}_R000_{type_code}_BH_524288H.ats"

    try:
        channel, raw_data = read_ats_direct(ats_file, xml_file)
    except FileNotFoundError:
        print(f"File not found for {comp_name}, skipping: {ats_file}")
        continue

    count = 0
    stacked_fft = None

    for start in range(0, len(raw_data), window_size):
        data = raw_data[start:start + window_size]
        if len(data) < window_size:
            break
        data_detrended = detrend(data)
        fft_data = np.fft.rfft(data_detrended * np.hanning(len(data_detrended)))
        fft_magnitude = np.abs(fft_data)

        if stacked_fft is None:
            stacked_fft = fft_magnitude
        else:
            stacked_fft += fft_magnitude
        count += 1

    stacked_fft /= count
    results[comp_name] = (channel['sample_rate'], stacked_fft)
    print(f"{comp_name}: {count} windows processed")

# ── plot all components in subplots ──────────────────────────────────────────
fmin = 1e3
fmax = 200e3

fig, axes = plt.subplots(len(results), 1, figsize=(20, 2 * len(results)), sharex=True)
if len(results) == 1:
    axes = [axes]

# Poster-friendly typography (readable but not oversized)
LABEL_FONTSIZE = 13
TICK_FONTSIZE = 11
LEGEND_FONTSIZE = 11
SAVE_DPI = 600

for ax, (comp_name, (sample_rate, stacked_fft)) in zip(axes, results.items()):
    frequencies = np.fft.rfftfreq(window_size, d=1 / sample_rate)
    ax.loglog(frequencies[1:] / 1000, stacked_fft[1:], label=comp_name, color='steelblue', linewidth=1.0)
    ax.set_xlim(fmin / 1000, fmax / 1000)
    ax.legend(fontsize=LEGEND_FONTSIZE, loc='upper left')
    ax.grid(True, which="both", ls="--")
    ax.tick_params(axis='both', labelsize=TICK_FONTSIZE)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}'))
    ax.xaxis.set_minor_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}'))

axes[-1].set_xlabel('Frequency (kHz)', fontsize=LABEL_FONTSIZE)
fig.supylabel('Amplitude', fontsize=LABEL_FONTSIZE)
plt.tight_layout(rect=[0.04, 0, 1, 1])

out_file = base_dir + 'FFT_all_components.png'
plt.savefig(
    out_file,
    dpi=SAVE_DPI,
    bbox_inches='tight',
    facecolor=fig.get_facecolor(),
    edgecolor='none',
)
print(f"Saved: {out_file}")
plt.show()
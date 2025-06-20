import os
import pickle
import re
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from collections import defaultdict
from astropy.table import Table, vstack, Column
from PyPDF2 import PdfMerger

def merge_pdfs_by_pattern(folder_path, pattern, output_filename):
    search_pattern = os.path.join(folder_path, pattern)
    pdf_files = sorted(glob.glob(search_pattern))

    if not pdf_files:
        print(f"No PDF files found for pattern: {pattern}")
        return

    merger = PdfMerger()
    for pdf in pdf_files:
        print(f"Merging: {pdf}")
        merger.append(pdf)

    output_path = os.path.join(folder_path, output_filename)
    merger.write(output_path)
    merger.close()
    print(f"Merged PDF saved to: {output_path}")

def extract_sightline():
    return os.path.basename(os.getcwd())

def extract_redshift(filename):
    match = re.search(r'z[_]?(\d+\.\d+)', filename)
    return float(match.group(1)) if match else None

def extract_instrument(filename):
    return filename.split('_')[0]

def plot_ion_panels(merged, zval, sightline):
    output_dir = f'results_z_{zval}'
    os.makedirs(output_dir, exist_ok=True)

    pdf_nav = PdfPages(os.path.join(output_dir, f'{sightline}_z_{zval}_Nav.pdf'))

    for ion, data in merged.items():
        if not isinstance(data, dict) or 'vel' not in data:
            continue

        vel = data['vel']
        EWlims = data['EWlims']
        v1, v2 = EWlims
        fnorm = data.get('flux', np.ones_like(vel)) / data.get('contin', np.ones_like(vel))
        fnorm_err = data.get('fnorm_err', np.zeros_like(vel))

        fig, axs = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
        axs[0].set_title(ion)
        axs[0].plot(vel, fnorm, color='blue', label='Normalized flux')
        axs[0].plot(vel, fnorm_err, color='black', label='Error')
        axs[0].axvline(v1, linestyle='--', color='blue')
        axs[0].axvline(v2, linestyle='--', color='blue')
        axs[0].axhline(1.0, linestyle=':', color='k')
        axs[0].axhline(0.0, color='black')
        axs[0].legend()
        axs[0].set_ylabel('Normalized flux')
        axs[0].set_xlim(-1000, 1000)

        if 'Nav_vel' in data and 'Nav' in data:
            axs[1].plot(data['Nav_vel'], data['Nav'], drawstyle='steps-mid', color='green', label='Na(v)')
            axs[1].axhline(0, linestyle=':', color='k')
            axs[1].axvline(v1, linestyle='--', color='blue')
            axs[1].axvline(v2, linestyle='--', color='blue')
            axs[1].set_ylabel(r'$N_a(v)$ [cm$^{-2}$ (km/s)$^{-1}$]')
            axs[1].legend()

        axs[1].set_xlabel('Velocity (km/s)')
        plt.tight_layout()
        pdf_nav.savefig(fig)
        plt.close(fig)

    pdf_nav.close()

# ---- Begin main merging logic ---- #
all_files = os.listdir('.')
pickle_files = [f for f in all_files if f.endswith('.p') and extract_redshift(f) is not None]
table_files = [f for f in all_files if f.endswith('.dat') and extract_redshift(f) is not None]

sightline = extract_sightline()
z_to_pickles = defaultdict(list)
z_to_tables = defaultdict(list)

for pf in pickle_files:
    z = extract_redshift(pf)
    z_to_pickles[z].append(pf)

for tf in table_files:
    z = extract_redshift(tf)
    z_to_tables[z].append(tf)

for z in sorted(z_to_pickles.keys()):
    merged = {}
    table_rows = []

    for pf in z_to_pickles[z]:
        with open(pf, 'rb') as f:
            data = pickle.load(f)
        instrument = extract_instrument(pf)

        for key in data:
            if key == 'Target':
                merged['Target'] = data['Target']
                merged['Sightline'] = sightline
                continue

            data[key]['instrument'] = instrument

            if key in merged and isinstance(merged[key], dict):
                merged[key].update(data[key])
            else:
                merged[key] = data[key]

    outdir = f'results_z_{z}'
    os.makedirs(outdir, exist_ok=True)

    out_pickle = os.path.join(outdir, f'{sightline}_z_{z}_Merged.pkl')
    with open(out_pickle, 'wb') as f:
        pickle.dump(merged, f)

    if z in z_to_tables:
        dat_tables = []
        for tf in z_to_tables[z]:
            try:
                t = Table.read(tf, format='ascii')
                instrument = extract_instrument(tf)
                t['instrument'] = Column([instrument] * len(t), dtype='str')
                dat_tables.append(t)
            except Exception as e:
                print(f"Could not read {tf}: {e}")

        if dat_tables:
            combined_table = vstack(dat_tables, join_type='exact', metadata_conflicts='silent')
            out_table = os.path.join(outdir, f'{sightline}_z_{z}_Merged_Table.dat')
            combined_table.write(out_table, format='ascii.commented_header', overwrite=True)

            out_csv = os.path.join(outdir, f'{sightline}_z_{z}_Merged_Table.csv')
            combined_table.write(out_csv, format='csv', overwrite=True)

            print(f"Saved merged table to: {out_table}")
            print(f"Saved CSV table to: {out_csv}")

    # Generate Na(v) and flux plots
    plot_ion_panels(merged, z, sightline)

    # Merge all *_z_<z>_Ions*.pdf for this redshift
    pattern = f"*z_{z}_Ions*.pdf"
    output_pdf = f"combined_{z}.pdf"
    merge_pdfs_by_pattern(folder_path=".", pattern=pattern, output_filename=output_pdf)

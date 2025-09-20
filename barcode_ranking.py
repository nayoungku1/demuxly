#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, gzip, collections, h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse

# -----------------------------------------------
# 1) R1 FASTQ 읽기
# -----------------------------------------------
def read_all_r1_fastq(fastq_folder, max_reads_per_file=None):
    barcode_umi_counts = collections.Counter()
    fastq_files = sorted([f for f in os.listdir(fastq_folder)
                          if f.endswith(".fastq.gz") and "_R1_" in f])
    if not fastq_files:
        raise FileNotFoundError(f"No R1 fastq.gz files found in {fastq_folder}")
    print(f"Found {len(fastq_files)} R1 files. Reading...")

    for f in fastq_files:
        fp = os.path.join(fastq_folder, f)
        print(f"  processing {f} ...")
        with gzip.open(fp, "rt") as fh:
            cnt = 0
            while True:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().strip()
                fh.readline()  # +
                fh.readline()  # qual
                cb = seq[:16]
                umi = seq[16:28]
                barcode_umi_counts[(cb, umi)] += 1
                cnt += 1
                if max_reads_per_file and cnt >= max_reads_per_file:
                    break
    print("Done reading FASTQ.")
    return barcode_umi_counts

# -----------------------------------------------
# 2) Cell Ranger h5에서 barcodes 읽기
# -----------------------------------------------
def load_cellranger_barcodes(h5_path):
    with h5py.File(h5_path, "r") as f:
        arr = None
        if "matrix" in f and "barcodes" in f["matrix"]:
            arr = f["matrix"]["barcodes"][:]
        elif "barcodes" in f:
            arr = f["barcodes"][:]
        else:
            def find_barcodes(group):
                for k in group:
                    obj = group[k]
                    if isinstance(obj, h5py.Dataset) and k == "barcodes":
                        return obj[:]
                    elif isinstance(obj, h5py.Group):
                        res = find_barcodes(obj)
                        if res is not None:
                            return res
                return None
            arr = find_barcodes(f)
        if arr is None:
            raise RuntimeError("Could not find barcodes dataset inside the h5 file.")
        barcodes = [b.decode("utf-8") if isinstance(b, (bytes, bytearray)) else str(b) for b in arr]
    return set(barcodes)

def tolerant_barcode_set(barcodes):
    s = set()
    for b in barcodes:
        s.add(b)
        if b.endswith("-1"):
            s.add(b[:-2])
        else:
            s.add(b + "-1")
    return s

# -----------------------------------------------
# 3) Plot & Save function
# -----------------------------------------------
def plot_manual_vs_cellranger(barcode_umi_counts, output_dir,
                              show_top_n_on_x=20000,
                              manual_topk=1000,
                              cellranger_h5_path=None,
                              expect_cells=1000,
                              plot_filename="barcode_rank_plot.png",
                              csv_filename="barcode_summary.csv"):

    os.makedirs(output_dir, exist_ok=True)

    barcode_counts = {}
    for (cb, umi), c in barcode_umi_counts.items():
        barcode_counts[cb] = barcode_counts.get(cb, 0) + c

    sorted_items = sorted(barcode_counts.items(), key=lambda x: x[1], reverse=True)
    sorted_barcodes = [b for b,_ in sorted_items]
    sorted_counts = [cnt for _,cnt in sorted_items]
    Nplot = min(show_top_n_on_x, len(sorted_counts))
    x = np.arange(1, Nplot+1)
    y = np.array(sorted_counts[:Nplot])

    manual_barcodes = set(sorted_barcodes[:manual_topk])
    manual_indices = [i for i,b in enumerate(sorted_barcodes[:Nplot]) if b in manual_barcodes]
    manual_counts = [sorted_counts[i] for i in manual_indices]

    cr_indices, cr_counts, cr_barcodes = [], [], set()
    if cellranger_h5_path:
        try:
            cr_barcodes = load_cellranger_barcodes(cellranger_h5_path)
            cr_barcodes = tolerant_barcode_set(cr_barcodes)
            cr_indices = [i for i,b in enumerate(sorted_barcodes[:Nplot]) if b in cr_barcodes]
            cr_counts = [sorted_counts[i] for i in cr_indices]
        except Exception as e:
            print("Warning: couldn't load cellranger barcodes:", e)

    # plot
    plt.figure(figsize=(10,6))
    plt.plot(x, y, marker='.', linestyle='none', alpha=0.4, label='All barcodes (manual raw)')
    if manual_indices:
        plt.scatter([i+1 for i in manual_indices], manual_counts, s=20, facecolors='none',
                    edgecolors='green', label=f'Manual top-{manual_topk}')
    if cr_indices:
        plt.scatter([i+1 for i in cr_indices], cr_counts, s=15, color='red', label=f'Cell Ranger selected ({len(cr_barcodes)})')
    plt.axvline(x=expect_cells, color='gray', linestyle='--', linewidth=1, label=f'--expect-cells={expect_cells}')
    if cr_indices:
        plt.axvline(x=len(cr_barcodes), color='red', linestyle=':', linewidth=1, label=f'CR detected = {len(cr_barcodes)}')
    plt.yscale('log')
    plt.xlabel('Barcode rank')
    plt.ylabel('UMI counts (log scale)')
    plt.title('Barcode rank: Manual top-K vs Cell Ranger selection')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.4)
    plt.tight_layout()

    plot_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_path)
    plt.close()
    print(f"Plot saved to {plot_path}")

    # save CSV summary
    df = pd.DataFrame({
        "barcode": sorted_barcodes[:Nplot],
        "UMI_count": sorted_counts[:Nplot],
        "manual_topk": [1 if i in manual_indices else 0 for i in range(Nplot)],
        "cellranger_selected": [1 if i in cr_indices else 0 for i in range(Nplot)]
    })
    csv_path = os.path.join(output_dir, csv_filename)
    df.to_csv(csv_path, index=False)
    print(f"CSV summary saved to {csv_path}")

# -----------------------------------------------
# 4) main
# -----------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Barcode rank plot & comparison with Cell Ranger")
    parser.add_argument("--fastq_folder", required=True, help="Folder containing R1 FASTQ files")
    parser.add_argument("--cellranger_h5", required=False, default=None, help="Cell Ranger raw_feature_bc_matrix.h5")
    parser.add_argument("--output_dir", required=True, help="Output directory to save plots and CSV")
    parser.add_argument("--manual_topk", type=int, default=1000, help="Number of top barcodes to manually select")
    parser.add_argument("--expect_cells", type=int, default=1000, help="--expect-cells value used in Cell Ranger")
    args = parser.parse_args()

    print("Reading FASTQ...")
    barcode_umi_counts = read_all_r1_fastq(args.fastq_folder)

    print("Generating plot and CSV...")
    plot_manual_vs_cellranger(barcode_umi_counts,
                              output_dir=args.output_dir,
                              manual_topk=args.manual_topk,
                              cellranger_h5_path=args.cellranger_h5,
                              expect_cells=args.expect_cells)

if __name__ == "__main__":
    main()

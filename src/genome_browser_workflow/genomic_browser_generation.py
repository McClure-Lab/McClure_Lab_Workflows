#!/usr/bin/env python3
"""
Generate smoothed per-chromosome BrdU genome browser plots.
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import ScalarFormatter

plt.rcParams["agg.path.chunksize"] = 10000

SUBTEL_SIZE = 30000

GENBANK_TO_CHR = {
    "CM007964.1": "1",
    "CM007965.1": "2",
    "CM007966.1": "3",
    "CM007967.1": "4",
    "CM007968.1": "5",
    "CM007969.1": "6",
    "CM007970.1": "7",
    "CM007971.1": "8",
    "CM007972.1": "9",
    "CM007973.1": "10",
    "CM007974.1": "11",
    "CM007975.1": "12",
    "CM007976.1": "13",
    "CM007977.1": "14",
    "CM007978.1": "15",
    "CM007979.1": "16",
    "CM007980.1": "p2-micron",
    "CM007981.1": "MT",
}

NCBI_REFSEQ_TO_CHR = {
    "NC_001133.9": "1",
    "NC_001134.8": "2",
    "NC_001135.5": "3",
    "NC_001136.10": "4",
    "NC_001137.3": "5",
    "NC_001138.5": "6",
    "NC_001139.9": "7",
    "NC_001140.6": "8",
    "NC_001141.2": "9",
    "NC_001142.9": "10",
    "NC_001143.9": "11",
    "NC_001144.5": "12",
    "NC_001145.3": "13",
    "NC_001146.8": "14",
    "NC_001147.6": "15",
    "NC_001148.4": "16",
    "NC_001224.1": "MT",
}

CHROM_LENGTHS = {
    "1": 230218,
    "2": 813184,
    "3": 316620,
    "4": 1531933,
    "5": 576874,
    "6": 270161,
    "7": 1090940,
    "8": 562643,
    "9": 439888,
    "10": 745751,
    "11": 666816,
    "12": 1078177,
    "13": 924431,
    "14": 784333,
    "15": 1091291,
    "16": 948066,
    "MT": 85779,
    "p2-micron": 6318,
}


def repo_root():
    return Path(__file__).resolve().parents[2]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate smoothed genome browser plots from BrdU bedgraphs."
    )
    add_common_args(parser)
    parser.add_argument(
        "--window",
        type=int,
        default=1000,
        help="Centered rolling window used for smoothed signal tracks.",
    )
    return parser.parse_args()


def add_common_args(parser):
    root = repo_root()
    bed_dir = root / "data" / "bed"
    parser.add_argument("--positive-bedgraph", required=True)
    parser.add_argument("--negative-bedgraph", required=True)
    parser.add_argument(
        "--output-dir",
        default=str(root / "results" / "genome_browser_results" / "smoothed"),
    )
    parser.add_argument(
        "--prefix",
        default="genome_browser",
        help="Prefix included in generated PNG filenames.",
    )
    parser.add_argument(
        "--g4-bed",
        default=str(bed_dir / "W303_g4_motifs.bed"),
        help="G4 motif BED file.",
    )
    parser.add_argument(
        "--trna-bed",
        default=str(bed_dir / "trna_coordinates.bed"),
        help="tRNA BED file.",
    )
    parser.add_argument(
        "--te-bed",
        default=str(bed_dir / "w303_te_and_ltrs.bed"),
        help="TE and LTR BED file.",
    )


def map_chromosome(chrom_id: str) -> Optional[str]:
    chrom = str(chrom_id)

    if chrom in GENBANK_TO_CHR:
        return GENBANK_TO_CHR[chrom]

    if chrom in NCBI_REFSEQ_TO_CHR:
        return NCBI_REFSEQ_TO_CHR[chrom]

    if chrom in CHROM_LENGTHS:
        return chrom

    if chrom.startswith("chr"):
        value = chrom[3:]
        if value in CHROM_LENGTHS:
            return value

    return None


def chrom_sort_key(value):
    value = str(value)
    if value.isdigit():
        return 0, int(value)
    return 1, value


def load_bedgraph(path):
    columns = ["chrom", "start", "end", "frac_mod", "Nmod", "Nvalid_cov"]
    df = pd.read_csv(path, sep="\t", header=None, names=columns)

    for column in columns[1:]:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    df = df.dropna(subset=columns)
    df["chrom"] = df["chrom"].map(map_chromosome)
    df = df.dropna(subset=["chrom"])
    df["chrom"] = df["chrom"].astype(str)
    return df


def prepare_dataframe(positive_bedgraph, negative_bedgraph):
    pos_df = load_bedgraph(positive_bedgraph)
    neg_df = load_bedgraph(negative_bedgraph)

    df = pd.concat([pos_df, neg_df], ignore_index=True)

    if df.empty:
        raise ValueError("No usable BrdU rows were found in the bedgraph files.")

    return df


def load_feature_bed(path, description):
    columns = ["chrom", "start", "end", "name", "score", "strand"]

    if not path:
        print(f"[WARN] No {description} BED provided; track will be empty.", file=sys.stderr)
        return pd.DataFrame(columns=columns)

    if not os.path.exists(path):
        print(f"[WARN] {description} BED not found: {path}; track will be empty.", file=sys.stderr)
        return pd.DataFrame(columns=columns)

    df = pd.read_csv(path, sep="\t", header=None, comment="#")

    if df.shape[1] < 3:
        print(f"[WARN] {description} BED has fewer than 3 columns: {path}", file=sys.stderr)
        return pd.DataFrame(columns=columns)

    df = df.iloc[:, : min(6, df.shape[1])].copy()
    df.columns = columns[: df.shape[1]]

    for missing in columns[df.shape[1] :]:
        df[missing] = "."

    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["chrom"] = df["chrom"].map(map_chromosome)
    df = df.dropna(subset=["chrom", "start", "end"])
    df["chrom"] = df["chrom"].astype(str)
    return df


def combine_strands(chrom_df):
    combined = (
        chrom_df.groupby(["chrom", "start", "end"], as_index=False)
        .agg(Nmod=("Nmod", "sum"), Nvalid_cov=("Nvalid_cov", "sum"))
        .sort_values("start", kind="mergesort")
    )
    combined["frac_mod"] = np.where(
        combined["Nvalid_cov"] > 0,
        combined["Nmod"] / combined["Nvalid_cov"],
        0.0,
    )
    combined["BrdU_pct"] = combined["frac_mod"] * 100.0
    return combined


def smooth_signals(chrom_df, window):
    chrom_df["Coverage_plot"] = (
        chrom_df["Nvalid_cov"].rolling(window=window, min_periods=1, center=True).mean()
    )
    chrom_df["BrdU_plot"] = (
        chrom_df["Nmod"].rolling(window=window, min_periods=1, center=True).mean()
    )
    nmod_window = (
        chrom_df["Nmod"].rolling(window=window, min_periods=1, center=True).sum()
    )
    cov_window = (
        chrom_df["Nvalid_cov"].rolling(window=window, min_periods=1, center=True).sum()
    )
    chrom_df["BrdU_pct_plot"] = np.where(
        cov_window > 0,
        100.0 * nmod_window / cov_window,
        0.0,
    )
    return chrom_df


def use_raw_signals(chrom_df):
    chrom_df["Coverage_plot"] = chrom_df["Nvalid_cov"]
    chrom_df["BrdU_plot"] = chrom_df["Nmod"]
    chrom_df["BrdU_pct_plot"] = chrom_df["BrdU_pct"]
    return chrom_df


def binned_signal(x, y, chr_length, bin_size):
    if bin_size <= 0:
        bin_size = 2500

    edges = np.arange(0, chr_length + bin_size, bin_size, dtype=np.int64)
    if edges.size < 2:
        edges = np.array([0, chr_length], dtype=np.int64)

    finite = np.isfinite(y)

    if not np.any(finite):
        centers = (edges[:-1] + edges[1:]) / 2.0
        return centers, np.zeros_like(centers, dtype=float)

    x = x[finite]
    y = y[finite]
    idx = np.clip(np.digitize(x, edges) - 1, 0, edges.size - 2)

    sums = np.zeros(edges.size - 1, dtype=float)
    counts = np.zeros(edges.size - 1, dtype=float)
    np.add.at(sums, idx, y)
    np.add.at(counts, idx, 1.0)

    output = np.zeros_like(sums)
    populated = counts > 0
    output[populated] = sums[populated] / counts[populated]
    centers = (edges[:-1] + edges[1:]) / 2.0
    return centers, output


def plot_stick_ruler(ax, chrom_df, x_col, y_col, chr_length, bin_size, max_height=1.0):
    x = chrom_df[x_col].to_numpy(dtype=np.float64)
    y = chrom_df[y_col].to_numpy(dtype=np.float64)
    centers, binned = binned_signal(x, y, chr_length=chr_length, bin_size=bin_size)
    bmax = float(np.nanmax(binned)) if binned.size else 0.0
    heights = (binned / bmax) * max_height if bmax > 0 else np.zeros_like(binned)

    ax.axhline(0, linewidth=0.8)
    ax.vlines(centers, 0, heights, linewidth=0.6)
    format_track_axis(ax, chr_length, max_height * 1.05)
    ax.text(
        0.995,
        0.82,
        f"max={bmax:.2f}",
        ha="right",
        va="center",
        transform=ax.transAxes,
        fontsize=8,
    )


def plot_feature_track(ax, features_df, chr_length, color):
    ax.axhline(0, linewidth=0.8)

    if len(features_df) > 0:
        mids = (features_df["start"].to_numpy(dtype=float) + features_df["end"].to_numpy(dtype=float)) / 2.0
        mids = np.sort(np.clip(mids, 0, chr_length))

        for index, position in enumerate(mids):
            y0, y1 = (0.15, 0.85) if index % 2 == 0 else (0.25, 0.75)
            ax.vlines(position, y0, y1, colors=color, linewidth=1.2, alpha=0.95)

    format_track_axis(ax, chr_length, 1.0, ymin=0.0)


def plot_subtel_track(ax, chr_length):
    ax.axhline(0, linewidth=0.8)
    left_end = min(SUBTEL_SIZE, chr_length)
    right_start = max(0, chr_length - SUBTEL_SIZE)
    ax.add_patch(
        plt.Rectangle((0, 0.15), left_end, 0.70, linewidth=0, facecolor="red", alpha=0.8)
    )
    ax.add_patch(
        plt.Rectangle(
            (right_start, 0.15),
            chr_length - right_start,
            0.70,
            linewidth=0,
            facecolor="red",
            alpha=0.8,
        )
    )
    format_track_axis(ax, chr_length, 1.0, ymin=0.0)


def format_track_axis(ax, chr_length, ymax, ymin=-0.02):
    ax.set_xlim(0, chr_length)
    ax.set_ylim(ymin, ymax)
    ax.set_yticks([])
    ax.grid(True, axis="x", alpha=0.08)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(axis="x", which="both", labelbottom=False, bottom=False)


def add_all_rulers_panel(fig, gs_slot, ax_sharex, chrom_df, chr_length, features, bin_size):
    sub = gs_slot.subgridspec(
        7,
        2,
        height_ratios=[0.75] * 7,
        width_ratios=[0.10, 0.90],
        hspace=0.10,
        wspace=0.0,
    )
    labels = ["Coverage", "BrdU", "BrdU %", "G4", "tRNA", "TE", "Subtel"]
    axes = []

    for index, label in enumerate(labels):
        ax_label = fig.add_subplot(sub[index, 0])
        ax_label.axis("off")
        ax_label.text(
            0.98,
            0.5,
            label,
            ha="right",
            va="center",
            fontsize=9,
            transform=ax_label.transAxes,
        )
        axes.append(fig.add_subplot(sub[index, 1], sharex=ax_sharex))

    ax_cov, ax_brd, ax_brd_pct, ax_g4, ax_trna, ax_te, ax_subtel = axes
    plot_stick_ruler(ax_cov, chrom_df, "start", "Coverage_plot", chr_length, bin_size)
    plot_stick_ruler(ax_brd, chrom_df, "start", "BrdU_plot", chr_length, bin_size)
    plot_stick_ruler(
        ax_brd_pct,
        chrom_df,
        "start",
        "BrdU_pct_plot",
        chr_length,
        bin_size,
        max_height=100.0,
    )
    plot_feature_track(ax_g4, features["g4"], chr_length, color="green")
    plot_feature_track(ax_trna, features["trna"], chr_length, color="purple")
    plot_feature_track(ax_te, features["te"], chr_length, color="orange")
    plot_subtel_track(ax_subtel, chr_length)
    ax_subtel.tick_params(axis="x", which="both", labelbottom=True, bottom=True)

    return {"coverage": ax_cov, "brdu": ax_brd, "brdu_pct": ax_brd_pct}


def replace_axis_label(ax, smooth, raw, label_smooth, label_raw):
    for text in list(ax.texts):
        if text.get_text().startswith("max="):
            text.remove()

    if label_raw:
        label_text = f"{label_smooth}={smooth} | {label_raw}={raw}"
    else:
        label_text = f"{label_smooth}={smooth}"

    ax.text(
        0.995,
        0.82,
        label_text,
        ha="right",
        va="center",
        transform=ax.transAxes,
        fontsize=8,
    )


def add_ruler_max_labels(ruler_axes, chrom_df, smoothed):
    if smoothed:
        replace_axis_label(
            ruler_axes["coverage"],
            f"{float(chrom_df['Coverage_plot'].max()):.2f}",
            f"{float(chrom_df['Nvalid_cov'].max()):.0f}",
            "smooth max",
            "raw max",
        )
        replace_axis_label(
            ruler_axes["brdu"],
            f"{float(chrom_df['BrdU_plot'].max()):.2f}",
            f"{float(chrom_df['Nmod'].max()):.0f}",
            "smooth max",
            "raw max",
        )
        replace_axis_label(
            ruler_axes["brdu_pct"],
            f"{float(chrom_df['BrdU_pct_plot'].max()):.2f}%",
            f"{float(chrom_df['BrdU_pct'].max()):.2f}%",
            "smooth max",
            "raw max",
        )
    else:
        replace_axis_label(
            ruler_axes["coverage"],
            f"{float(chrom_df['Nvalid_cov'].max()):.0f}",
            "",
            "max",
            "",
        )
        replace_axis_label(
            ruler_axes["brdu"],
            f"{float(chrom_df['Nmod'].max()):.0f}",
            "",
            "max",
            "",
        )
        replace_axis_label(
            ruler_axes["brdu_pct"],
            f"{float(chrom_df['BrdU_pct'].max()):.2f}%",
            "",
            "max",
            "",
        )


def downsampled_main_signal(chrom_df, chr_length, smoothed):
    if smoothed:
        return (
            chrom_df["start"].to_numpy(dtype=float),
            chrom_df["Coverage_plot"].to_numpy(dtype=float),
            chrom_df["BrdU_plot"].to_numpy(dtype=float),
            "line",
        )

    max_points = 200000
    stride = max(1, int(np.ceil(len(chrom_df) / max_points)))
    xs = chrom_df["start"].to_numpy(dtype=float)[::stride]
    coverage = chrom_df["Coverage_plot"].to_numpy(dtype=float)[::stride]
    brdu = chrom_df["BrdU_plot"].to_numpy(dtype=float)[::stride]
    display_bin_size = max(1, int(chr_length // 2000))
    centers, coverage_binned = binned_signal(xs, coverage, chr_length, display_bin_size)
    _, brdu_binned = binned_signal(xs, brdu, chr_length, display_bin_size)
    return centers, coverage_binned, brdu_binned, "step"


def plot_chromosome(chrom, chrom_df, features_by_chrom, output_dir, prefix, smoothed):
    chr_length = int(CHROM_LENGTHS.get(chrom, chrom_df["end"].max()))
    mode_label = "smoothed" if smoothed else "unsmoothed"

    fig = plt.figure(figsize=(15, 7))
    gs = fig.add_gridspec(3, 1, height_ratios=[3.5, 0.40, 2.5], hspace=0.15)
    ax = fig.add_subplot(gs[0, 0])

    ax.axvspan(0, SUBTEL_SIZE, alpha=0.2, color="red", label="Subtelomeric region")
    ax.axvspan(chr_length - SUBTEL_SIZE, chr_length, alpha=0.2, color="red")

    x, coverage, brdu, draw_style = downsampled_main_signal(chrom_df, chr_length, smoothed)

    if draw_style == "step":
        ax.fill_between(
            x,
            0,
            coverage,
            color="lightgrey",
            alpha=0.6,
            label="Coverage (Nvalid_cov, unsmoothed)",
            step="pre",
            rasterized=True,
        )
        ax.step(
            x,
            brdu,
            color="blue",
            linewidth=0.8,
            label="BrdU count (Nmod, unsmoothed)",
            where="pre",
            rasterized=True,
        )
    else:
        ax.fill_between(
            x,
            0,
            coverage,
            color="lightgrey",
            alpha=0.6,
            label="Coverage (Nvalid_cov, smoothed)",
        )
        ax.plot(x, brdu, color="blue", linewidth=1.2, label="BrdU count (Nmod, smoothed)")

    signal_max = max(
        float(np.nanmax(coverage)) if np.any(np.isfinite(coverage)) else 0.0,
        float(np.nanmax(brdu)) if np.any(np.isfinite(brdu)) else 0.0,
    )
    ax.set_ylim(0, signal_max * 1.05 if signal_max > 0 else 1.0)
    ax.set_xlim(0, chr_length)
    ax.set_ylabel("Read count")
    ax.set_xlabel("Genomic position (bp)")
    ax.set_title(f"BrdU pileup along chromosome {chrom} ({mode_label}; subtelomeres highlighted)")
    ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    ax.ticklabel_format(style="plain", axis="x")
    ax.grid(True, axis="x", alpha=0.15)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0)

    bin_size = max(50, int(chr_length // 2000))
    features = {
        "g4": features_by_chrom["g4"].get(chrom, empty_feature_df()),
        "trna": features_by_chrom["trna"].get(chrom, empty_feature_df()),
        "te": features_by_chrom["te"].get(chrom, empty_feature_df()),
    }
    ruler_axes = add_all_rulers_panel(fig, gs[2, 0], ax, chrom_df, chr_length, features, bin_size)
    add_ruler_max_labels(ruler_axes, chrom_df, smoothed=smoothed)

    plt.tight_layout()
    save_path = output_dir / f"{prefix}.chromosome_{chrom}.{mode_label}.png"
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    return save_path


def empty_feature_df():
    return pd.DataFrame(columns=["chrom", "start", "end", "name", "score", "strand"])


def group_features(features):
    grouped = {}
    for key, df in features.items():
        grouped[key] = {chrom: sub for chrom, sub in df.groupby("chrom", sort=False)}
    return grouped


def generate_genome_browser(
    positive_bedgraph,
    negative_bedgraph,
    output_dir,
    prefix,
    g4_bed=None,
    trna_bed=None,
    te_bed=None,
    smoothed=True,
    window=1000,
):
    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    df = prepare_dataframe(positive_bedgraph, negative_bedgraph)
    features = {
        "g4": load_feature_bed(g4_bed, "G4 motif"),
        "trna": load_feature_bed(trna_bed, "tRNA motif"),
        "te": load_feature_bed(te_bed, "TE motif"),
    }
    features_by_chrom = group_features(features)

    generated = []

    for chrom in sorted(df["chrom"].unique(), key=chrom_sort_key):
        chrom_df = combine_strands(df[df["chrom"] == chrom].copy())

        if smoothed:
            chrom_df = smooth_signals(chrom_df, window=window)
        else:
            chrom_df = use_raw_signals(chrom_df)

        generated.append(
            plot_chromosome(
                chrom=chrom,
                chrom_df=chrom_df,
                features_by_chrom=features_by_chrom,
                output_dir=output_dir,
                prefix=prefix,
                smoothed=smoothed,
            )
        )

    print(f"[INFO] Generated {len(generated)} plots in {output_dir}", file=sys.stderr)
    return generated


def main():
    args = parse_args()
    generate_genome_browser(
        positive_bedgraph=args.positive_bedgraph,
        negative_bedgraph=args.negative_bedgraph,
        output_dir=args.output_dir,
        prefix=args.prefix,
        g4_bed=args.g4_bed,
        trna_bed=args.trna_bed,
        te_bed=args.te_bed,
        smoothed=True,
        window=args.window,
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
spike_overview.py — one figure summarising intrahost polymorphism across every
sample, along the whole spike CDS.

Each sample is mapped to its OWN spike CDS reference, so the per-sample
coordinates do not line up with each other. This script uses a codon-aware
multiple sequence alignment of those references to put every sample on one
common amino-acid axis, then draws:

  A  the within-sample non-reference fraction of every position, coloured and
     shaped by sample. Every frequency is computed inside one animal against
     that animal's own consensus reference; samples are never pooled to make a
     frequency, only placed on a shared axis. The sequencing-error floor stays
     visible as the grey cloud, and positions where the non-reference base is
     the majority (i.e. the reference is not that sample's consensus there —
     a between-animal difference) are excluded, and by default not drawn.
  B  a heat map of the same values, one row per animal.
  C  how many of the animals carry an intrahost variant in each window, over a
     pale bar giving how many had usable coverage there. Counting animals, not
     sites, keeps one noisy sample from inventing a peak.

The furin cleavage site is marked across all three panels.

Outputs
  <prefix>.pdf / .png          the figure
  <prefix>.sites.tsv           every position x sample: depth, counts, fraction
  <prefix>.by_codon.tsv        per alignment codon: how many samples, max fraction

Example
  python3 spike_overview.py \
      --samples samples.tsv \
      --msa refs/FECV_refs_all_spike_cds_linsi.fas \
      --fcs-sample fecv_NARC4 --fcs-aa 791-795 \
      --font /Users/skawakubo/python/fonts/Arial.ttf \
      --out spike_overview

samples.tsv is a headerless TSV of:  label <TAB> bam <TAB> reference_fasta <TAB> contig
The contig must match a record name in the MSA.

Requires: pysam, matplotlib, numpy
"""

from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np
import pysam
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

BASES = ["A", "C", "G", "T"]
INK = "#1F2328"
INK_SOFT = "#5A6169"
GRID = "#E3E5E8"
MASK = "#DDE0E4"          # not enough coverage to say anything
NOISE = "#B9BEC4"         # the everything-included cloud in panel A
CALL = "#2166AC"          # sites at or above the display level
FCS = "#B2182B"
SEQ_RAMP = [              # single-hue sequential ramp, light -> dark
    "#DCE6F1", "#B9D0E5", "#8FB4D9", "#6394C8", "#3B72B0", "#245C96", "#123F6D",
]

plt.rcParams.update({
    "font.size": 8,
    "axes.edgecolor": INK_SOFT,
    "axes.labelcolor": INK,
    "axes.linewidth": 0.7,
    "xtick.color": INK_SOFT,
    "ytick.color": INK_SOFT,
    "text.color": INK,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def setup_font(path=None):
    matplotlib.rcParams["svg.fonttype"] = "none"
    cands = [path] if path else [
        "/Users/skawakubo/python/fonts/Arial.ttf",
        "/Library/Fonts/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Arial.ttf",
    ]
    for c in cands:
        if c and os.path.exists(c):
            from matplotlib import font_manager
            font_manager.fontManager.addfont(c)
            name = font_manager.FontProperties(fname=c).get_name()
            matplotlib.rcParams["font.family"] = name
            print(f"[info] font: {name}  ({c})", file=sys.stderr)
            return name
    if path:
        print(f"[warning] font not found: {path} — using the default",
              file=sys.stderr)
    return None


# --------------------------------------------------------------------------
# alignment
# --------------------------------------------------------------------------
def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


def check_codon_aware(aln):
    """A codon axis is only meaningful if every gap is a whole, in-frame codon."""
    length = len(next(iter(aln.values())))
    problems = []
    if length % 3:
        problems.append(f"alignment length {length} is not a multiple of 3")
    for name, seq in aln.items():
        for m in re.finditer(r"-+", seq):
            if (m.end() - m.start()) % 3 or m.start() % 3:
                problems.append(f"{name}: gap at {m.start() + 1}-{m.end()} "
                                f"is not an in-frame whole codon")
    return length, problems


def column_map(seq):
    """1-based ungapped position -> 1-based alignment column."""
    out = {}
    p = 0
    for col, ch in enumerate(seq, start=1):
        if ch != "-":
            p += 1
            out[p] = col
    return out


# --------------------------------------------------------------------------
# counting
# --------------------------------------------------------------------------
def scan_sample(bam, ref_fa, chrom, min_bq, min_mq, max_depth):
    """Per-position depth and non-reference fraction over the whole contig."""
    fa = pysam.FastaFile(ref_fa)
    ref = fa.fetch(chrom).upper()
    fa.close()
    n = len(ref)
    counts = {b: np.zeros(n, dtype=np.int32) for b in BASES}
    dele = np.zeros(n, dtype=np.int32)

    af = pysam.AlignmentFile(bam, "rb")
    for col in af.pileup(chrom, 0, n, truncate=True, stepper="samtools",
                         min_base_quality=min_bq, min_mapping_quality=min_mq,
                         max_depth=max_depth, ignore_overlaps=True,
                         ignore_orphans=False,     # must match isnv_pileup.py:
                         # pysam defaults this to True, which silently discards
                         # reads whose mate did not map to this small reference
                         # and made the depths here disagree with the per-sample
                         # figures
                         compute_baq=False):
        i = col.reference_pos
        if not (0 <= i < n):
            continue
        for pr in col.pileups:
            if pr.is_refskip:
                continue
            if pr.is_del:
                dele[i] += 1
                continue
            q = pr.query_position
            if q is None:
                continue
            b = pr.alignment.query_sequence[q].upper()
            if b in counts:
                counts[b][i] += 1
    af.close()

    depth = sum(counts[b] for b in BASES) + dele
    refcount = np.zeros(n, dtype=np.int32)
    for i, rb in enumerate(ref):
        if rb in counts:
            refcount[i] = counts[rb][i]
    alt = depth - refcount
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(depth > 0, alt / np.maximum(depth, 1), 0.0)
    return dict(ref=ref, counts=counts, dele=dele, depth=depth,
                alt=alt, frac=frac)



# --------------------------------------------------------------------------
# optional panel C: a per-base pileup for one representative animal, drawn by
# reusing isnv_pileup.py so the two figures cannot drift apart
# --------------------------------------------------------------------------
def draw_pileup_panel(fig, slot, args, letter=None):
    """Render one sample's FCS pileup into `slot` of an existing figure."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import isnv_pileup as ip

    a1, a2 = (int(x) for x in args.pileup_aa.split("-"))
    orf = args.pileup_orf_start
    region = f"{args.pileup_chrom}:{orf + (a1 - 1) * 3}-{orf + a2 * 3 - 1}"
    argv = ["--bam", f"{args.pileup_bam}={args.pileup_label}",
            "--ref", args.pileup_ref, "--region", region,
            "--orf-start", str(orf), "--out", "/tmp/_panelc"]
    if args.pileup_highlight_aa:
        h1, h2 = (int(x) for x in args.pileup_highlight_aa.split("-"))
        argv += ["--highlight", f"{orf + (h1 - 1) * 3}-{orf + h2 * 3 - 1}"]
    argv += ["--min-depth", str(args.pileup_min_depth),
             "--min-aa-depth", str(args.pileup_min_depth)]
    pargs = ip.build_parser().parse_args(argv)

    chrom, start, end = ip.parse_region(region)
    ref_fa = pysam.FastaFile(args.pileup_ref)
    c = ip.count_window(args.pileup_bam, chrom, start, end, ref_fa,
                        pargs.min_bq, pargs.min_mq, pargs.max_depth)
    codons = ip.codon_map(start, end, orf, None)
    cc = ip.count_codon_haplotypes(args.pileup_bam, chrom, codons,
                                   pargs.min_bq, pargs.min_mq)
    prof = ip.aa_profile(codons, cc, ref_fa, chrom, pargs.min_aa_depth,
                         pargs.min_af, pargs.aa_noise_floor, pargs.min_alt_reads)
    rows = ip.site_table(c, pargs.min_depth, pargs.min_af,
                         pargs.max_strand_bias_p, pargs.metric,
                         pargs.min_alt_reads)
    rows = ip.annotate_aa(rows, ref_fa, chrom, orf, None)
    aa_map = {r["pos"]: r["aa_change"] for r in rows if r.get("aa_change")}
    ref_fa.close()

    highlight = ip.parse_span(pargs.highlight) if pargs.highlight else None
    label = args.pileup_label
    if letter:
        label = (f"{letter}   per-base composition around the furin cleavage site "
                 f"in one representative animal\n{label}")
    ip.draw_pileup_block(fig, slot, c, label, pargs,
                         highlight=highlight, aa_map=aa_map, aa_prof=prof)
    return ip


# --------------------------------------------------------------------------
# figure
# --------------------------------------------------------------------------
def sample_style(i):
    """Composite encoding: 5 validated hues x 2 markers gives 10 distinct samples
    without inventing a 9th or 10th hue."""
    hues = ["#2166AC", "#1B7F5A", "#B8860B", "#B2182B", "#762A83"]
    marks = ["o", "^"]
    return hues[i % len(hues)], marks[(i // len(hues)) % len(marks)]


def build_figure(samples, data, colmap, aln_len, args, fcs_cols):
    """Every value plotted is a WITHIN-sample frequency.

    Each BAM is counted against that sample's own consensus reference, so a
    point is the fraction of that one animal's viral population carrying a
    non-reference base. Nothing is pooled to make a frequency; samples are only
    brought together on a shared coordinate axis. Positions where the
    non-reference base is at or above --max-af are NOT intrahost variation --
    there the reference simply is not that sample's consensus -- so they are
    drawn apart and excluded from the counts.
    """
    from matplotlib.gridspec import GridSpec

    n_codon = aln_len // 3
    n_s = len(samples)
    floor = args.floor
    lo_af, hi_af = args.level, args.max_af

    mat = np.full((n_s, n_codon), np.nan)
    per_sample_pts = []                       # (x, y) per sample, iSNVs only
    noise_x, noise_y = [], []                 # every position
    fix_x, fix_y = [], []                     # near-fixed = reference mismatch
    has_data = np.zeros((n_s, n_codon), dtype=bool)
    has_isnv = np.zeros((n_s, n_codon), dtype=bool)

    for si, s in enumerate(samples):
        d = data[s["label"]]
        cm = colmap[s["contig"]]
        xs, ys = [], []
        for p1 in range(1, len(d["ref"]) + 1):
            col = cm.get(p1)
            if col is None:
                continue
            dp = int(d["depth"][p1 - 1])
            if dp < args.min_depth:
                continue
            ci = (col - 1) // 3
            aa = ci + 1
            fr = float(d["frac"][p1 - 1])
            ac = int(d["alt"][p1 - 1])
            has_data[si, ci] = True
            if fr > hi_af and not args.keep_mismatch:
                # The reference base is not this animal's consensus here, so the
                # value is a difference between animals, not within one. It is
                # already excluded from every count; by default it is not drawn
                # either, which keeps panel A to the 0-50% range it measures.
                # --keep-mismatch puts the crosses back.
                continue
            noise_x.append(aa)
            noise_y.append(max(fr, floor * 0.35) if fr > 0 else floor * 0.35)
            if fr > hi_af:
                fix_x.append(aa)
                fix_y.append(fr)
            elif fr >= lo_af and ac >= args.min_alt_reads:
                xs.append(aa)
                ys.append(fr)
                has_isnv[si, ci] = True
                prev = mat[si, ci]
                mat[si, ci] = fr if np.isnan(prev) else max(prev, fr)
        per_sample_pts.append((xs, ys))
    # covered but no iSNV -> just under the colour scale, so it renders white
    for si in range(n_s):
        for ci in range(n_codon):
            if has_data[si, ci] and np.isnan(mat[si, ci]):
                mat[si, ci] = floor * 0.1

    mode = args.figure
    want_pileup = mode == "main" and args.pileup_bam
    if mode == "supp":
        rows_h, keys = [0.34 * n_s + 0.5], ["heat"]
    elif mode == "main":
        rows_h, keys = [2.7, 1.3], ["scatter", "density"]
        if want_pileup:
            rows_h.append(3.1)
            keys.append("pileup")
    else:
        rows_h = [2.7, 0.34 * n_s + 0.5, 1.3]
        keys = ["scatter", "heat", "density"]

    fig = plt.figure(figsize=(13.5, sum(rows_h) + 1.0))
    gs = GridSpec(len(rows_h), 1, figure=fig, height_ratios=rows_h,
                  hspace=0.55 if want_pileup else 0.30)
    slot = {k: gs[i] for i, k in enumerate(keys)}
    axA = fig.add_subplot(slot["scatter"]) if "scatter" in slot else None
    first = axA
    axB = None
    if "heat" in slot:
        axB = fig.add_subplot(slot["heat"], sharex=first) if first \
            else fig.add_subplot(slot["heat"])
        first = first or axB
    axC = fig.add_subplot(slot["density"], sharex=first) if "density" in slot else None
    letters = {}
    for i, k in enumerate(keys):
        letters[k] = "ABC"[i]
    for ax in (a for a in (axA, axB, axC) if a is not None):
        ax.set_xlim(0.5, n_codon + 0.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    fcs_aa = None
    if fcs_cols:
        fcs_aa = ((fcs_cols[0] - 1) // 3 + 1, (fcs_cols[1] - 1) // 3 + 1)
        for ax in (a for a in (axA, axB, axC) if a is not None):
            ax.axvspan(fcs_aa[0] - 0.5, fcs_aa[1] + 0.5, color=FCS, alpha=0.13,
                       linewidth=0, zorder=1)

    # ---------------- A ----------------
    if axA is not None:
        axA.set_yscale("log")
        axA.set_ylim(floor * 0.3, 1.4 if args.keep_mismatch else hi_af)
        axA.grid(axis="y", color=GRID, linewidth=0.6, zorder=0)
        axA.set_axisbelow(True)
        axA.scatter(noise_x, noise_y, s=1.2, color=NOISE, alpha=0.3, linewidths=0,
                    zorder=2, rasterized=True)
        if fix_x:
            axA.scatter(fix_x, fix_y, s=14, color="#7C848C", linewidths=0.8,
                        marker="x", zorder=3)
        for si, s in enumerate(samples):
            xs, ys = per_sample_pts[si]
            if args.colour_by_sample:
                colr, mark = sample_style(si)
            else:
                colr, mark = CALL, "o"
            axA.scatter(xs, ys, s=14, color=colr, marker=mark, alpha=0.75,
                        linewidths=0.3, edgecolors="white", zorder=4)
        for y, lab in ((lo_af, f"{lo_af * 100:g}%"), (hi_af, f"{hi_af * 100:g}%")):
            axA.axhline(y, color=INK_SOFT, linestyle=(0, (3, 2)), linewidth=0.7,
                        zorder=3)
            axA.annotate(lab, (n_codon, y), xytext=(4, 0), textcoords="offset points",
                         va="center", fontsize=6, color=INK_SOFT,
                         annotation_clip=False)
        axA.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: f"{v * 100:g}"))
        axA.set_ylabel("within-sample\nnon-reference fraction (%)", fontsize=8)
        axA.tick_params(axis="x", labelbottom=False, length=0)
        sub = ("colour and shape identify the animal" if args.colour_by_sample
               else "each point is one position in one animal")
        axA.set_title(f"{letters['scatter']}   intrahost frequency of every "
                      f"position, computed separately within each animal — {sub}",
                      fontsize=9, loc="left", pad=6)

    # ---------------- B ----------------
    if axB is not None:
        cmap = ListedColormap(SEQ_RAMP)
        cmap.set_bad(MASK)
        cmap.set_under("#FFFFFF")
        im = axB.imshow(np.ma.masked_invalid(mat), aspect="auto",
                        interpolation="nearest",
                        extent=(0.5, n_codon + 0.5, n_s - 0.5, -0.5),
                        cmap=cmap, norm=LogNorm(vmin=floor, vmax=hi_af), zorder=2)
        med = {s["label"]: int(np.median(data[s["label"]]["depth"]
                                        [data[s["label"]]["depth"] > 0]))
               for s in samples}
        axB.set_yticks(range(n_s))
        axB.set_yticklabels([f"{s['label']}  ({med[s['label']]}x)" for s in samples],
                            fontsize=7)
        if args.colour_by_sample:
            for tick, si in zip(axB.get_yticklabels(), range(n_s)):
                tick.set_color(sample_style(si)[0])
        axB.tick_params(axis="y", length=0)
        if mode == "supp":
            axB.set_xlabel("amino acid position in the spike CDS alignment",
                           fontsize=8)
        else:
            axB.tick_params(axis="x", labelbottom=False, length=0)
        axB.set_ylim(n_s - 0.5, -0.5)
        for sp in ("left", "bottom"):
            axB.spines[sp].set_visible(False)
        axB.set_title(f"{letters['heat']}   one row per animal — white: covered, "
                      "no intrahost variant · grey: coverage too low to judge",
                      fontsize=9, loc="left", pad=6)
        cax = axB.inset_axes([1.008, 0.0, 0.010, 1.0])
        cb = fig.colorbar(im, cax=cax)
        cb.set_label("within-sample fraction (%)", fontsize=7)
        cb.ax.tick_params(labelsize=6, length=2)
        cb.ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: f"{v * 100:g}"))

    # ---------------- C ----------------
    if axC is not None:
        b = args.bin_size
        nb = int(np.ceil(n_codon / b))
        n_with_isnv = np.zeros(nb)
        n_with_data = np.zeros(nb)
        for k in range(nb):
            sl = slice(k * b, min((k + 1) * b, n_codon))
            n_with_isnv[k] = has_isnv[:, sl].any(axis=1).sum()
            n_with_data[k] = has_data[:, sl].any(axis=1).sum()
        centres = (np.arange(nb) + 0.5) * b
        axC.bar(centres, n_with_data, width=b * 0.95, color="#EDEFF1",
                linewidth=0, zorder=2)
        axC.bar(centres, n_with_isnv, width=b * 0.95, color=CALL, linewidth=0,
                zorder=3)
        axC.grid(axis="y", color=GRID, linewidth=0.6, zorder=0)
        axC.set_axisbelow(True)
        axC.set_ylim(0, n_s + 0.5)
        axC.yaxis.set_major_locator(MaxNLocator(5, integer=True))
        axC.set_ylabel(f"animals with an\nintrahost variant\nper {b} codons",
                       fontsize=8)
        axC.set_xlabel("amino acid position in the spike CDS alignment", fontsize=8)
        axC.set_title(f"{letters['density']}   how many of the {n_s} animals carry "
                      f"an intrahost variant in each {b}-codon window "
                      f"(pale bar = animals with usable coverage there)",
                      fontsize=9, loc="left", pad=6)

    if want_pileup:
        ip = draw_pileup_panel(fig, slot["pileup"], args,
                               letter=letters["pileup"])
        fig.text(0.0, 0.0, "", alpha=0)   # keep tight bbox stable

    if fcs_aa and axA is not None:
        mid = (fcs_aa[0] + fcs_aa[1]) / 2
        axA.annotate("FCS", (mid, 0.965), xycoords=("data", "axes fraction"),
                     ha="center", va="top", fontsize=8, color=FCS,
                     fontweight="bold")
        for ax in (a for a in (axA, axC) if a is not None):
            ax.axvline(mid, color=FCS, linewidth=0.6, alpha=0.5,
                       linestyle=(0, (2, 2)), zorder=2)

    if mode == "supp":
        handles = [Patch(facecolor=SEQ_RAMP[3], label="intrahost variant"),
                   Patch(facecolor="#FFFFFF", edgecolor="#C9CDD2",
                         label="covered, no variant"),
                   Patch(facecolor=MASK, label=f"depth < {args.min_depth}"),
                   Patch(facecolor=FCS, alpha=0.18, label="furin cleavage site")]
        fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.0),
                   ncol=4, frameon=False, fontsize=7)
        if args.title:
            fig.suptitle(args.title, fontsize=11, x=0.0, ha="left", y=1.06,
                         va="bottom")
        for ext in ("pdf", "png", "svg"):
            fig.savefig(f"{args.out}.{ext}", dpi=300)
        plt.close(fig)
        return mat, has_isnv, has_data, fcs_aa

    if args.colour_by_sample:
        handles = [plt.Line2D([], [], marker=sample_style(i)[1], linestyle="none",
                              color=sample_style(i)[0], markersize=5,
                              label=s["label"]) for i, s in enumerate(samples)]
    else:
        handles = [plt.Line2D([], [], marker="o", linestyle="none", color=CALL,
                              markersize=5,
                              label="intrahost variant (one position in one animal)")]
    handles += [
        plt.Line2D([], [], marker="o", linestyle="none", color=NOISE,
                   markersize=4, label="all positions (error floor)"),
    ]
    if args.keep_mismatch:
        handles.append(
            plt.Line2D([], [], marker="x", linestyle="none", color="#8A9099",
                       markersize=5,
                       label=f"> {hi_af * 100:g}%: reference mismatch, "
                             f"not intrahost"))
    handles += [
        Patch(facecolor=FCS, alpha=0.18, label="furin cleavage site"),
    ]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.0),
               ncol=7 if args.colour_by_sample else 4, frameon=False,
               fontsize=7, columnspacing=1.1, handletextpad=0.4)
    if args.title:
        fig.suptitle(args.title, fontsize=11, x=0.0, ha="left", y=1.075,
                     va="bottom")

    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{args.out}.{ext}", dpi=300)
    plt.close(fig)
    return mat, has_isnv, has_data, fcs_aa


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", required=True,
                    help="TSV: label, bam, reference fasta, contig (no header)")
    ap.add_argument("--msa", required=True,
                    help="codon-aware alignment of the per-sample references")
    ap.add_argument("--out", default="spike_overview")
    ap.add_argument("--fcs-sample", help="contig whose coordinates --fcs-aa uses")
    ap.add_argument("--fcs-aa", help="FCS amino-acid range in that sample, e.g. 791-795")
    ap.add_argument("--min-depth", type=int, default=20,
                    help="positions below this depth are shown as 'no data' "
                         "rather than as 0%% (default 20). This is a validity "
                         "mask, not a frequency threshold.")
    ap.add_argument("--min-alt-reads", type=int, default=5, metavar="N",
                    help="a position counts as polymorphic only if at least this "
                         "many reads carry the non-reference base (default 5). "
                         "Coverage here spans 18x to 2500x, so a frequency alone "
                         "cannot separate signal from error: at 20x a single "
                         "misread base is 5%%. This is a read-count floor, not a "
                         "frequency threshold — panel A still shows every "
                         "position and the TSV still contains all of them.")
    ap.add_argument("--floor", type=float, default=0.001,
                    help="bottom of the log axis and colour scale (default 0.001)")
    ap.add_argument("--level", type=float, default=0.01,
                    help="lower bound of the intrahost window (default 0.01). "
                         "Panel A still shows every position as a grey dot and "
                         "the TSV still lists all of them.")
    ap.add_argument("--max-af", type=float, default=0.5, metavar="F",
                    help="upper bound of the intrahost window (default 0.5). "
                         "Above this the non-reference base is the majority, "
                         "which means the reference is not that sample's own "
                         "consensus at that position — a between-sample "
                         "difference, not variation within the animal. Such "
                         "positions are excluded from every count and, by "
                         "default, are not drawn at all; see --keep-mismatch.")
    ap.add_argument("--keep-mismatch", action="store_true",
                    help="draw the positions above --max-af in panel A as grey "
                         "crosses, and extend the y axis above the cut to make "
                         "room for them. Off by default, so panel A shows only "
                         "the 0 to --max-af range that the analysis measures.")
    ap.add_argument("--colour-by-sample", "--color-by-sample",
                    dest="colour_by_sample", action="store_true",
                    help="give each animal its own colour and marker in panel A. "
                         "Off by default: the points are one colour, which stays "
                         "readable as the number of animals grows. Panel B "
                         "already carries per-animal identity.")
    ap.add_argument("--bin-size", type=int, default=3,
                    help="codons per bin in panel C (default 3)")
    ap.add_argument("--min-bq", type=int, default=20)
    ap.add_argument("--min-mq", type=int, default=20)
    ap.add_argument("--max-depth", type=int, default=1_000_000)
    ap.add_argument("--figure", choices=["full", "main", "supp"], default="full",
                    help="'full' (default): scatter + heat map + density in one "
                         "figure. 'main': scatter + density, plus a pileup panel "
                         "for one representative animal if --pileup-bam is given. "
                         "'supp': the heat map alone.")
    g = ap.add_argument_group("panel C of --figure main (one representative animal)")
    g.add_argument("--pileup-bam")
    g.add_argument("--pileup-ref")
    g.add_argument("--pileup-chrom")
    g.add_argument("--pileup-label", default="representative sample")
    g.add_argument("--pileup-aa", help="window in that sample's aa coordinates, e.g. 786-800")
    g.add_argument("--pileup-highlight-aa", help="e.g. 791-795")
    g.add_argument("--pileup-orf-start", type=int, default=1)
    g.add_argument("--pileup-min-depth", type=int, default=50)
    ap.add_argument("--font", default=None)
    ap.add_argument("--title", default=None)
    args = ap.parse_args(argv)

    setup_font(args.font)

    samples = []
    with open(args.samples) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 4:
                raise SystemExit(f"[error] samples file needs 4 columns: {line!r}")
            samples.append(dict(label=f[0], bam=f[1], ref=f[2], contig=f[3]))

    aln = read_fasta(args.msa)
    aln_len, problems = check_codon_aware(aln)
    if problems:
        print("[warning] the alignment is not strictly codon-aware; the amino-acid "
              "axis may be off:", file=sys.stderr)
        for p in problems[:10]:
            print("          " + p, file=sys.stderr)
    missing = [s["contig"] for s in samples if s["contig"] not in aln]
    if missing:
        raise SystemExit(f"[error] not in the alignment: {missing}\n"
                         f"        records present: {sorted(aln)}")
    colmap = {name: column_map(seq) for name, seq in aln.items()}

    fcs_cols = None
    if args.fcs_sample and args.fcs_aa:
        a1, a2 = (int(x) for x in args.fcs_aa.split("-"))
        cm = colmap[args.fcs_sample]
        fcs_cols = (cm[(a1 - 1) * 3 + 1], cm[a2 * 3])
        print(f"[info] FCS {args.fcs_sample} aa {a1}-{a2} -> alignment columns "
              f"{fcs_cols[0]}-{fcs_cols[1]} (alignment aa "
              f"{(fcs_cols[0] - 1) // 3 + 1}-{(fcs_cols[1] - 1) // 3 + 1})",
              file=sys.stderr)

    data = {}
    for s in samples:
        print(f"[scan] {s['label']} ...", file=sys.stderr)
        data[s["label"]] = scan_sample(s["bam"], s["ref"], s["contig"],
                                       args.min_bq, args.min_mq, args.max_depth)

    # ---- per-position table (everything, no filtering) ----
    path = f"{args.out}.sites.tsv"
    with open(path, "w") as fh:
        fh.write("\t".join([
            "sample", "contig", "pos_sample", "aln_col", "aln_aa", "codon_pos",
            "ref_base", "depth", "A", "C", "G", "T", "del",
            "alt_count", "alt_fraction", "class"]) + "\n")
        for s in samples:
            d = data[s["label"]]
            cm = colmap[s["contig"]]
            for p1 in range(1, len(d["ref"]) + 1):
                dp = int(d["depth"][p1 - 1])
                if dp < args.min_depth:
                    continue
                col = cm.get(p1)
                if col is None:
                    continue
                fh.write("\t".join(map(str, [
                    s["label"], s["contig"], p1, col, (col - 1) // 3 + 1,
                    (col - 1) % 3 + 1, d["ref"][p1 - 1], dp,
                    int(d["counts"]["A"][p1 - 1]), int(d["counts"]["C"][p1 - 1]),
                    int(d["counts"]["G"][p1 - 1]), int(d["counts"]["T"][p1 - 1]),
                    int(d["dele"][p1 - 1]), int(d["alt"][p1 - 1]),
                    f"{float(d['frac'][p1 - 1]):.6f}",
                    ("reference_mismatch"
                     if float(d["frac"][p1 - 1]) > args.max_af else
                     "intrahost_variant"
                     if (float(d["frac"][p1 - 1]) >= args.level
                         and int(d["alt"][p1 - 1]) >= args.min_alt_reads)
                     else "no_variant")])) + "\n")
    print(f"[written] {path}", file=sys.stderr)

    mat, has_isnv, has_data, fcs_aa = build_figure(
        samples, data, colmap, aln_len, args, fcs_cols)

    # ---- per-codon summary ----
    path = f"{args.out}.by_codon.tsv"
    n_codon = aln_len // 3
    with open(path, "w") as fh:
        fh.write("\t".join([
            "aln_aa", "in_fcs", "n_samples_with_coverage",
            "n_samples_with_intrahost_variant", "fraction_of_covered_samples",
            "max_within_sample_fraction", "samples_with_variant"]) + "\n")
        for ci in range(n_codon):
            nd = int(has_data[:, ci].sum())
            if not nd:
                continue
            names = [samples[i]["label"] for i in range(len(samples))
                     if has_isnv[i, ci]]
            colv = mat[:, ci]
            best = np.nanmax(colv) if np.isfinite(colv).any() else 0.0
            fh.write("\t".join(map(str, [
                ci + 1,
                "yes" if (fcs_aa and fcs_aa[0] <= ci + 1 <= fcs_aa[1]) else "no",
                nd, len(names), f"{len(names) / nd:.4f}",
                f"{max(best, 0.0):.6f}", ";".join(names) or "-"])) + "\n")
    print(f"[written] {path}", file=sys.stderr)
    print(f"[written] {args.out}.pdf / {args.out}.png", file=sys.stderr)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
isnv_pileup.py — visualise intrahost polymorphic sites from mapped BAM files.

Designed for viral amplicon / shotgun NGS (e.g. FCoV spike, furin cleavage site),
but works for any reference.

Two modes (can be combined):

  1. REGION mode  (--region CHR:START-END)
     Draw a per-base stacked-bar "pileup" figure for an explicitly given window.

  2. SCAN mode    (--scan CHR:START-END)
     Scan a larger window (e.g. the whole spike ORF), auto-detect sites whose
     minor-allele frequency passes the filters, then draw
       (a) an overview MAF-vs-position plot for the scanned window, and
       (b) a zoom panel (+/- --flank bp) around each detected site.

With --orf-start the figures also carry an amino-acid ribbon under the x-axis:
one box per codon labelled with the one-letter code and the AA position, with
non-synonymous changes, synonymous-only variation and deletions called out.
Amino acids are counted from reads spanning the whole codon, so two variants in
the same codon are resolved correctly rather than being combined per position.

Outputs
  <prefix>.counts.tsv   per-position base counts / fractions (all plotted positions)
  <prefix>.codons.tsv   per-codon amino-acid frequencies (needs --orf-start)
  <prefix>.sites.tsv    detected polymorphic sites with stats (scan mode)
  <prefix>.region.pdf   region-mode figure (+ .png)
  <prefix>.scan.pdf     scan overview figure (+ .png)
  <prefix>.sites.pdf    per-site zoom panels (+ .png)

Examples
  # explicit window around the FCS
  python isnv_pileup.py --bam cat1.bam --ref FCoV_ref.fasta \
      --region NC_002306.3:23600-23700 --orf-start 20400 --highlight 23640-23655 \
      --out fcs_cat1

  # scan the whole spike, auto-detect polymorphic sites
  python isnv_pileup.py --bam cat1.bam=Cat1 --bam cat2.bam=Cat2 --ref FCoV_ref.fasta \
      --scan NC_002306.3:20400-24800 --min-af 0.03 --min-depth 200 --out spike_scan

Requires: pysam, matplotlib, numpy
"""

from __future__ import annotations

import argparse
import bisect
import math
import os
import re
import sys
from dataclasses import dataclass, field

import numpy as np
import pysam
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

# --------------------------------------------------------------------------
# Style — colourblind-safe base palette that still follows the IGV convention
# (A green, C blue, G orange/brown, T red). Validated: all adjacent pairs pass
# deutan/protan/tritan separation and >= 3:1 contrast on a white surface.
# --------------------------------------------------------------------------
BASES = ["A", "C", "G", "T"]
BASE_COLORS = {
    "A": "#1B7F5A",
    "C": "#2166AC",
    "G": "#B8860B",
    "T": "#B2182B",
}
DEL_COLOR = "#8C8C8C"      # deletion (read supports a gap at this position)
INS_COLOR = "#762A83"      # insertion marker
COV_COLOR = "#C9CDD2"      # coverage fill
GRID = "#E3E5E8"
INK = "#1F2328"
INK_SOFT = "#5A6169"

DEFAULT_FONT_CANDIDATES = [
    "/Users/skawakubo/python/fonts/Arial.ttf",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "Arial.ttf"),
    "/Library/Fonts/Arial.ttf",
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    "C:/Windows/Fonts/arial.ttf",
]


def setup_font(path=None, quiet=False):
    """Register a TrueType font file and make it the figure font.

    Text stays as text: svg.fonttype='none' keeps SVG text editable and
    pdf.fonttype=42 embeds the TrueType outlines, so both formats open cleanly
    in Illustrator with the font still selectable.
    """
    matplotlib.rcParams["svg.fonttype"] = "none"
    candidates = [path] if path else DEFAULT_FONT_CANDIDATES
    for cand in candidates:
        if cand and os.path.exists(cand):
            try:
                from matplotlib import font_manager
                font_manager.fontManager.addfont(cand)
                name = font_manager.FontProperties(fname=cand).get_name()
            except Exception as exc:                       # unreadable / not a font
                print(f"[warning] could not load font {cand}: {exc}", file=sys.stderr)
                continue
            matplotlib.rcParams["font.family"] = name
            if not quiet:
                print(f"[info] font: {name}  ({cand})", file=sys.stderr)
            return name
    if path:
        print(f"[warning] font file not found: {path}\n"
              f"          falling back to {matplotlib.rcParams['font.family']}. "
              f"Pass --font /path/to/Arial.ttf, or put Arial.ttf next to this script.",
              file=sys.stderr)
    return None


plt.rcParams.update(
    {
        "font.size": 8,
        "font.family": "DejaVu Sans",
        "axes.edgecolor": INK_SOFT,
        "axes.labelcolor": INK,
        "axes.linewidth": 0.7,
        "xtick.color": INK_SOFT,
        "ytick.color": INK_SOFT,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "text.color": INK,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,      # editable text in Illustrator
        "ps.fonttype": 42,
    }
)

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


# --------------------------------------------------------------------------
# Data structures
# --------------------------------------------------------------------------
@dataclass
class Sample:
    label: str
    path: str


@dataclass
class Counts:
    """Per-position counts for one sample over one window (1-based inclusive)."""

    chrom: str
    start: int
    end: int
    positions: np.ndarray                      # 1-based genomic positions
    base: dict = field(default_factory=dict)   # base -> np.ndarray of counts
    fwd: dict = field(default_factory=dict)    # base -> forward-strand counts
    rev: dict = field(default_factory=dict)    # base -> reverse-strand counts
    dele: np.ndarray = None                    # deletion counts
    ins: np.ndarray = None                     # insertion-start counts
    depth: np.ndarray = None                   # A+C+G+T+del (plotted denominator)
    ref: list = field(default_factory=list)    # reference base per position ('N' if none)


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def parse_region(s: str):
    m = re.match(r"^([^:]+):([\d,]+)-([\d,]+)$", s.strip())
    if not m:
        raise ValueError(f"region must look like CHR:START-END, got {s!r}")
    chrom = m.group(1)
    start = int(m.group(2).replace(",", ""))
    end = int(m.group(3).replace(",", ""))
    if end < start:
        raise ValueError(f"region end < start: {s!r}")
    return chrom, start, end


def parse_span(s: str):
    m = re.match(r"^([\d,]+)-([\d,]+)$", s.strip())
    if not m:
        raise ValueError(f"span must look like START-END, got {s!r}")
    return int(m.group(1).replace(",", "")), int(m.group(2).replace(",", ""))


def wilson_ci(k: int, n: int, z: float = 1.96):
    """Wilson score interval for a proportion — no scipy needed."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = (z / d) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, centre - half), min(1.0, centre + half))


def _lnfact(n: int) -> float:
    return math.lgamma(n + 1)


def fisher_exact_2x2(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p-value via lgamma (dependency-free)."""
    n = a + b + c + d
    if n == 0:
        return 1.0
    r1, r2 = a + b, c + d
    c1, c2 = a + c, b + d
    const = _lnfact(r1) + _lnfact(r2) + _lnfact(c1) + _lnfact(c2) - _lnfact(n)

    def logp(x):
        y = r1 - x
        z = c1 - x
        w = d + (a - x)
        if min(x, y, z, w) < 0:
            return None
        return const - (_lnfact(x) + _lnfact(y) + _lnfact(z) + _lnfact(w))

    obs = logp(a)
    if obs is None:
        return 1.0
    tot = 0.0
    for x in range(max(0, c1 - r2), min(r1, c1) + 1):
        lp = logp(x)
        if lp is not None and lp <= obs + 1e-9:
            tot += math.exp(lp)
    return min(1.0, tot)


def strand_bias_p(alt_fwd, alt_rev, ref_fwd, ref_rev) -> float:
    return fisher_exact_2x2(alt_fwd, alt_rev, ref_fwd, ref_rev)


# --------------------------------------------------------------------------
# Pileup counting
# --------------------------------------------------------------------------
def count_window(bam_path, chrom, start, end, ref_fa=None, min_bq=20, min_mq=20,
                 max_depth=1_000_000, ignore_orphans=False) -> Counts:
    """Count A/C/G/T, deletions and insertion starts per position (1-based, inclusive)."""
    n = end - start + 1
    positions = np.arange(start, end + 1, dtype=int)
    c = Counts(chrom=chrom, start=start, end=end, positions=positions)
    for b in BASES:
        c.base[b] = np.zeros(n, dtype=int)
        c.fwd[b] = np.zeros(n, dtype=int)
        c.rev[b] = np.zeros(n, dtype=int)
    c.dele = np.zeros(n, dtype=int)
    c.ins = np.zeros(n, dtype=int)

    af = pysam.AlignmentFile(bam_path, "rb")
    if chrom not in af.references:
        af.close()
        raise SystemExit(
            f"[error] contig {chrom!r} not in {bam_path}.\n"
            f"        available: {', '.join(af.references[:20])}"
        )

    for col in af.pileup(
        chrom,
        start - 1,
        end,
        truncate=True,
        stepper="samtools",
        min_base_quality=min_bq,
        min_mapping_quality=min_mq,
        max_depth=max_depth,
        ignore_overlaps=True,
        ignore_orphans=ignore_orphans,
        compute_baq=False,
    ):
        i = col.reference_pos + 1 - start
        if i < 0 or i >= n:
            continue
        for pr in col.pileups:
            if pr.is_refskip:
                continue
            if pr.is_del:
                c.dele[i] += 1
                continue
            qpos = pr.query_position
            if qpos is None:
                continue
            b = pr.alignment.query_sequence[qpos].upper()
            if b not in BASE_COLORS:
                continue
            c.base[b][i] += 1
            if pr.alignment.is_reverse:
                c.rev[b][i] += 1
            else:
                c.fwd[b][i] += 1
            if pr.indel > 0:
                c.ins[i] += 1
    af.close()

    c.depth = sum(c.base[b] for b in BASES) + c.dele

    if ref_fa is not None:
        seq = ref_fa.fetch(chrom, start - 1, end).upper()
        c.ref = list(seq) + ["N"] * (n - len(seq))
    else:
        # fall back to the observed majority base
        stack = np.vstack([c.base[b] for b in BASES])
        c.ref = [BASES[k] if c.depth[i] > 0 else "N"
                 for i, k in enumerate(stack.argmax(axis=0))]
    return c


# --------------------------------------------------------------------------
# Codon / amino-acid level counting
#
# Amino acids are counted from whole reads, not by combining per-position base
# frequencies: a read must cover all three bases of a codon to contribute. That
# is the only way to get the right answer when one codon carries two variants,
# because it preserves the linkage between them.
# --------------------------------------------------------------------------
def codon_map(window_start, window_end, orf_start, orf_end=None):
    """Codons overlapping a window -> [(aa_pos, (p1, p2, p3)), ...] (1-based, + strand)."""
    if not orf_start:
        return []
    first = max(0, (window_start - orf_start) // 3)
    last = (window_end - orf_start) // 3
    out = []
    for ci in range(first, last + 1):
        p1 = orf_start + ci * 3
        if p1 < orf_start:
            continue
        if orf_end and p1 + 2 > orf_end:
            break
        if p1 + 2 < window_start or p1 > window_end:
            continue
        out.append((ci + 1, (p1, p1 + 1, p1 + 2)))
    return out


def count_codon_haplotypes(bam_path, chrom, codons, min_bq=20, min_mq=20):
    """{aa_pos: {codon_string: n_fragments}} from reads spanning the whole codon.

    Counting is per FRAGMENT, not per read: with paired-end data whose mates
    overlap, both mates often cover the same codon and would otherwise be counted
    twice, inflating both the depth and any variant carried by that fragment.
    Mates that disagree at a codon are discarded for that codon.
    """
    if not codons:
        return {}
    wanted = set()
    for _, ps in codons:
        wanted.update(ps)
    lo, hi = min(wanted), max(wanted)
    codons = sorted(codons, key=lambda t: t[1][0])
    starts = [ps[0] for _, ps in codons]
    out = {aa: {} for aa, _ in codons}
    frag = {}          # (read name, aa_pos) -> codon string, or None if mates disagree
    discordant = 0

    af = pysam.AlignmentFile(bam_path, "rb")
    for read in af.fetch(chrom, lo - 1, hi):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.is_duplicate or read.is_qcfail:
            continue
        if read.mapping_quality < min_mq:
            continue
        seq = read.query_sequence
        quals = read.query_qualities
        if seq is None:
            continue
        pos2base, aligned = {}, set()
        for qpos, rpos in enumerate(read.get_reference_positions(full_length=True)):
            if rpos is None:
                continue
            rp1 = rpos + 1
            if rp1 in wanted:
                aligned.add(rp1)
                if quals is not None and quals[qpos] < min_bq:
                    continue
                pos2base[rp1] = seq[qpos].upper()
        r_start, r_end = read.reference_start + 1, read.reference_end
        i0 = bisect.bisect_left(starts, r_start)
        i1 = bisect.bisect_right(starts, r_end)
        for aa_pos, ps in codons[i0:i1]:
            trio = [pos2base.get(p) for p in ps]
            if all(b is not None and b in BASE_COLORS for b in trio):
                cod = "".join(trio)
            elif (r_start <= ps[0] and r_end >= ps[2]
                  and any(p not in aligned for p in ps)):
                # the read covers the codon but one or more bases are deleted
                cod = "del"
            else:
                continue  # partial coverage or low base quality — not counted
            key = (read.query_name, aa_pos)
            if key in frag:
                if frag[key] is not None and frag[key] != cod:
                    frag[key] = None          # mates disagree — drop this fragment
                    discordant += 1
                continue                      # already counted this fragment
            frag[key] = cod
    af.close()

    for (_, aa_pos), cod in frag.items():
        if cod is None:
            continue
        out[aa_pos][cod] = out[aa_pos].get(cod, 0) + 1
    if discordant:
        print(f"[note] {os.path.basename(bam_path)}: {discordant} codon(s) dropped "
              f"because the two mates of a read pair disagreed", file=sys.stderr)
    return out


def aa_profile(codons, codon_counts, ref_fa, chrom, min_depth, min_af,
               noise_floor=0.002, min_alt_reads=2, sub_frac=0.5, sub_floor=0.005):
    """Per codon: reference AA, observed AA frequencies, and the variant AAs.

    The non-synonymous call uses the TOTAL fraction of reads coding a different
    amino acid, not just the most frequent alternative — three different alt
    amino acids at 1.5% each are 4.5% of the population, not 1.5%.
    """
    prof = {}
    for aa_pos, ps in codons:
        counts = codon_counts.get(aa_pos, {})
        total = sum(counts.values())
        ref_cod = ""
        if ref_fa is not None:
            try:
                ref_cod = ref_fa.fetch(chrom, ps[0] - 1, ps[2]).upper()
            except (ValueError, KeyError):
                ref_cod = ""
        if len(ref_cod) < 3:
            subs = {k: v for k, v in counts.items() if k != "del"}
            ref_cod = max(subs, key=subs.get) if subs else ""
        aa_ref = CODON_TABLE.get(ref_cod, "X") if ref_cod else "?"

        aa_counts = {}
        for cod, n in counts.items():
            aa = "-" if cod == "del" else CODON_TABLE.get(cod, "X")
            aa_counts[aa] = aa_counts.get(aa, 0) + n
        del_freq = counts.get("del", 0) / total if total else 0.0
        alts = {a: n for a, n in aa_counts.items() if a not in (aa_ref, "-")}
        # total non-reference amino-acid fraction — this is what the call uses.
        # Amino acids seen at less than the noise floor are treated as sequencing
        # error and excluded, so the total is not inflated by 19 rare alternatives.
        if total:
            alts = {a: n for a, n in alts.items() if n / total >= noise_floor}
        alt_total_n = sum(alts.values())
        alt_total_f = (alt_total_n / total) if total else 0.0
        alt_aa, alt_n, alt_f = "", 0, 0.0
        if alts and total:
            alt_aa = max(alts, key=alts.get)
            alt_n = alts[alt_aa]
            alt_f = alt_n / total
        deep_enough = total >= min_depth
        enough_reads = alt_total_n >= min_alt_reads
        nonsyn = bool(alts) and deep_enough and enough_reads and alt_total_f >= min_af
        # present but below the calling threshold — shown faintly rather than hidden
        nonsyn_sub = (bool(alts) and not nonsyn and deep_enough and enough_reads
                      and alt_total_f >= max(min_af * sub_frac, sub_floor))
        # synonymous variation: codon differs from the reference codon but the AA does not
        syn_n = sum(n for cod, n in counts.items()
                    if cod not in (ref_cod, "del")
                    and CODON_TABLE.get(cod, "X") == aa_ref)
        prof[aa_pos] = dict(
            aa_pos=aa_pos, positions=ps, ref_codon=ref_cod, aa_ref=aa_ref,
            total=total, aa_counts=aa_counts, alt_aa=alt_aa, alt_n=alt_n,
            alt_freq=alt_f, alt_total_n=alt_total_n, alt_total_freq=alt_total_f,
            n_alt_aa=len(alts), nonsyn=nonsyn, nonsyn_sub=nonsyn_sub,
            del_freq=del_freq,
            has_del=(del_freq >= min_af and total >= min_depth
                     and counts.get("del", 0) >= min_alt_reads),
            syn_freq=(syn_n / total if total else 0.0),
            codon_counts=counts,
        )
    return prof


def _suggest_depth(median):
    """A round threshold comfortably below the observed depth."""
    for v in (500, 200, 100, 50, 30, 20, 10, 5):
        if v <= median * 0.8:
            return v
    return max(1, int(median * 0.8))


def warn_if_thresholds_exceed_depth(counts_by_label, profs_by_label, args):
    """Loudly flag the case where the filters are above the data's actual depth.

    Without this the figure comes back with nothing marked and no explanation,
    which looks like 'the variant is missing' rather than 'the filter removed it'.
    """
    aa_depth = args.min_aa_depth if args.min_aa_depth is not None else args.min_depth
    msgs = []
    for label, c in counts_by_label.items():
        dp = [int(x) for x in c.depth if x > 0]
        med_nt = int(np.median(dp)) if dp else 0
        prof = (profs_by_label or {}).get(label) or {}
        frags = [p["total"] for p in prof.values() if p["total"] > 0]
        med_aa = int(np.median(frags)) if frags else None

        nt_bad = med_nt and med_nt < args.min_depth
        aa_bad = med_aa is not None and med_aa < aa_depth
        if not (nt_bad or aa_bad):
            continue
        m = [f"[warning] {label}: the depth filters are ABOVE this sample's actual "
             f"depth, so nothing can be called."]
        m.append(f"          median per-base depth      : {med_nt}"
                 f"{'   <-- below --min-depth ' + str(args.min_depth) if nt_bad else ''}")
        if med_aa is not None:
            m.append(f"          median fragments per codon : {med_aa}"
                     f"{'   <-- below --min-aa-depth ' + str(aa_depth) if aa_bad else ''}")
        sug_nt = _suggest_depth(med_nt) if med_nt else args.min_depth
        sug_aa = _suggest_depth(med_aa) if med_aa else sug_nt
        m.append(f"          try:  --min-depth {sug_nt} --min-aa-depth {sug_aa}")
        m.append(f"          At this depth a 10% variant is only ~{max(1, round(0.1 * (med_aa or med_nt)))} "
                 f"reads; check ci_low/ci_high in the TSV before trusting a frequency.")
        msgs.append("\n".join(m))
    if msgs:
        print("\n" + "\n\n".join(msgs) + "\n", file=sys.stderr)


def explain_codons(prof, counts_by_label, args, chrom, wanted):
    """Say in plain language why each named codon was or was not called."""
    aa_depth = args.min_aa_depth if args.min_aa_depth is not None else args.min_depth
    aa_af = args.min_aa_af if args.min_aa_af is not None else args.min_af
    out = []
    for label, prof_l in prof.items():
        c = counts_by_label.get(label)
        for aa_pos in wanted:
            p = prof_l.get(aa_pos)
            out.append("")
            head = f"[explain] {label}  codon {aa_pos}"
            if p is None:
                out.append(head)
                out.append("  not in the plotted window — widen --region/--aa, "
                           "or check --orf-start is the ORF's first base")
                continue
            ps = p["positions"]
            out.append(f"{head}  (nt {ps[0]}-{ps[2]} on {chrom})")
            out.append(f"  reference codon         : {p['ref_codon']} -> {p['aa_ref']}")
            dist = "; ".join(f"{a}:{n}" for a, n in
                             sorted(p["aa_counts"].items(), key=lambda kv: -kv[1]))
            out.append(f"  amino acids observed    : {dist or '(none)'}")
            out.append(f"  fragments spanning codon: {p['total']}  (read pairs counted once)")
            if c is not None:
                dps = []
                for q in ps:
                    idx = q - int(c.positions[0])
                    if 0 <= idx < len(c.depth):
                        dps.append(int(c.depth[idx]))
                if dps:
                    out.append(f"  per-base depth here     : {min(dps)}-{max(dps)}")
                    if p["total"] < 0.5 * max(dps):
                        out.append("    ^ NOTE codon-level depth is less than half the "
                                   "per-base depth. A fragment must cover all three bases "
                                   "at BQ >= --min-bq; short reads, an amplicon boundary "
                                   "or the window edge cause this. Consider a lower "
                                   "--min-aa-depth.")
            out.append(f"  non-reference AA total  : {p['alt_total_n']}/{p['total']} "
                       f"= {p['alt_total_freq'] * 100:.3f}%"
                       f"   (top: {p['alt_aa'] or 'none'} {p['alt_freq'] * 100:.3f}%)")
            if p["syn_freq"] > 0:
                out.append(f"  synonymous-only         : {p['syn_freq'] * 100:.3f}%")
            if p["del_freq"] > 0:
                out.append(f"  deletion                : {p['del_freq'] * 100:.3f}%")
            ok_d = p["total"] >= aa_depth
            ok_f = p["alt_total_freq"] >= aa_af
            out.append(f"  test depth >= {aa_depth:<10}: "
                       f"{'PASS' if ok_d else 'FAIL'}")
            out.append(f"  test AA fraction >= {aa_af * 100:g}%    : "
                       f"{'PASS' if ok_f else 'FAIL'}")
            if p["nonsyn"]:
                out.append("  => called NON-SYNONYMOUS (drawn red)")
            elif not p["alt_total_n"]:
                if p["del_freq"] >= aa_af:
                    out.append("  => NOT a substitution: the variation here is a "
                               "DELETION (drawn hatched, not red)")
                elif p["syn_freq"] >= aa_af:
                    out.append("  => NOT called: the variation here is SYNONYMOUS "
                               "(drawn grey, not red)")
                else:
                    out.append("  => no non-reference amino acid above the noise "
                               f"floor (--aa-noise-floor {args.aa_noise_floor:g})")
            elif not ok_d:
                out.append(f"  => NOT called: only {p['total']} fragments span the codon; "
                           f"raise coverage or lower --min-aa-depth (now {aa_depth})")
            else:
                out.append(f"  => NOT called: {p['alt_total_freq'] * 100:.3f}% is below "
                           f"--min-aa-af ({aa_af * 100:g}%). Lower it, e.g. "
                           f"--min-aa-af {max(0.001, p['alt_total_freq'] * 0.9):.3f}")
    print("\n".join(out), file=sys.stderr)


def site_table(c: Counts, min_depth: int, min_af: float, max_sb_p: float,
               ref_mode: str = "ref", min_alt_reads: int = 2):
    """Return list of dicts for positions passing the polymorphism filters."""
    rows = []
    n = len(c.positions)
    for i in range(n):
        dp = int(c.depth[i])
        if dp < min_depth:
            continue
        counts = {b: int(c.base[b][i]) for b in BASES}
        counts["del"] = int(c.dele[i])
        ref_b = c.ref[i]
        order = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)
        major, major_n = order[0]
        minor, minor_n = order[1]
        if minor_n == 0:
            continue
        if ref_mode == "ref" and ref_b in counts:
            alt_n = dp - counts[ref_b]
            alt_best = max(((k, v) for k, v in counts.items() if k != ref_b),
                           key=lambda kv: kv[1])
        else:
            alt_n = dp - major_n
            alt_best = (minor, minor_n)
        af = alt_n / dp if dp else 0.0
        maf = minor_n / dp if dp else 0.0
        metric = af if ref_mode == "ref" else maf
        metric_n = alt_n if ref_mode == "ref" else minor_n
        if metric < min_af or metric_n < min_alt_reads:
            continue
        alt_b = alt_best[0]
        if alt_b in BASES:
            a_f, a_r = int(c.fwd[alt_b][i]), int(c.rev[alt_b][i])
        else:
            a_f = a_r = 0
        if ref_b in BASES:
            r_f, r_r = int(c.fwd[ref_b][i]), int(c.rev[ref_b][i])
        else:
            r_f = r_r = 0
        sb = strand_bias_p(a_f, a_r, r_f, r_r) if (a_f + a_r) else 1.0
        if max_sb_p is not None and (a_f + a_r) > 0 and sb < max_sb_p:
            continue
        lo, hi = wilson_ci(int(round(metric * dp)), dp)
        rows.append(
            dict(
                chrom=c.chrom, pos=int(c.positions[i]), ref=ref_b,
                depth=dp, major=major, minor=minor, alt=alt_b,
                alt_count=alt_n, af=metric, ci_low=lo, ci_high=hi,
                A=counts["A"], C=counts["C"], G=counts["G"], T=counts["T"],
                deletion=counts["del"], insertion=int(c.ins[i]),
                alt_fwd=a_f, alt_rev=a_r, ref_fwd=r_f, ref_rev=r_r,
                strand_bias_p=sb,
            )
        )
    return rows


def annotate_aa(rows, ref_fa, chrom, orf_start, orf_end=None):
    """Add codon / amino-acid change columns for a + strand ORF (1-based orf_start)."""
    if ref_fa is None or orf_start is None:
        return rows
    for r in rows:
        pos = r["pos"]
        if pos < orf_start or (orf_end and pos > orf_end):
            r["codon"] = ""
            r["aa_change"] = ""
            r["effect"] = ""
            continue
        off = pos - orf_start
        codon_i = off // 3
        within = off % 3
        cstart = orf_start + codon_i * 3
        codon = ref_fa.fetch(chrom, cstart - 1, cstart + 2).upper()
        if len(codon) < 3 or r["alt"] not in BASES:
            r["codon"] = codon
            r["aa_change"] = ""
            r["effect"] = ""
            continue
        mut = list(codon)
        mut[within] = r["alt"]
        aa_ref = CODON_TABLE.get(codon, "X")
        aa_alt = CODON_TABLE.get("".join(mut), "X")
        r["codon"] = f"{codon}>{''.join(mut)}"
        r["aa_change"] = f"{aa_ref}{codon_i + 1}{aa_alt}"
        r["effect"] = "syn" if aa_ref == aa_alt else "nonsyn"
    return rows


# --------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------
def _fractions(c: Counts):
    dp = np.where(c.depth > 0, c.depth, 1)
    frac = {b: c.base[b] / dp for b in BASES}
    frac["del"] = c.dele / dp
    return frac


def _stack(ax, c: Counts, frac, keys, colors, bar_w=0.82, gap=True):
    bottom = np.zeros(len(c.positions))
    drawn = []
    for k in keys:
        v = frac[k]
        if not np.any(v > 0):
            continue
        ax.bar(c.positions, v, bottom=bottom, width=bar_w,
               color=colors[k], linewidth=0.0, zorder=3)
        bottom = bottom + v
        drawn.append(k)
    return bottom, drawn


def _thin_labels(values, depth, args, npos):
    """Indices to label: strongest first, never within `label_gap` of a kept one."""
    gap = args.label_gap if args.label_gap is not None else max(2, round(npos / 14))
    kept = []
    for i in np.argsort(values)[::-1]:
        if len(kept) >= args.max_labels:
            break
        if values[i] < args.min_af or depth[i] < args.min_depth:
            continue
        if any(abs(int(i) - j) < gap for j in kept):
            continue
        kept.append(int(i))
    return kept


def draw_aa_track(ax, c: Counts, aa_prof, args, compact=False, show_labels=True):
    """A codon ribbon under the x-axis: one-letter AA + AA position, variants called out."""
    LBL = 5.5 if compact else 6.5
    aa_depth = args.min_aa_depth if args.min_aa_depth is not None else args.min_depth
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.patch.set_visible(False)   # let the nucleotide labels above show through
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="both", length=0, labelsize=LBL)
    ax.set_xticks([])
    if not aa_prof:
        return

    # a deletion spanning several codons is one event — annotate it once, not per codon
    ordered = sorted(aa_prof)
    del_group, run = {}, []
    for aa_pos in ordered + [None]:
        p = aa_prof.get(aa_pos) if aa_pos is not None else None
        same_run = (p and p.get("has_del") and run
                    and abs(p["del_freq"] - aa_prof[run[-1]]["del_freq"]) < 0.005
                    and aa_pos == run[-1] + 1)
        if same_run:
            run.append(aa_pos)
            continue
        if run:
            for q in run:
                del_group[q] = (run[0], run[-1])
        run = [aa_pos] if (p and p.get("has_del")) else []

    prev_noted = False
    note_rows = 0
    for k, aa_pos in enumerate(sorted(aa_prof)):
        p = aa_prof[aa_pos]
        lo, hi = p["positions"][0], p["positions"][2]
        nonsyn = p["nonsyn"]
        nonsyn_sub = p.get("nonsyn_sub") and not args.hide_below_threshold
        has_del = p.get("has_del")
        syn = (not nonsyn and not has_del
               and p["syn_freq"] >= args.min_af and p["total"] >= aa_depth)
        style = dict(linewidth=0.9, linestyle="solid")
        if nonsyn:
            face, edge = "#FDE9E7", "#B2182B"
        elif has_del:
            face, edge = "#EDEDED", "#5A6169"
        elif nonsyn_sub:
            face, edge = "#FFFFFF", "#E39A94"
            style = dict(linewidth=0.8, linestyle=(0, (2, 1.5)))
        elif syn:
            face, edge = "#EFEFEF", "#8C8C8C"
        else:
            face, edge = ("#F6F7F8" if k % 2 == 0 else "#FFFFFF"), "#DCDFE3"
            style = dict(linewidth=0.5, linestyle="solid")
        ax.add_patch(
            plt.Rectangle((lo - 0.5, 0.30), (hi - lo) + 1.0, 0.66,
                          facecolor=face, edgecolor=edge,
                          hatch="///" if has_del else None, zorder=2, **style)
        )
        if not show_labels:
            continue
        centre = (lo + hi) / 2.0
        ax.text(centre, 0.63, f"{p['aa_ref']}{aa_pos}", ha="center", va="center",
                fontsize=LBL, color=INK, zorder=4,
                fontweight="bold" if (nonsyn or has_del) else "normal",
                alpha=1.0)
        notes = []
        if nonsyn or nonsyn_sub:
            extra = p.get("n_alt_aa", 1) - 1
            name = f"{p['aa_ref']}{aa_pos}{p['alt_aa']}"
            if extra > 0:
                name += f" +{extra}aa"
            notes.append((name,
                          "#B2182B" if nonsyn else "#C39490", centre))
        if has_del:
            g0, g1 = del_group.get(aa_pos, (aa_pos, aa_pos))
            if aa_pos == g0:      # only the first codon of the run carries the note
                name = (f"{p['aa_ref']}{g0}del" if g0 == g1
                        else f"{aa_prof[g0]['aa_ref']}{g0}-"
                             f"{aa_prof[g1]['aa_ref']}{g1}del")
                gx = (aa_prof[g0]["positions"][0] + aa_prof[g1]["positions"][2]) / 2.0
                notes.append((name, "#5A6169", gx))
        drop = 0.80 if prev_noted else 0.0     # stagger neighbours
        if notes:
            note_rows = max(note_rows, (1 if drop else 0) + len(notes))
        for j, (txt, col, nx) in enumerate(notes):
            ax.annotate(txt, (nx, 0.24 - drop - j * 0.80),
                        xycoords=("data", "data"), ha="center", va="top",
                        fontsize=LBL, color=col, zorder=5, annotation_clip=False)
        prev_noted = bool(notes) and not prev_noted
    if not compact:      # in compact panels the x-label already says "amino acid"
        ax.set_ylabel("amino acid", fontsize=LBL, rotation=0,
                      ha="right", va="center", labelpad=8)
    return note_rows


def draw_pileup_block(fig, gs_row, c: Counts, label, args, highlight=None,
                      aa_map=None, aa_prof=None, show_xlabels=True, compact=False):
    """Draw coverage + 100% composition + alt-only panels (+ AA track) for one sample."""
    from matplotlib.gridspec import GridSpecFromSubplotSpec

    show_aa = (bool(aa_prof) and args.xaxis in ("aa", "both")
               and len(c.positions) <= args.max_base_labels)
    show_nt = args.xaxis in ("nt", "both") or not show_aa

    panels = []
    if not args.no_coverage:
        panels.append(("cov", 0.55))
    if args.layout in ("full", "both"):
        panels.append(("full", 1.35))
    if args.layout in ("alt", "both"):
        panels.append(("alt", 1.15))
    if show_aa:
        panels.append(("aa", 0.42 if compact else 0.5))

    inner = GridSpecFromSubplotSpec(
        len(panels), 1, subplot_spec=gs_row,
        height_ratios=[h for _, h in panels], hspace=0.12
    )
    LBL = 6 if compact else 7
    frac = _fractions(c)
    colors = dict(BASE_COLORS)
    colors["del"] = DEL_COLOR
    ref_seq = c.ref
    x0, x1 = c.positions[0] - 0.6, c.positions[-1] + 0.6
    axes = {}
    aa_note_rows = 0

    for k, (name, _) in enumerate(panels):
        ax = fig.add_subplot(inner[k])
        axes[name] = ax
        ax.set_xlim(x0, x1)
        if name == "aa":
            aa_note_rows = draw_aa_track(ax, c, aa_prof, args, compact=compact)
            continue
        ax.set_axisbelow(True)
        ax.grid(axis="y", color=GRID, linewidth=0.6, zorder=0)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(labelsize=LBL)
        if highlight:
            ax.axvspan(highlight[0] - 0.5, highlight[1] + 0.5,
                       color="#F2C14E", alpha=0.18, zorder=1, linewidth=0)

        if name == "cov":
            ax.fill_between(c.positions, c.depth, step="mid",
                            color=COV_COLOR, zorder=2, linewidth=0)
            ax.plot(c.positions, c.depth, drawstyle="steps-mid",
                    color="#8A9099", linewidth=0.9, zorder=3)
            if args.min_depth > 0:
                ax.axhline(args.min_depth, color=INK_SOFT, linestyle=(0, (3, 2)),
                           linewidth=0.7, zorder=4)
            ax.set_ylabel("depth", fontsize=LBL)
            ax.yaxis.set_major_locator(MaxNLocator(3))
            if args.log_coverage:
                ax.set_yscale("log")

        elif name == "full":
            _stack(ax, c, frac, BASES + ["del"], colors)
            ax.set_ylim(0, 1)
            ax.set_yticks([0, 0.5, 1.0])
            ax.set_yticklabels(["0", "50", "100"])
            ax.set_ylabel("comp. (%)" if compact else "base\ncomposition (%)",
                          fontsize=LBL)

        elif name == "alt":
            altfrac = {}
            for b in BASES:
                v = frac[b].copy()
                for i, rb in enumerate(ref_seq):
                    if rb == b:
                        v[i] = 0.0
                altfrac[b] = v
            altfrac["del"] = frac["del"]
            top, _ = _stack(ax, c, altfrac, BASES + ["del"], colors)
            ymax = max(float(top.max()) * 1.55, args.min_af * 3.0, 0.02)
            ymax = min(ymax, 1.0)
            ax.set_ylim(0, ymax)
            ax.axhline(args.min_af, color=INK_SOFT, linestyle=(0, (3, 2)),
                       linewidth=0.7, zorder=4)
            ax.set_ylabel("non-ref (%)" if compact else "non-reference\nfraction (%)",
                          fontsize=LBL)
            ax.yaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_formatter(
                matplotlib.ticker.FuncFormatter(lambda v, _: f"{v * 100:g}")
            )
            # insertion markers
            insf = np.divide(c.ins, np.where(c.depth > 0, c.depth, 1))
            hits = np.where(insf >= args.min_af)[0]
            for i in hits:
                ax.plot(c.positions[i], ymax * 0.93, marker="v", markersize=4,
                        color=INS_COLOR, zorder=6, clip_on=False)
            # label the strongest sites, thinning collisions
            npos_ = len(c.positions)
            for i in _thin_labels(top, c.depth, args, len(c.positions)):
                txt = f"{c.positions[i]}\n{top[i] * 100:.1f}%"
                if aa_map and c.positions[i] in aa_map and aa_map[c.positions[i]]:
                    txt = f"{c.positions[i]} {aa_map[c.positions[i]]}\n{top[i] * 100:.1f}%"
                # keep the label inside the axes at the window edges
                frac = (c.positions[i] - c.positions[0]) / max(1, npos_ - 1)
                ha = "left" if frac < 0.06 else ("right" if frac > 0.94 else "center")
                dx = 2 if ha == "left" else (-2 if ha == "right" else 0)
                ax.annotate(txt, (c.positions[i], top[i]),
                            textcoords="offset points", xytext=(dx, 4),
                            ha=ha, fontsize=6, color=INK, zorder=7)

        ax.set_xticklabels([])
        ax.set_xlabel("")
        ax.tick_params(axis="x", length=0)

    # the panel that carries the nucleotide labels is the last non-AA panel
    nt_names = [n for n, _ in panels if n != "aa"]
    nt_ax = axes[nt_names[-1]]
    npos = len(c.positions)
    ruler_added = False
    if show_aa:
        # a coordinate ruler along the top of the block; the AA ribbon owns the
        # bottom. In compact panels the base letters are dropped for room.
        if show_nt:
            ruler = axes[panels[0][0]].secondary_xaxis("top")
            ruler.xaxis.set_major_locator(MaxNLocator(6, integer=True))
            ruler.tick_params(labelsize=LBL, length=2, pad=1, colors=INK_SOFT)
            ruler_added = True
            if not compact:
                nt_ax.set_xticks(c.positions)
                nt_ax.set_xticklabels([f"{b}" for b in ref_seq], fontsize=6)
                nt_ax.tick_params(axis="x", length=0, pad=1)
    elif show_xlabels and show_nt:
        if npos <= args.max_base_labels:
            nt_ax.set_xticks(c.positions)
            nt_ax.set_xticklabels([f"{b}" for b in ref_seq], fontsize=6)
            nt_ax.tick_params(axis="x", length=0, pad=1)
            step = max(1, npos // 10)
            for i in range(0, npos, step):
                nt_ax.annotate(
                    str(c.positions[i]), (c.positions[i], 0),
                    xycoords=("data", "axes fraction"), xytext=(0, -16),
                    textcoords="offset points", ha="center", fontsize=6,
                    color=INK_SOFT, rotation=90, annotation_clip=False,
                )
            nt_ax.set_xlabel("reference base" if compact
                             else f"{c.chrom} (reference base)",
                             fontsize=LBL, labelpad=18)
        else:
            nt_ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            nt_ax.tick_params(axis="x", length=3, labelbottom=True, labelsize=LBL)
            nt_ax.set_xlabel(f"position on {c.chrom}", fontsize=LBL)

    if show_aa and show_xlabels:
        axes["aa"].set_xlabel(
            "amino acid" if compact else "amino acid position (one-letter code)",
            fontsize=LBL,
            labelpad=(16 if show_nt else 8) + 15 * max(0, aa_note_rows - 1),
        )

    top_ax = axes[panels[0][0]]
    if ruler_added and "\n" not in label:
        label = f"{label}   ({c.chrom})"
    top_ax.set_title(label, fontsize=8 if compact else 9, loc="left",
                     pad=15 if ruler_added else 4, color=INK)
    return axes


def legend_handles(include_del=True, include_ins=True, include_aa=False):
    h = [Patch(facecolor=BASE_COLORS[b], label=b) for b in BASES]
    if include_del:
        h.append(Patch(facecolor=DEL_COLOR, label="deletion"))
    if include_ins:
        h.append(plt.Line2D([], [], marker="v", linestyle="none",
                            color=INS_COLOR, markersize=5, label="insertion"))
    if include_aa:
        h += [
            Patch(facecolor="#FDE9E7", edgecolor="#B2182B", label="aa: non-synonymous"),
            Patch(facecolor="#FFFFFF", edgecolor="#E39A94", linestyle=(0, (2, 1.5)),
                  label="aa: below threshold"),
            Patch(facecolor="#EFEFEF", edgecolor="#8C8C8C", label="aa: synonymous only"),
            Patch(facecolor="#EDEDED", edgecolor="#5A6169", hatch="///",
                  label="aa: deletion"),
        ]
    return h


def _header(fig, title, aa_key, include_ins=True):
    """Legend and title sit entirely ABOVE the axes, so nothing can collide."""
    fig.legend(handles=legend_handles(include_ins=include_ins, include_aa=aa_key),
               loc="lower right", bbox_to_anchor=(1.0, 1.0),
               ncol=5 if aa_key else 6, frameon=False, fontsize=7)
    if title:
        fig.suptitle(title, fontsize=10, x=0.0, ha="left", y=1.0, va="bottom")


def figure_region(counts_by_sample, args, highlight, aa_map, out_stem, title,
                  aa_prof_by_sample=None):
    from matplotlib.gridspec import GridSpec

    nrow = len(counts_by_sample)
    npos = len(next(iter(counts_by_sample.values())).positions)
    width = max(6.0, min(16.0, 0.14 * npos + 2.4))
    per = 2.5 if args.layout == "both" else 1.9
    if aa_prof_by_sample and args.xaxis in ("aa", "both"):
        per += 0.45
    fig = plt.figure(figsize=(width, per * nrow + 0.9))
    gs = GridSpec(nrow, 1, figure=fig, hspace=0.62)
    for i, (label, c) in enumerate(counts_by_sample.items()):
        draw_pileup_block(fig, gs[i], c, label, args, highlight=highlight,
                          show_xlabels=(i == nrow - 1), aa_map=aa_map,
                          aa_prof=(aa_prof_by_sample or {}).get(label))
    aa_key = bool(aa_prof_by_sample) and args.xaxis in ("aa", "both")
    _header(fig, title, aa_key)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{out_stem}.{ext}", dpi=300)
    plt.close(fig)
    return f"{out_stem}.pdf"


def figure_scan(counts_by_sample, sites_by_sample, args, out_stem, title):
    """Overview: alt fraction vs position across the scanned window + coverage."""
    from matplotlib.gridspec import GridSpec

    nrow = len(counts_by_sample)
    fig = plt.figure(figsize=(9.5, 2.3 * nrow + 0.8))
    gs = GridSpec(nrow, 1, figure=fig, hspace=0.55)
    for i, (label, c) in enumerate(counts_by_sample.items()):
        from matplotlib.gridspec import GridSpecFromSubplotSpec

        inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[i],
                                        height_ratios=[2.2, 0.6], hspace=0.1)
        ax = fig.add_subplot(inner[0])
        axc = fig.add_subplot(inner[1], sharex=ax)
        for a in (ax, axc):
            a.set_axisbelow(True)
            a.grid(axis="y", color=GRID, linewidth=0.6)
            for s in ("top", "right"):
                a.spines[s].set_visible(False)
        rows = sites_by_sample[label]
        if rows:
            for r in rows:
                col = BASE_COLORS.get(r["alt"], DEL_COLOR)
                ax.vlines(r["pos"], 0, r["af"], color=col, linewidth=1.0, zorder=3)
                ax.plot(r["pos"], r["af"], marker="o", markersize=4.2, color=col,
                        markeredgecolor="white", markeredgewidth=0.6, zorder=4)
            span = max(1, c.positions[-1] - c.positions[0])
            mingap = span / 18.0
            kept = []
            for r in sorted(rows, key=lambda r: -r["af"]):
                if len(kept) >= args.max_labels:
                    break
                if any(abs(r["pos"] - k["pos"]) < mingap for k in kept):
                    continue
                kept.append(r)
            for r in kept:
                lab = f"{r['pos']}"
                if r.get("aa_change"):
                    lab += f"\n{r['aa_change']}"
                ax.annotate(lab, (r["pos"], r["af"]), textcoords="offset points",
                            xytext=(0, 5), ha="center", fontsize=6, color=INK)
            ax.set_ylim(0, max(0.06, max(r["af"] for r in rows) * 1.55))
        else:
            ax.set_ylim(0, 0.1)
            ax.text(0.5, 0.5, "no site passed the filters",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=8, color=INK_SOFT)
        ax.axhline(args.min_af, color=INK_SOFT, linestyle=(0, (3, 2)), linewidth=0.7)
        ax.set_ylabel("non-reference\nfraction (%)", fontsize=7)
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: f"{v * 100:g}")
        )
        ax.set_title(f"{label}  ({len(rows)} site{'s' if len(rows) != 1 else ''})",
                     fontsize=9, loc="left", pad=4)
        ax.tick_params(axis="x", labelbottom=False, length=0)

        axc.fill_between(c.positions, c.depth, color=COV_COLOR, linewidth=0)
        axc.plot(c.positions, c.depth, color="#8A9099", linewidth=0.7)
        axc.axhline(args.min_depth, color=INK_SOFT, linestyle=(0, (3, 2)),
                    linewidth=0.7)
        axc.set_ylabel("depth", fontsize=7)
        axc.yaxis.set_major_locator(MaxNLocator(3))
        if args.log_coverage:
            axc.set_yscale("log")
        if i == nrow - 1:
            axc.set_xlabel(f"position on {c.chrom}", fontsize=7)
    _header(fig, title, False, include_ins=False)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{out_stem}.{ext}", dpi=300)
    plt.close(fig)
    return f"{out_stem}.pdf"


def figure_site_panels(sample_paths, ref_fa, site_list, args, out_stem, title):
    """One zoom panel per detected site, samples stacked inside each panel."""
    from matplotlib.gridspec import GridSpec

    if not site_list:
        return None

    # cluster sites that would share the same zoom window into one panel
    merge_dist = args.merge_dist if args.merge_dist is not None else args.flank
    by_pos = sorted(site_list, key=lambda r: r["pos"])
    clusters = []
    for r in by_pos:
        if clusters and r["pos"] - clusters[-1][-1]["pos"] <= merge_dist:
            clusters[-1].append(r)
        else:
            clusters.append([r])
    clusters.sort(key=lambda cl: -max(x["af"] for x in cl))
    clusters = clusters[: args.max_sites]

    ncol = min(args.panel_cols, len(clusters))
    nrow = math.ceil(len(clusters) / ncol)
    nsamp = len(sample_paths)
    per_h = (2.0 if args.layout == "both" else 1.5) * nsamp + 0.5
    if args.orf_start and args.xaxis in ("aa", "both"):
        per_h += 0.4 * nsamp
    fig = plt.figure(figsize=(4.6 * ncol, per_h * nrow))
    gs = GridSpec(nrow, ncol, figure=fig, hspace=0.75, wspace=0.36)

    for k, cl in enumerate(clusters):
        r, cidx = divmod(k, ncol)
        from matplotlib.gridspec import GridSpecFromSubplotSpec

        inner = GridSpecFromSubplotSpec(nsamp, 1, subplot_spec=gs[r, cidx], hspace=0.45)
        lo = max(1, min(x["pos"] for x in cl) - args.flank)
        hi = max(x["pos"] for x in cl) + args.flank
        chrom = cl[0]["chrom"]
        parts = []
        for x in sorted(cl, key=lambda x: -x["af"])[:3]:
            t = f"{x['pos']} {x['ref']}>{x['alt']}"
            if x.get("aa_change"):
                t += f" ({x['aa_change']})"
            parts.append(t)
        head = f"{chrom}: " + "; ".join(parts)
        if len(cl) > 3:
            head += f"; +{len(cl) - 3} more"
        codons = codon_map(lo, hi, args.orf_start, args.orf_end)
        for si, s in enumerate(sample_paths):
            c = count_window(s.path, chrom, lo, hi, ref_fa,
                             args.min_bq, args.min_mq, args.max_depth)
            prof = None
            if codons and args.xaxis in ("aa", "both"):
                cc = count_codon_haplotypes(s.path, chrom, codons,
                                            args.min_bq, args.min_mq)
                prof = aa_profile(codons, cc, ref_fa, chrom,
                                  args.min_aa_depth if args.min_aa_depth is not None
                                  else args.min_depth,
                                  args.min_aa_af if args.min_aa_af is not None
                                  else args.min_af, args.aa_noise_floor,
                                  args.min_alt_reads)
            lab = f"{head}\n{s.label}" if si == 0 else s.label
            draw_pileup_block(fig, inner[si], c, lab, args,
                              highlight=(min(x["pos"] for x in cl),
                                         max(x["pos"] for x in cl)),
                              show_xlabels=(si == nsamp - 1), aa_prof=prof,
                              compact=True)
    aa_key = bool(args.orf_start) and args.xaxis in ("aa", "both")
    _header(fig, title, aa_key)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{out_stem}.{ext}", dpi=300)
    plt.close(fig)
    return f"{out_stem}.pdf"


# --------------------------------------------------------------------------
# TSV writers
# --------------------------------------------------------------------------
def write_counts_tsv(path, counts_by_sample):
    cols = ["sample", "chrom", "pos", "ref", "depth", "A", "C", "G", "T",
            "del", "ins", "alt_fraction", "minor_allele_fraction"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for label, c in counts_by_sample.items():
            for i, p in enumerate(c.positions):
                dp = int(c.depth[i])
                vals = {b: int(c.base[b][i]) for b in BASES}
                vals["del"] = int(c.dele[i])
                rb = c.ref[i]
                altf = (dp - vals.get(rb, 0)) / dp if dp else 0.0
                srt = sorted(vals.values(), reverse=True)
                maf = (srt[1] / dp) if dp and len(srt) > 1 else 0.0
                fh.write("\t".join(map(str, [
                    label, c.chrom, int(p), rb, dp,
                    vals["A"], vals["C"], vals["G"], vals["T"], vals["del"],
                    int(c.ins[i]), f"{altf:.5f}", f"{maf:.5f}",
                ])) + "\n")


def write_codons_tsv(path, profs_by_sample, chrom):
    """Per-codon amino-acid frequencies, counted from reads spanning the whole codon."""
    cols = ["sample", "chrom", "aa_pos", "nt_start", "nt_end", "ref_codon", "aa_ref",
            "fragments_spanning", "alt_aa", "alt_aa_count", "alt_aa_freq",
            "alt_aa_total_count", "alt_aa_total_freq", "n_alt_aa",
            "syn_freq", "del_freq", "nonsyn", "below_threshold",
            "aa_distribution", "codon_distribution"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for label, prof in profs_by_sample.items():
            for aa_pos in sorted(prof):
                p = prof[aa_pos]
                aas = ";".join(f"{a}:{n}" for a, n in
                               sorted(p["aa_counts"].items(), key=lambda kv: -kv[1]))
                cods = ";".join(f"{cd}:{n}" for cd, n in
                                sorted(p["codon_counts"].items(), key=lambda kv: -kv[1]))
                fh.write("\t".join(map(str, [
                    label, chrom, aa_pos, p["positions"][0], p["positions"][2],
                    p["ref_codon"], p["aa_ref"], p["total"], p["alt_aa"], p["alt_n"],
                    f"{p['alt_freq']:.5f}",
                    p["alt_total_n"], f"{p['alt_total_freq']:.5f}", p["n_alt_aa"],
                    f"{p['syn_freq']:.5f}", f"{p['del_freq']:.5f}",
                    "yes" if p["nonsyn"] else "no",
                    "yes" if p.get("nonsyn_sub") else "no", aas, cods,
                ])) + "\n")


def write_sites_tsv(path, sites_by_sample):
    cols = ["sample", "chrom", "pos", "ref", "alt", "depth", "alt_count", "af",
            "ci_low", "ci_high", "A", "C", "G", "T", "deletion", "insertion",
            "alt_fwd", "alt_rev", "ref_fwd", "ref_rev", "strand_bias_p",
            "codon", "aa_change", "effect"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for label, rows in sites_by_sample.items():
            for r in rows:
                fh.write("\t".join(
                    [label] + [str(r.get(c, "")) if not isinstance(r.get(c), float)
                               else f"{r.get(c):.5g}" for c in cols[1:]]
                ) + "\n")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        description="Plot per-base composition (pileup) at intrahost polymorphic sites from BAM.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Examples")[-1],
    )
    p.add_argument("--bam", action="append", required=True, metavar="PATH[=LABEL]",
                   help="mapped, sorted, indexed BAM. Repeat for several samples. "
                        "Append '=Label' to set the display name.")
    p.add_argument("--ref", help="reference FASTA (indexed or indexable). "
                                 "Strongly recommended; without it the majority "
                                 "base is treated as the reference.")
    p.add_argument("--region", help="explicit window to draw, CHR:START-END (1-based, inclusive)")
    p.add_argument("--scan", help="larger window to search for polymorphic sites, CHR:START-END")
    p.add_argument("--highlight", help="span to shade in region mode, e.g. 23640-23655 "
                                       "(the FCS motif)")
    p.add_argument("--aa", metavar="FIRST-LAST",
                   help="give the drawn window as an amino-acid range instead of "
                        "nucleotides, e.g. --aa 780-800 (needs --orf-start and, "
                        "unless --region is also given, --chrom)")
    p.add_argument("--highlight-aa", metavar="FIRST-LAST",
                   help="shade an amino-acid range, e.g. --highlight-aa 790-793 "
                        "for the FCS motif (needs --orf-start)")
    p.add_argument("--chrom", help="contig name, used with --aa when --region is absent")
    p.add_argument("--orf-start", type=int,
                   help="1-based genomic coordinate of the first base of the ORF "
                        "(+ strand) used for codon / amino-acid annotation, e.g. the spike ATG")
    p.add_argument("--orf-end", type=int, help="1-based end of that ORF (optional)")
    p.add_argument("--explain", metavar="AA_POS[,AA_POS...]",
                   help="print a plain-language diagnosis for these amino-acid "
                        "positions: codon-level depth, the amino acids seen, and "
                        "exactly which threshold decided the call. Use this when a "
                        "codon you expected to be red is not.")
    p.add_argument("--out", default="isnv", help="output prefix (default: isnv)")

    f = p.add_argument_group("filters")
    f.add_argument("--min-af", type=float, default=0.01,
                   help="minimum non-reference fraction to call a site (default 0.01)")
    f.add_argument("--min-alt-reads", type=int, default=2, metavar="N",
                   help="a variant must also be supported by at least this many "
                        "reads/fragments (default 2). This is what makes a 1%% "
                        "frequency threshold safe across uneven coverage: at 5000x "
                        "1%% is 50 reads, but at 60x it is well under one, so "
                        "without a read floor a single sequencing error clears the "
                        "threshold. Set 1 to disable.")
    f.add_argument("--min-depth", type=int, default=100,
                   help="minimum depth at a site (default 100)")
    f.add_argument("--min-bq", type=int, default=20, help="minimum base quality (default 20)")
    f.add_argument("--min-mq", type=int, default=20, help="minimum mapping quality (default 20)")
    f.add_argument("--max-depth", type=int, default=1_000_000,
                   help="pysam pileup depth cap (default 1e6 = effectively off)")
    f.add_argument("--max-strand-bias-p", type=float, default=None, metavar="P",
                   help="drop sites whose Fisher strand-bias p is BELOW this "
                        "(e.g. 0.001). Default: no strand filter, p is reported.")
    f.add_argument("--min-aa-af", type=float, default=None, metavar="F",
                   help="amino-acid level threshold: minimum TOTAL fraction of "
                        "reads coding a non-reference amino acid for a codon to "
                        "be called non-synonymous (red). Default: same as --min-af")
    f.add_argument("--min-aa-depth", type=int, default=None, metavar="N",
                   help="minimum number of reads spanning the whole codon. "
                        "Default: same as --min-depth. Note this is always lower "
                        "than the per-base depth, since a fragment must cover "
                        "all three bases at sufficient quality")
    f.add_argument("--aa-noise-floor", type=float, default=0.002, metavar="F",
                   help="amino acids seen below this fraction of codon-spanning "
                        "fragments are treated as sequencing error and excluded from "
                        "the codon's non-reference total (default 0.002)")
    f.add_argument("--hide-below-threshold", action="store_true",
                   help="do not mark codons whose non-reference amino acid is "
                        "present but below the calling threshold (by default "
                        "these get a dashed pink outline so they are not "
                        "silently indistinguishable from invariant codons)")
    f.add_argument("--metric", choices=["alt", "minor"], default="alt",
                   help="'alt' = fraction of reads differing from the reference "
                        "(default); 'minor' = 1 - majority fraction")

    d = p.add_argument_group("drawing")
    d.add_argument("--layout", choices=["both", "full", "alt"], default="both",
                   help="'full' = 100%% stacked composition only, 'alt' = "
                        "non-reference fraction only, 'both' (default)")
    d.add_argument("--no-coverage", action="store_true", help="hide the depth panel")
    d.add_argument("--log-coverage", action="store_true", help="log-scale the depth panel")
    d.add_argument("--flank", type=int, default=12,
                   help="bases either side of a detected site in the zoom panels (default 12)")
    d.add_argument("--max-sites", type=int, default=12,
                   help="max number of zoom panels to draw (default 12)")
    d.add_argument("--merge-dist", type=int, default=None,
                   help="sites closer than this are drawn in one shared panel "
                        "(default: --flank)")
    d.add_argument("--panel-cols", type=int, default=2, help="columns of zoom panels (default 2)")
    d.add_argument("--max-labels", type=int, default=6,
                   help="max number of sites to label per panel (default 6)")
    d.add_argument("--label-gap", type=int, default=None,
                   help="minimum spacing (in positions) between labels in region "
                        "mode; default is auto from the window width")
    d.add_argument("--max-base-labels", type=int, default=120,
                   help="write the reference base under each position when the "
                        "window is at most this wide (default 120)")
    d.add_argument("--xaxis", choices=["both", "nt", "aa"], default="both",
                   help="x-axis annotation: 'nt' = reference base + nucleotide "
                        "position only, 'aa' = amino-acid one-letter code + AA "
                        "position track only, 'both' (default). The AA track "
                        "needs --orf-start.")
    d.add_argument("--font", default=None, metavar="PATH",
                   help="TrueType font file to use for all text (e.g. "
                        "/Users/you/python/fonts/Arial.ttf). If omitted, a few "
                        "usual Arial locations are tried, then matplotlib's "
                        "default. SVG text is kept as text and PDF embeds the "
                        "font, so both stay editable in Illustrator.")
    d.add_argument("--title", default=None, help="figure title")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    setup_font(args.font)
    if not (args.region or args.scan or args.aa):
        raise SystemExit("[error] give at least one of --region, --scan or --aa")

    samples = []
    for spec in args.bam:
        if "=" in spec:
            path, label = spec.rsplit("=", 1)
        else:
            path, label = spec, os.path.basename(spec).replace(".bam", "")
        if not os.path.exists(path):
            raise SystemExit(f"[error] BAM not found: {path}")
        if not (os.path.exists(path + ".bai") or os.path.exists(path[:-4] + ".bai")):
            raise SystemExit(f"[error] BAM index missing — run: samtools index {path}")
        samples.append(Sample(label=label, path=path))

    ref_fa = None
    if args.ref:
        if not os.path.exists(args.ref + ".fai"):
            pysam.faidx(args.ref)
        ref_fa = pysam.FastaFile(args.ref)

    # amino-acid convenience: convert AA ranges to nucleotide coordinates
    if args.aa or args.highlight_aa:
        if not args.orf_start:
            raise SystemExit("[error] --aa / --highlight-aa need --orf-start")
    if args.aa:
        a1, a2 = parse_span(args.aa)
        nt1 = args.orf_start + (a1 - 1) * 3
        nt2 = args.orf_start + a2 * 3 - 1
        chrom = args.chrom or (parse_region(args.region)[0] if args.region
                               else (parse_region(args.scan)[0] if args.scan else None))
        if chrom is None:
            raise SystemExit("[error] --aa needs --chrom (or --region / --scan)")
        args.region = f"{chrom}:{nt1}-{nt2}"
        print(f"[info] --aa {args.aa} -> --region {args.region}", file=sys.stderr)
    if args.highlight_aa:
        h1, h2 = parse_span(args.highlight_aa)
        args.highlight = (f"{args.orf_start + (h1 - 1) * 3}-"
                          f"{args.orf_start + h2 * 3 - 1}")
        print(f"[info] --highlight-aa {args.highlight_aa} -> --highlight {args.highlight}",
              file=sys.stderr)

    highlight = parse_span(args.highlight) if args.highlight else None
    written = []

    # ---------------- region mode ----------------
    if args.region:
        chrom, start, end = parse_region(args.region)
        counts, profs = {}, {}
        aa_map = {}
        codons = codon_map(start, end, args.orf_start, args.orf_end)
        for s in samples:
            c = count_window(s.path, chrom, start, end, ref_fa,
                             args.min_bq, args.min_mq, args.max_depth)
            counts[s.label] = c
            rows = site_table(c, args.min_depth, args.min_af,
                              args.max_strand_bias_p, args.metric,
                              args.min_alt_reads)
            rows = annotate_aa(rows, ref_fa, chrom, args.orf_start, args.orf_end)
            for r in rows:
                if r.get("aa_change"):
                    aa_map[r["pos"]] = r["aa_change"]
            if codons and args.xaxis in ("aa", "both"):
                cc = count_codon_haplotypes(s.path, chrom, codons,
                                            args.min_bq, args.min_mq)
                profs[s.label] = aa_profile(codons, cc, ref_fa, chrom,
                                            args.min_aa_depth if args.min_aa_depth
                                            is not None else args.min_depth,
                                            args.min_aa_af if args.min_aa_af
                                            is not None else args.min_af,
                                            args.aa_noise_floor,
                                            args.min_alt_reads)
        write_counts_tsv(f"{args.out}.region.counts.tsv", counts)
        written.append(f"{args.out}.region.counts.tsv")
        if profs:
            write_codons_tsv(f"{args.out}.region.codons.tsv", profs, chrom)
            written.append(f"{args.out}.region.codons.tsv")
        warn_if_thresholds_exceed_depth(counts, profs, args)
        if args.explain and profs:
            explain_codons(profs, counts, args, chrom,
                           [int(x) for x in args.explain.replace(" ", "").split(",") if x])
        written.append(
            figure_region(counts, args, highlight, aa_map, f"{args.out}.region",
                          args.title or f"{chrom}:{start:,}-{end:,}",
                          aa_prof_by_sample=profs)
        )

    # ---------------- scan mode ----------------
    if args.scan:
        chrom, start, end = parse_region(args.scan)
        counts, sites = {}, {}
        for s in samples:
            c = count_window(s.path, chrom, start, end, ref_fa,
                             args.min_bq, args.min_mq, args.max_depth)
            counts[s.label] = c
            rows = site_table(c, args.min_depth, args.min_af,
                              args.max_strand_bias_p, args.metric,
                              args.min_alt_reads)
            rows = annotate_aa(rows, ref_fa, chrom, args.orf_start, args.orf_end)
            sites[s.label] = rows
            print(f"[scan] {s.label}: {len(rows)} polymorphic site(s) "
                  f"(depth>={args.min_depth}, AF>={args.min_af})", file=sys.stderr)
        warn_if_thresholds_exceed_depth(counts, None, args)
        write_counts_tsv(f"{args.out}.scan.counts.tsv", counts)
        write_sites_tsv(f"{args.out}.sites.tsv", sites)
        written += [f"{args.out}.scan.counts.tsv", f"{args.out}.sites.tsv"]
        written.append(
            figure_scan(counts, sites, args, f"{args.out}.scan",
                        args.title or f"{chrom}:{start:,}-{end:,}")
        )
        # union of sites across samples, strongest first
        merged, seen = [], set()
        for rows in sites.values():
            for r in rows:
                if r["pos"] not in seen:
                    seen.add(r["pos"])
                    merged.append(r)
        merged.sort(key=lambda r: -r["af"])
        if args.orf_start and merged:
            windows = set()
            for r in merged[: args.max_sites * 4]:
                for cd in codon_map(r["pos"] - args.flank, r["pos"] + args.flank,
                                    args.orf_start, args.orf_end):
                    windows.add(cd)
            site_codons = sorted(windows, key=lambda t: t[0])
            profs = {}
            for s in samples:
                cc = count_codon_haplotypes(s.path, chrom, site_codons,
                                            args.min_bq, args.min_mq)
                profs[s.label] = aa_profile(site_codons, cc, ref_fa, chrom,
                                            args.min_aa_depth if args.min_aa_depth
                                            is not None else args.min_depth,
                                            args.min_aa_af if args.min_aa_af
                                            is not None else args.min_af,
                                            args.aa_noise_floor,
                                            args.min_alt_reads)
            write_codons_tsv(f"{args.out}.sites.codons.tsv", profs, chrom)
            written.append(f"{args.out}.sites.codons.tsv")
        out = figure_site_panels(samples, ref_fa, merged, args,
                                 f"{args.out}.sites", args.title)
        if out:
            written.append(out)

    print("\n".join("[written] " + w for w in written if w), file=sys.stderr)


if __name__ == "__main__":
    main()

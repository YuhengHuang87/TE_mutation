#!/usr/bin/env python3
"""
MA_mutation_TE.py
=================

Python port of three Perl scripts that analyze de novo mutations around TEs in
the MA-line WGS Strelka calls:

  1. MAlines_callable_sites_joint_consider_each_TEs.pl
       -> per-sample callable totals + per-TE callable summed across samples,
          using ancestor + per-sample Strelka gVCFs (genome.S{N}.vcf).
  2. MAlines_mutation_joint_consider_sample_ID.pl
       -> per-sample-unique de novo SNPs from the joint variants.vcf in TE
          flanking windows, requiring >=5 supporting alt reads.
  3. Summary
       -> per-TE de novo mutation events (mutations within 1 kb in the same
          sample collapse into one event) and per-TE callable sums.

Inputs (all required):

  --group                   Act / Control / Ubi
  --reference               A4 / mc_initial
  --joint-vcf               .../strelka/{group}_{ref}/results/variants/variants.vcf
  --gvcf-dir                .../strelka/{group}_{ref}/results/variants
                            (contains genome.S1.vcf .. genome.S{N+1}.vcf)
  --repeatmasker            RepeatMasker .out file (used to build the repeat mask)
  --te-annotation           TE annotation .txt (cols 1-3 = chr, left, right)
  --euchromatin-bed         BED restricting analysis to the same euchromatin
                            regions Strelka was given as --callRegions
  --num-ma-samples          number of MA samples (e.g. Act=48, Control=47, Ubi=52)
  --output-prefix           prefix for output files

Optional:

  --flank-outer 5000   --flank-inner 500
  --collapse-distance 1000
  --min-alt-count 5
  --chroms 2L 2R 3L 3R X        # filtered against (allows trailing _RagTag)

The ancestor is assumed to occupy column 9 of the joint VCF (i.e. the first
sample column) and to correspond to genome.S1.vcf. MA samples occupy columns
10..(9+N-1) and correspond to genome.S2.vcf .. genome.S{N+1}.vcf, matching the
order in which BAMs were passed to configureStrelkaGermlineWorkflow.py.

Sample names are read from the joint VCF #CHROM header.

Outputs (written next to --output-prefix):

  {prefix}.per_sample_callable.tsv
        sample_id  total_callable  te_flanking_callable
        n_mutations_total_callable  n_mutations_te_flanking
  {prefix}.de_novo_mutations.tsv
        sample_id  chrom  pos  ref  alt  anc_gt  sample_gt  alt_count  ref_count
        in_te_flanking
  {prefix}.per_te_summary.tsv
        te_id  chrom  left  right
        n_events_collapsed_1kb  n_samples_with_events
        callable_summed_across_samples
"""

from __future__ import annotations

import argparse
import csv
import sys
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

# ---------------------------------------------------------------------------
# Sorted-interval utilities. Each interval list is sorted by start, with
# [start, end] inclusive on both ends and non-overlapping (after merge).
# ---------------------------------------------------------------------------

Interval = Tuple[int, int]


def merge_intervals(ivs: List[Interval]) -> List[Interval]:
    if not ivs:
        return []
    ivs = sorted(ivs)
    merged = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= merged[-1][1] + 1:
            if e > merged[-1][1]:
                merged[-1][1] = e
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def intersect_intervals(a: List[Interval], b: List[Interval]) -> List[Interval]:
    out: List[Interval] = []
    i = j = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0])
        e = min(a[i][1], b[j][1])
        if s <= e:
            out.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def subtract_intervals(a: List[Interval], b: List[Interval]) -> List[Interval]:
    """a minus b. Both sorted, merged, non-overlapping."""
    out: List[Interval] = []
    j = 0
    for s, e in a:
        cur = s
        while j < len(b) and b[j][1] < cur:
            j += 1
        k = j
        while k < len(b) and b[k][0] <= e:
            bs, be = b[k]
            if bs > cur:
                out.append((cur, min(bs - 1, e)))
            cur = max(cur, be + 1)
            if cur > e:
                break
            k += 1
        if cur <= e:
            out.append((cur, e))
    return out


def total_length(ivs: List[Interval]) -> int:
    return sum(e - s + 1 for s, e in ivs)


def position_in_intervals(p: int, ivs: List[Interval]) -> bool:
    if not ivs:
        return False
    starts = [s for s, _ in ivs]
    idx = bisect_right(starts, p) - 1
    return idx >= 0 and ivs[idx][1] >= p


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def load_repeatmasker(path: str, chroms: set) -> Dict[str, List[Interval]]:
    """RepeatMasker .out: skip the 3-line header; columns are whitespace-split.
    Use cols 5,6,7 (1-based) = chrom, start, end."""
    by_chrom: Dict[str, List[Interval]] = defaultdict(list)
    with open(path) as f:
        for i, line in enumerate(f):
            if i < 3:
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            chrom = parts[4]
            if chrom not in chroms:
                continue
            try:
                s = int(parts[5])
                e = int(parts[6])
            except ValueError:
                continue
            by_chrom[chrom].append((s, e))
    return {c: merge_intervals(v) for c, v in by_chrom.items()}


def load_te_annotation(path: str, chroms: set) -> Dict[str, List[Tuple[int, int, str]]]:
    """TE annotation TSV: cols 0,1,2 = chrom, left, right. te_id = chrom:left-right.
    Returns per-chrom list sorted by left, with all 3 fields."""
    by_chrom: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            chrom = parts[0]
            if chrom not in chroms:
                continue
            try:
                left = int(parts[1])
                right = int(parts[2])
            except ValueError:
                continue
            te_id = f"{chrom}:{left}-{right}"
            by_chrom[chrom].append((left, right, te_id))
    for c in by_chrom:
        by_chrom[c].sort()
    return dict(by_chrom)


def load_bed(path: str, chroms: set) -> Dict[str, List[Interval]]:
    """Plain BED -> per-chrom sorted+merged intervals (BED is half-open;
    we convert to inclusive [start+1, end])."""
    by_chrom: Dict[str, List[Interval]] = defaultdict(list)
    opener = open
    if path.endswith('.gz'):
        import gzip
        opener = gzip.open
    with opener(path, 'rt') as f:
        for line in f:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            chrom = parts[0]
            if chrom not in chroms:
                continue
            s = int(parts[1]) + 1  # BED is 0-based half-open
            e = int(parts[2])
            if e >= s:
                by_chrom[chrom].append((s, e))
    return {c: merge_intervals(v) for c, v in by_chrom.items()}


def build_te_flanking(
    tes_by_chrom: Dict[str, List[Tuple[int, int, str]]],
    flank_outer: int,
    flank_inner: int,
) -> Dict[str, List[Tuple[int, int, str]]]:
    """For each TE, flanking = [left-outer, left-inner-1] U [right+inner+1, right+outer].
    Returned as per-chrom sorted list of (start, end, te_id). May overlap across TEs."""
    out: Dict[str, List[Tuple[int, int, str]]] = {}
    for chrom, tes in tes_by_chrom.items():
        ivs: List[Tuple[int, int, str]] = []
        for left, right, te_id in tes:
            ls = max(1, left - flank_outer)
            le = left - flank_inner - 1
            if le >= ls:
                ivs.append((ls, le, te_id))
            rs = right + flank_inner + 1
            re_ = right + flank_outer
            if re_ >= rs:
                ivs.append((rs, re_, te_id))
        ivs.sort()
        out[chrom] = ivs
    return out


# ---------------------------------------------------------------------------
# gVCF parsing -> per-chrom PASS intervals
# ---------------------------------------------------------------------------

def parse_gvcf_pass_intervals(path: str, chroms: set) -> Dict[str, List[Interval]]:
    by_chrom: Dict[str, List[Interval]] = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not line or line[0] == '#':
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 10:
                continue
            chrom = parts[0]
            if chrom not in chroms:
                continue
            if parts[6] != 'PASS':
                continue
            sample_field = parts[9]
            if '|' in sample_field:  # phased; skip per the original Perl
                continue
            try:
                start = int(parts[1])
            except ValueError:
                continue
            end = start
            info = parts[7]
            if 'END=' in info:
                for kv in info.split(';'):
                    if kv.startswith('END='):
                        try:
                            end = int(kv[4:])
                        except ValueError:
                            end = start
                        break
            by_chrom[chrom].append((start, end))
    return {c: merge_intervals(v) for c, v in by_chrom.items()}


# ---------------------------------------------------------------------------
# Per-TE callable accumulation
# ---------------------------------------------------------------------------

def accumulate_te_callable(
    callable_by_chrom: Dict[str, List[Interval]],
    te_flanking: Dict[str, List[Tuple[int, int, str]]],
) -> Tuple[Dict[str, int], int]:
    """Return ({te_id: bp_callable_in_flanking}, total_te_flanking_callable)
    for the given sample's callable intervals."""
    per_te: Dict[str, int] = defaultdict(int)
    total = 0
    for chrom, ivs in callable_by_chrom.items():
        flanks = te_flanking.get(chrom, [])
        if not flanks or not ivs:
            continue
        flank_starts = [s for s, _, _ in flanks]
        for cs, ce in ivs:
            # candidate flanking intervals overlap [cs, ce]
            i = bisect_right(flank_starts, ce) - 1
            j = bisect_left(flank_starts, cs - 1)  # earliest possibly overlapping
            j = max(0, min(j, len(flanks) - 1))
            # Walk back while preceding interval might extend into [cs, ce]
            while j > 0 and flanks[j - 1][1] >= cs:
                j -= 1
            for k in range(j, i + 1):
                fs, fe, te_id = flanks[k]
                if fe < cs:
                    continue
                if fs > ce:
                    break
                ov = min(ce, fe) - max(cs, fs) + 1
                if ov > 0:
                    per_te[te_id] += ov
                    total += ov
    return dict(per_te), total


# ---------------------------------------------------------------------------
# Joint VCF parsing for de novo SNPs
# ---------------------------------------------------------------------------

def is_snv(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and ref in 'ACGT' and alt in 'ACGT'


def call_de_novo(
    joint_vcf: str,
    chroms: set,
    repeats_by_chrom: Dict[str, List[Interval]],
    euchromatin_by_chrom: Dict[str, List[Interval]],
    ancestor_col: int,
    num_ma_samples: int,
    min_alt_count: int,
    sample_names_override: List[str] | None = None,
):
    """Yield dicts of de novo mutation rows."""
    sample_names: List[str] = []
    with open(joint_vcf) as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                cols = line.rstrip('\n').split('\t')
                sample_names = cols[9:]
                break

        if not sample_names:
            raise RuntimeError(f"no #CHROM header in {joint_vcf}")
        if len(sample_names) < 1 + num_ma_samples:
            raise RuntimeError(
                f"{joint_vcf} has {len(sample_names)} sample columns; "
                f"expected >= {1 + num_ma_samples}"
            )

        if sample_names_override is not None:
            sample_names = list(sample_names_override) + sample_names[len(sample_names_override):]

        ma_cols = list(range(ancestor_col + 1, ancestor_col + 1 + num_ma_samples))

        for line in f:
            if not line or line[0] == '#':
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= ma_cols[-1]:
                continue
            chrom = parts[0]
            if chrom not in chroms:
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue

            if not is_snv(parts[3], parts[4]):
                continue

            # restrict to euchromatin
            if not position_in_intervals(pos, euchromatin_by_chrom.get(chrom, [])):
                continue
            # exclude repeat-masked positions
            if position_in_intervals(pos, repeats_by_chrom.get(chrom, [])):
                continue

            anc = parts[ancestor_col]
            if 'PASS' not in anc:
                continue
            anc_fields = anc.split(':')
            anc_gt = anc_fields[0]
            if anc_gt not in ('0/0', '1/1'):
                continue

            # walk MA samples; require exactly one differing call
            hits = []
            for col in ma_cols:
                samp = parts[col]
                samp_fields = samp.split(':')
                samp_gt = samp_fields[0]
                if anc_gt == '0/0' and samp_gt in ('1/1', '0/1'):
                    hits.append((col, samp, samp_gt, samp_fields))
                elif anc_gt == '1/1' and samp_gt in ('0/0', '0/1'):
                    hits.append((col, samp, samp_gt, samp_fields))
                if len(hits) > 1:
                    break

            if len(hits) != 1:
                continue
            col, samp, samp_gt, samp_fields = hits[0]
            if 'PASS' not in samp:
                continue

            # AD index in FORMAT for Strelka germline is 5
            # (GT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL)
            if len(samp_fields) <= 5:
                continue
            ad_parts = samp_fields[5].split(',')
            if len(ad_parts) < 2:
                continue
            try:
                ref_count = int(ad_parts[0])
                alt_count = int(ad_parts[1])
            except ValueError:
                continue

            ok = (anc_gt == '0/0' and alt_count >= min_alt_count) or \
                 (anc_gt == '1/1' and ref_count >= min_alt_count)
            if not ok:
                continue

            yield {
                'sample_idx': col,                  # 0-based VCF column index
                'sample_name': sample_names[col - 9],
                'chrom': chrom,
                'pos': pos,
                'ref': parts[3],
                'alt': parts[4],
                'anc_gt': anc_gt,
                'sample_gt': samp_gt,
                'ref_count': ref_count,
                'alt_count': alt_count,
            }


# ---------------------------------------------------------------------------
# Per-TE summary: collapse mutations within --collapse-distance per sample,
# attribute each event to a TE flanking window.
# ---------------------------------------------------------------------------

def assign_to_te(
    chrom: str,
    pos: int,
    te_flanking: Dict[str, List[Tuple[int, int, str]]],
    tes_by_chrom: Dict[str, List[Tuple[int, int, str]]],
) -> str | None:
    """Assign a position to the TE whose flanking interval contains it.
    On overlap, pick the TE whose body is closest to pos (ties: smaller left)."""
    flanks = te_flanking.get(chrom, [])
    if not flanks:
        return None
    starts = [s for s, _, _ in flanks]
    # candidate flanks contain pos: their start <= pos
    idx = bisect_right(starts, pos) - 1
    if idx < 0:
        return None
    candidates: List[str] = []
    k = idx
    while k >= 0 and flanks[k][0] >= pos - 5_000_000:  # safety bound
        s, e, tid = flanks[k]
        if e < pos:
            k -= 1
            continue
        if s <= pos <= e:
            candidates.append(tid)
        k -= 1
        # continue walking back: an earlier (smaller start) interval may still
        # contain pos if its end >= pos
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    # tie-break: closest TE body
    body_by_id = {tid: (l, r) for l, r, tid in tes_by_chrom.get(chrom, [])}
    def dist(tid: str) -> Tuple[int, int]:
        l, r = body_by_id[tid]
        if pos < l:
            return (l - pos, l)
        if pos > r:
            return (pos - r, l)
        return (0, l)
    candidates.sort(key=dist)
    return candidates[0]


def collapse_events(
    mutations: List[dict],
    distance: int,
) -> List[dict]:
    """Within each sample, sort mutations by (chrom, pos) and collapse
    consecutive ones whose gap (per chrom) is <= distance into one event.
    The event is represented by the first mutation of the cluster."""
    by_sample: Dict[str, List[dict]] = defaultdict(list)
    for m in mutations:
        by_sample[m['sample_name']].append(m)

    out: List[dict] = []
    for sample, muts in by_sample.items():
        muts.sort(key=lambda m: (m['chrom'], m['pos']))
        last_chrom = None
        last_pos = None
        for m in muts:
            if last_chrom == m['chrom'] and last_pos is not None and m['pos'] - last_pos <= distance:
                # merged into previous event; advance the trailing edge so
                # subsequent positions chain off the most recent SNP, not the
                # cluster head.
                last_pos = m['pos']
                continue
            out.append(m)
            last_chrom = m['chrom']
            last_pos = m['pos']
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--group', required=True, choices=['Act', 'Control', 'Ubi'])
    p.add_argument('--reference', required=True)
    p.add_argument('--joint-vcf', required=True)
    p.add_argument('--gvcf-dir', required=True)
    p.add_argument('--repeatmasker', required=True)
    p.add_argument('--te-annotation', required=True)
    p.add_argument('--euchromatin-bed', required=True)
    p.add_argument('--num-ma-samples', type=int, required=True)
    p.add_argument('--output-prefix', required=True)
    p.add_argument('--ancestor-col', type=int, default=9,
                   help='0-based column index of the ancestor in the joint VCF (default: 9)')
    p.add_argument('--sample-names', nargs='+',
                   help='Override sample column names from the joint VCF #CHROM header. '
                        'Pass 1+N names in BAM order: ancestor first, then MA samples. '
                        'Use this when BAMs were RG-tagged with a non-unique SM '
                        '(e.g. all "Sample1") so the VCF header lost the per-sample '
                        'identity. The column-index <-> position mapping is unambiguous '
                        '(column 9 = BAM 1, column 10 = BAM 2, ...).')
    p.add_argument('--flank-outer', type=int, default=5000)
    p.add_argument('--flank-inner', type=int, default=500)
    p.add_argument('--collapse-distance', type=int, default=1000)
    p.add_argument('--min-alt-count', type=int, default=5)
    p.add_argument('--chroms', nargs='+',
                   default=['2L', '2R', '3L', '3R', 'X',
                            '2L_RagTag', '2R_RagTag', '3L_RagTag', '3R_RagTag', 'X_RagTag'])
    return p.parse_args()


def main() -> int:
    args = parse_args()
    chroms = set(args.chroms)
    log = lambda msg: print(f"[MA_mutation_TE] {msg}", file=sys.stderr, flush=True)

    log(f"loading repeatmasker: {args.repeatmasker}")
    repeats = load_repeatmasker(args.repeatmasker, chroms)

    log(f"loading TE annotation: {args.te_annotation}")
    tes = load_te_annotation(args.te_annotation, chroms)

    # also mask the TE bodies themselves
    te_bodies: Dict[str, List[Interval]] = defaultdict(list)
    for c, lst in tes.items():
        for l, r, _ in lst:
            te_bodies[c].append((l, r))
    repeats_full: Dict[str, List[Interval]] = {}
    for c in set(list(repeats.keys()) + list(te_bodies.keys())):
        merged = merge_intervals(repeats.get(c, []) + te_bodies.get(c, []))
        repeats_full[c] = merged

    log(f"loading euchromatin BED: {args.euchromatin_bed}")
    euchromatin = load_bed(args.euchromatin_bed, chroms)

    te_flanking = build_te_flanking(tes, args.flank_outer, args.flank_inner)

    # Restrict TE flanking + good-callable scope to euchromatin minus repeats
    # by intersecting per-sample callable with euchromatin then subtracting
    # repeats; the per-TE accumulation step then naturally only counts
    # in-flanking positions.

    gvcf_dir = Path(args.gvcf_dir)
    ancestor_gvcf = gvcf_dir / 'genome.S1.vcf'
    log(f"parsing ancestor gVCF: {ancestor_gvcf}")
    anc_pass = parse_gvcf_pass_intervals(str(ancestor_gvcf), chroms)
    anc_callable: Dict[str, List[Interval]] = {}
    for c in set(list(anc_pass.keys()) + list(euchromatin.keys())):
        x = intersect_intervals(anc_pass.get(c, []), euchromatin.get(c, []))
        x = subtract_intervals(x, repeats_full.get(c, []))
        anc_callable[c] = x

    # ----- per-sample callable -----
    per_sample_rows: List[dict] = []
    per_te_callable_total: Dict[str, int] = defaultdict(int)

    # Read sample names from joint VCF header (also sanity-check num_ma_samples).
    with open(args.joint_vcf) as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                sample_names = line.rstrip('\n').split('\t')[9:]
                break
    if len(sample_names) < 1 + args.num_ma_samples:
        log(f"WARNING: joint VCF has {len(sample_names)} sample columns "
            f"but --num-ma-samples={args.num_ma_samples} expects "
            f">= {1 + args.num_ma_samples}")

    # --sample-names override: replace the (possibly non-unique) #CHROM names
    # with caller-supplied labels, in BAM order. Column 9 = override[0]
    # (ancestor), column 10 = override[1] (first MA), and so on.
    if args.sample_names is not None:
        expected = 1 + args.num_ma_samples
        if len(args.sample_names) != expected:
            raise SystemExit(
                f"--sample-names got {len(args.sample_names)} entries; "
                f"expected {expected} (1 ancestor + {args.num_ma_samples} MA samples)"
            )
        sample_names = args.sample_names + sample_names[len(args.sample_names):]
        log(f"overrode sample names with --sample-names "
            f"(ancestor={args.sample_names[0]}, first MA={args.sample_names[1]})")

    for i in range(args.num_ma_samples):
        col = args.ancestor_col + 1 + i        # joint VCF column
        s_id = i + 2                            # genome.S{N}.vcf index
        gvcf = gvcf_dir / f'genome.S{s_id}.vcf'
        sample_name = sample_names[col - 9] if col - 9 < len(sample_names) else f'col{col}'
        log(f"[{i+1}/{args.num_ma_samples}] {sample_name}: parsing {gvcf.name}")

        sa_pass = parse_gvcf_pass_intervals(str(gvcf), chroms)
        shared: Dict[str, List[Interval]] = {}
        for c in set(list(sa_pass.keys()) + list(anc_callable.keys())):
            x = intersect_intervals(sa_pass.get(c, []), anc_callable.get(c, []))
            shared[c] = x

        total_callable = sum(total_length(v) for v in shared.values())
        per_te, te_total = accumulate_te_callable(shared, te_flanking)
        for tid, bp in per_te.items():
            per_te_callable_total[tid] += bp

        per_sample_rows.append({
            'sample_id': sample_name,
            'genome_S_index': s_id,
            'joint_vcf_col': col,
            'total_callable': total_callable,
            'te_flanking_callable': te_total,
        })

    # ----- de novo mutations -----
    log(f"scanning joint VCF for de novo SNPs: {args.joint_vcf}")
    mutations: List[dict] = list(call_de_novo(
        args.joint_vcf, chroms, repeats_full, euchromatin,
        ancestor_col=args.ancestor_col,
        num_ma_samples=args.num_ma_samples,
        min_alt_count=args.min_alt_count,
        sample_names_override=args.sample_names,
    ))
    log(f"  found {len(mutations)} de novo SNPs")

    # ----- annotate mutations with TE-flanking membership and count per sample -----
    mut_total_by_sample: Dict[str, int] = defaultdict(int)
    mut_te_by_sample: Dict[str, int] = defaultdict(int)
    for m in mutations:
        m['in_te_flanking'] = assign_to_te(m['chrom'], m['pos'], te_flanking, tes) is not None
        sname = m['sample_name']
        mut_total_by_sample[sname] += 1
        if m['in_te_flanking']:
            mut_te_by_sample[sname] += 1
    for r in per_sample_rows:
        r['n_mutations_total_callable'] = mut_total_by_sample.get(r['sample_id'], 0)
        r['n_mutations_te_flanking'] = mut_te_by_sample.get(r['sample_id'], 0)

    # ----- per-TE summary -----
    events = collapse_events(mutations, args.collapse_distance)
    log(f"  collapsed to {len(events)} events at <= {args.collapse_distance} bp")

    te_event_counts: Dict[str, int] = defaultdict(int)
    te_event_samples: Dict[str, set] = defaultdict(set)
    for ev in events:
        tid = assign_to_te(ev['chrom'], ev['pos'], te_flanking, tes)
        if tid is None:
            continue
        te_event_counts[tid] += 1
        te_event_samples[tid].add(ev['sample_name'])

    # ----- write outputs -----
    out_prefix = args.output_prefix
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)

    with open(f'{out_prefix}.per_sample_callable.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['sample_id', 'genome_S_index', 'joint_vcf_col',
                                            'total_callable', 'te_flanking_callable',
                                            'n_mutations_total_callable',
                                            'n_mutations_te_flanking'],
                           delimiter='\t')
        w.writeheader()
        for r in per_sample_rows:
            w.writerow(r)

    with open(f'{out_prefix}.de_novo_mutations.tsv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['sample_id', 'chrom', 'pos', 'ref', 'alt',
                                            'anc_gt', 'sample_gt', 'ref_count', 'alt_count',
                                            'in_te_flanking'],
                           delimiter='\t')
        w.writeheader()
        for m in mutations:
            w.writerow({
                'sample_id': m['sample_name'],
                'chrom': m['chrom'],
                'pos': m['pos'],
                'ref': m['ref'],
                'alt': m['alt'],
                'anc_gt': m['anc_gt'],
                'sample_gt': m['sample_gt'],
                'ref_count': m['ref_count'],
                'alt_count': m['alt_count'],
                'in_te_flanking': int(m['in_te_flanking']),
            })

    with open(f'{out_prefix}.per_te_summary.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['te_id', 'chrom', 'left', 'right',
                    'n_events_collapsed_1kb', 'n_samples_with_events',
                    'callable_summed_across_samples'])
        for chrom, lst in tes.items():
            for left, right, tid in lst:
                w.writerow([
                    tid, chrom, left, right,
                    te_event_counts.get(tid, 0),
                    len(te_event_samples.get(tid, ())),
                    per_te_callable_total.get(tid, 0),
                ])

    log(f"wrote outputs with prefix {out_prefix}.*")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

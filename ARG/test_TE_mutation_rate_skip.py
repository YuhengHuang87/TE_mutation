#!/usr/bin/env python3
"""
Test whether focal mutations (TEs, strong SNPs, or neutral SNPs) affect
the local SNP mutation rate using a paired-lineage comparison across
multiple posterior tree samples.

For each focal mutation at position P carried by haplotypes C:
  Partition branches into with-focal (all descendants carry the mutation)
  and without-focal (no descendants carry it), counting only
  post-insertion branch length and mutations.

Results from each tree sample are placed in separate columns.
Rows represent different focal mutations.
"""

import argparse
import csv
import os
import sys
import time

import numpy as np
import tskit


# ── VCF reader ───────────────────────────────────────────────────────────

def read_focal_entries(vcf_path, chrom, focal="TE"):
    """Read focal mutation entries from VCF.

    focal: "TE" for transposable elements,
           "strong" for HIGH-impact SNPs (snpEff annotation),
           "neutral" for LOW-impact SNPs (snpEff annotation).
    """
    entries = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if cols[0] != chrom:
                continue

            is_snp = "SNP_" in cols[2]
            if focal == "TE":
                if is_snp:
                    continue
            elif focal == "strong":
                if not is_snp or "|HIGH|" not in cols[7]:
                    continue
            elif focal == "neutral":
                if not is_snp or "|LOW|" not in cols[7]:
                    continue

            pos = float(cols[1])
            genotypes = cols[9:]
            carriers = []
            for i, gt in enumerate(genotypes):
                a0, a1 = gt.split("|")
                if a0 == "1":
                    carriers.append(2 * i)
                if a1 == "1":
                    carriers.append(2 * i + 1)
            entries.append({"name": cols[2], "position": pos, "carriers": carriers})
    return entries


def filter_nearby_neutrals(entries, min_dist=5000):
    """If two neutral mutations are within min_dist bp, keep only the one
    with more carriers (higher frequency).  Greedy: sort by num_carriers
    descending, mark neighbours as removed."""
    sorted_by_pos = sorted(entries, key=lambda t: t["position"])
    n = len(sorted_by_pos)
    removed = [False] * n

    # Build index sorted by carrier count (descending) to greedily keep
    # the highest-frequency mutation first.
    order = sorted(range(n), key=lambda i: len(sorted_by_pos[i]["carriers"]))

    for idx in order:
        if removed[idx]:
            continue
        pos_i = sorted_by_pos[idx]["position"]
        # Remove lower-frequency neighbours within min_dist
        for j in range(idx - 1, -1, -1):
            if pos_i - sorted_by_pos[j]["position"] > min_dist:
                break
            if not removed[j]:
                removed[j] = True
        for j in range(idx + 1, n):
            if sorted_by_pos[j]["position"] - pos_i > min_dist:
                break
            if not removed[j]:
                removed[j] = True

    return [t for t, rm in zip(sorted_by_pos, removed) if not rm]


def filter_isolated_tes(te_list, window_half):
    """Keep only TEs with no different-carrier neighbor within ±window_half bp."""
    te_list_sorted = sorted(te_list, key=lambda t: t["position"])
    carrier_sets = [frozenset(t["carriers"]) for t in te_list_sorted]
    n = len(te_list_sorted)
    keep = [True] * n
    for i in range(n):
        if not keep[i]:
            continue
        pos_i = te_list_sorted[i]["position"]
        cs_i = carrier_sets[i]
        for j in range(i + 1, n):
            if te_list_sorted[j]["position"] - pos_i > window_half:
                break
            if carrier_sets[j] != cs_i:
                keep[i] = False
                keep[j] = False
        for j in range(i - 1, -1, -1):
            if pos_i - te_list_sorted[j]["position"] > window_half:
                break
            if carrier_sets[j] != cs_i:
                keep[i] = False
                keep[j] = False
    return [t for t, k in zip(te_list_sorted, keep) if k]


# ── Neighbor-contamination mask ───────────────────────────────────────────

def _compute_contaminated_nodes(tree, nearby_carrier_sets, all_samples):
    """Return the set of node IDs whose subtree carries at least one nearby TE.

    Uses bottom-up bitmask propagation: each sample gets a bitmask of
    which nearby TEs it carries, and each internal node gets the bitwise
    AND of its children's masks.  A non-zero mask means every descendant
    sample is a carrier of that nearby TE, so the branch is
    "contaminated" and should be excluded from the with/without comparison.
    """
    if not nearby_carrier_sets:
        return set()

    n_nearby = len(nearby_carrier_sets)
    all_mask = (1 << n_nearby) - 1

    sample_mask = {}
    for s in all_samples:
        m = 0
        for j, cs in enumerate(nearby_carrier_sets):
            if s in cs:
                m |= (1 << j)
        sample_mask[s] = m

    node_mask = {}
    contaminated = set()

    for node in tree.postorder():
        if tree.is_sample(node):
            m = sample_mask.get(node, 0)
        else:
            m = all_mask
            for child in tree.children(node):
                m &= node_mask.get(child, 0)
                if m == 0:
                    break
        node_mask[node] = m
        if m != 0:
            contaminated.add(node)

    return contaminated


def _precompute_nearby(te_list, window_half):
    """For each TE, collect frozensets of carriers of nearby TEs on
    different haplotypes (within ±window_half bp).  Deduplicated."""
    te_sorted = sorted(te_list, key=lambda t: t["position"])
    carrier_fs = [frozenset(t["carriers"]) for t in te_sorted]
    n = len(te_sorted)

    for i in range(n):
        nearby = set()
        cs_i = carrier_fs[i]
        pos_i = te_sorted[i]["position"]
        for j in range(i - 1, -1, -1):
            if pos_i - te_sorted[j]["position"] > window_half:
                break
            if carrier_fs[j] != cs_i:
                nearby.add(carrier_fs[j])
        for j in range(i + 1, n):
            if te_sorted[j]["position"] - pos_i > window_half:
                break
            if carrier_fs[j] != cs_i:
                nearby.add(carrier_fs[j])
        te_sorted[i]["nearby_carrier_sets"] = list(nearby)

    return te_sorted


# ── Per-TE analysis ──────────────────────────────────────────────────────

def analyze_te(ts, te, node_times, window_half, exclude_half=0,
               nearby_carrier_sets=None, bin_edges=None):
    """Return a result dict or None if the TE cannot be analysed."""

    carriers = te["carriers"]
    pos = te["position"]
    n_samples = ts.num_samples

    if len(carriers) == 0 or len(carriers) == n_samples:
        return None
    if pos < 0 or pos >= ts.sequence_length:
        return None

    ref = ts.at(pos)
    mrca = carriers[0]
    for c in carriers[1:]:
        mrca = ref.mrca(mrca, c)
        if mrca == tskit.NULL:
            return None

    mrca_time = node_times[mrca]
    mrca_parent = ref.parent(mrca)
    if mrca_parent == tskit.NULL:
        return None
    mrca_parent_time = node_times[mrca_parent]
    te_age_midpoint = (mrca_time + mrca_parent_time) / 2.0
    te_age = mrca_time
    if te_age <= 0:
        return None
    is_mono = ref.num_samples(mrca) == len(carriers)

    seq_len = ts.sequence_length

    if bin_edges is not None:
        n_bins = len(bin_edges) - 1
        intervals = []
        for b in range(n_bins):
            inner, outer = bin_edges[b], bin_edges[b + 1]
            intervals.append((max(pos - outer, 0.0), max(pos - inner, 0.0), b))
            intervals.append((min(pos + inner, seq_len), min(pos + outer, seq_len), b))
    else:
        n_bins = 1
        intervals = [
            (max(pos - window_half, 0.0), max(pos - exclude_half, 0.0), 0),
            (min(pos + exclude_half, seq_len), min(pos + window_half, seq_len), 0),
        ]

    tree = tskit.Tree(ts, tracked_samples=carriers)
    nearby = nearby_carrier_sets or []
    all_samples = set(ts.samples()) if nearby else set()

    with_bl = [0.0] * n_bins
    without_bl = [0.0] * n_bins
    with_muts = [0] * n_bins
    without_muts = [0] * n_bins

    for wl, wr, b in intervals:
        if wr <= wl:
            continue
        tree.seek(wl)

        while tree.interval.left < wr:
            tl = max(tree.interval.left, wl)
            tr = min(tree.interval.right, wr)
            span = tr - tl

            contaminated = _compute_contaminated_nodes(
                tree, nearby, all_samples) if nearby else set()

            for node in tree.preorder():
                par = tree.parent(node)
                if par == tskit.NULL:
                    continue
                if node in contaminated:
                    continue
                ct = node_times[node]
                if ct >= te_age:
                    continue
                pt = node_times[par]
                # Skip branches that straddle te_age entirely
                if pt >= te_age:
                    continue
                bl_val = (pt - ct) * span

                nt = tree.num_tracked_samples(node)
                ns = tree.num_samples(node)
                if nt == ns:
                    with_bl[b] += bl_val
                elif nt == 0:
                    without_bl[b] += bl_val

            for site in tree.sites():
                if site.position < wl or site.position >= wr:
                    continue
                for mut in site.mutations:
                    if mut.node in contaminated:
                        continue
                    mp = tree.parent(mut.node)
                    if mp == tskit.NULL:
                        continue
                    mpt = node_times[mp]
                    if mpt >= te_age:
                        continue

                    nt = tree.num_tracked_samples(mut.node)
                    ns = tree.num_samples(mut.node)
                    if nt == ns:
                        with_muts[b] += 1
                    elif nt == 0:
                        without_muts[b] += 1

            if not tree.next():
                break

    tot_w_bl = sum(with_bl)
    tot_wo_bl = sum(without_bl)
    tot_w_m = sum(with_muts)
    tot_wo_m = sum(without_muts)

    r_w = tot_w_m / tot_w_bl if tot_w_bl > 0 else np.nan
    r_wo = tot_wo_m / tot_wo_bl if tot_wo_bl > 0 else np.nan
    ratio = r_w / r_wo if (tot_wo_bl > 0 and tot_w_bl > 0 and r_wo > 0) else np.nan

    result = {
        "te_name": te["name"],
        "position": int(pos),
        "te_age": te_age,
        "te_age_midpoint": te_age_midpoint,
        "num_carriers": len(carriers),
        "is_monophyletic": is_mono,
        "with_te_bl": tot_w_bl,
        "without_te_bl": tot_wo_bl,
        "with_te_muts": tot_w_m,
        "without_te_muts": tot_wo_m,
        "with_te_rate": r_w,
        "without_te_rate": r_wo,
        "rate_ratio": ratio,
    }

    if bin_edges is not None:
        for bi in range(n_bins):
            lo, hi = int(bin_edges[bi]), int(bin_edges[bi + 1])
            tag = f"bin{bi}_{lo}_{hi}_"
            bw = with_muts[bi] / with_bl[bi] if with_bl[bi] > 0 else np.nan
            bwo = without_muts[bi] / without_bl[bi] if without_bl[bi] > 0 else np.nan
            br = bw / bwo if (without_bl[bi] > 0 and with_bl[bi] > 0 and bwo > 0) else np.nan
            result[f"{tag}with_te_bl"] = with_bl[bi]
            result[f"{tag}without_te_bl"] = without_bl[bi]
            result[f"{tag}with_te_muts"] = with_muts[bi]
            result[f"{tag}without_te_muts"] = without_muts[bi]
            result[f"{tag}with_te_rate"] = bw
            result[f"{tag}without_te_rate"] = bwo
            result[f"{tag}rate_ratio"] = br

    return result


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Paired-lineage mutation-rate test around TE insertions "
                    "across multiple posterior tree samples."
    )
    ap.add_argument("-trees", nargs="+", required=True,
                    help="Redated .trees files (one per posterior sample)")
    ap.add_argument("-vcf", required=True,
                    help="Annotated VCF with TEs and snpEff annotations")
    ap.add_argument("-focal", choices=["TE", "strong", "neutral"], default="TE",
                    help="Focal mutation type: TE, strong (HIGH-impact SNPs), "
                         "or neutral (LOW-impact SNPs) (default: TE)")
    ap.add_argument("-chrom", default="X")
    ap.add_argument("-window", type=float, default=5000,
                    help="Flanking window size in bp on each side (default: 5000)")
    ap.add_argument("-exclude", type=float, default=0,
                    help="Exclude +/-N bp around the TE to avoid mapping bias "
                         "(default: 0; set to read length, e.g. 150)")
    ap.add_argument("-output", default=None, help="Per-TE results TSV")
    ap.add_argument("-isolated", action="store_true",
                    help="Only analyse TEs with no different-carrier "
                         "neighbor within the window")
    ap.add_argument("-mask_neighbors", action="store_true",
                    help="Exclude branches that carry any other nearby TE "
                         "(different haplotypes within window) from the "
                         "with/without comparison")
    ap.add_argument("-bins", type=str, default=None,
                    help="Comma-separated distance-bin edges in bp "
                         "(e.g. 500,1000,2000,5000). First edge = exclusion "
                         "zone, last = window size. Overrides -exclude "
                         "and -window.")
    ap.add_argument("-min_carriers", type=int, default=4,
                    help="Only analyse TEs with at least this many carriers "
                         "(default: 4)")
    ap.add_argument("-max_te", type=int, default=None,
                    help="Process at most N TEs (for testing)")
    args = ap.parse_args()

    log = lambda msg: print(msg, file=sys.stderr)

    # ── Parse bins / window ─────────────────────────────────────────────
    if args.bins:
        bin_edges = [float(x) for x in args.bins.split(",")]
        excl = bin_edges[0]
        window = bin_edges[-1]
        n_bins = len(bin_edges) - 1
    else:
        bin_edges = None
        excl = args.exclude
        window = args.window
        n_bins = 0

    # ── Read and filter focal mutations (independent of trees) ─────────
    log(f"Reading focal mutations (type={args.focal}) from VCF (chrom={args.chrom})")
    te_list = read_focal_entries(args.vcf, args.chrom, focal=args.focal)
    log(f"  {len(te_list)} {args.focal} entries read from VCF")
    if args.min_carriers > 0:
        before = len(te_list)
        te_list = [t for t in te_list if len(t["carriers"]) >= args.min_carriers]
        log(f"  Filtered to {len(te_list)} TEs with >= {args.min_carriers} carriers "
            f"(dropped {before - len(te_list)})")
    if args.focal == "neutral":
        before = len(te_list)
        te_list = filter_nearby_neutrals(te_list, min_dist=5000)
        log(f"  Thinned nearby neutrals (5 kb): {len(te_list)} kept "
            f"(dropped {before - len(te_list)} lower-frequency neighbours)")
    if args.isolated:
        before = len(te_list)
        te_list = filter_isolated_tes(te_list, window)
        log(f"  Filtered to {len(te_list)} isolated TEs "
            f"(dropped {before - len(te_list)} with different-carrier "
            f"neighbors within +/-{window:.0f} bp)")
    if args.mask_neighbors:
        te_list = _precompute_nearby(te_list, window)
        n_with = sum(1 for t in te_list if t.get("nearby_carrier_sets"))
        log(f"  Mask-neighbors mode: {n_with}/{len(te_list)} TEs have "
            f"nearby different-carrier TEs to mask")
    if args.max_te:
        te_list = te_list[:args.max_te]
    log(f"  {len(te_list)} TEs to process")

    # ── Extract tree IDs from filenames ─────────────────────────────────
    tree_ids = []
    for path in args.trees:
        basename = os.path.basename(path)
        name = basename.rsplit(".", 1)[0]       # strip .trees
        tid = name.rsplit("_", 1)[-1]           # last segment after _
        tree_ids.append(tid)
    log(f"  {len(args.trees)} tree files to process: "
        f"sample IDs {tree_ids[0]}..{tree_ids[-1]}")

    mask_note = ", masking neighbor-TE branches" if args.mask_neighbors else ""
    if bin_edges:
        bin_strs = [f"{int(bin_edges[i])}-{int(bin_edges[i+1])}"
                    for i in range(n_bins)]
        log(f"Analysing bins (bp from TE): {', '.join(bin_strs)}"
            f"{mask_note} ...")
    else:
        log(f"Analysing (window = +/-{window:.0f} bp"
            f"{f', excluding +/-{excl:.0f} bp' if excl > 0 else ''}"
            f"{mask_note}) ...")

    # ── Initialize result rows (one per TE) ─────────────────────────────
    te_results = {}
    for te in te_list:
        key = (te["name"], int(te["position"]))
        te_results[key] = {
            "te_name": te["name"],
            "position": int(te["position"]),
            "num_carriers": len(te["carriers"]),
        }

    per_tree_keys = [
        "te_age", "te_age_midpoint", "is_monophyletic",
        "with_te_bl", "without_te_bl",
        "with_te_muts", "without_te_muts",
        "with_te_rate", "without_te_rate", "rate_ratio",
    ]
    bin_keys = [
        "with_te_bl", "without_te_bl",
        "with_te_muts", "without_te_muts",
        "with_te_rate", "without_te_rate", "rate_ratio",
    ]

    # ── Process each tree file ──────────────────────────────────────────
    t0 = time.time()
    for ti, trees_path in enumerate(args.trees):
        tid = tree_ids[ti]
        log(f"  [{ti+1}/{len(args.trees)}] Loading tree sample {tid}: "
            f"{os.path.basename(trees_path)}")
        ts = tskit.load(trees_path)
        node_times = ts.nodes_time

        for te in te_list:
            nearby = te.get("nearby_carrier_sets") if args.mask_neighbors else None
            r = analyze_te(ts, te, node_times, window, excl,
                           nearby_carrier_sets=nearby, bin_edges=bin_edges)

            prefix = f"t{tid}_"
            key = (te["name"], int(te["position"]))
            row = te_results[key]
            if r is not None:
                for key in per_tree_keys:
                    row[f"{prefix}{key}"] = r[key]
                if bin_edges is not None:
                    for bi in range(n_bins):
                        lo, hi = int(bin_edges[bi]), int(bin_edges[bi + 1])
                        btag = f"bin{bi}_{lo}_{hi}_"
                        for s in bin_keys:
                            row[f"{prefix}{btag}{s}"] = r.get(f"{btag}{s}", np.nan)
            else:
                for key in per_tree_keys:
                    row[f"{prefix}{key}"] = np.nan
                if bin_edges is not None:
                    for bi in range(n_bins):
                        lo, hi = int(bin_edges[bi]), int(bin_edges[bi + 1])
                        btag = f"bin{bi}_{lo}_{hi}_"
                        for s in bin_keys:
                            row[f"{prefix}{btag}{s}"] = np.nan

        log(f"    Done ({time.time() - t0:.0f}s elapsed)")

    elapsed = time.time() - t0
    log(f"  All {len(args.trees)} tree samples processed in {elapsed:.0f}s")

    # ── Build column list ───────────────────────────────────────────────
    fixed_fields = ["te_name", "position", "num_carriers"]
    fields = list(fixed_fields)
    for tid in tree_ids:
        prefix = f"t{tid}_"
        for key in per_tree_keys:
            fields.append(f"{prefix}{key}")
        if bin_edges is not None:
            for bi in range(n_bins):
                lo, hi = int(bin_edges[bi]), int(bin_edges[bi + 1])
                btag = f"bin{bi}_{lo}_{hi}_"
                for s in bin_keys:
                    fields.append(f"{prefix}{btag}{s}")

    # ── Write output ────────────────────────────────────────────────────
    out = open(args.output, "w", newline="") if args.output else sys.stdout
    w = csv.DictWriter(out, fieldnames=fields, delimiter="\t",
                       extrasaction="ignore")
    w.writeheader()
    for te in te_list:
        key = (te["name"], int(te["position"]))
        w.writerow(te_results[key])
    if args.output:
        out.close()
        log(f"  Per-TE results -> {args.output}")


if __name__ == "__main__":
    main()

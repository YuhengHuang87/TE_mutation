#!/usr/bin/env python3
"""TE annotation post-processing for the MA-line nanopore pipeline.

Python port of the post-RepeatMasker Perl scripts:

  step1_extract_euchromatic
      <- extract_TE_repeatmasker_nano.pl
      Filter RepeatMasker .out hits to euchromatic regions and selected
      TE classes; emit a per-TE FASTA library cut from the scaffold.

  step2_compute_distances
      <- TE_nearby_distance_nanopore.pl
      Annotate each kept hit with the distance to its left/right
      neighbour ("NA" across chromosome boundaries).

  step3_merge_nearby_family
      <- TE_merge_exclude_nearby_family_combined_nanopore.pl
      Merge runs of adjacent same-family TEs falling within
      --family-merge-dist bp of each other.

  step4_filter_length
      <- filter_TE_length_INE.pl
      Drop hits whose family matches --exclude-family-pattern (e.g. INE-1)
      and hits shorter than --min-length.

  step5_merge_by_class
      <- merge_te_script.py
      Merge remaining hits by DNA-vs-RNA class within --class-merge-dist.

Inputs:
  - RepeatMasker .out file  (--rm-out)
  - Scaffolded assembly FASTA  (--scaffold)
  - BED of euchromatic regions  (--euchromatin)

Outputs (in --outdir):
  filtered.out         step 1: euchromatic-restricted RM rows
  library.fasta        step 1: TE sequences pulled from the scaffold
  distance.txt         step 2: filtered.out + dist_prev + dist_next
  merged_family.txt    step 3: 8-col table with merged/unmerged status
  length_filtered.txt  step 4: 5-col table (chr, start, end, fam, super)
  merged_class.txt     step 5: 5-col table merged by RNA/DNA class
"""

import argparse
import os
import re
import sys


def read_euchromatin(path):
    regions = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            regions.setdefault(chrom, []).append((start, end))
    return regions


def in_euchromatin(chrom, left, right, regions):
    for s, e in regions.get(chrom, ()):
        if left > s and right < e:
            return True
    return False


def step1_extract_euchromatic(rm_out, scaffold, regions, te_classes,
                              out_filtered, out_fasta):
    pos_to_name = {}
    n_kept = 0

    with open(rm_out) as fin, open(out_filtered, "w") as fout:
        for line in fin:
            parts = line.split()
            if len(parts) < 11:
                continue
            try:
                left = int(parts[5])
                right = int(parts[6])
            except ValueError:
                continue  # RM .out header rows
            chrom = parts[4]
            class_family = parts[10]

            if not in_euchromatin(chrom, left, right, regions):
                continue
            if not any(c in class_family for c in te_classes):
                continue

            name = f"{class_family}:{chrom}:{left}:{right}"
            fout.write(line)
            for i in range(left, right + 1):
                pos_to_name[(chrom, i)] = name
            n_kept += 1

    sequences = {}
    chrom = None
    pos = 0
    with open(scaffold) as fin:
        for line in fin:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                chrom = line[1:].split()[0]
                pos = 0
                continue
            if chrom is None:
                continue
            for ch in line:
                pos += 1
                name = pos_to_name.get((chrom, pos))
                if name is not None:
                    sequences.setdefault(name, []).append(ch)

    with open(out_fasta, "w") as fout:
        for name, chars in sequences.items():
            fout.write(f">{name}\n{''.join(chars)}\n")

    print(f"step1: {n_kept} euchromatic TE hits kept", file=sys.stderr)


def step2_compute_distances(in_filtered, out_distances):
    rows = []
    with open(in_filtered) as fin:
        for line in fin:
            parts = line.split()
            if len(parts) < 11:
                continue
            try:
                left = int(parts[5])
                right = int(parts[6])
            except ValueError:
                continue
            rows.append((parts[4], left, right, line.rstrip("\n")))

    with open(out_distances, "w") as fout:
        for i, (chrom, left, right, orig) in enumerate(rows):
            if i > 0 and rows[i - 1][0] == chrom:
                dist_prev = str(left - rows[i - 1][2])
            else:
                dist_prev = "NA"
            if i + 1 < len(rows) and rows[i + 1][0] == chrom:
                dist_next = str(rows[i + 1][1] - right)
            else:
                dist_next = "NA"
            fout.write(f"{orig} {dist_prev} {dist_next}\n")


def _close_flag(value, threshold):
    if value == "NA":
        return False
    return int(value) < threshold


def step3_merge_nearby_family(in_distances, out_merged, dist_threshold):
    g_chr = ""
    g_start = 0
    g_family = ""
    g_super = ""
    n_out = 0

    with open(in_distances) as fin, open(out_merged, "w") as fout:
        for line in fin:
            b = line.split()
            if len(b) < 13:
                continue
            chrom = b[4]
            left = int(b[5])
            right = int(b[6])
            family = b[9]
            super_family = b[10]
            dist_prev = b[-2]
            dist_next = b[-1]
            close_prev = _close_flag(dist_prev, dist_threshold)
            close_next = _close_flag(dist_next, dist_threshold)

            if not close_prev and not close_next:
                fout.write(
                    f"{chrom}\t{left}\t{right}\t{family}\t{super_family}"
                    f"\t{dist_prev}\t{dist_next}\tunmerged\n"
                )
                n_out += 1
            elif not close_prev and close_next:
                g_chr = chrom
                g_start = left
                g_family = family
                g_super = super_family
            elif close_prev and not close_next:
                if family == g_family:
                    out_family = family
                    out_super = super_family
                else:
                    out_family = f"{g_family}_{family}"
                    out_super = f"{g_super}_{super_family}"
                fout.write(
                    f"{g_chr}\t{g_start}\t{right}\t{out_family}\t{out_super}"
                    f"\t{dist_prev}\t{dist_next}\tmerged\n"
                )
                n_out += 1
            else:
                if family != g_family:
                    g_family = f"{g_family}_{family}"
                    g_super = f"{g_super}_{super_family}"

    print(f"step3: {n_out} merged/unmerged records written", file=sys.stderr)


def step4_filter_length(in_merged, out_filtered, min_length, exclude_pattern):
    pat = re.compile(exclude_pattern)
    n_out = 0
    with open(in_merged) as fin, open(out_filtered, "w") as fout:
        for line in fin:
            b = line.rstrip("\n").split("\t")
            if len(b) < 5:
                continue
            chrom, left, right, family, super_family = b[0], int(b[1]), int(b[2]), b[3], b[4]
            if pat.search(family):
                continue
            if (right - left) < min_length:
                continue
            fout.write(f"{chrom}\t{left}\t{right}\t{family}\t{super_family}\n")
            n_out += 1
    print(f"step4: {n_out} TEs after length/family filter", file=sys.stderr)


def _te_class(super_family):
    is_dna = "DNA" in super_family
    is_rna = "LTR" in super_family or "LINE" in super_family
    if is_dna and not is_rna:
        return "DNA"
    if is_rna and not is_dna:
        return "RNA"
    if is_dna and is_rna:
        return "Ambiguous"
    return "Other"


def step5_merge_by_class(in_filtered, out_merged, dist_threshold):
    tes = []
    with open(in_filtered) as fin:
        for line in fin:
            parts = line.strip().split()
            if len(parts) < 5:
                continue
            try:
                tes.append({
                    "chrom": parts[0],
                    "start": int(parts[1]),
                    "end": int(parts[2]),
                    "name": parts[3],
                    "super": parts[4],
                    "class": _te_class(parts[4]),
                })
            except ValueError:
                continue

    tes.sort(key=lambda x: (x["chrom"], x["start"]))

    if not tes:
        open(out_merged, "w").close()
        print("step5: no TEs to merge", file=sys.stderr)
        return

    groups = [[tes[0]]]
    for te in tes[1:]:
        prev = groups[-1][-1]
        gclass = groups[-1][0]["class"]
        same_chrom = prev["chrom"] == te["chrom"]
        within = (te["start"] - prev["end"]) <= dist_threshold
        match = gclass == te["class"] and gclass in ("DNA", "RNA")
        if same_chrom and within and match:
            groups[-1].append(te)
        else:
            groups.append([te])

    with open(out_merged, "w") as fout:
        for g in groups:
            chrom = g[0]["chrom"]
            start = g[0]["start"]
            end = max(t["end"] for t in g)
            names = ",".join(t["name"] for t in g)
            supers = ",".join(t["super"] for t in g)
            fout.write(f"{chrom}\t{start}\t{end}\t{names}\t{supers}\n")

    print(f"step5: {len(tes)} -> {len(groups)} TEs after class merge",
          file=sys.stderr)


def main():
    p = argparse.ArgumentParser(
        description="TE annotation post-processing for MA-line nanopore assemblies.",
    )
    p.add_argument("--rm-out", required=True, help="RepeatMasker .out file")
    p.add_argument("--scaffold", required=True, help="Scaffolded assembly FASTA")
    p.add_argument("--euchromatin", required=True,
                   help="BED of euchromatic regions: chrom<TAB>start<TAB>end")
    p.add_argument("--outdir", required=True)
    p.add_argument("--te-classes", default="LINE,LTR,DNA,Unknown",
                   help="Comma-separated substrings to retain in RM class/family")
    p.add_argument("--family-merge-dist", type=int, default=200,
                   help="step 3: max bp gap to merge same-family neighbours")
    p.add_argument("--min-length", type=int, default=500,
                   help="step 4: minimum TE length after family merge")
    p.add_argument("--exclude-family-pattern", default="INE-1",
                   help="step 4: regex of family names to drop")
    p.add_argument("--class-merge-dist", type=int, default=5000,
                   help="step 5: max bp gap to merge by RNA/DNA class")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    classes = [c.strip() for c in args.te_classes.split(",") if c.strip()]
    regions = read_euchromatin(args.euchromatin)

    out1_filt = os.path.join(args.outdir, "filtered.out")
    out1_fa = os.path.join(args.outdir, "library.fasta")
    out2 = os.path.join(args.outdir, "distance.txt")
    out3 = os.path.join(args.outdir, "merged_family.txt")
    out4 = os.path.join(args.outdir, "length_filtered.txt")
    out5 = os.path.join(args.outdir, "merged_class.txt")

    step1_extract_euchromatic(args.rm_out, args.scaffold, regions, classes,
                              out1_filt, out1_fa)
    step2_compute_distances(out1_filt, out2)
    step3_merge_nearby_family(out2, out3, args.family_merge_dist)
    step4_filter_length(out3, out4, args.min_length, args.exclude_family_pattern)
    step5_merge_by_class(out4, out5, args.class_merge_dist)


if __name__ == "__main__":
    main()

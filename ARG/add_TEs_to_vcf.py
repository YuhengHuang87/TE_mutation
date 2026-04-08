#!/usr/bin/env python3
"""
Add TE annotations as pseudo-SNP entries in a VCF file.

For each TE in each strain's BED file:
  - POS = midpoint of TE coordinates: (col2 + col3) // 2
  - ID  = last column of the BED file (TE family name)
  - Genotype = 1|1 for strains carrying the TE, 0|0 for others

TE entries are merged with existing SNPs and sorted by chromosome and position.
"""

import os
import sys
from collections import defaultdict

# ── Paths ──────────────────────────────────────────────────────────────────
VCF_FILE = "./PhasedSNPsFitTEs.vcf"
BED_DIR  = "../TE_annotations/ReferenceCoordinates"
OUTPUT   = "./PhasedSNPsFitTEs_with_TEs.vcf"

CHROM_ORDER = ["2L", "2R", "3L", "3R", "X"]
chrom_rank  = {c: i for i, c in enumerate(CHROM_ORDER)}

# ── 1. Read VCF: first header block + all data lines ──────────────────────
print("Reading VCF …")
header_lines = []
sample_names = []
# Store data as list of (chrom, pos, line)
snp_rows = []

with open(VCF_FILE) as fh:
    in_first_header = True
    header_seen = False
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith("#"):
            # Only keep the first contiguous header block
            if in_first_header:
                header_lines.append(line)
                if line.startswith("#CHROM"):
                    fields = line.split("\t")
                    sample_names = fields[9:]
                    header_seen = True
            continue
        # First non-header line → end of first header block
        in_first_header = False
        fields = line.split("\t", 3)  # only split enough to get chrom + pos
        chrom = fields[0]
        pos   = int(fields[1])
        snp_rows.append((chrom, pos, line))

n_samples = len(sample_names)
sample_idx = {name: i for i, name in enumerate(sample_names)}
print(f"  {n_samples} samples, {len(snp_rows)} existing SNP records")

# ── 2. Read TE BED files ──────────────────────────────────────────────────
# Key: (chrom, start, end, family) → set of strains
print("Reading TE BED files …")
te_dict = defaultdict(set)   # (chrom, start, end, family) → {strain, …}
te_family = {}               # (chrom, start, end, family) → family  (redundant but clear)

bed_files = sorted(f for f in os.listdir(BED_DIR) if f.endswith("_Ref_Coord.bed"))
for bed_file in bed_files:
    strain = bed_file.replace("_Ref_Coord.bed", "")
    if strain not in sample_idx:
        print(f"  Skipping {strain} (not in VCF samples)")
        continue

    bed_path = os.path.join(BED_DIR, bed_file)
    with open(bed_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom  = parts[0]
            start  = int(parts[1])
            end    = int(parts[2])
            family = parts[-1]           # last column = TE family

            key = (chrom, start, end, family)
            te_dict[key].add(strain)

print(f"  {len(te_dict)} unique TE insertions across all strains")

# ── 3. Build TE VCF rows ──────────────────────────────────────────────────
print("Building TE VCF rows …")
te_rows = []
for (chrom, start, end, family), strains in te_dict.items():
    pos = (start + end) // 2            # midpoint (integer)
    te_id = family                      # ID = TE family name

    # Genotype columns
    gts = []
    for s in sample_names:
        gts.append("1|1" if s in strains else "0|0")

    ac = sum(2 for s in sample_names if s in strains)
    af = ac / (2 * n_samples) if n_samples else 0

    row = "\t".join([
        chrom,
        str(pos),
        te_id,
        "A",                            # REF  (not a real SNP)
        "C",                         # ALT
        ".",                            # QUAL
        ".",                            # FILTER
        f"AC={ac};AF={af:.6g};CM=0",         # INFO
        "GT",                           # FORMAT
    ] + gts)

    te_rows.append((chrom, pos, row))

print(f"  {len(te_rows)} TE rows created")

# ── 4. Merge & sort ───────────────────────────────────────────────────────
print("Merging and sorting …")

def sort_key(item):
    chrom, pos, _ = item
    return (chrom_rank.get(chrom, 999), pos)

all_rows = snp_rows + te_rows
all_rows.sort(key=sort_key)

# ── 5. Write output VCF ───────────────────────────────────────────────────
print(f"Writing {len(all_rows)} records to {OUTPUT} …")
with open(OUTPUT, "w") as out:
    for h in header_lines:
        out.write(h + "\n")
    for _, _, row in all_rows:
        out.write(row + "\n")

print("Done.")

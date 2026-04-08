import argparse
import glob
import os
import tskit

parser = argparse.ArgumentParser(
    description="Remove focal mutation sites from tree sequences. "
                "Always removes TEs; optionally also removes strong or neutral SNPs."
)
parser.add_argument("--chr", required=True, help="Chromosome arm (e.g. 2L, 3R, X)")
parser.add_argument("--focal", required=True, choices=["TE", "strong", "neutral"],
                    help="TE: remove only TEs; strong: remove TEs + HIGH-impact SNPs; "
                         "neutral: remove TEs + LOW-impact SNPs")
parser.add_argument("--trees_dir", required=True,
                    help="Directory containing input singer_{chr}_*.trees files")
parser.add_argument("--output_dir", required=True,
                    help="Output directory for filtered trees")
parser.add_argument("--te_vcf", required=True,
                    help="VCF with TE entries (non-SNP rows)")
parser.add_argument("--ann_vcf", default=None,
                    help="Annotated VCF with snpEff annotations "
                         "(required for --focal strong or neutral)")
args = parser.parse_args()

CHROMOSOME = args.chr
os.makedirs(args.output_dir, exist_ok=True)

# 1. Always collect TE positions (non-SNP entries)
positions_to_remove = set()
with open(args.te_vcf) as f:
    for line in f:
        if line.startswith('#'):
            continue
        cols = line.strip().split('\t')
        if "SNP_" not in cols[2]:
            positions_to_remove.add(float(cols[1]))

print(f"Found {len(positions_to_remove)} TE positions to remove")

# 2. For strong/neutral, also collect the corresponding SNP positions
if args.focal in ("strong", "neutral"):
    if not args.ann_vcf:
        raise ValueError(f"--ann_vcf is required for --focal {args.focal}")
    impact = "|HIGH|" if args.focal == "strong" else "|LOW|"
    count_before = len(positions_to_remove)
    with open(args.ann_vcf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if impact in cols[7]:
                positions_to_remove.add(float(cols[1]))
    print(f"Found {len(positions_to_remove) - count_before} "
          f"{args.focal}-impact SNP positions to remove")

print(f"Total positions to remove: {len(positions_to_remove)}")

# 3. Process every singer_{chr}_*.trees file
tree_files = sorted(glob.glob(
    os.path.join(args.trees_dir, f"singer_{CHROMOSOME}_*.trees")))
print(f"Found {len(tree_files)} trees files to process")

for tree_file in tree_files:
    basename = os.path.basename(tree_file)
    ts = tskit.load(tree_file)

    site_ids_to_remove = [s.id for s in ts.sites()
                          if s.position in positions_to_remove]
    ts_filtered = ts.delete_sites(site_ids_to_remove)

    out_path = os.path.join(args.output_dir, basename)
    ts_filtered.dump(out_path)

    print(f"{basename}: removed {len(site_ids_to_remove)} sites "
          f"({ts.num_sites} -> {ts_filtered.num_sites} sites)")

print(f"\nDone. Filtered files saved to: {args.output_dir}")

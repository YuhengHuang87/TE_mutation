# MA_mutation_TE — de novo SNPs in TE flanking windows

This pipeline ports three legacy Perl scripts to a single Python tool plus a
Snakemake driver, covering all six (group × reference) MA-line combinations.

- Original Perl:
  - [MAlines_callable_sites_joint_consider_each_TEs.pl](MAlines_callable_sites_joint_consider_each_TEs.pl)
  - [MAlines_mutation_joint_consider_sample_ID.pl](MAlines_mutation_joint_consider_sample_ID.pl)
- Python rewrite: [MA_mutation_TE.py](MA_mutation_TE.py)
- Snakemake driver: [MA_mutation_TE.smk](MA_mutation_TE.smk)
- Spec the rewrite is based on: [MA_lines_mutation_around_TE.md](MA_lines_mutation_around_TE.md)

## What it does

For each MA-line group (`Act` / `Control` / `Ubi`) and each reference (`A4` /
`mc_initial`):

1. **Per-sample callable sites in TE flanking windows.**
   For every MA sample `i`, intersect that sample's PASS regions in
   `genome.S{i+1}.vcf` (Strelka per-sample gVCF) with the ancestor's PASS
   regions in `genome.S1.vcf`, restrict to the euchromatin BED, and subtract
   RepeatMasker repeats and the TE bodies themselves. Then sum bp falling in
   any TE's `[left-5000, left-500)` ∪ `(right+500, right+5000]` flanking
   window.

2. **De novo SNPs from the joint Strelka VCF.**
   Walk the joint `variants.vcf`. For each biallelic single-base SNV at a
   non-repeat euchromatin position whose ancestor genotype is `0/0` or `1/1`
   and `PASS`, find MA samples with a differing genotype. Emit only sites
   where exactly one MA sample differs, that sample is `PASS`, and the new
   allele's read depth (AD₂ if ancestor `0/0`; AD₁ if ancestor `1/1`) is
   `≥ 5`.

3. **Per-TE summary.**
   Collapse mutations within `1000 bp` per sample into single events, attribute
   each event to the TE whose flanking window contains it (closest TE body on
   ties), and report per-TE event counts, contributing-sample counts, and
   the per-TE callable-site sum.

The 6 reference parameter sets:

| Reference   | RepeatMasker                                               | TE annotation                                                        | Euchromatin BED                                  | Chromosomes                                                  |
|-------------|------------------------------------------------------------|----------------------------------------------------------------------|--------------------------------------------------|--------------------------------------------------------------|
| A4          | `A4.Hifi.Scaf.Repeat.MT.fasta.out`                          | `A4_A7.UU_MU.Zscore.SigFlag.txt`                                     | `A4_hifi_euchromatin_region.bed.gz`              | 2L, 2R, 3L, 3R, X                                            |
| mc_initial  | `nanopore.../repeatmasker/mc_TE_eu_ini_scaffold.out`        | `nanopore.../te_annotation/merged_class.txt`                          | `mc_nanopore_euchromatin_region.bed.gz`           | 2L_RagTag, 2R_RagTag, 3L_RagTag, 3R_RagTag, X_RagTag         |

The MA-line counts are pulled from the `GROUPS` dict in
[repository/TE_mutation/MA_line_analysis/TE_anno_mapping.smk](repository/TE_mutation/MA_line_analysis/TE_anno_mapping.smk):
Act = 48, Control = 47, Ubi = 52 (each plus one ancestor at column 9 of the
joint VCF).

## Inputs (per `group_ref`)

The Snakemake rule consumes the upstream Strelka outputs produced by
[TE_anno_mapping.smk](repository/TE_mutation/MA_line_analysis/TE_anno_mapping.smk):

```
{STRELKA_DIR}/{group}_{ref}/results/variants/
├── variants.vcf            # joint VCF (--joint-vcf)
├── genome.S1.vcf           # ancestor gVCF (genome.S1)
├── genome.S2.vcf           # MA sample 1 gVCF
├── ...
└── genome.S{1+N}.vcf       # MA sample N gVCF (--gvcf-dir is the directory)
```

The mapping `column index ↔ genome.S*.vcf index` follows directly from the
order BAMs were given to `configureStrelkaGermlineWorkflow.py`. In the new
pipeline, that order is `[ancestor] + GROUPS[group]["samples"]` (see
[TE_anno_mapping.smk:222-225](repository/TE_mutation/MA_line_analysis/TE_anno_mapping.smk#L222-L225)),
so column 9 of `variants.vcf` is the ancestor and corresponds to
`genome.S1.vcf`; column 10 is the first MA sample and corresponds to
`genome.S2.vcf`; and so on.

## Outputs

Each `(group, ref)` writes three TSVs to
`{OUT_ROOT}/{group}_{ref}/`. Filenames are prefixed with `{group}_{ref}.` so
they remain self-identifying when copied or aggregated outside the run dir:

- `{group}_{ref}.per_sample_callable.tsv` — `sample_id, genome_S_index,
  joint_vcf_col, total_callable, te_flanking_callable,
  n_mutations_total_callable, n_mutations_te_flanking`
- `{group}_{ref}.de_novo_mutations.tsv` — `sample_id, chrom, pos, ref, alt,
  anc_gt, sample_gt, ref_count, alt_count, in_te_flanking` (`in_te_flanking`
  is `1` if the SNP falls in any TE's `[outer..inner)` flanking window, else `0`)
- `{group}_{ref}.per_te_summary.tsv` — `te_id, chrom, left, right,
  n_events_collapsed_1kb, n_samples_with_events,
  callable_summed_across_samples` (one row per TE; supersedes the
  per-TE callable totals)

The `te_id` is `chrom:left-right`, taken directly from the TE-annotation file
the run was given.

## Running

```bash
# All six combinations
snakemake -s MA_mutation_TE.smk --cores 6

# Just one combination
snakemake -s MA_mutation_TE.smk --cores 1 Act_A4

# SLURM (mirrors TE_anno_mapping.smk's launch style)
snakemake -s MA_mutation_TE.smk \
  --executor slurm \
  --jobs 6 \
  --default-resources slurm_account=grylee_lab slurm_partition=standard runtime=$((2*24*60))
```

`mutation_te_analysis` runs single-threaded; the wall-clock cost is dominated
by streaming the per-sample gVCFs once each. Memory is bounded by the
PASS-interval lists (kilobytes per sample after merging), not by per-base
hashes — the original Perl stored every base as a hash key, which here is
replaced with sorted-interval intersection / subtraction.

## Differences from the original Perl

- **Interval algebra instead of per-base hashes.** Repeat masking, ancestor
  callable lookup, TE flanking attribution, and per-sample intersections are
  all done with sorted `[start, end]` intervals plus binary search. The Perl
  version's `$repeat_site{$chr."\t".$j}=1` loop is gone. This is the main
  reason the script can run on a laptop.
- **Euchromatin filter is data-driven.** The Perl hard-coded
  `2L_RagTag>97997 && <21400000`, etc. The Python script reads
  `--euchromatin-bed` (the same BED Strelka was given as `--callRegions`) so
  changing the euchromatin definition does not require editing source.
- **Fixed the `0/0` flank exclusion.** The Perl callable script enumerates
  positions in `[left-5000, left-1]` and `[right+1, right+5000]` — i.e. it
  does *not* exclude the inner 500 bp from code; it relies on the input TE
  file already having the inner buffer baked in. Here, `--flank-inner` (default
  500) is applied explicitly, so any TE annotation can be used safely.
- **Multi-mutation collapse.** New: per-sample mutations within
  `--collapse-distance` (default 1000 bp) collapse to a single event before
  per-TE attribution. This is step 3 from
  [MA_lines_mutation_around_TE.md](MA_lines_mutation_around_TE.md).
- **One script for both steps + summary.** The three Perl scripts (callable,
  mutation, summary) plus their per-group / per-reference duplication collapse
  into one `MA_mutation_TE.py` invocation per `(group, ref)` driven by the
  Snakefile.
- **Output files are TSV with headers.** The original `OUT1`/`OUT2` produced
  unnamed-column whitespace files; outputs here are tab-separated with named
  columns to drop straight into pandas / R.
- **No genotype-string `|` heuristic for MA samples in the joint VCF.** The
  Perl skipped any line whose ancestor field contained `|` (a phased
  genotype); Strelka germline does not emit phased genotypes, so the
  Python version omits this filter for the joint VCF. It is retained when
  parsing the per-sample gVCFs (where it was acting as a lightweight
  "skip indel-block" proxy).
- **Tie-break for overlapping flanking windows.** When a position falls in
  more than one TE's flanking window, the Python version assigns it to the
  TE whose body is closest. The Perl let "last-write-wins" via hash overwrite;
  the new behavior is order-independent. With the current TE annotations
  (`...merged_5kb...txt` for mc_initial), overlap is rare in practice.

## Open items / things to double-check before scaling up

- The TE-annotation file for `A4` (`A4_A7.UU_MU.Zscore.SigFlag.txt`) carries
  Z-score / "SigFlag" columns past column 3. Only columns 1–3 (chr, left,
  right) are used. If you want to subset to "significant" TEs only, filter
  the file before passing it to `--te-annotation`, or extend the script to
  honor a `SigFlag` column.
- `--num-ma-samples` is asserted against the `#CHROM` header of the joint
  VCF; if Strelka was run with a different BAM count than `NUM_MA[group]`,
  the script logs a warning and proceeds with whichever number is smaller.
  Verify the joint VCF header matches `GROUPS[group]["samples"]` in
  [TE_anno_mapping.smk](repository/TE_mutation/MA_line_analysis/TE_anno_mapping.smk).
- `--sample-names` overrides the names parsed from the joint VCF `#CHROM`
  header. Use it when the BAMs going into Strelka shared a generic
  `RGSM` (e.g. `Sample1`) so the joint VCF can no longer distinguish samples
  by name. The Snakefile populates this from the `SAMPLES` dict in BAM order
  (`[ancestor] + samples`, with a `_Sample` suffix to mirror the
  `add_read_groups` rule's RGSM convention). The column-index → BAM-order
  mapping is fixed (column 9 = BAM 1, column 10 = BAM 2, …), so the override
  is unambiguous as long as `SAMPLES[group]` is kept in sync with
  `GROUPS[group]` in TE_anno_mapping.smk.
- The script intentionally **does not** treat the joint VCF site-level FILTER
  column (`parts[6]`) as a hard filter, matching the original Perl. The
  per-sample FT sub-field is gated for `PASS` on both the ancestor and the
  focal MA sample (see `'PASS' not in anc` /
  `'PASS' not in samp` in `call_de_novo`). If you want to additionally
  require `parts[6] == 'PASS'`, add a one-liner in `call_de_novo`.

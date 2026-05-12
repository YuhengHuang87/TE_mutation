# ARG.smk — ARG-Based Mutation Rate Workflow

Snakemake (v9) workflow for estimating per-site mutation rates using Ancestral Recombination Graphs (ARGs). For each chromosome arm, the pipeline annotates a phased SNP VCF with transposable element (TE) positions, haploidizes it, builds ARGs with SINGER per 5 Mb window and merges the windows, strips TE sites from the merged trees, redates them with Polegon, and measures mutation rates for three focal classes (TE, strong SNPs, neutral SNPs) on the TE-filtered trees.

## Usage

```bash
# Full run on SLURM
~/.pixi/bin/snakemake -s ARG.smk --executor slurm -j 5 --latency-wait 60 &

# Dry run
snakemake -s ARG.smk --executor dryrun -j 10 --latency-wait 60

# Target specific outputs
snakemake -s ARG.smk --executor slurm -j 10 \
  .../results/TE_mutation_rate_weight_on_TE_trees_2R.tsv
```

To restrict focal classes without editing the file, override `focal=FOCAL_ANALYSIS_TYPES` in `rule all` (e.g. `focal=["TE"]`), or request specific output files on the CLI.

> **Note:** Before running, update the hardcoded paths inside `code/add_TEs_to_vcf.py` to match the current HPC layout.

## Configuration

| Variable | Purpose |
|---|---|
| `BASE_DIR` | Root of ARG analysis: `/dfs7/grylee/yuhenh3/mutation_accumulation/ARG` |
| `VCF_DIR` | `BASE_DIR/Rech_2022_data/SNPs_vcf` — input/output VCFs and results |
| `BED_DIR` | `BASE_DIR/Rech_2022_data/TE_annotations/ReferenceCoordinates` |
| `RAW_VCF` | `VCF_DIR/PhasedSNPsFitTEs.vcf` — input phased SNP VCF |
| `CODE_DIR` | `BASE_DIR/code` — helper Python scripts |
| `SINGER_DIR` | `BASE_DIR/singer-0.1.9-beta-linux-x86_64` — Singer 0.1.9 release dir |
| `SINGER` | `SINGER_DIR/singer_master` — per-window ARG inference binary |
| `MERGE_ARG` | `SINGER_DIR/merge_ARG.py` — stitches per-window outputs into one `.trees` |
| `POLEGON` | Path to `polegon_master` binary |
| `OUT_DIR` | `VCF_DIR/singer_<chr>/` — per-window txt outputs (under `windows/`) and merged `.trees` |
| `REDATE_DIR` | `VCF_DIR/redate_no_TEs_<chr>/` — Polegon-redated trees |
| `TSKIT_PY` | Python in the `tskit` conda env |
| `SNPEFF_BIN` | `snpEff` binary in the `snpeff` conda env |
| `CHROMS` | Chromosome arms to process (default `["2L"]`) |

Environments are invoked via absolute binary paths after `module load mamba/24.3.0`, avoiding `conda activate` issues on SLURM.

### Euchromatic Regions

Only SNPs within euchromatic boundaries are analyzed:

| Arm | Start | End |
|---|---|---|
| 2L | 500,000 | 21,501,009 |
| 2R | 5,898,184 | 24,786,936 |
| 3L | 500,000 | 22,462,476 |
| 3R | 5,052,934 | 31,579,331 |
| X  | 500,000 | 21,659,299 |

### SINGER Parameters

```python
Ne = 1e6, mu = 5.0e-9, ratio = 2.0
N_REPLICATES = 100       # Singer -n, MCMC samples per window
THIN = 20                # Singer -thin
WINDOW_SIZE = 2_000_000  # 2 Mb per singer_master invocation
```

Windows are generated per chromosome by `singer_windows(chrom)` — contiguous 2 Mb blocks covering the euchromatic range (the final window may be shorter).

Polegon redates only post-burn-in samples: `POLEGON_SAMPLES = range(50, 100)` (the first 50 are burn-in and discarded).

### Focal Analysis Classes

Analysis runs against the TE-filtered redated trees for three focal classes:

```python
FOCAL_ANALYSIS_TYPES = ["TE", "strong", "neutral"]
```

The filter/redate pipeline has a **single** branch (TE-only removal); the focal class is only used downstream to select which sites drive the mutation-rate test.

Helper: `te_trees_input(wildcards)` returns the 50 post-burn-in TE-filtered redated tree files for the given chromosome.

## Target (`rule all`)

For each (chromosome, focal) pair, one result TSV is produced:

```
VCF_DIR/results/{focal}_mutation_rate_weight_on_TE_trees_{chr}.tsv
```

## Workflow Stages

### Step I — VCF Preprocessing

| Rule | In | Out | Purpose |
|---|---|---|---|
| [add_TEs_to_vcf](ARG.smk#L89) | `RAW_VCF`, `BED_DIR` | `PhasedSNPsFitTEs_with_TEs.vcf` | Inject TE entries into the VCF via `add_TEs_to_vcf.py`. |
| [extract_euchromatic](ARG.smk#L109) | `RAW_VCF` | `PhasedSNPsFitTEs_{chr}_euchromatic.vcf` | Keep header + euchromatic SNPs on `{chr}`. |
| [extract_euchromatic_with_TEs](ARG.smk#L134) | `..._with_TEs.vcf` | `PhasedSNPsFitTEs_{chr}_euchromatic_with_TEs.vcf` | Same slicing on the TE-augmented VCF. |
| [snpeff](ARG.smk#L159) | `..._euchromatic.vcf` | `..._euchromatic.ann.vcf` | Annotate SNPs with snpEff (`BDGP6.115`). |
| [annotate_te_vcf](ARG.smk#L178) | `_with_TEs.vcf` + `.ann.vcf` | `..._with_TEs.ann.vcf` | Copy snpEff INFO onto the TE-bearing VCF by `chr:pos`. |
| [haploidize_vcf](ARG.smk#L201) | `..._with_TEs.ann.vcf` | `..._with_TEs.ann.haploid.vcf` | Randomly pick one haplotype per sample (`seed=42`) via `haploidize_vcf.py`. |

### Step II — ARG Inference (windowed SINGER + merge)

The haploid VCF feeds into Singer 0.1.9, which runs per 2 Mb window; a separate merge step stitches the per-window outputs into one `.trees` per MCMC replicate.

- **[singer_master](ARG.smk#L224)** — per `(chr, w)`, runs `singer_master -Ne 1e6 -m 5e-9 -ratio 2.0 -n 100 -thin 20 -ploidy 1 -start <w_start> -end <w_end>` on `...haploid.vcf`. Writes `OUT_DIR{chr}/windows/singer_{chr}_{w}_{nodes,branches,muts}_{i}.txt` for `i` in `0..99` (middle `{w}` = window, trailing `{i}` = replicate).
- **[merge_ARG](ARG.smk#L267)** — per `(chr, i)`, builds a `sub_file_table` listing the per-window files for replicate `i` along with their block-start coordinates, then calls `merge_ARG.py --file_table ... --output singer_{chr}_{i}.trees`.

### Step III — TE Filter + Redate

The filter/redate chain is TE-only; SNP mutations (strong, neutral) are kept on the tree and scored later in Step IV.

- **[filter_trees](ARG.smk#L301)** — calls `remove_TE_from_trees.py --focal TE` to strip TE sites from all 100 merged tree sequences, writing to `OUT_DIR{chr}/filtered_no_TEs/`.
- **[polegon](ARG.smk#L331)** — per-sample rule that runs `polegon_master -m 5.0e-9 -num_samples 100 -thin 20` on each TE-filtered tree, producing redated trees under `REDATE_DIR{chr}/`. Wildcard `{i}` is scattered over `POLEGON_SAMPLES` (50..99).

### Step IV — Mutation Rate Analysis

A single rule ([mutation_rate_weight_on_TE_trees](ARG.smk#L360)) consumes the 50 post-burn-in TE-filtered redated trees plus the annotated haploid VCF, binning by recombination distance (`-bins 500,5000`), requiring `min_carriers=2`, and applying `-mask_neighbors`.

| Rule                               | Script                                           | Output                                                | Trees Source          |
|------------------------------------|--------------------------------------------------|-------------------------------------------------------|-----------------------|
| `mutation_rate_weight_on_TE_trees` | `test_SNP_mutation_rate_weight_on_TE_trees.py`   | `{focal}_mutation_rate_weight_on_TE_trees_{chr}.tsv`  | TE-filtered (always)  |

The rule expands over `focal ∈ {TE, strong, neutral}` to produce three TSVs per chromosome. SLURM resources: `grylee_lab` account, `standard` partition, 5-day runtime.

## DAG Summary

```
RAW_VCF ──► add_TEs_to_vcf ──► extract_euchromatic_with_TEs ──► annotate_te_vcf ──► haploidize_vcf
   │                                                                ▲                    │
   └─► extract_euchromatic ──► snpeff ─────────────────────────────┘                    │
                                                                                         ▼
                                                         (per window w)           singer_master
                                                                                         │
                                                         (per replicate i, across w)     ▼
                                                                                    merge_ARG ──► singer_{chr}_{i}.trees
                                                                                         │
                                                                                         ▼
                                                                                 filter_trees (--focal TE)
                                                                                         │
                                                                                         ▼
                                                                                    polegon_{i}
                                                                                         │
                                                                                         ▼
                                             mutation_rate_weight_on_TE_trees  (focal ∈ {TE, strong, neutral})
```

## Helper Scripts (under `CODE_DIR`)

- `add_TEs_to_vcf.py` — injects TE records from BED into the VCF.
- `haploidize_vcf.py` — emits a random haplotype per sample.
- `remove_TE_from_trees.py` — strips focal sites from tree sequences (called with `--focal TE` here).
- `test_SNP_mutation_rate_weight_on_TE_trees.py` — mutation-rate estimator on TE-filtered trees.

## Logs

All rules write to `VCF_DIR/logs/<rule>_{chr}[_{i|w}].log`.

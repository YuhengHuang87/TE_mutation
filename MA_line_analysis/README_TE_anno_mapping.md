<!-- markdownlint-disable MD024 -->

# MA-Line Analysis: Snakemake Pipelines

Two Snakemake workflows used by the mutation-accumulation (MA) line project:

- **[Pipeline A: TE_annotation.smk](#pipeline-a-te_annotationsmk)** — basecalled nanopore BAM → hifiasm assembly → RagTag scaffold → RepeatMasker → merged, length-filtered TE table.
- **[Pipeline B: TE_mapping.smk](#pipeline-b-te_mutationsmk)** — Illumina short reads from three MA groups (Act / Control / Ubi) and their ancestors → bwa mapping against the A4 HiFi reference and the nanopore `mc_initial` scaffold → joint Strelka germline VCFs per group × reference.

The two pipelines share the directory but are independent: Pipeline A's `mc_initial` scaffold feeds Pipeline B's reference set, but B can be run on its own once that scaffold exists. Both pipelines consume an externally-defined euchromatin BED (e.g. derived from CUT&Tag data on the `mc_initial` genome).

---

## Files

| File | Purpose |
|---|---|
| [TE_annotation.smk](TE_annotation.smk) | Pipeline A rules, DAG, and per-rule SLURM resources |
| [TE_mutation.smk](TE_mutation.smk) | Pipeline B rules: short-read mapping + joint Strelka calling for the MA / ancestor lines |
| [te_annotate.py](te_annotate.py) | Post-RepeatMasker TE annotation: euchromatic filter, distance, family merge, length/INE filter, class merge |
| [config.yaml](config.yaml) | Pipeline A: sample, paths, threads, mamba env prefixes, and TE annotation parameters. Pipeline B is self-contained — see the constants block at the top of [TE_mutation.smk](TE_mutation.smk). |
| [cluster.yaml](cluster.yaml) | **Legacy.** Snakemake-7 `--cluster-config` file. No longer used; per-rule resources now live in `resources:` blocks inside [TE_annotation.smk](TE_annotation.smk). Kept for reference only. |
| [hifiasm_nano.sub](hifiasm_nano.sub) | Original SLURM script (assembly path) — kept for reference |
| [repeatmasker_nano.sub](repeatmasker_nano.sub) | Original SLURM script (RepeatMasker path) — kept for reference |
| `SNP_calling_MAlines_new.sub`, `SNP_calling_strelka_MAlines_{Act,Control,Ubi}{,_nanopore}.sub` | Original SLURM scripts ported by Pipeline B — kept for reference (in the parent project directory) |

---

## Pipeline A: TE_annotation.smk

A Snakemake workflow that wraps the logic of [hifiasm_nano.sub](hifiasm_nano.sub), [repeatmasker_nano.sub](repeatmasker_nano.sub), and the post-RepeatMasker Perl chain into a single, reproducible DAG. It takes a basecalled nanopore BAM, produces a barcode-separated *de novo* assembly with [hifiasm](https://github.com/chhylp123/hifiasm), scaffolds against a reference with [RagTag](https://github.com/malonge/RagTag), annotates repeats with [RepeatMasker](https://www.repeatmasker.org/), and post-processes the RepeatMasker output into a merged, length-filtered TE table restricted to a user-supplied euchromatin BED.

### Pipeline DAG

```
  raw BAM
     │
     ▼
bam_to_fastq                 (samtools; filters by QV, emits fastq)
     │
     ▼
separate_barcode             (awk; pulls reads whose header matches {barcode})
     │
     ▼
hifiasm_assemble             (hifiasm --ont; primary contigs GFA)
     │
     ▼
gfa_to_fasta                 (awk; GFA S-lines → FASTA)
     │
     ▼
ragtag_correct               (RagTag; misassembly correction vs. reference)
     │
     ▼
ragtag_scaffold              (RagTag; reference-guided scaffolding)
     │
     ├──► assembly_stats         (assembly-stats; N50 etc.)
     │
     └──► repeatmasker ──► te_annotate
                          (te_annotate.py; euchromatic filter →
                           distance → family merge → length/INE
                           filter → DNA/RNA class merge.
                           Euchromatin BED supplied via
                           config.euchromatin_bed.)
```

Final targets declared in `rule all`:

- `results/assembly_stats.txt`
- `results/repeatmasker/ragtag.scaffold.fasta.out`
- `results/te_annotation/merged_class.txt`

### Rules

#### `bam_to_fastq`

Filters the run BAM by mean quality and converts to fastq.

- **Input:** `config.raw_bam`
- **Output:** `results/reads/suvar_nanopore.QVMin{QV_MIN}.barcode_all.fastq`
- **Tools:** `samtools view -e '[qs]>={qv}'` piped to `samtools fastq -T qs,RG`
- **Env:** `conda_envs.hifiasm`

#### `separate_barcode`

Extracts reads whose header line contains the target barcode. Uses a 4-line window (fastq record = 4 lines).

- **Input:** the QV-filtered fastq above
- **Output:** `results/reads/{barcode}_reads.fastq`
- **Wildcard:** `{barcode}` — e.g. `barcode06`
- **Env:** none (pure awk)

#### `hifiasm_assemble`

Runs hifiasm in ONT mode with purging disabled (`-l 0`), matching the original script.

- **Input:** per-barcode fastq
- **Output:** `results/hifiasm/{barcode}.asm.bp.p_ctg.gfa` (plus hifiasm's other side-files)
- **Threads:** `config.hifiasm_threads` (default 32)
- **Env:** `conda_envs.hifiasm`

#### `gfa_to_fasta`

Converts hifiasm's primary-contig GFA to FASTA.

- **Input:** the primary GFA
- **Output:** `results/hifiasm/{barcode}.asm.p_ctg.fa`
- **Env:** none (pure awk)

#### `ragtag_correct`

Corrects misassemblies in the hifiasm FASTA using the reference.

- **Input:** reference + hifiasm FASTA
- **Output:** `results/ragtag_correct/ragtag.correct.fasta`
- **Env:** `conda_envs.ragtag`

#### `ragtag_scaffold`

Scaffolds the corrected assembly against the same reference.

- **Input:** reference + corrected FASTA
- **Output:** `results/ragtag_scaffold/ragtag.scaffold.fasta`
- **Env:** `conda_envs.ragtag`

#### `assembly_stats`

Summary statistics (N50, longest contig, total length, etc.).

- **Input:** scaffolded FASTA
- **Output:** `results/assembly_stats.txt`
- **Env:** `conda_envs.assembly_stats`

#### `repeatmasker`

Soft-masks and annotates repeats using a custom TE library.

- **Input:** scaffolded FASTA + TE library FASTA
- **Output:** `results/repeatmasker/ragtag.scaffold.fasta.out` (and `.gff`, `.tbl`, `.masked`)
- **Threads:** `config.repeatmasker_threads` (default 32) passed to `-pa`
- **Engine:** `rmblast`
- **Env:** `conda_envs.repeatmasker`

#### `te_annotate`

Post-RepeatMasker TE annotation. Single Python entry point ([te_annotate.py](te_annotate.py)) that ports the original Perl chain into five sequential steps:

1. **Extract euchromatic** — keep only RepeatMasker hits inside the regions listed in `config.euchromatin_bed` and whose class/family contains one of `te_annotation.classes` (default `LINE`/`LTR`/`DNA`/`Unknown`). Cuts each kept TE's sequence out of the scaffold FASTA.
2. **Compute distances** — annotate each kept hit with the bp gap to its left/right neighbour (`NA` across chromosome boundaries).
3. **Merge nearby family** — collapse runs of adjacent same-family TEs whose gap is below `te_annotation.family_merge_dist` (default 200 bp). Mixed-family runs get a `_`-joined family/superfamily name.
4. **Filter length/INE** — drop any record whose family matches `te_annotation.exclude_family_pattern` (default regex `INE-1`) or is shorter than `te_annotation.min_length` (default 500 bp).
5. **Merge by class** — second-pass merge of nearby DNA-only or RNA-only (`LTR`/`LINE`) records within `te_annotation.class_merge_dist` (default 5000 bp). Ambiguous and Other records are left unmerged.

Replaces these scripts:

| Step | Original |
| --- | --- |
| 1 | `extract_TE_repeatmasker_nano.pl` |
| 2 | `TE_nearby_distance_nanopore.pl` |
| 3 | `TE_merge_exclude_nearby_family_combined_nanopore.pl` |
| 4 | `filter_TE_length_INE.pl` |
| 5 | `merge_te_script.py` |

- **Input:** RepeatMasker `.out`, scaffolded FASTA, euchromatic BED
- **Output:** `results/te_annotation/{filtered.out, library.fasta, distance.txt, merged_family.txt, length_filtered.txt, merged_class.txt}`
- **Env:** none (Python stdlib only)

### Configuration ([config.yaml](config.yaml))

```yaml
barcode: barcode06                 # which barcode to extract and assemble
raw_bam: suvar_nanopore_run1.bam   # basecalled reads (methylated-aware BAM)
qv_min: 10                         # minimum mean-quality filter on the BAM

reference_fasta: /dfs7/.../A4.Hifi.MT.fasta
te_library:      /dfs7/.../D_mel_transposon_sequence_set.fa

hifiasm_threads: 32
repeatmasker_threads: 32

# Per-assembly euchromatin BED on scaffold chromosome names
# (e.g. "2L_RagTag\t97997\t21400000"). Defined externally — e.g. from
# CUT&Tag data on the mc_initial genome.
euchromatin_bed: barcode06_euchromatin.bed

# Optional. Defaults match the original Perl/Python pipeline.
te_annotation:
  classes: "LINE,LTR,DNA,Unknown"      # RM class/family substrings to keep
  family_merge_dist: 200               # bp gap for step-3 same-family merge
  min_length: 500                      # bp minimum length after step 3
  exclude_family_pattern: "INE-1"      # regex of families to drop in step 4
  class_merge_dist: 5000               # bp gap for step-5 DNA/RNA-class merge

conda_envs:
  hifiasm:        /data/homezvol2/yuhenh3/.conda/envs/Hifiasm
  ragtag:         /data/homezvol2/yuhenh3/.conda/envs/RagTag
  assembly_stats: /data/homezvol2/yuhenh3/.conda/envs/assembly-stats
  repeatmasker:   /data/homezvol2/yuhenh3/.conda/envs/repeatmasker
```

**Switching barcode:** set `barcode:` to `barcode05`, `barcode07`, etc. The wildcard propagates through `separate_barcode` → `hifiasm_assemble` → `gfa_to_fasta`. Downstream RagTag/RepeatMasker rules depend on the active `{barcode}` indirectly via `rule ragtag_correct`, which picks up the current `config.barcode`.

**Conda envs:** values are *prefix paths*, not env names. Snakemake detects paths and uses the env as-is (no solve, no build). Regenerate the list with `mamba env list` if paths change.

### How to run

#### 1. Dry run (login node, always do this first)

```bash
snakemake -s TE_annotation.smk -n --use-conda
```

Prints the DAG and confirms input files resolve. Nothing runs.

#### 2. SLURM submission (recommended for real runs)

The hifiasm step needs ~32 cores for days; RepeatMasker even longer. Never run those on a login node.

Snakemake 8 dropped `--cluster` / `--cluster-config` in favor of executor plugins. This workflow uses [`snakemake-executor-plugin-slurm`](https://github.com/snakemake/snakemake-executor-plugin-slurm); install once with `pip install snakemake-executor-plugin-slurm`.

```bash
module load mamba/24.3.0
. ~/.mymambainit-24.3.0

snakemake -s TE_annotation.smk \
  --executor slurm \
  --use-conda \
  --jobs 20 \
  --default-resources \
    slurm_account=grylee_lab \
    slurm_partition=standard \
    runtime=1440
```

What changed vs. the v7 command:

- `--executor slurm` replaces `--cluster "sbatch ..."`. The plugin builds the `sbatch` invocation from `resources:` keys.
- `--default-resources` sets cluster-wide defaults (account, partition, default 1-day runtime). These were the `__default__` block in `cluster.yaml`.
- Per-rule overrides (longer runtimes for `hifiasm_assemble` and `repeatmasker`) now live in `resources:` blocks inside [TE_annotation.smk](TE_annotation.smk), not `cluster.yaml`.
- `threads:` is auto-mapped to `cpus_per_task` by the plugin, so no separate `cpus` field is needed.

Run the `snakemake` controller inside `tmux` or `screen` — it needs to stay alive to submit downstream jobs as each upstream job finishes. If SSH drops and the controller dies, already-running `sbatch` jobs keep going, but no new ones get queued.

#### 3. Local run (small test only)

```bash
snakemake --use-conda --cores 32
```

Uses whatever host you're logged into. Fine for a few-read test fastq; not for a real run.

### Per-rule SLURM resources

Defined inside [TE_annotation.smk](TE_annotation.smk) via `threads:` (→ `cpus_per_task`) and `resources: runtime=<minutes>`. Defaults (account, partition, baseline runtime) come from `--default-resources` on the command line. Mirrors the `#SBATCH` headers from the original `.sub` scripts:

| Rule | `threads:` (cpus_per_task) | `resources.runtime` (minutes) | Walltime |
| --- | --- | --- | --- |
| (default) | 1 | 1440 | 1d |
| `bam_to_fastq` | 4 | 1440 | 1d |
| `hifiasm_assemble` | 32 | 4080 | 2d20h |
| `ragtag_correct` | 8 | 1440 | 1d |
| `ragtag_scaffold` | 8 | 1440 | 1d |
| `repeatmasker` | 32 | 15600 | 10d20h |
| `te_annotate` | 1 | 1440 | 1d |

Edit the `resources:` blocks in [TE_annotation.smk](TE_annotation.smk) if your cluster limits differ from UCI HPC3. Account/partition defaults can be overridden per-run by changing the `--default-resources` flag.

### Outputs

```
results/
├── reads/
│   ├── suvar_nanopore.QVMin10.barcode_all.fastq
│   └── {barcode}_reads.fastq
├── hifiasm/
│   ├── {barcode}.asm.bp.p_ctg.gfa   (+ hifiasm side files)
│   └── {barcode}.asm.p_ctg.fa
├── ragtag_correct/
│   └── ragtag.correct.fasta
├── ragtag_scaffold/
│   └── ragtag.scaffold.fasta
├── repeatmasker/
│   ├── ragtag.scaffold.fasta.out
│   ├── ragtag.scaffold.fasta.masked
│   ├── ragtag.scaffold.fasta.tbl
│   └── ragtag.scaffold.fasta.out.gff
├── te_annotation/
│   ├── filtered.out            (step 1: euchromatic-restricted RM rows)
│   ├── library.fasta           (step 1: TE sequences cut from scaffold)
│   ├── distance.txt            (step 2: rows + dist_prev + dist_next)
│   ├── merged_family.txt       (step 3: same-family merge, 8 cols)
│   ├── length_filtered.txt     (step 4: INE-1 dropped, ≥ min_length)
│   └── merged_class.txt        (step 5: DNA/RNA-class merge, 5 cols)
└── assembly_stats.txt
```

### Differences from the original `.sub` scripts

- **Intermediate compression removed logically but preserved.** The original `samtools fastq | gzip > f.gz && gunzip f.gz` dance is kept verbatim since it doesn't affect correctness. Safe to simplify to a single uncompressed stream later.
- **Barcode is a wildcard, not a hardcoded constant.** The original awk pattern `/barcode06/` is now `$0 ~ bc` with `bc` bound from the rule wildcard. Change `config.barcode` to retarget.
- **RagTag output directories are explicit.** The original script relied on default `ragtag_output/` subfolders being created in the CWD; this Snakefile passes `-o` so outputs land under `results/`.
- **RepeatMasker output is redirected to `results/repeatmasker/` via `-dir`.** Original script wrote next to the input FASTA.
- **Module loads (`module load samtools/1.15.1`) are replaced by conda envs.** The `hifiasm` env bundles both `hifiasm` and `samtools`, so no `module load` inside Snakemake rules.
- **Post-RepeatMasker Perl chain is now one Python script.** `extract_TE_repeatmasker_nano.pl` → `TE_nearby_distance_nanopore.pl` → `TE_merge_exclude_nearby_family_combined_nanopore.pl` → `filter_TE_length_INE.pl` → `merge_te_script.py` are folded into [te_annotate.py](te_annotate.py), driven by the `te_annotate` rule. The hard-coded `mc`-strain euchromatic boundaries are replaced by a user-supplied BED (`config.euchromatin_bed`) so the rule generalizes across barcodes. The "Harsh approach" merge variant was not ported.

### Troubleshooting

- **`Missing input files for rule bam_to_fastq: suvar_nanopore_run1.bam`** — the BAM path in `config.raw_bam` is relative; run Snakemake from the directory containing it, or set an absolute path.
- **`RepeatMasker: no TE library found at ...`** — `config.te_library` path is a cluster path (`/dfs7/...`); confirm it's mounted on whatever node runs the `repeatmasker` rule.
- **Conda env not activating under `--use-conda`** — Snakemake expects `conda` on `PATH`. Your HPC has it, so this should just work. If you switch to a machine where only `mamba` is available, drop `--use-conda` and source the env inside each rule's shell.
- **RagTag rule re-runs every invocation** — the output directory already exists. RagTag refuses by default; the `-o` flag handles this, but if you manually seeded the dir, delete it first.
- **`te_annotate` produces an empty `library.fasta` or `merged_class.txt`** — most often the chromosome names in `config.euchromatin_bed` don't match the scaffold's headers (e.g. you wrote `2L` but RagTag emits `2L_RagTag`). Confirm with `grep '^>' results/ragtag_scaffold/ragtag.scaffold.fasta`.

---

## Pipeline B: TE_mutation.smk

A Snakemake port of [TE_mutation_calling.md](../../TE_mutation_calling.md) and the five SLURM scripts it documents (`SNP_calling_MAlines_new.sub`, `SNP_calling_strelka_MAlines_{Act,Control,Ubi}{,_nanopore}.sub`). Each MA-line and ancestor sample is mapped against both the A4 HiFi reference and the nanopore `mc_initial` scaffold (gunzip → trim_galore → bwa mem | samtools sort → MarkDuplicates → AddOrReplaceReadGroups), then a single Strelka germline workflow is run per (group, reference) pair, including the group's ancestor BAM and restricted to the matching euchromatin BED.

### Pipeline DAG

Per sample × reference (mapping path):

```
{sample}-READ{1,2}-Sequences.txt.gz   (raw paired-end fastq.gz)
     │
     ▼
gunzip_reads
     │
     ▼
trim_galore                           (paired, --fastqc)
     │
     ▼
bwa_index                             (one-shot per reference: bwa index +
     │                                 samtools faidx + picard dict)
     ▼
map_and_dedup                         (bwa mem | samtools sort piped into
     │                                 gatk MarkDuplicates; sorted BAM lives
     │                                 in $TMPDIR, never declared)
     ▼
add_read_groups                       (gatk AddOrReplaceReadGroups,
     │                                 RGSM={sample}_Sample, CREATE_INDEX=True)
     ▼
bam/{sample}_{ref}_reads_with_RG.bam (+ .bai)
```

Per group × reference (joint-call path):

```
ancestor BAM + N MA-line BAMs (all _{ref}_reads_with_RG.bam)
     │
     ▼
strelka_configure                     (configureStrelkaGermlineWorkflow.py
     │                                 --bam ... --referenceFasta ...
     │                                 --runDir ... --callRegions {bed}.gz)
     ▼
strelka_run                           (./runWorkflow.py -m local -j {threads},
     │                                 then gunzip -c variants.vcf.gz > variants.vcf)
     ▼
strelka/{group}_{ref}/results/variants/variants.vcf
```

Default `rule all` produces six decompressed joint-call VCFs (3 groups × 2 references):

| Group | Reference | VCF |
| --- | --- | --- |
| Act | A4 | `strelka/Act_A4/results/variants/variants.vcf` |
| Act | mc_initial | `strelka/Act_mc_initial/results/variants/variants.vcf` |
| Control | A4 | `strelka/Control_A4/results/variants/variants.vcf` |
| Control | mc_initial | `strelka/Control_mc_initial/results/variants/variants.vcf` |
| Ubi | A4 | `strelka/Ubi_A4/results/variants/variants.vcf` |
| Ubi | mc_initial | `strelka/Ubi_mc_initial/results/variants/variants.vcf` |

The bgzipped `variants.vcf.gz` is preserved as an upstream output of `strelka_run`; `decompress_vcf` adds a plain text copy for direct inspection.

Two convenience targets are also defined:

- `mapping_all` — every `bam/{sample}_{ref}_reads_with_RG.bam` (step 1 only, no calling)
- `strelka_all` — same set as `all`

### Rules

#### `gunzip_reads`

Decompresses the paired-end fastq.gz files into the local `reads/` directory. Mirrors the commented-out `gunzip ${id}-READ{1,2}-Sequences.txt.gz` lines in `SNP_calling_MAlines_new.sub`.

- **Input:** `{SAMPLE_DIR[sample]}/{sample}-READ{1,2}-Sequences.txt.gz`
- **Output:** `reads/{sample}-READ{1,2}-Sequences.txt`

#### `trim_galore`

Adapter / quality trimming with Trim Galore in paired mode plus FastQC.

- **Input:** the two uncompressed fastq files
- **Output:** `trim/{sample}-READ{1,2}-Sequences.txt_val_{1,2}.fq` (Trim Galore's default suffix; matches the original `.sub`)
- **Threads:** 4

#### `bwa_index`

One-time indexing for each reference. Touch flag `ref/{ref}.indexed` gates the mapping rule.

- **Input:** `REFERENCES[ref].fasta`
- **Output:** `ref/{ref}.indexed`
- **Tools:** `bwa index`, `samtools faidx`, `picard CreateSequenceDictionary`

#### `map_and_dedup`

Maps trimmed reads with `bwa mem`, coordinate-sorts the alignment with `samtools sort`, and removes duplicates with `gatk MarkDuplicates` — all inside one rule. The intermediate sorted BAM is held in a `mktemp` file under `$TMPDIR` and deleted on exit, so it never appears in the workflow DAG. This is what lets re-runs short-circuit when `marked_duplicates.bam` already exists on disk: there is no missing `output_*_sorted.bam` to walk back through.

- **Input:** trimmed fastq + `ref/{ref}.indexed`
- **Output:** `bam/{sample}_{ref}_marked_duplicates.bam` (declared `temp()`) + `bam/{sample}_{ref}_marked_dup_metrics.txt` (kept)
- **Threads:** 16
- **Command:** `bwa mem | samtools sort -o $tmp; gatk MarkDuplicates -I $tmp -O ...`

#### `add_read_groups`

Adds a coordinate-sorted, indexed BAM with read-group metadata. RGSM is `{sample}_Sample` (matches the original).

- **Input:** dedup BAM
- **Output:** `bam/{sample}_{ref}_reads_with_RG.bam` (+ `.bai`)

#### `strelka_configure`

Builds a Strelka germline run directory. Pulls every BAM in the group (ancestor + MA lines) into a single configure call, restricted to the reference's euchromatin BED.

- **Input:** every `bam/{member}_{ref}_reads_with_RG.bam` for `member ∈ [ancestor] + GROUPS[group].samples`
- **Output:** `strelka/{group}_{ref}/runWorkflow.py`
- **Notes:** the run directory is removed (`rm -rf`) before configure, since `configureStrelkaGermlineWorkflow.py` refuses to reuse one.

#### `strelka_run`

Executes the Strelka workflow locally with `-j {threads}`, then decompresses the bgzipped output into a plain VCF. Only the decompressed `variants.vcf` is declared as the rule's output, so once it exists Snakemake treats the rule as satisfied even if `variants.vcf.gz` was later deleted.

- **Input:** the configured `runWorkflow.py`
- **Output:** `strelka/{group}_{ref}/results/variants/variants.vcf` (the bgzipped `variants.vcf.gz` is left on disk by Strelka, but is not a declared output)
- **Threads:** 32

### Configuration

Pipeline B is self-contained — no `config.yaml`. Edit the constants block at the top of [TE_mutation.smk](TE_mutation.smk):

| Variable | Purpose |
| --- | --- |
| `WGS_ROOT` | MA-line WGS root (per-group `final_samples/{group}/raw_reads`) |
| `ANC_ROOT` | Ancestor Illumina-WGS directory |
| `NANO_REF_DIR` | Nanopore reference (`mc_initial_ragtag.scaffold.fasta`) directory |
| `OUTPUT_ROOT` | Where the workflow reads and writes per-sample / per-group outputs. Defaults to `{WGS_ROOT}/final_samples`, which is where the existing `bam/` and `strelka/` directories already live, so re-running picks up prior results without copying. |
| `BAM_DIR`, `STRELKA_DIR`, `READS_DIR`, `TRIM_DIR`, `REF_FLAG_DIR` | Subdirs of `OUTPUT_ROOT`: `{OUTPUT_ROOT}/bam`, `{OUTPUT_ROOT}/strelka`, etc. Override individually if a particular stage needs to land elsewhere. |
| `REFERENCES` | dict `name → {fasta, call_regions}`. Default: `A4` (HiFi) + `mc_initial` (nanopore scaffold). The `call_regions` BED is the bgzip-tabix-indexed euchromatin region that Strelka restricts calling to. |
| `GROUPS` | dict `name → {reads_dir, ancestor, samples[]}`. Default: `Act`, `Control`, `Ubi`. Sample lists were extracted verbatim from the four `SNP_calling_strelka_MAlines_*.sub` files (47–52 lines per group). |

To restrict a run, either edit `REFERENCES` / `GROUPS` to drop entries, or invoke a specific target on the command line:

```bash
# only Act vs A4 (path is OUTPUT_ROOT-relative; substitute the actual value)
snakemake -s TE_mutation.smk --cores 32 \
  /dfs7/grylee/yuhenh3/mutation_accumulation/MA_lines_WGS/final_samples/strelka/Act_A4/results/variants/variants.vcf

# all six VCFs (default)
snakemake -s TE_mutation.smk --cores 32

# stop after mapping (skip Strelka)
snakemake -s TE_mutation.smk --cores 32 mapping_all
```

### How to run

Tools that must be on `PATH` when Snakemake runs:

- `gunzip` (system)
- `trim_galore`, `bwa`, `samtools`, `picard`, `gatk` (4.x) — UCI HPC3: `module load trimgalore/0.6.6 bwa/0.7.17 samtools/1.15.1 gatk/4.2.6.1 picard-tools/2.27.1`
- `configureStrelkaGermlineWorkflow.py`, `runWorkflow.py` — installed via the `Strelka` mamba env in the original scripts (`mamba activate Strelka`)

#### 1. Dry run

```bash
snakemake -s TE_mutation.smk -n
```

Confirms the DAG and that every input fastq.gz resolves to a real path under `WGS_ROOT` / `ANC_ROOT`.

#### 2. Local run

```bash
module load trimgalore/0.6.6 bwa/0.7.17 samtools/1.15.1 gatk/4.2.6.1 picard-tools/2.27.1
module load mamba/24.3.0 && . ~/.mymambainit-24.3.0 && mamba activate Strelka

snakemake -s TE_mutation.smk --cores 32
```

Each Strelka job uses 32 threads and the per-sample mapping rule uses 16, so `--cores 32` lets one Strelka or two mapping jobs run at a time.

#### 3. SLURM submission

The Snakemake-8 `--executor slurm` plugin builds `sbatch` invocations from each rule's `threads:` and any `resources:`. Pipeline B does not (yet) declare per-rule `resources: runtime=...`, so cap the runtime via `--default-resources`:

```bash
snakemake -s TE_mutation.smk \
  --executor slurm \
  --jobs 30 \
  --default-resources \
    slurm_account=grylee_lab \
    slurm_partition=standard \
    runtime=$((5*24*60))     # 5 days; matches the original .sub headers
```

Run the controller in `tmux` / `screen` — the longest job (`bwa_mem_sort` for a deep sample, or a 32-way `strelka_run`) takes hours to days.

### Outputs

```text
reads/                                        # temp(): deleted once trim_galore consumes them
└── {sample}-READ{1,2}-Sequences.txt          (decompressed fastq)
trim/                                         # temp(): deleted once bwa_mem_sort consumes them
├── {sample}-READ1-Sequences.txt_val_1.fq     (Trim Galore output)
├── {sample}-READ2-Sequences.txt_val_2.fq
└── {sample}-READ{1,2}-Sequences.txt_trimming_report.txt
ref/
└── {ref}.indexed                             (touch flag; bwa/.fai/.dict live next to the source FASTA)
bam/
├── {sample}_{ref}_marked_duplicates.bam      # temp(): deleted once add_read_groups consumes it
├── {sample}_{ref}_marked_dup_metrics.txt     (kept; MarkDuplicates QC report)
├── {sample}_{ref}_reads_with_RG.bam          (FINAL per-sample BAM)
└── {sample}_{ref}_reads_with_RG.bai
# Note: the bwa-mem-sorted BAM lives in $TMPDIR for the duration of map_and_dedup
# and is deleted by the rule's shell trap; it is never declared as a workflow file.
strelka/{group}_{ref}/
├── runWorkflow.py                            (per-group Strelka run dir)
├── runWorkflow.py.config.pickle
├── results/
│   ├── variants/
│   │   ├── variants.vcf                      (FINAL — declared output of strelka_run)
│   │   ├── variants.vcf.gz                   (Strelka's native bgzipped output; not declared, safe to delete)
│   │   └── variants.vcf.gz.tbi
│   └── stats/
└── workspace/                                (Strelka scratch)
```

**Re-runs do not regenerate intermediates.** `reads/`, `trim/`, the sorted BAM, and the dedup BAM are wrapped in [`temp()`](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#protected-and-temporary-files), so Snakemake auto-deletes them once the consuming job completes. On a subsequent invocation, if the per-sample `*_reads_with_RG.bam` and the per-group `variants.vcf` are already up-to-date, Snakemake walks the DAG, sees the finals are satisfied, and short-circuits — the temp intermediates are never recreated. To force a full rebuild for a sample, delete its `*_reads_with_RG.bam` (and downstream VCFs); to invalidate just the calling step, delete `strelka/{group}_{ref}/` instead.

### Differences from the original `.sub` scripts

- **Per-sample loops are replaced by Snakemake wildcards.** The hard-coded `for index in 23` driving `SNP_calling_MAlines_new.sub` is gone; every MA-line and ancestor sample fans out as a `{sample}` wildcard against both `{ref}` values. `--cores` controls parallelism.
- **The four hand-maintained Strelka sub-scripts collapse to one rule.** `SNP_calling_strelka_MAlines_{Act,Control,Ubi}{,_nanopore}.sub` — six files, ~600 lines total — become `strelka_configure` + `strelka_run` with the BAM lists driven from `GROUPS[group].samples`. To add or drop a sample, edit one Python list in [TE_mutation.smk](TE_mutation.smk).
- **Both references are run from the same workflow.** The original split A4 mapping and `mc_initial` (nanopore) mapping into separate scripts; `REFERENCES` here lets the same rules cover both, distinguished only by the `{ref}` wildcard.
- **`module load` lines are not embedded in rules.** Load the modules in the launching shell (or wrap the rules in `conda:` envs later). Mirrors how Pipeline A originally evolved before `conda_envs` were introduced.
- **Strelka run-dir collisions are handled.** The configure rule runs `rm -rf {run_dir}` before each invocation; `configureStrelkaGermlineWorkflow.py` otherwise refuses to overwrite.
- **Optional/commented-out steps from the original script are not ported.** `samtools depth`, `gatk HaplotypeCaller`, the `bcftools isec` ancestor/sample intersection, and the `--SORTING_COLLECTION_SIZE_RATIO 0.1` MarkDuplicates tweak (`SNP_calling_MAlines_new.sub` lines 61–68 and 65) are absent. Re-add to the corresponding rule's shell if a sample needs them.

### Troubleshooting

- **`Missing input files: {sample}-READ1-Sequences.txt.gz`** — the sample isn't in any `GROUPS[*].samples` list, so it has no entry in `SAMPLE_DIR` and no resolved reads dir. Add it, or use one of the existing samples.
- **`bwa_index` re-runs every invocation** — `ref/{ref}.indexed` is a touch flag; deleting the reference's `.bwt` does *not* invalidate it. `rm ref/{ref}.indexed` to force a re-index.
- **`picard CreateSequenceDictionary` errors that the `.dict` already exists** — picard refuses to overwrite. Either delete the `.dict` first, or wrap the rule's call with `[ -f {dict} ] || picard ...` if you re-run often.
- **Strelka `configureStrelkaGermlineWorkflow.py: error: Run directory already exists`** — the `rm -rf {run_dir}` fallback in `strelka_configure` should prevent this; if it still fires (e.g. permissions), delete `strelka/{group}_{ref}/` manually.
- **`mark_duplicates` runs out of memory on a deep sample** — re-add `--SORTING_COLLECTION_SIZE_RATIO 0.1` (used in `SNP_calling_MAlines_new.sub` for the Mahul-mapped path) to the rule's shell, or bump `resources: mem_mb=...` and the SLURM allocation.
- **`callRegions` BED has no contigs in common with the reference** — `A4_hifi_euchromatin_region.bed.gz` uses the A4 chromosome names, `mc_nanopore_euchromatin_region.bed.gz` uses the `*_RagTag` scaffold names. Crossing them silently produces an empty Strelka output. Confirm with `zcat {bed.gz} | head` against `samtools view -H` of any RG BAM.
- **An ancestor or sample appears under both groups by mistake** — `SAMPLE_DIR` is a flat dict, so a duplicate sample ID across groups overwrites silently and the reads will be looked up in whichever group came last. The default lists are disjoint; keep them that way.

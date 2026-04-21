# ARG analysis workflow
#

# NOTE: update the paths in add_TEs_to_vcf.py to match HPC paths before running.

# ── Configuration ─────────────────────────────────────────────────────────────

BASE_DIR   = "/dfs7/grylee/yuhenh3/mutation_accumulation/ARG"
VCF_DIR    = BASE_DIR + "/Rech_2022_data/SNPs_vcf"
BED_DIR    = BASE_DIR + "/Rech_2022_data/TE_annotations/ReferenceCoordinates"
RAW_VCF    = VCF_DIR  + "/PhasedSNPsFitTEs.vcf"
CODE_DIR   = BASE_DIR + "/code"

# Singer 0.1.9: singer_master runs per window; merge_ARG.py stitches windows into a .trees file
SINGER_DIR = BASE_DIR + "/singer-0.1.9-beta-linux-x86_64"
SINGER     = SINGER_DIR + "/singer_master"
MERGE_ARG  = SINGER_DIR + "/merge_ARG.py"

POLEGON  = BASE_DIR + "/polegon-0.1.3-alpha-linux-x86_64/polegon_master"

# Absolute paths so remove_TE_from_trees.py's hardcoded TREES_DIR stays consistent
OUT_DIR    = VCF_DIR + "/singer_"                    # e.g. .../SNPs_vcf/singer_2L/
REDATE_DIR = VCF_DIR + "/redate_no_TEs_"             # e.g. .../redate_no_TEs_2L/

# Use Python/binaries directly from the env — avoids conda/mamba activate issues
# $MAMBA_ROOT_PREFIX is set automatically after: module load mamba/24.3.0
TSKIT_PY   = "$HOME/.conda/envs/tskit/bin/python3"
SNPEFF_BIN = "$HOME/.conda/envs/snpeff/bin/snpEff"

CHROMS = ["2L"]   # extend as needed: ["2L", "2R", "3L", "3R", "X"]

# Euchromatic boundaries per chromosome arm
EUCHROMATIC = {
    "2L": (500000,  21501009),
    "2R": (5898184, 24786936),
    "3L": (500000,  22462476),
    "3R": (5052934, 31579331),
    "X":  (500000,  21659299),
}

SINGER_PARAMS = dict(
    Ne    = "1e6",
    mu    = "5.0e-9",
    ratio = "2.0",
)
N_REPLICATES = 100                 # MCMC samples per window (Singer's -n)
THIN         = 20                  # Singer's -thin
WINDOW_SIZE  = 5_000_000           # singer_master runs per 5 Mb segment

def singer_windows(chrom):
    """Non-overlapping 5 Mb windows covering the euchromatic range of `chrom`.

    Returns a list of (window_index, start, end) tuples.
    """
    lo, hi = EUCHROMATIC[chrom]
    out, idx, start = [], 0, lo
    while start < hi:
        end = min(start + WINDOW_SIZE, hi)
        out.append((idx, start, end))
        idx += 1
        start = end
    return out

WINDOWS = {chrom: singer_windows(chrom) for chrom in CHROMS}

wildcard_constraints:
    chr = "|".join(["2L", "2R", "3L", "3R", "X"]),
    i   = r"\d+",
    w   = r"\d+",

# Polegon only redates post-burn-in samples (first 50 discarded as burn-in)
POLEGON_SAMPLES = range(50, 100)

# Focal mutation types for mutation-rate analysis (analysis-only; filter path is TE-only)
FOCAL_ANALYSIS_TYPES = ["TE", "strong", "neutral"]

def te_trees_input(wildcards):
    """TE-filtered Polegon tree samples — used by mutation_rate_weight_on_TE_trees."""
    return [
        f"{REDATE_DIR}{wildcards.chr}/redate_singer_{wildcards.chr}_{i}.trees"
        for i in POLEGON_SAMPLES
    ]

# ── Target rule ───────────────────────────────────────────────────────────────

rule all:
    input:
        expand(
            VCF_DIR + "/results/{focal}_mutation_rate_weight_on_TE_trees_{chr}.tsv",
            chr=CHROMS,
            focal=FOCAL_ANALYSIS_TYPES,
        ),

# ── Step I-a: add TE annotations to VCF ──────────────────────────────────────

rule add_TEs_to_vcf:
    input:
        vcf = RAW_VCF,
        bed = BED_DIR,
    output:
        vcf = VCF_DIR + "/PhasedSNPsFitTEs_with_TEs.vcf",
    resources:
        slurm_account   = "grylee_lab",
        slurm_partition = "standard",
        runtime         = 120,
    log:
        VCF_DIR + "/logs/add_TEs_to_vcf.log"
    shell:
        """
        module load mamba/24.3.0
        {TSKIT_PY} {CODE_DIR}/add_TEs_to_vcf.py > {log} 2>&1
        """

# ── Step I-b: extract euchromatic SNPs per chromosome ─────────────────────────

rule extract_euchromatic:
    input:
        vcf = RAW_VCF,
    output:
        vcf = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic.vcf",
    params:
        lo = lambda wc: EUCHROMATIC[wc.chr][0],
        hi = lambda wc: EUCHROMATIC[wc.chr][1],
    resources:
        slurm_account   = "grylee_lab",
        slurm_partition = "standard",
        runtime         = 30,
    log:
        VCF_DIR + "/logs/extract_euchromatic_{chr}.log"
    shell:
        """
        awk '/^[^#]/{{exit}} {{print}}' {input.vcf} \
          | awk -v chrom="{wildcards.chr}" \
              '/^##contig=/ {{ if ($0 ~ "ID=" chrom ">") print; next }} {{ print }}' \
          > {output.vcf}

        awk -v chrom="{wildcards.chr}" -v lo={params.lo} -v hi={params.hi} \
            '$1==chrom && $2+0>=lo && $2+0<=hi' {input.vcf} >> {output.vcf}
        """

rule extract_euchromatic_with_TEs:
    input:
        vcf = VCF_DIR + "/PhasedSNPsFitTEs_with_TEs.vcf",
    output:
        vcf = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic_with_TEs.vcf",
    params:
        lo = lambda wc: EUCHROMATIC[wc.chr][0],
        hi = lambda wc: EUCHROMATIC[wc.chr][1],
    resources:
        slurm_account   = "grylee_lab",
        slurm_partition = "standard",
        runtime         = 30,
    log:
        VCF_DIR + "/logs/extract_euchromatic_with_TEs_{chr}.log"
    shell:
        """
        awk '/^[^#]/{{exit}} {{print}}' {input.vcf} \
          | awk -v chrom="{wildcards.chr}" \
              '/^##contig=/ {{ if ($0 ~ "ID=" chrom ">") print; next }} {{ print }}' \
          > {output.vcf}

        awk -v chrom="{wildcards.chr}" -v lo={params.lo} -v hi={params.hi} \
            '$1==chrom && $2+0>=lo && $2+0<=hi' {input.vcf} >> {output.vcf}
        """

# ── Step I-c: snpEff annotation per chromosome ────────────────────────────────

rule snpeff:
    input:
        vcf = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic.vcf",
    output:
        vcf = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic.ann.vcf",
    resources:
        slurm_account   = "grylee_lab",
        slurm_partition = "standard",
        runtime         = 7200,   # minutes (= 5 days)
    log:
        VCF_DIR + "/logs/snpeff_{chr}.log"
    shell:
        """
        module load mamba/24.3.0
        {SNPEFF_BIN} BDGP6.115 {input.vcf} > {output.vcf} 2> {log}
        """

# ── Step I-d: add snpEff annotations to the with-TEs VCF ─────────────────────

rule annotate_te_vcf:
    input:
        te_vcf  = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic_with_TEs.vcf",
        ann_vcf = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic.ann.vcf",
    output:
        vcf = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic_with_TEs.ann.vcf",
    resources:
        slurm_account   = "grylee_lab",
        slurm_partition = "standard",
        runtime         = 30,
    log:
        VCF_DIR + "/logs/annotate_te_vcf_{chr}.log"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}} \
            NR==FNR {{ if ($0 !~ /^#/) info[$1":"$2] = $8; next }} \
            /^#/ {{ print; next }} \
            {{ key = $1":"$2; if (key in info) $8 = info[key]; print }}' \
            {input.ann_vcf} {input.te_vcf} > {output.vcf} 2> {log}
        """

# ── Step I-e: haploidize (one randomly chosen haplotype per sample) ──────────

rule haploidize_vcf:
    input:
        vcf = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic_with_TEs.ann.vcf",
    output:
        vcf = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic_with_TEs.ann.haploid.vcf",
    params:
        seed = 42,
    resources:
        slurm_account   = "grylee_lab",
        slurm_partition = "standard",
        runtime         = 30,
    log:
        VCF_DIR + "/logs/haploidize_vcf_{chr}.log"
    shell:
        """
        {TSKIT_PY} {CODE_DIR}/haploidize_vcf.py \
            --input {input.vcf} \
            --output {output.vcf} \
            --seed {params.seed} 2> {log}
        """

# ── Step II-a: singer_master per 5 Mb window ─────────────────────────────────
#
# Each invocation writes <outprefix>_nodes_<i>.txt (and _branches_, _muts_) for
# i in 0..N_REPLICATES-1. With outprefix = "singer_{chr}_{w}", the full filename
# is singer_{chr}_{w}_nodes_{i}.txt where w = window index, i = MCMC replicate.

rule singer_master:
    input:
        vcf = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic_with_TEs.ann.haploid.vcf",
    output:
        files = expand(
            OUT_DIR + "{{chr}}/windows/singer_{{chr}}_{{w}}_{kind}_{i}.txt",
            kind=["nodes", "branches", "muts"],
            i=range(N_REPLICATES),
        ),
    params:
        outprefix  = lambda wc: f"{OUT_DIR}{wc.chr}/windows/singer_{wc.chr}_{wc.w}",
        vcf_prefix = lambda wc, input: input.vcf[:-4],   # singer_master -vcf expects prefix (no .vcf)
        start      = lambda wc: WINDOWS[wc.chr][int(wc.w)][1],
        end        = lambda wc: WINDOWS[wc.chr][int(wc.w)][2],
        Ne         = SINGER_PARAMS["Ne"],
        mu         = SINGER_PARAMS["mu"],
        ratio      = SINGER_PARAMS["ratio"],
        n          = N_REPLICATES,
        thin       = THIN,
        ploidy     = 1,
    resources:
        slurm_account   = "grylee_lab",
        slurm_partition = "standard",
        nodes           = 1,
        tasks           = 1,
        mem_mb          = 32000,  # 32 GB — dense 2 Mb windows can need >default
        runtime         = 2880,   # minutes (= 2 days)
    log:
        VCF_DIR + "/logs/singer_master_{chr}_{w}.log"
    shell:
        """
        module load mamba/24.3.0
        mkdir -p $(dirname {params.outprefix})
        {TSKIT_PY} {SINGER} \
            -Ne {params.Ne} \
            -m  {params.mu} \
            -ratio {params.ratio} \
            -vcf {params.vcf_prefix} \
            -output {params.outprefix} \
            -start {params.start} -end {params.end} \
            -n {params.n} \
            -thin {params.thin} \
            -ploidy {params.ploidy} \
        > {log} 2>&1
        """

# ── Step II-b: merge windows into one .trees per replicate ───────────────────
#
# For replicate i, collect the per-window files and build a sub_file_table with
# columns: <nodes_file> <branches_file> <muts_file> <block_start>.

rule merge_ARG:
    input:
        nodes    = lambda wc: [
            f"{OUT_DIR}{wc.chr}/windows/singer_{wc.chr}_{w}_nodes_{wc.i}.txt"
            for w, _, _ in WINDOWS[wc.chr]
        ],
        branches = lambda wc: [
            f"{OUT_DIR}{wc.chr}/windows/singer_{wc.chr}_{w}_branches_{wc.i}.txt"
            for w, _, _ in WINDOWS[wc.chr]
        ],
        muts     = lambda wc: [
            f"{OUT_DIR}{wc.chr}/windows/singer_{wc.chr}_{w}_muts_{wc.i}.txt"
            for w, _, _ in WINDOWS[wc.chr]
        ],
    output:
        trees      = OUT_DIR + "{chr}/singer_{chr}_{i}.trees",
        file_table = temp(OUT_DIR + "{chr}/file_tables/merge_{chr}_{i}.txt"),
    params:
        table_content = lambda wc: "".join(
            f"{OUT_DIR}{wc.chr}/windows/singer_{wc.chr}_{w}_nodes_{wc.i}.txt "
            f"{OUT_DIR}{wc.chr}/windows/singer_{wc.chr}_{w}_branches_{wc.i}.txt "
            f"{OUT_DIR}{wc.chr}/windows/singer_{wc.chr}_{w}_muts_{wc.i}.txt "
            f"{start}\n"
            for w, start, _ in WINDOWS[wc.chr]
        ),
    resources:
        slurm_account   = "grylee_lab",
        slurm_partition = "standard",
        runtime         = 60,
    log:
        VCF_DIR + "/logs/merge_ARG_{chr}_{i}.log"
    shell:
        """
        module load mamba/24.3.0
        mkdir -p $(dirname {output.file_table})
        printf '%s' "{params.table_content}" > {output.file_table}
        {TSKIT_PY} {MERGE_ARG} \
            --file_table {output.file_table} \
            --output {output.trees} \
        > {log} 2>&1
        """

# ── Step III-a: strip TE sites from each merged tree sequence ────────────────

rule filter_trees:
    input:
        trees   = expand(
            OUT_DIR + "{{chr}}/singer_{{chr}}_{i}.trees",
            i=range(N_REPLICATES),
        ),
        vcf_tes = VCF_DIR + "/PhasedSNPsFitTEs_with_TEs.vcf",
    output:
        trees = expand(
            OUT_DIR + "{{chr}}/filtered_no_TEs/singer_{{chr}}_{i}.trees",
            i=range(N_REPLICATES),
        ),
    params:
        output_dir = OUT_DIR + "{chr}/filtered_no_TEs",
    resources:
        slurm_account   = "grylee_lab",
        slurm_partition = "standard",
        runtime         = 240,
    log:
        VCF_DIR + "/logs/filter_trees_{chr}.log"
    shell:
        """
        module load mamba/24.3.0
        {TSKIT_PY} {CODE_DIR}/remove_TE_from_trees.py \
            --chr {wildcards.chr} \
            --focal TE \
            --trees_dir {OUT_DIR}{wildcards.chr} \
            --output_dir {params.output_dir} \
            --te_vcf {input.vcf_tes} \
        > {log} 2>&1
        """

# ── Step III-b: redate TE-filtered trees with Polegon ────────────────────────

rule polegon:
    input:
        trees = OUT_DIR + "{chr}/filtered_no_TEs/singer_{chr}_{i}.trees",
    output:
        trees = REDATE_DIR + "{chr}/redate_singer_{chr}_{i}.trees",
    params:
        inprefix  = OUT_DIR + "{chr}/filtered_no_TEs/singer_{chr}_{i}",
        outprefix = REDATE_DIR + "{chr}/redate_singer_{chr}_{i}",
    threads: 2
    resources:
        slurm_account   = "grylee_lab",
        slurm_partition = "standard",
        runtime         = 4320,   # minutes (= 3 days)
    log:
        VCF_DIR + "/logs/polegon_{chr}_{i}.log"
    shell:
        """
        module load mamba/24.3.0
        mkdir -p {REDATE_DIR}{wildcards.chr}
        {TSKIT_PY} {POLEGON} \
            -m 5.0e-9 \
            -input {params.inprefix} \
            -output {params.outprefix} \
            -num_samples 100 \
            -thin 20 \
        > {log} 2>&1
        """

# ── Step IV: mutation rate analysis on TE-filtered trees ─────────────────────
#
# A single rule, run per focal class (TE / strong / neutral) against the
# TE-filtered redated trees.

rule mutation_rate_weight_on_TE_trees:
    input:
        trees = te_trees_input,
        vcf   = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic_with_TEs.ann.haploid.vcf",
    output:
        tsv = VCF_DIR + "/results/{focal}_mutation_rate_weight_on_TE_trees_{chr}.tsv",
    params:
        chrom = "{chr}",
        focal = "{focal}",
    resources:
        slurm_account   = "grylee_lab",
        slurm_partition = "standard",
        runtime         = 7200,   # minutes (= 5 days)
    log:
        VCF_DIR + "/logs/{focal}_mutation_rate_weight_on_TE_trees_{chr}.log"
    shell:
        """
        module load mamba/24.3.0
        mkdir -p $(dirname {output.tsv})
        {TSKIT_PY} {CODE_DIR}/test_SNP_mutation_rate_weight_on_TE_trees.py \
            -trees {input.trees} \
            -vcf {input.vcf} \
            -chrom {params.chrom} \
            -focal {params.focal} \
            -bins 500,5000 \
            -min_carriers 2 \
            -mask_neighbors \
            -output {output.tsv} \
        > {log} 2>&1
        """

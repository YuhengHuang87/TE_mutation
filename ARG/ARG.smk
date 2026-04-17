# ARG analysis workflow
#
# Usage (Snakemake v9):
#   ~/.pixi/bin/snakemake -s ARG.smk --executor slurm -j 25 --latency-wait 60 &
#snakemake -s ARG.smk --executor dryrun -j 10 --latency-wait 60

#Snakemake will submit all 3 jobs (potentially in parallel, depending on your -j limit). 
#VCF_DIR/results/TE_mutation_rate_2R.tsv
#VCF_DIR/results/strong_mutation_rate_2R.tsv
#VCF_DIR/results/neutral_mutation_rate_2R.tsv

#If you want to run only a subset, you can change line 98:
#focal=FOCAL_TYPES.keys(), to e.g.: focal=["TE"],
#or target specific outputs on the command line without changing the file:
#snakemake -s ARG.smk --executor slurm -j 10 \
#  .../results/TE_mutation_rate_2R.tsv
#snakemake -s ARG.smk --executor slurm -j 10 \
#  .../results/TE_mutation_rate_weight_2L.tsv

# NOTE: update the paths in add_TEs_to_vcf.py to match HPC paths before running.
#Right now with CHROMS = ["2R"], rule all will request all 3 focal types for chromosome 2R simultaneously

# ── Configuration ─────────────────────────────────────────────────────────────

BASE_DIR   = "/dfs7/grylee/yuhenh3/mutation_accumulation/ARG"
VCF_DIR    = BASE_DIR + "/Rech_2022_data/SNPs_vcf"
BED_DIR    = BASE_DIR + "/Rech_2022_data/TE_annotations/ReferenceCoordinates"
RAW_VCF    = VCF_DIR  + "/PhasedSNPsFitTEs.vcf"
CODE_DIR   = BASE_DIR + "/code"

SINGER   = BASE_DIR + "/SINGER-0.1.8-beta/SINGER/SINGER/parallel_singer_modified"
POLEGON  = BASE_DIR + "/polegon-0.1.3-alpha-linux-x86_64/polegon_master"

# Absolute paths so remove_TE_from_trees.py's hardcoded TREES_DIR stays consistent
OUT_DIR    = VCF_DIR + "/singer_"                    # e.g. .../SNPs_vcf/singer_2L/
REDATE_DIR = VCF_DIR + "/redate_no_TEs_"  # e.g. .../redate_..._2L/

# Use Python/binaries directly from the env — avoids conda/mamba activate issues
# $MAMBA_ROOT_PREFIX is set automatically after: module load mamba/24.3.0
TSKIT_PY   = "$HOME/.conda/envs/tskit/bin/python3"
SNPEFF_BIN = "$HOME/.conda/envs/snpeff/bin/snpEff"

CHROMS = ["2L"]   # extend as needed: ["2L", "2R", "3L", "3R", "X"]

wildcard_constraints:
    chr = "|".join(["2L", "2R", "3L", "3R", "X"]),

# Euchromatic boundaries per chromosome arm
EUCHROMATIC = {
    "2L": (500000,  21501009),
    "2R": (5898184, 24786936),
    "3L": (500000,  22462476),
    "3R": (5052934, 31579331),
    "X":  (500000,  21659299),
}

SINGER_PARAMS = dict(
    Ne        = "1e6",
    mu        = "5.0e-9",
    ratio     = "2.0",
    n_samples = 100,
    thin      = 20,
    num_cores = 20,
)

# Polegon only redates post-burn-in samples (first 50 discarded as burn-in)
POLEGON_SAMPLES = range(50, 100)

# Focal mutation types for mutation-rate analysis
FOCAL_TYPES = {
    "TE": {
        "filter_subdir": "filtered_no_TEs",
        "trees_dir": VCF_DIR + "/redate_no_TEs_",
        "trees_name": "redate_singer_",
    },
    "strong": {
        "filter_subdir": "filtered_no_TEs_no_strong_SNPs",
        "trees_dir": VCF_DIR + "/redate_no_TEs_no_strong_SNPs_",
        "trees_name": "redate_singer_for_strong_SNPs_",
    },
    "neutral": {
        "filter_subdir": "filtered_no_TEs_no_neutral_SNPs",
        "trees_dir": VCF_DIR + "/redate_no_TEs_no_neutral_SNPs_",
        "trees_name": "redate_singer_for_neutral_SNPs_",
    },
}

def focal_trees_input(wildcards):
    cfg = FOCAL_TYPES[wildcards.focal]
    return [
        f"{cfg['trees_dir']}{wildcards.chr}/{cfg['trees_name']}{wildcards.chr}_{i}.trees"
        for i in POLEGON_SAMPLES
    ]

def te_trees_input(wildcards):
    """TE-based Polegon tree samples (used regardless of focal type)."""
    cfg = FOCAL_TYPES["TE"]
    return [
        f"{cfg['trees_dir']}{wildcards.chr}/{cfg['trees_name']}{wildcards.chr}_{i}.trees"
        for i in POLEGON_SAMPLES
    ]

# ── Target rule ───────────────────────────────────────────────────────────────

rule all:
    input:
        expand(
            VCF_DIR + "/results/{focal}_mutation_rate_{chr}.tsv",
            chr=CHROMS,
            focal=FOCAL_TYPES.keys(),
        ),
        expand(
            VCF_DIR + "/results/{focal}_mutation_rate_weight_{chr}.tsv",
            chr=CHROMS,
            focal=FOCAL_TYPES.keys(),
        ),
        expand(
            VCF_DIR + "/results/{focal}_mutation_rate_skip_{chr}.tsv",
            chr=CHROMS,
            focal=FOCAL_TYPES.keys(),
        ),
        expand(
            VCF_DIR + "/results/{focal}_mutation_rate_weight_on_TE_trees_{chr}.tsv",
            chr=CHROMS,
            focal=FOCAL_TYPES.keys(),
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

# ── Step II: run SINGER per chromosome ────────────────────────────────────────

rule singer:
    input:
        vcf = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic_with_TEs.vcf",
    output:
        trees = expand(
            OUT_DIR + "{{chr}}/singer_{{chr}}_{i}.trees",
            i=range(SINGER_PARAMS["n_samples"]),
        )
    params:
        outprefix = lambda wc: f"{OUT_DIR}{wc.chr}/singer_{wc.chr}",
        Ne        = SINGER_PARAMS["Ne"],
        mu        = SINGER_PARAMS["mu"],
        ratio     = SINGER_PARAMS["ratio"],
        n         = SINGER_PARAMS["n_samples"],
        thin      = SINGER_PARAMS["thin"],
        num_cores = SINGER_PARAMS["num_cores"],
    threads: SINGER_PARAMS["num_cores"]
    resources:
        slurm_account   = "grylee_lab",
        slurm_partition = "standard",
        nodes           = 1,
        tasks           = 1,
        runtime         = 2880,   # minutes (= 2 days)
    log:
        VCF_DIR + "/logs/singer_{chr}.log"
    shell:
        """
        module load mamba/24.3.0
        mkdir -p {OUT_DIR}{wildcards.chr}
        {SINGER} \
            -Ne {params.Ne} \
            -m  {params.mu} \
            -ratio {params.ratio} \
            -vcf {input.vcf} \
            -output {params.outprefix} \
            -n {params.n} \
            -thin {params.thin} \
            -num_cores {params.num_cores} \
        > {log} 2>&1
        """

# ── Step III: filter trees and redate with Polegon (per focal type) ───────────
#
# For each focal type, two rules are generated:
#   filter_trees_{focal} — remove focal sites from singer trees
#   polegon_{focal}      — redate filtered trees with Polegon

for _focal, _cfg in FOCAL_TYPES.items():

    rule:
        name:
            f"filter_trees_{_focal}"
        input:
            trees   = expand(
                OUT_DIR + "{{chr}}/singer_{{chr}}_{i}.trees",
                i=range(SINGER_PARAMS["n_samples"]),
            ),
            vcf_tes = VCF_DIR + "/PhasedSNPsFitTEs_with_TEs.vcf",
            ann_vcf = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic.ann.vcf",
        output:
            trees = expand(
                OUT_DIR + "{{chr}}/" + _cfg["filter_subdir"] + "/singer_{{chr}}_{i}.trees",
                i=range(SINGER_PARAMS["n_samples"]),
            ),
        params:
            focal      = _focal,
            output_dir = OUT_DIR + "{chr}/" + _cfg["filter_subdir"],
        resources:
            slurm_account   = "grylee_lab",
            slurm_partition = "standard",
            runtime         = 240,
        log:
            VCF_DIR + f"/logs/filter_trees_{_focal}_{{chr}}.log"
        shell:
            """
            module load mamba/24.3.0
            {TSKIT_PY} {CODE_DIR}/remove_TE_from_trees.py \
                --chr {wildcards.chr} \
                --focal {params.focal} \
                --trees_dir {OUT_DIR}{wildcards.chr} \
                --output_dir {params.output_dir} \
                --te_vcf {input.vcf_tes} \
                --ann_vcf {input.ann_vcf} \
            > {log} 2>&1
            """

    rule:
        name:
            f"polegon_{_focal}"
        input:
            trees = OUT_DIR + "{chr}/" + _cfg["filter_subdir"] + "/singer_{chr}_{i}.trees",
        output:
            trees = _cfg["trees_dir"] + "{chr}/" + _cfg["trees_name"] + "{chr}_{i}.trees",
        params:
            inprefix  = OUT_DIR + "{chr}/" + _cfg["filter_subdir"] + "/singer_{chr}_{i}",
            outprefix = _cfg["trees_dir"] + "{chr}/" + _cfg["trees_name"] + "{chr}_{i}",
            trees_dir = _cfg["trees_dir"],
        threads: 2
        resources:
            slurm_account   = "grylee_lab",
            slurm_partition = "standard",
            runtime         = 4320,   # minutes (= 3 days)
        log:
            VCF_DIR + f"/logs/polegon_{_focal}_{{chr}}_{{i}}.log"
        shell:
            """
            module load mamba/24.3.0
            mkdir -p {params.trees_dir}{wildcards.chr}
            {TSKIT_PY} {POLEGON} \
                -m 5.0e-9 \
                -input {params.inprefix} \
                -output {params.outprefix} \
                -num_samples 100 \
                -thin 20 \
            > {log} 2>&1
            """

# ── Step IV: mutation rate analysis (all posterior samples at once) ───────

rule mutation_rate:
    input:
        trees = focal_trees_input,
        vcf = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic_with_TEs.ann.vcf",
    output:
        tsv = VCF_DIR + "/results/{focal}_mutation_rate_{chr}.tsv",
    params:
        chrom = "{chr}",
        focal = "{focal}",
    resources:
        slurm_account   = "grylee_lab",
        slurm_partition = "standard",
        runtime         = 7200,   # minutes (= 5 days)
    log:
        VCF_DIR + "/logs/{focal}_mutation_rate_{chr}.log"
    shell:
        """
        module load mamba/24.3.0
        mkdir -p $(dirname {output.tsv})
        {TSKIT_PY} {CODE_DIR}/test_TE_mutation_rate.py \
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

rule mutation_rate_weight:
    input:
        trees = focal_trees_input,
        vcf = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic_with_TEs.ann.vcf",
    output:
        tsv = VCF_DIR + "/results/{focal}_mutation_rate_weight_{chr}.tsv",
    params:
        chrom = "{chr}",
        focal = "{focal}",
    resources:
        slurm_account   = "grylee_lab",
        slurm_partition = "standard",
        runtime         = 7200,   # minutes (= 5 days)
    log:
        VCF_DIR + "/logs/{focal}_mutation_rate_weight_{chr}.log"
    shell:
        """
        module load mamba/24.3.0
        mkdir -p $(dirname {output.tsv})
        {TSKIT_PY} {CODE_DIR}/test_TE_mutation_rate_weight.py \
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

rule mutation_rate_weight_on_TE_trees:
    input:
        trees = te_trees_input,
        vcf = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic_with_TEs.ann.vcf",
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

rule mutation_rate_skip:
    input:
        trees = focal_trees_input,
        vcf = VCF_DIR + "/PhasedSNPsFitTEs_{chr}_euchromatic_with_TEs.ann.vcf",
    output:
        tsv = VCF_DIR + "/results/{focal}_mutation_rate_skip_{chr}.tsv",
    params:
        chrom = "{chr}",
        focal = "{focal}",
    resources:
        slurm_account   = "grylee_lab",
        slurm_partition = "standard",
        runtime         = 7200,   # minutes (= 5 days)
    log:
        VCF_DIR + "/logs/{focal}_mutation_rate_skip_{chr}.log"
    shell:
        """
        module load mamba/24.3.0
        mkdir -p $(dirname {output.tsv})
        {TSKIT_PY} {CODE_DIR}/test_TE_mutation_rate_skip.py \
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

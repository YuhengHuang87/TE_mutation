# TE_mutation.smk
# Snakemake port of TE_mutation_calling.md:
#   1. Map MA-line + ancestor short reads against the A4 and nanopore references
#      (gunzip -> trim_galore -> bwa mem | samtools sort -> MarkDuplicates ->
#       AddOrReplaceReadGroups)
#   2. Joint Strelka germline call per MA group (Act / Control / Ubi) together
#      with its ancestor, restricted to the matching euchromatin region, against
#      both references.
#
# Run:  snakemake -s TE_mutation.smk --cores 16
#       snakemake -s TE_mutation.smk --cores 16 strelka_all      # only step 2

from pathlib import Path

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
WGS_ROOT  = "/dfs7/grylee/yuhenh3/mutation_accumulation/MA_lines_WGS"
ANC_ROOT  = "/dfs7/grylee/yuhenh3/mutation_accumulation/SuVar3_9_lines_WGS_nanopore/Illumina_seq_ancestor_lines"
NANO_REF_DIR = "/dfs7/grylee/yuhenh3/mutation_accumulation/SuVar3_9_lines_WGS_nanopore/nanopore_seq_ancestor_lines"

# Where final/intermediate outputs are written. Existing per-sample BAMs already
# live under {OUTPUT_ROOT}/bam/ and the Strelka run dirs under {OUTPUT_ROOT}/strelka/,
# so re-running the workflow against the same OUTPUT_ROOT picks them up.
OUTPUT_ROOT = f"{WGS_ROOT}/final_samples"
BAM_DIR     = f"{OUTPUT_ROOT}/bam"
STRELKA_DIR = f"{OUTPUT_ROOT}/strelka"
READS_DIR   = f"{OUTPUT_ROOT}/reads"
TRIM_DIR    = f"{OUTPUT_ROOT}/trim"
REF_FLAG_DIR = f"{OUTPUT_ROOT}/ref"

REFERENCES = {
    "A4": {
        "fasta":          f"{WGS_ROOT}/A4.Hifi.MT.fasta",
        "call_regions":   "/dfs7/grylee/yuhenh3/mutation_accumulation/A4_hifi_euchromatin_region.bed.gz",
    },
    "mc_initial": {
        "fasta":          f"{NANO_REF_DIR}/mc_initial_ragtag.scaffold.fasta",
        "call_regions":   "/dfs7/grylee/yuhenh3/mutation_accumulation/mc_nanopore_euchromatin_region.bed.gz",
    },
}

# MA-line groups: reads_dir + ancestor + sample IDs.
GROUPS = {
    "Act": {
        "reads_dir": f"{WGS_ROOT}/final_samples/Act/raw_reads",
        "ancestor":  "anc_Act_S",
        "samples": [
            "Act_1","Act_3","Act_4","Act_6","Act_7","Act_8","Act_11","Act_12",
            "Act_16","Act_17","Act_19","Act_20","Act_21","Act_23","Act_25",
            "Act_27","Act_28","Act_29","Act_30","Act_32","Act_33","Act_34",
            "Act_36","Act_37","Act_38","Act_40","Act_42","Act_43","Act_45",
            "Act_46","Act_49","Act_52","Act_53","Act_54","Act_55","Act_58",
            "Act_59","Act_60","Act_61","Act_63","Act_64","Act_65","Act_66",
            "Act_68","Act_70","Act_71","Act_73","Act_75",
        ],
    },
    "Control": {
        "reads_dir": f"{WGS_ROOT}/final_samples/Control/raw_reads",
        "ancestor":  "anc_new_mc8_control",
        "samples": [
            "Control_1","Control_3","Control_5","Control_13","Control_14",
            "Control_15","Control_16","Control_19","Control_21","Control_22",
            "Control_23","Control_24","Control_25","Control_27","Control_28",
            "Control_29","Control_31","Control_32","Control_34","Control_36",
            "Control_37","Control_40","Control_44","Control_46","Control_47",
            "Control_50","Control_52","Control_53","Control_54","Control_55",
            "Control_56","Control_57","Control_58","Control_59","Control_60",
            "Control_61","Control_62","Control_63","Control_64","Control_66",
            "Control_67","Control_72","Control_73","Control_74","Control_75",
            "Control_76","Control_77",
        ],
    },
    "Ubi": {
        "reads_dir": f"{WGS_ROOT}/final_samples/Ubi/raw_reads",
        "ancestor":  "anc_Ubi_S",
        "samples": [
            "Ubi_1","Ubi_2","Ubi_3","Ubi_5","Ubi_7","Ubi_8","Ubi_9","Ubi_10",
            "Ubi_13","Ubi_14","Ubi_15","Ubi_16","Ubi_17","Ubi_21","Ubi_25",
            "Ubi_27","Ubi_28","Ubi_30","Ubi_31","Ubi_32","Ubi_33","Ubi_34",
            "Ubi_35","Ubi_36","Ubi_37","Ubi_38","Ubi_39","Ubi_40","Ubi_41",
            "Ubi_43","Ubi_44","Ubi_45","Ubi_46","Ubi_47","Ubi_48","Ubi_49",
            "Ubi_52","Ubi_54","Ubi_55","Ubi_56","Ubi_57","Ubi_58","Ubi_60",
            "Ubi_61","Ubi_64","Ubi_65","Ubi_66","Ubi_69","Ubi_70","Ubi_72",
            "Ubi_73","Ubi_74",
        ],
    },
}

# Flat per-sample registry: sample_id -> directory holding raw fastq.gz reads.
SAMPLE_DIR = {}
for grp, info in GROUPS.items():
    for s in info["samples"]:
        SAMPLE_DIR[s] = info["reads_dir"]
    SAMPLE_DIR[info["ancestor"]] = ANC_ROOT

# -----------------------------------------------------------------------------
# Targets
# -----------------------------------------------------------------------------
def all_rg_bams():
    targets = []
    for ref in REFERENCES:
        for s in SAMPLE_DIR:
            targets.append(f"{BAM_DIR}/{s}_{ref}_reads_with_RG.bam")
    return targets

def all_strelka_vcfs():
    return [
        f"{STRELKA_DIR}/{grp}_{ref}/results/variants/variants.vcf"
        for grp in GROUPS
        for ref in REFERENCES
    ]

rule all:
    input:
        all_strelka_vcfs()

rule mapping_all:
    input:
        all_rg_bams()

rule strelka_all:
    input:
        all_strelka_vcfs()

# -----------------------------------------------------------------------------
# Step 1: per-sample read processing & mapping
# -----------------------------------------------------------------------------
rule gunzip_reads:
    input:
        r1 = lambda wc: f"{SAMPLE_DIR[wc.sample]}/{wc.sample}-READ1-Sequences.txt.gz",
        r2 = lambda wc: f"{SAMPLE_DIR[wc.sample]}/{wc.sample}-READ2-Sequences.txt.gz",
    output:
        r1 = temp(READS_DIR + "/{sample}-READ1-Sequences.txt"),
        r2 = temp(READS_DIR + "/{sample}-READ2-Sequences.txt"),
    shell:
        r"""
        gunzip -c {input.r1} > {output.r1}
        gunzip -c {input.r2} > {output.r2}
        """

rule trim_galore:
    input:
        r1 = READS_DIR + "/{sample}-READ1-Sequences.txt",
        r2 = READS_DIR + "/{sample}-READ2-Sequences.txt",
    output:
        r1 = temp(TRIM_DIR + "/{sample}-READ1-Sequences.txt_val_1.fq"),
        r2 = temp(TRIM_DIR + "/{sample}-READ2-Sequences.txt_val_2.fq"),
    params:
        outdir = TRIM_DIR,
    threads: 4
    shell:
        r"""
        trim_galore --paired --fastqc -o {params.outdir} {input.r1} {input.r2}
        """

rule bwa_index:
    input:
        fa = lambda wc: REFERENCES[wc.ref]["fasta"],
    output:
        touch(REF_FLAG_DIR + "/{ref}.indexed"),
    shell:
        r"""
        bwa index {input.fa}
        samtools faidx {input.fa}
        picard CreateSequenceDictionary -R {input.fa}
        """

rule map_and_dedup:
    """bwa mem | samtools sort | gatk MarkDuplicates fused into one rule so the
    intermediate sorted BAM is never declared as a Snakemake output. Existing
    final_samples/bam/{sample}_{ref}_marked_duplicates.bam files are then enough
    to satisfy the DAG; the ad-hoc sorted BAM lives in $TMPDIR for the duration
    of the job."""
    input:
        r1   = TRIM_DIR + "/{sample}-READ1-Sequences.txt_val_1.fq",
        r2   = TRIM_DIR + "/{sample}-READ2-Sequences.txt_val_2.fq",
        idx  = REF_FLAG_DIR + "/{ref}.indexed",
    output:
        bam     = temp(BAM_DIR + "/{sample}_{ref}_marked_duplicates.bam"),
        metrics = BAM_DIR + "/{sample}_{ref}_marked_dup_metrics.txt",
    params:
        fa = lambda wc: REFERENCES[wc.ref]["fasta"],
    threads: 16
    shell:
        r"""
        sorted_bam=$(mktemp --suffix=.bam)
        trap 'rm -f "$sorted_bam"' EXIT

        bwa mem -t {threads} {params.fa} {input.r1} {input.r2} \
            | samtools sort -@ {threads} -o "$sorted_bam"

        gatk MarkDuplicates \
            -I "$sorted_bam" \
            -O {output.bam} \
            -M {output.metrics} \
            -REMOVE_DUPLICATES True \
            -ASSUME_SORTED True \
            -VALIDATION_STRINGENCY LENIENT
        """

rule add_read_groups:
    input:
        bam = BAM_DIR + "/{sample}_{ref}_marked_duplicates.bam",
    output:
        bam = BAM_DIR + "/{sample}_{ref}_reads_with_RG.bam",
        bai = BAM_DIR + "/{sample}_{ref}_reads_with_RG.bai",
    shell:
        r"""
        gatk AddOrReplaceReadGroups \
            -I {input.bam} \
            -O {output.bam} \
            -SORT_ORDER coordinate \
            -RGID foo -RGLB bar -RGPL illumina -RGPU unit1 \
            -RGSM {wildcards.sample}_Sample \
            -CREATE_INDEX True
        """

# -----------------------------------------------------------------------------
# Step 2: joint Strelka germline calling per (group, reference)
# -----------------------------------------------------------------------------
def strelka_input_bams(wc):
    grp = GROUPS[wc.group]
    members = [grp["ancestor"]] + grp["samples"]
    return [f"{BAM_DIR}/{s}_{wc.ref}_reads_with_RG.bam" for s in members]

rule strelka_configure:
    input:
        bams = strelka_input_bams,
    output:
        runner = STRELKA_DIR + "/{group}_{ref}/runWorkflow.py",
    params:
        fa          = lambda wc: REFERENCES[wc.ref]["fasta"],
        regions     = lambda wc: REFERENCES[wc.ref]["call_regions"],
        run_dir     = lambda wc: f"{STRELKA_DIR}/{wc.group}_{wc.ref}",
        bam_args    = lambda wc, input: " ".join(f"--bam {b}" for b in input.bams),
    shell:
        r"""
        rm -rf {params.run_dir}
        configureStrelkaGermlineWorkflow.py \
            {params.bam_args} \
            --referenceFasta {params.fa} \
            --runDir {params.run_dir} \
            --callRegions {params.regions}
        """

rule strelka_run:
    input:
        runner = STRELKA_DIR + "/{group}_{ref}/runWorkflow.py",
    output:
        vcf = STRELKA_DIR + "/{group}_{ref}/results/variants/variants.vcf",
    threads: 32
    shell:
        r"""
        {input.runner} -m local -j {threads}
        gunzip -c {output.vcf}.gz > {output.vcf}
        """

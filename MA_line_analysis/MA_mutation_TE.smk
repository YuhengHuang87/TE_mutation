# MA_mutation_TE.smk
#
# Run MA_mutation_TE.py for every (group, reference) combination produced by
# the upstream TE_anno_mapping.smk pipeline (Strelka joint germline calls).
# Six combinations total: {Act, Control, Ubi} x {A4, mc_initial}.
#
#
# Each (group, ref) job emits four TSVs under {OUT_ROOT}/{group}_{ref}/:
#   per_sample_callable.tsv
#   per_te_callable.tsv
#   de_novo_mutations.tsv
#   per_te_summary.tsv

from pathlib import Path

# -----------------------------------------------------------------------------
# Configuration (mirrors TE_anno_mapping.smk so paths line up automatically)
# -----------------------------------------------------------------------------
WGS_ROOT     = "/dfs7/grylee/yuhenh3/mutation_accumulation/MA_lines_WGS"
NANO_REF_DIR = "/dfs7/grylee/yuhenh3/mutation_accumulation/SuVar3_9_lines_WGS_nanopore/nanopore_seq_ancestor_lines/real_run/results"

OUTPUT_ROOT  = f"{WGS_ROOT}/final_samples"
STRELKA_DIR  = f"{OUTPUT_ROOT}/strelka"
OUT_ROOT     = f"{OUTPUT_ROOT}/MA_mutation_TE"

SCRIPT = str(Path(workflow.basedir) / "MA_mutation_TE.py")

REFERENCES = {
    "A4": {
        "repeatmasker":     f"{WGS_ROOT}/A4.Hifi.Scaf.Repeat.MT.fasta.out",
        "te_annotation":    f"{WGS_ROOT}/A4_A7.UU_MU.Zscore.SigFlag.txt",
        "euchromatin_bed":  "/dfs7/grylee/yuhenh3/mutation_accumulation/A4_hifi_euchromatin_region.bed.gz",
        "chroms":           ["2L", "2R", "3L", "3R", "X"],
    },
    "mc_initial": {
        "repeatmasker":     f"{NANO_REF_DIR}/repeatmasker/mc_TE_eu_ini_scaffold.out",
        "te_annotation":    f"{NANO_REF_DIR}/te_annotation/merged_class.txt",
        "euchromatin_bed":  "/dfs7/grylee/yuhenh3/mutation_accumulation/mc_nanopore_euchromatin_region.bed.gz",
        "chroms":           ["2L_RagTag", "2R_RagTag", "3L_RagTag", "3R_RagTag", "X_RagTag"],
    },
}

# Sample lists per group. Must stay in sync with GROUPS[group] in
# TE_anno_mapping.smk -- the BAM order Strelka was given. Column 9 of the
# joint VCF is the ancestor; columns 10..(9+N-1) are the MA samples in this
# exact order. genome.S1.vcf is the ancestor, genome.S{i+1}.vcf is samples[i].
#
# These names get a "_Sample" suffix appended to match the RGSM that
# `add_read_groups` writes (see TE_anno_mapping.smk -RGSM {wildcards.sample}_Sample),
# which is what Strelka would have used as the column name had the BAMs been
# re-RG'd. They are passed via --sample-names to override the (often
# non-unique, e.g. all "Sample1") names that survived in older joint VCFs.
SAMPLES = {
    "Act": {
        "ancestor": "anc_Act_S",
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
        "ancestor": "anc_new_mc8_control",
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
        "ancestor": "anc_Ubi_S",
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

NUM_MA = {g: len(SAMPLES[g]["samples"]) for g in SAMPLES}

def sample_name_args(group):
    """Build the --sample-names argv list (ancestor first, then MA samples in
    BAM order). The "_Sample" suffix matches the RGSM convention in
    TE_anno_mapping.smk's add_read_groups rule."""
    s = SAMPLES[group]
    return [s["ancestor"] + "_Sample"] + [m + "_Sample" for m in s["samples"]]

GROUPS = list(SAMPLES.keys())
REFS   = list(REFERENCES.keys())

# Pin wildcard regexes so the `{group}_{ref}` join in output paths is parsed
# unambiguously. Without this, `Act_mc_initial` greedily matches as
# group=Act_mc, ref=initial.
wildcard_constraints:
    group = "|".join(GROUPS),
    ref   = "|".join(REFS),

# -----------------------------------------------------------------------------
# Targets
# -----------------------------------------------------------------------------
OUT_SUFFIXES = [
    "per_sample_callable.tsv",
    "de_novo_mutations.tsv",
    "per_te_summary.tsv",
]

def outputs_for(group, ref):
    return [
        f"{OUT_ROOT}/{group}_{ref}/{group}_{ref}.{suffix}"
        for suffix in OUT_SUFFIXES
    ]

def all_outputs():
    return [
        path
        for group in GROUPS
        for ref in REFS
        for path in outputs_for(group, ref)
    ]

rule all:
    input:
        all_outputs()

# Convenience aliases: `snakemake -s MA_mutation_TE.smk Act_A4` etc.
for group in GROUPS:
    for ref in REFS:
        rule:
            name: f"{group}_{ref}"
            input:
                outputs_for(group, ref)

# -----------------------------------------------------------------------------
# Per (group, ref) analysis
# -----------------------------------------------------------------------------
rule mutation_te_analysis:
    input:
        joint_vcf      = lambda wc: f"{STRELKA_DIR}/{wc.group}_{wc.ref}/results/variants/variants.vcf",
        repeatmasker   = lambda wc: REFERENCES[wc.ref]["repeatmasker"],
        te_annotation  = lambda wc: REFERENCES[wc.ref]["te_annotation"],
        euchromatin    = lambda wc: REFERENCES[wc.ref]["euchromatin_bed"],
        script         = SCRIPT,
    output:
        per_sample = OUT_ROOT + "/{group}_{ref}/{group}_{ref}.per_sample_callable.tsv",
        muts       = OUT_ROOT + "/{group}_{ref}/{group}_{ref}.de_novo_mutations.tsv",
        summary    = OUT_ROOT + "/{group}_{ref}/{group}_{ref}.per_te_summary.tsv",
    params:
        gvcf_dir     = lambda wc: f"{STRELKA_DIR}/{wc.group}_{wc.ref}/results/variants",
        prefix       = lambda wc: f"{OUT_ROOT}/{wc.group}_{wc.ref}/{wc.group}_{wc.ref}",
        out_dir      = lambda wc: f"{OUT_ROOT}/{wc.group}_{wc.ref}",
        num_ma       = lambda wc: NUM_MA[wc.group],
        chroms       = lambda wc: " ".join(REFERENCES[wc.ref]["chroms"]),
        sample_names = lambda wc: " ".join(sample_name_args(wc.group)),
    threads: 1
    resources:
        mem_mb  = 16000,
        runtime = 24 * 60,
    shell:
        r"""
        mkdir -p {params.out_dir}
        python {input.script} \
            --group {wildcards.group} \
            --reference {wildcards.ref} \
            --joint-vcf {input.joint_vcf} \
            --gvcf-dir {params.gvcf_dir} \
            --repeatmasker {input.repeatmasker} \
            --te-annotation {input.te_annotation} \
            --euchromatin-bed {input.euchromatin} \
            --num-ma-samples {params.num_ma} \
            --output-prefix {params.prefix} \
            --chroms {params.chroms} \
            --sample-names {params.sample_names}
        """

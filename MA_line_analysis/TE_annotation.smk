### Snakemake pipeline: nanopore assembly + scaffolding + RepeatMasker
### Translates hifiasm_nano.sub and repeatmasker_nano.sub into a single workflow.
###
### Uses pre-existing mamba envs (not built by Snakemake). Env paths come from
### config.yaml -> conda_envs. `conda:` values that look like a directory path
### are treated by Snakemake as existing env prefixes — no solving, no build.
###
### Local run:    snakemake -s TE_annotation.smk --use-conda --cores 32
### SLURM run:    snakemake -s TE_annotation.smk --executor slurm --use-conda --jobs 20 \
###                 --default-resources slurm_account=grylee_lab \
###                                     slurm_partition=standard runtime=1440
###
### Per-rule runtime/cpu overrides live in `resources:` blocks below.
### `threads:` is auto-mapped to SLURM `cpus_per_task` by the slurm executor.

configfile: "config.yaml"

BARCODE         = config["barcode"]              # e.g. "barcode06"
RAW_BAM         = config["raw_bam"]              # suvar_nanopore_run1.bam
QV_MIN          = config.get("qv_min", 10)
REF_FASTA       = config["reference_fasta"]      # A4.Hifi.MT.fasta
TELIB           = config["te_library"]           # D_mel_transposon_sequence_set.fa
HIFIASM_THREADS = config.get("hifiasm_threads", 32)
RM_THREADS      = config.get("repeatmasker_threads", 32)

ENVS = config["conda_envs"]                      # mapping name -> env prefix path

# TE annotation post-processing (te_annotate.py).
# `euchromatin_bed` is a per-assembly BED on the scaffold's chromosome
# names (e.g. "2L_RagTag\t97997\t21400000"). Defined externally — e.g.
# from CUT&Tag data on the mc_initial genome.
TE_OUTDIR        = "results/te_annotation"
EUCH_BED         = config.get("euchromatin_bed")
if not EUCH_BED:
    raise ValueError(
        "config must set 'euchromatin_bed' to a per-assembly euchromatin "
        "BED on scaffold chromosome names."
    )

TE_ANN           = config.get("te_annotation", {})
TE_CLASSES       = TE_ANN.get("classes", "LINE,LTR,DNA,Unknown")
FAMILY_MERGE_DIS = TE_ANN.get("family_merge_dist", 200)
MIN_LENGTH       = TE_ANN.get("min_length", 500)
EXCL_FAMILY      = TE_ANN.get("exclude_family_pattern", "INE-1")
CLASS_MERGE_DIS  = TE_ANN.get("class_merge_dist", 5000)


rule all:
    input:
        "results/assembly_stats.txt",
        "results/repeatmasker/ragtag.scaffold.fasta.out",
        f"{TE_OUTDIR}/merged_class.txt"


rule bam_to_fastq:
    input:
        bam = RAW_BAM
    output:
        fq = f"results/reads/suvar_nanopore.QVMin{QV_MIN}.barcode_all.fastq"
    params:
        qv = QV_MIN
    threads: 4
    conda: ENVS["hifiasm"]
    shell:
        r"""
        samtools view -h -e '[qs]>={params.qv}' {input.bam} \
          | samtools fastq -T qs,RG - \
          | gzip --stdout > {output.fq}.gz
        gunzip -f {output.fq}.gz
        """


rule separate_barcode:
    input:
        fq = f"results/reads/suvar_nanopore.QVMin{QV_MIN}.barcode_all.fastq"
    output:
        fq = "results/reads/{barcode}_reads.fastq"
    shell:
        r"""
        awk -v bc={wildcards.barcode} '
          BEGIN {{ n=0 }}
          $0 ~ bc && n==0 {{ n=4 }}
          n>0 {{ print; n-- }}
        ' {input.fq} > {output.fq}
        """


rule hifiasm_assemble:
    input:
        fq = "results/reads/{barcode}_reads.fastq"
    output:
        gfa = "results/hifiasm/{barcode}.asm.bp.p_ctg.gfa"
    params:
        prefix = "results/hifiasm/{barcode}.asm"
    threads: HIFIASM_THREADS
    resources:
        runtime = 4080  # 2-20:00:00
    conda: ENVS["hifiasm"]
    shell:
        r"""
        mkdir -p results/hifiasm
        hifiasm --ont -o {params.prefix} -t {threads} -l 0 {input.fq}
        """


rule gfa_to_fasta:
    input:
        gfa = "results/hifiasm/{barcode}.asm.bp.p_ctg.gfa"
    output:
        fa = "results/hifiasm/{barcode}.asm.p_ctg.fa"
    shell:
        r"""
        awk '/^S/{{print ">"$2; print $3}}' {input.gfa} > {output.fa}
        """


rule ragtag_correct:
    input:
        ref = REF_FASTA,
        query = f"results/hifiasm/{BARCODE}.asm.p_ctg.fa"
    output:
        fa = "results/ragtag_correct/ragtag.correct.fasta"
    params:
        outdir = "results/ragtag_correct"
    threads: 8
    conda: ENVS["ragtag"]
    shell:
        r"""
        ragtag.py correct -o {params.outdir} -t {threads} {input.ref} {input.query}
        """


rule ragtag_scaffold:
    input:
        ref = REF_FASTA,
        corrected = "results/ragtag_correct/ragtag.correct.fasta"
    output:
        fa = "results/ragtag_scaffold/ragtag.scaffold.fasta"
    params:
        outdir = "results/ragtag_scaffold"
    threads: 8
    conda: ENVS["ragtag"]
    shell:
        r"""
        ragtag.py scaffold -o {params.outdir} -t {threads} {input.ref} {input.corrected}
        """


rule assembly_stats:
    input:
        fa = "results/ragtag_scaffold/ragtag.scaffold.fasta"
    output:
        stats = "results/assembly_stats.txt"
    conda: ENVS["assembly_stats"]
    shell:
        r"""
        assembly-stats {input.fa} > {output.stats}
        """


rule repeatmasker:
    input:
        fa = "results/ragtag_scaffold/ragtag.scaffold.fasta",
        lib = TELIB
    output:
        out = "results/repeatmasker/ragtag.scaffold.fasta.out"
    params:
        outdir = "results/repeatmasker"
    threads: RM_THREADS
    resources:
        runtime = 15600  # 10-20:00:00
    conda: ENVS["repeatmasker"]
    shell:
        r"""
        mkdir -p {params.outdir}
        RepeatMasker -pa {threads} -s -e rmblast \
            -lib {input.lib} \
            -dir {params.outdir} \
            {input.fa}
        """


rule te_annotate:
    """Post-RepeatMasker TE annotation: euchromatic filter, distance
    annotation, family merge, length/INE filter, and DNA/RNA class merge.
    Ports the post-RepeatMasker Perl chain into te_annotate.py."""
    input:
        rm_out      = "results/repeatmasker/ragtag.scaffold.fasta.out",
        scaffold    = "results/ragtag_scaffold/ragtag.scaffold.fasta",
        euchromatin = EUCH_BED
    output:
        filtered      = f"{TE_OUTDIR}/filtered.out",
        library       = f"{TE_OUTDIR}/library.fasta",
        distance      = f"{TE_OUTDIR}/distance.txt",
        merged_family = f"{TE_OUTDIR}/merged_family.txt",
        length_filt   = f"{TE_OUTDIR}/length_filtered.txt",
        merged_class  = f"{TE_OUTDIR}/merged_class.txt"
    params:
        outdir     = TE_OUTDIR,
        classes    = TE_CLASSES,
        fam_dist   = FAMILY_MERGE_DIS,
        min_len    = MIN_LENGTH,
        excl       = EXCL_FAMILY,
        class_dist = CLASS_MERGE_DIS
    shell:
        r"""
        python3 te_annotate.py \
            --rm-out {input.rm_out} \
            --scaffold {input.scaffold} \
            --euchromatin {input.euchromatin} \
            --outdir {params.outdir} \
            --te-classes {params.classes} \
            --family-merge-dist {params.fam_dist} \
            --min-length {params.min_len} \
            --exclude-family-pattern {params.excl} \
            --class-merge-dist {params.class_dist}
        """

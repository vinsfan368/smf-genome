# SMF bismark pipeline
# trim_galore → bismark (non-directional PE) → deduplicate_bismark → samtools sort/index
#
# Run locally:   snakemake -c8
# Run on SLURM:  snakemake --executor cluster-generic \
#                  --cluster-generic-submit-cmd "sbatch --cpus-per-task={threads} --mem=32G -t 4:00:00"

import pandas as pd

configfile: "config.yaml"

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
SAMPLES = samples.index.tolist()

# Bismark names the PE BAM after the R1 trimmed stem
_S = "{sample}_val_1_bismark_bt2_pe"


rule all:
    input:
        expand(f"aligns/{_S}.deduplicated.sorted.bam",     sample=SAMPLES),
        expand(f"aligns/{_S}.deduplicated.sorted.bam.bai", sample=SAMPLES),


rule trim:
    """Adapter trimming: stringency=3, 9 bp 5'-clip, paired FastQC."""
    input:
        r1=lambda wc: samples.loc[wc.sample, "R1"],
        r2=lambda wc: samples.loc[wc.sample, "R2"],
    output:
        r1="trimgalore/{sample}_val_1.fq.gz",
        r2="trimgalore/{sample}_val_2.fq.gz",
    log:
        "logs/trim/{sample}.log",
    threads: config["trim_cores"]
    shell:
        """
        trim_galore \
            --fastqc \
            --stringency 3 \
            --cores {threads} \
            --clip_R1 9 --clip_R2 9 \
            --paired \
            --basename {wildcards.sample} \
            -o trimgalore/ \
            {input.r1} {input.r2} \
            > {log} 2>&1
        """


rule align:
    """Bismark PE alignment, non-directional, with strandID."""
    input:
        r1="trimgalore/{sample}_val_1.fq.gz",
        r2="trimgalore/{sample}_val_2.fq.gz",
    output:
        # temp() removes pre-dedup BAM once dedup completes; comment out to keep
        bam=temp(f"aligns/{_S}.bam"),
    params:
        idx=config["bismark_index"],
    log:
        "logs/align/{sample}.log",
    threads: config["bismark_parallel"] * 4   # bismark uses ~4 threads per --parallel instance
    shell:
        """
        bismark \
            {params.idx} \
            -1 {input.r1} -2 {input.r2} \
            --parallel {config[bismark_parallel]} \
            --non_directional \
            --minins 76 --maxins 1000 \
            --strandID \
            -o aligns/ \
            > {log} 2>&1
        """


rule dedup:
    """Deduplicate Bismark PE BAM."""
    input:
        bam=f"aligns/{_S}.bam",
    output:
        bam=f"aligns/{_S}.deduplicated.bam",
    params:
        bin=config["dedup_bin"],
    log:
        "logs/dedup/{sample}.log",
    shell:
        "{params.bin} --bam {input.bam} > {log} 2>&1"


rule sort_index:
    """Sort and index deduplicated BAM (required for pysam random access)."""
    input:
        bam=f"aligns/{_S}.deduplicated.bam",
    output:
        bam=f"aligns/{_S}.deduplicated.sorted.bam",
        bai=f"aligns/{_S}.deduplicated.sorted.bam.bai",
    log:
        "logs/sort_index/{sample}.log",
    shell:
        """
        samtools sort -o {output.bam} {input.bam} > {log} 2>&1
        samtools index {output.bam} >> {log} 2>&1
        """

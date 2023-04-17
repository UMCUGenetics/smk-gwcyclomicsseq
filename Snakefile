# Auxiliary functions
import pandas as pd
from os.path import join as opj

def read_samples(file):
    df = pd.read_csv(file)
    return df

def get_no_bb_bams_per_sample(wildcards):
    run_names = SAMPLES[SAMPLES['sample_name'] == wildcards.sample_name]['run_name'].values
    bam_files = [f"{wildcards}_{run_name}_adapter_cleaned_noBB.bam" for run_name in run_names]
    bam_files = [opj(out_dir,bamfile) for bamfile in bam_files]

    return bam_files


def get_w_bb_bams_per_sample(wildcards):
    run_names = SAMPLES[SAMPLES['sample_name'] == wildcards.sample_name]['run_name'].values
    bam_files = [f"{wildcards}_{run_name}_adapter_cleaned_deduplicated_with_backbone.bam" for run_name in run_names]
    bam_files = [opj(out_dir,bamfile) for bamfile in bam_files]
    return bam_files


def get_all_bams_per_sample(wildcards):
    run_names = SAMPLES[SAMPLES['sample_name'] == wildcards.sample_name]['run_name'].values
    bam_files = [f"{wildcards}_{run_name}_deduplicated_all.bam" for run_name in run_names]
    bam_files = [opj(out_dir,bamfile) for bamfile in bam_files]

    return bam_files


SAMPLES = read_samples(config["sample_tsv"])
OUTPUT_DIR = config["output_dir"]
BB = config["backbones"] # "/Users/liting/01_data/20230407/backbones/BB41.fasta"
REF = config["references"] # "references/GRCh37_decoy/references_hs37d5_hs37d5.fa"

input_dir = "/Users/liting/01_data/20230407/"
out_dir = "/Users/liting/00_project/gwCyclomicsSeq/output"

input_bams = []
for sample in SAMPLES["sample_name"].unique():
    input_bams += list(SAMPLES[SAMPLES["sample_name"] == sample]['bam_name'].values)

def get_run_name(wildcards):
    run_names = SAMPLES.loc[wildcards.sample_name]['run_name']
    return run_names


SAMPLES['combined_prefix'] = SAMPLES.apply(lambda x: x['sample_name'] + '_'+ x['run_name'], axis = 1)
combine_dedup = [x+"_deduplicated_all.bam" for x in SAMPLES['combined_prefix']]
rule all:
    input:
        input_bams,
        combine_dedup,
        expand("{sample_name}_stats.txt", sample_name=SAMPLES["sample_name"]),


rule split_by_backbone:
    input:
        bam=opj(input_dir, "bams/{sample_name}_{run_name}.YM_gt_3.bam"),
        bb=opj(input_dir, "backbones/BB41.fasta"),
        ref=opj(input_dir, "references/GRCh37_decoy/references_hs37d5_hs37d5.fa"),
    output:
        with_bb="{sample_name}_{run_name}_reads_with_backbone.bam"Z,
        no_bb="{sample_name}_{run_name}_reads_without_backbone.bam",
    shell:
        """
        scripts/split_by_backbone.py --file_bam {input.bam} \\
                                 --file-out-insert-with-backbone {output.with_bb} \\
                                 --file-out-insert-no-backbone {output.no_bb}
        """

rule cutadapt_remove_bb:
    input:
        split_by_backbone=opj(out_dir, "{sample_name}_{run_name}_reads_with_backbone.bam"),
        bb=opj(input_dir, "backbones/BB41.fasta"),
    output:
        adapter_cleaned_with_bb=opj(out_dir, "{sample_name}_{run_name}_reads_with_backbone_removed_adapter.bam")
    shell:
        """
        bin/cutadapt_remove_bb.sh --input {input.split_by_backbone} \\
                                  --output {output.adapter_cleaned_with_bb} \\
                                  --backbone {input.bb}
        """

rule deduplication_bam_with_bb:
    input:
        adapter_cleaned_with_bb=opj(out_dir, "{sample_name}_{run_name}_reads_with_backbone_removed_adapter.bam"),
        ref=opj(input_dir, "references/GRCh37_decoy/references_hs37d5_hs37d5.fa"),
    output:
        dedup_with_bb=opj(out_dir, "{sample_name}_{run_name}_adapter_cleaned_deduplicated_with_backbone.bam"),
    shell:
        """
        bin/deduplication.py {input.adapter_cleaned_with_bb} {wildcards.sample_name} {wildcards.run_name} {input.ref}
        """

rule deduplication_without_bb:
    input:
        no_bb=opj(out_dir, "{sample_name}_{run_name}_reads_without_backbone.bam"),
        ref=opj(input_dir, "references/GRCh37_decoy/references_hs37d5_hs37d5.fa"),
    output:
        dedup_no_bb="{sample_name}_{run_name}_adapter_cleaned_noBB.bam"
    shell:
        """
        bin/deduplication.py {input.no_bb} {wildcards.sample_name} {wildcards.run_name} {input.ref}
        """

rule combine_dedup:
    input:
        dedup_with_bb=opj(out_dir, "{sample_name}_{run_name}_adapter_cleaned_deduplicated_with_backbone.bam"),
        dedup_no_bb=opj(out_dir, "{sample_name}_{run_name}_adapter_cleaned_noBB.bam"),
    output:
        dedup_combined=opj( out_dir, "{sample_name}_{run_name}_deduplicated_all.bam"),
    shell:
        """
        samtools merge -f {output.dedup_combined} {input.dedup_no_bb} {input.dedup_with_bb}
        """



rule merge_no_bb:
    input:
        bams=get_no_bb_bams_per_sample
    output:
        bam_per_sample_no_bb=opj(out_dir, "{sample_name}_all_noBB.bam")
    shell:
        """
        samtools merge -f {output.bam_per_sample_no_bb} {input.bams}
        """

rule merge_with_bb:
    input:
        bams=get_w_bb_bams_per_sample
    output:
        bam_per_sample_no_bb=opj(out_dir, "{sample_name}_all_withBB.bam")
    shell:
        """
        samtools merge -f {output.bam_per_sample_no_bb} {input.bams}
        """


rule merge_all:
    input:
        bams=get_all_bams_per_sample
    output:
        bam_per_sample_no_bb=opj(out_dir, "{sample_name}_allRuns_worwo_BB.bam")
    shell:
        """
        samtools merge -f {output.bam_per_sample_no_bb} {input.bams}
        """

rule get_sample_stats:
    input:
        bam = rules.merge_all.output.bam_per_sample_no_bb
    output:
        stats = "{sample_name}_stats.txt"
    shell:
        """
        bin/get_stats.sh --bam {input.bam}
        """



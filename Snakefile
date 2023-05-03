# Auxiliary functions
import pandas as pd
from os.path import join as opj
wildcard_constraints:
    sample_name = '[^_\W]+'


complement_translate = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')
def reverse_complement(seq):
    """Obtain reverse complement of seq
    returns:
        reverse complement (str)
    """
    return seq.translate(complement_translate)[::-1]


def get_backbone_sequence(backbone_file):
    with open(backbone_file) as f:
        x = f.readlines()
    dict_barcode = {}
    for line in x:
        if line.startswith(">"):
            seq_name = line.split(">")[-1].strip()
        else:
            sequence = line.strip()
            dict_barcode[seq_name] = sequence
    return dict_barcode


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
    bam_files = [f"{wildcards}_{run_name}_clean.bam" for run_name in run_names]
    bam_files = [opj(out_dir,bamfile) for bamfile in bam_files]
    return bam_files


SAMPLES = read_samples(config["sample_tsv"])
BB = config["backbones"] # "/Users/liting/01_data/20230407/backbones/BB41.fasta"

input_dir = config["input_dir"]
out_dir = config["out_dir"]

input_bams = []
for sample in SAMPLES["sample_name"].unique():
    input_bams += list(SAMPLES[SAMPLES["sample_name"] == sample]['bam_name'].values)
input_bams = [opj(input_dir, input1) for input1 in input_bams]

def get_run_name(wildcards):
    run_names = SAMPLES.loc[wildcards.sample_name]['run_name']
    return run_names


SAMPLES['combined_prefix'] = SAMPLES.apply(lambda x: x['sample_name'] + '_'+ x['run_name'], axis = 1)
combine_dedup = [opj(out_dir, x+"_clean.bam")  for x in SAMPLES['combined_prefix']]
rule all:
    input:
        input_bams,
        combine_dedup,
        expand(opj(out_dir, "{sample_name}_stats.txt"), sample_name=SAMPLES["sample_name"]),


rule split_by_backbone:
    input:
        bam=opj(input_dir, "{sample_name}_{run_name}.YM_gt_3.bam"),
        bb=BB,
        ref=config['reference'],
    output:
        with_bb=opj(out_dir, "{sample_name}_{run_name}_reads_with_backbone.bam"),
        no_bb=opj(out_dir, "{sample_name}_{run_name}_reads_without_backbone.bam"),
    resources:
        runtime_min=15,
        cpus = 1,
        mem_mb = 5000,
    priority: 50
    conda:
        "envs/align.yaml"
    shell:
        """
        python scripts/split_by_backbone.py -i {input.bam} \\
                                 -b {output.with_bb} \\
                                 -n {output.no_bb}
        """


rule deduplication_bam_with_bb:
    input:
        bam=opj(out_dir, "{sample_name}_{run_name}_reads_with_backbone.bam"),
        ref=config['reference'],
    output:
        out_bam=opj(out_dir, "{sample_name}_{run_name}_reads_with_backbone.dedup.bam"),
    conda:
        "envs/dedup-pip.yaml"
    priority: 2
    resources:
        runtime_min=400,
        cpus = 1,
        mem_mb = 10000,
    shell:
        """
        dedup --read_bam {input.bam} \\
              --out_bam {output.out_bam} \\
              --ref {input.ref}
        """

# backbone strings
bb_type = "BB41C"
backbone_forward = get_backbone_sequence(config["backbones"])[bb_type]
backbone_reverse = reverse_complement(backbone_forward)
rule cutadapt_remove_bb:
    input:
        split_by_backbone = rules.deduplication_bam_with_bb.output.out_bam,
    params:
        backbone_forward = backbone_forward,
        backbone_reverse = backbone_reverse,
    output:
        fastq = temp( opj(out_dir,"{sample_name}_{run_name}-extracted.fastq")),
        adapter_cleaned_with_bb_fastq = opj(out_dir,"{sample_name}_{run_name}_reads_with_backbone_removed_adapter.fastq"),
        info =  opj(out_dir,"{sample_name}_{run_name}-cutadapt.info"),
        summary = opj(out_dir,"{sample_name}_{run_name}-cutadapt.summary.txt"),
    priority: 50
    resources:
        runtime_min=30,
        cpus = 2,
        mem_mb = 5000,
    conda:
        "envs/align.yaml"
    shell:
        """
        bedtools bamtofastq -i {input.split_by_backbone} -fq {output.fastq};
        cutadapt -e 0.00064 -a {params.backbone_forward} -a {params.backbone_reverse} \\
                            -g {params.backbone_forward} -g {params.backbone_reverse} \\
                            --info-file {output.info} -o {output.adapter_cleaned_with_bb_fastq} {output.fastq} > {output.summary}
        """


rule map_after_cutadapt:
    input:
        bam = rules.cutadapt_remove_bb.output.adapter_cleaned_with_bb_fastq,
    params:
        ref=config['reference']
    threads:
        16
    output:
        sam=temp(opj(out_dir, "{sample_name}_{run_name}_reads_with_backbone_removed_adapter.sam")),
        adapter_cleaned_with_bb=opj(out_dir, "{sample_name}_{run_name}_adapter_cleaned_deduplicated_with_backbone.bam"),
    conda:
        "envs/align.yaml"
    priority: 50
    resources:
        runtime_min=120,
        cpus = 32,
        mem_mb = 16000,
    shell:
        """
        bwa mem -t {threads} {params.ref} {input.bam} > {output.sam};
        samtools sort {output.sam} > {output.adapter_cleaned_with_bb};
        samtools index {output.adapter_cleaned_with_bb};
        """



rule deduplication_without_bb:
    input:
        bam=opj(out_dir, "{sample_name}_{run_name}_reads_without_backbone.bam"),
        ref=config['reference'],
    output:
        out_bam=opj(out_dir, "{sample_name}_{run_name}_adapter_cleaned_noBB.bam"),
    conda:
        "envs/dedup-pip.yaml"
    priority: 2
    resources:
        runtime_min=400,
        cpus = 1,
        mem_mb = 10000,
    shell:
        """
        dedup --read_bam {input.bam} \\
              --out_bam {output.out_bam} \\
              --ref {input.ref}
        """


rule combine_dedup:
    input:
        dedup_with_bb=opj(out_dir, "{sample_name}_{run_name}_adapter_cleaned_deduplicated_with_backbone.bam"),
        dedup_no_bb=opj(out_dir, "{sample_name}_{run_name}_adapter_cleaned_noBB.bam"),
    output:
        dedup_combined=opj(out_dir, "{sample_name}_{run_name}_clean.bam"),
    conda:
        "envs/align.yaml"
    priority: 50
    resources:
        runtime_min=180,
        cpus = 2,
        mem_mb = 16000,
    shell:
        """
        samtools merge -f {output.dedup_combined} {input.dedup_no_bb} {input.dedup_with_bb};
        samtools index {output.dedup_combined}
        """


rule merge_no_bb:
    input:
        bams=get_no_bb_bams_per_sample
    output:
        bam_per_sample_no_bb=opj(out_dir, "{sample_name}_all_noBB.bam")
    conda:
        "envs/align.yaml"
    priority: 50
    resources:
        runtime_min=180,
        cpus = 2,
        mem_mb = 16000,
    shell:
        """
        samtools merge -f {output.bam_per_sample_no_bb} {input.bams};
        samtools index {output.bam_per_sample_no_bb};
        """


rule merge_with_bb:
    input:
        bams=get_w_bb_bams_per_sample
    output:
        bam_per_sample_no_bb=opj(out_dir, "{sample_name}_all_withBB.bam")
    conda:
        "envs/align.yaml"
    priority: 49
    resources:
        runtime_min=180,
        cpus = 2,
        mem_mb = 16000,
    shell:
        """
        samtools merge -f {output.bam_per_sample_no_bb} {input.bams};
        samtools index {output.bam_per_sample_no_bb};
        """


rule merge_all:
    input:
        bams=get_all_bams_per_sample
    output:
        bam_per_sample_no_bb=opj(out_dir, "{sample_name}_clean_merge_runs.bam")
    conda:
        "envs/align.yaml"
    priority: 48
    resources:
        runtime_min=180,
        cpus = 2,
        mem_mb = 16000,
    shell:
        """
        samtools merge -f {output.bam_per_sample_no_bb} {input.bams};
        samtools index {output.bam_per_sample_no_bb}
        """


rule get_sample_stats:
    input:
        bam = rules.merge_all.output.bam_per_sample_no_bb
    output:
        stats = opj(out_dir, "{sample_name}_stats.txt")
    conda:
        "envs/align.yaml"
    priority: 50
    resources:
        runtime_min=180,
        cpus = 2,
        mem_mb = 16000,
    shell:
        """
        scripts/get_stats.sh {input.bam} {output.stats}
        """

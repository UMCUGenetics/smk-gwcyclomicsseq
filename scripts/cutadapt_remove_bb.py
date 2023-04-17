import argparse
import os
import subprocess
import shlex


complement_translate = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')


def reverse_complement(seq):
    """Obtain reverse complement of seq
    returns:
        reverse complement (str)
    """
    return seq.translate(complement_translate)[::-1]


def complement(seq):
    """Obtain complement of seq
    returns:
        complement (str)
    """
    return seq.translate(complement_translate)


def get_basename(file_name):
    return os.path.splitext(file_name)[0]

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



def main(in_bam, out_bam, backbone_file, backbone_type, reference, info, summary):
    # define parameters
    name = get_basename(in_bam)
    print("basename: ", name)
    backbone_forward = get_backbone_sequence(backbone_file)[backbone_type]
    backbone_reverse = reverse_complement(backbone_forward)
    trimmed = f"{name}-trimmed.fastq"
    sam = f"{name}.sam"
    # defines commands
    cmd1 = f"bedtools bamtofastq -i {in_bam} -fq {name}.fastq"
    # cmd2 = f"samtools view -c -F 260 {name}.fastq"
    file3 = open(summary, "w")
    cmd3 = f"cutadapt -e 0.00064 -a {backbone_forward} -a {backbone_reverse} -g {backbone_forward} -g {backbone_reverse} --info-file {info} -o {trimmed} {name}.fastq"
    file4 = open(sam, "w")
    cmd4 = f"bwa mem {reference} {trimmed}"
    file5 = open(out_bam, "wb")
    cmd5 = f"samtools sort {sam}"
    cmd6 = f"samtools index {out_bam}"
    cmd7 = f"rm {sam}"
    # execute commands
    subprocess.run(shlex.split(cmd1))
    # subprocess.run(shlex.split(cmd2))
    subprocess.run(shlex.split(cmd3), stdout=file3)
    subprocess.run(shlex.split(cmd4), stdout=file4)
    subprocess.run(shlex.split(cmd5), stdout=file5)
    subprocess.run(shlex.split(cmd6))
    subprocess.run(shlex.split(cmd7))

    print("done with backbone removal + map to reference again.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Remove adapter for reads containing backbone - aggressive"
    )
    parser.add_argument("-i", "--in-bam", type=str)
    parser.add_argument("-o", "--out-bam", type=str)
    parser.add_argument("-b", "--backbone-file", type=str)
    parser.add_argument("-t", "--backbone-type", type=str, default="BB41C")
    parser.add_argument("-r", "--reference-genome", type=str)
    parser.add_argument("-s", "--summary", type=str)
    parser.add_argument("-s2", "--info", type=str)

    args = parser.parse_args()

    main(args.in_bam,
         args.out_bam,
         args.backbone_file,
         args.backbone_type,
         args.reference_genome,
         args.info,
         args.summary,
         )

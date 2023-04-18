import argparse
import pysam
from pathlib import Path
from datetime import datetime
from tqdm import tqdm
import os
import subprocess
import shlex




def backbone_in_mates(read,
                      tag_mate="YP",
                      backbone_contig_name="BB",
                      ):
    """
    Check if at least X backbone sequence is in mate of an insert read.
    tag_mate: value for tag containing mate information
    read: a pysam.AlignedSegment
    """
    YP_contains_BB = False
    # if "BB" in read.query_name:
    #     print(read.reference_name)
    if read.has_tag(tag_mate):
        if backbone_contig_name in read.get_tag(tag_mate):
            YP_contains_BB = True
        else:
            YP_contains_BB = False
    else:
        YP_contains_BB = False

    return YP_contains_BB

# To check if adapter trimming is necessary:

def get_backbone_contig_length():
    pass


def backbone_mate_mapping(read,
                          backbone_contig_length,
                          tag_mate="YP",
                          backbone_contig_name="BB",
                          tag_for_trimming="YD",
                      ):
    """
    Find the location of all subreads for a given raw read
    """
    if read.get_tag(tag_mate).contains(backbone_contig_name):
        tags = read.get_tag(tag_mate).split("|")
        for tag in tags:
            if tag.startswith("BB"):
                BB_start_pos = tag.split(":")[-1].split("-")[0]
                BB_end_pos = tag.split(":")[-1].split("-")[1]
                if BB_start_pos != 0 & BB_end_pos != backbone_contig_length:
                    trimming = True
                    # Add how many base to trim to this parameter: tag_for_trimming
    else:
        trimming = False
    #return trimming
    return 0


def main(in_bam_path, out_bam_path_with_backbone, out_bam_path_no_backbone):
    """
    Look up the Y tags in the metadata for all reads in a bam.
    """

    start_time = datetime.now()

    in_bam = pysam.AlignmentFile(in_bam_path, "rb")
    out_bam_with_backbone = pysam.AlignmentFile(out_bam_path_with_backbone, "wb", header=in_bam.header)
    out_bam_no_backbone = pysam.AlignmentFile(out_bam_path_no_backbone, "wb", header=in_bam.header)

    aln_count = 0
    wbb_count = 0
    nbb_count = 0
    for aln in tqdm(in_bam.fetch()):
        aln_count += 1
        if backbone_in_mates(aln):
            out_bam_with_backbone.write(aln)
            wbb_count += 1
        else:
            out_bam_no_backbone.write(aln)
            nbb_count += 1

    out_bam_with_backbone.close()
    out_bam_no_backbone.close()

    print(
        f"{in_bam_path} Split {aln_count} alignments with or without backbone in {datetime.now() - start_time}. {round(wbb_count/aln_count, 3)*100}% has backbone."
    )




if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process the information in the metadata and add it to the bam."
    )

    parser.add_argument("-i", "--file_bam", type=Path)
    parser.add_argument("-b", "--file_out_insert_with_backbone", type=Path)
    parser.add_argument("-n", "--file_out_insert_no_backbone", type=Path)
    args = parser.parse_args()
    intermediate_bam_w_bb = args.file_out_insert_with_backbone.with_suffix(".unsorted.bam")
    intermediate_bam_no_bb = args.file_out_insert_no_backbone.with_suffix(".unsorted.bam")


    main(args.file_bam,
         intermediate_bam_w_bb,
         intermediate_bam_no_bb)

    # pysam requires pure strings
    pysam.sort("-o", str(args.file_out_insert_with_backbone), str(intermediate_bam_w_bb))
    pysam.sort("-o", str(args.file_out_insert_no_backbone), str(intermediate_bam_no_bb))
    pysam.index(str(args.file_out_insert_with_backbone))
    pysam.index(str(args.file_out_insert_no_backbone))
    os.remove(intermediate_bam_no_bb)
    os.remove(intermediate_bam_w_bb)


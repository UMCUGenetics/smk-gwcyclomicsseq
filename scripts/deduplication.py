import pysam
import pandas as pd
from collections import Counter
import singlecellmultiomics
from singlecellmultiomics.molecule import MoleculeIterator
from singlecellmultiomics.bamProcessing import sorted_bam_file
from singlecellmultiomics.fragment import CHICFragment
from singlecellmultiomics.molecule.chic import CHICMolecule
from singlecellmultiomics.utils.prefetch import UnitialisedClass
from singlecellmultiomics.fastaProcessing import CachedFastaNoHandle
from singlecellmultiomics.bamProcessing.bamFunctions import get_reference_from_pysam_alignmentFile
import argparse
import time
import os
from tqdm import tqdm



def write_sm_tag_to_bam(input_bam, SM_bam, SM_tag_value):
    with pysam.AlignmentFile(input_bam, "r", ignore_truncation=True) as g:
        with sorted_bam_file(SM_bam, origin_bam=g) as f:
            for i, read in enumerate(g):
                read.set_tag("SM", SM_tag_value)
                f.write(read)
    return SM_bam


def write_deduplicate(input_bam_path, target_bam_path, reference, SM_tag_value, merge_max, rca_tag="YM"):
    unique_molecules_copies = Counter()
    with pysam.AlignmentFile(input_bam_path) as f:
        with sorted_bam_file(target_bam_path, origin_bam=f, ) as target_bam:
            for i, m in tqdm(enumerate(MoleculeIterator(
                    alignments=f,
                    moleculeClass=CHICMolecule,
                    fragmentClass=CHICFragment,
                    every_fragment_as_molecule=False,
                    perform_qflag=False,
                    molecule_class_args={"reference": reference, "max_associated_fragments": 100},
                    fragment_class_args={"assignment_radius": merge_max, "rca_tag": rca_tag, "sample": SM_tag_value},
                    max_buffer_size=1000000,
                    yield_overflow=False,
            ))):
                read_name = m.read_names[0]
                m.write_pysam(target_bam,
                              consensus=True,
                              consensus_name=read_name,
                              no_source_reads=True,
                              )
                unique_molecules_copies[len(m)] += 1
    print(sorted(unique_molecules_copies.items()))
    ### Write deduplication counter
    name = target_bam_path.split("/")[-1].split(".")[0]
    with open(out_bam+"unique_molecule_count.tsv", "w") as f:
        f.write( f"rep\t{name}\n" )
        for i in range(1, 101):
            appears_n_times = i
            k_molecules = unique_molecules_copies[i]
            f.write( f"{appears_n_times}\t{k_molecules}\n" )
    # pd.to_pickle(unique_molecules_copies, out_bam+".pickle.gz")


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Deduplicate base on molecular location.')
    parser.add_argument('--read_bam', type=str,
                    help='bam file path')
    parser.add_argument('--out_bam', type=str,
                    help='out bam including path')
    parser.add_argument('--SM', type=str, default=None, help="value to put in SM tag. Can skip is SM tag existed.")
    parser.add_argument('--ref', type=str, help='reference which bam file is mapped to. (str). Auto-detect is possible.')
    parser.add_argument('--merge', action='store_true', help='merge reads from different nanopore reads but covering the same start,'
                                                  'end sites within max X range. Direction of the read is set by the first added read.')
    parser.add_argument('--merge_max', type= int, default=5,
                        help='Merge nanopore reads within max X bp range. Strand of the read is the first added read.')
    args = parser.parse_args()

    timeA = time.time()

    out_bam = os.path.splitext(args.out_bam)[0]
    out_path = out_bam.split("/")[0]
    prefix = out_bam.split("/")[-1]
    # SM_bam
    SM_bam = f"{out_path}.SMtagged.sorted.bam"
    # deduplicated output bam
    t_bam = args.out_bam

    # autodetect reference:
    reference = None
    if args.ref is None:
        args.ref = get_reference_from_pysam_alignmentFile(args.read_bam)

    if args.ref is not None:
        try:
            reference = UnitialisedClass(CachedFastaNoHandle, args.ref)
            print(f'Loaded reference from {args.ref}')
        except Exception as e:
            print("Error when loading the reference file, continuing without a reference")
            reference = None

    # Check if SM tag addition is necessary (based on manual parameter)
    # TODO: convert the check to automatic check with SM tag in args.read_bam file.
    if args.SM is not None:
        input_bam = write_sm_tag_to_bam(args.read_bam, SM_bam, args.SM)
    else:
        input_bam = args.read_bam

    write_deduplicate(input_bam, t_bam, reference, args.SM, args.merge_max, rca_tag="YM")

    pysam.index(str(t_bam))

    # cleanup_intermediate_files
    if args.SM is not None:
        os.remove(input_bam)
        os.remove(f'{input_bam}.bai')

print((time.time() - timeA)/60, 'min')

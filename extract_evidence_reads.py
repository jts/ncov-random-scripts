#! /usr/bin/python

import argparse
import pysam
import sys
import os

parser = argparse.ArgumentParser()
parser.add_argument('--bam', type=str, default="", required=True)
parser.add_argument('--long-reads', dest="long_reads", action="store_true")
parser.add_argument('--position', type=int, required=True)
args = parser.parse_args()

#  https://github.com/pysam-developers/pysam/issues/939
save = pysam.set_verbosity(0)
bam_file = pysam.AlignmentFile(args.bam, "rb")
pysam.set_verbosity(save)

sample_name = "sample_" + os.path.basename(args.bam).split(".")[0].split("_")[0][-4:]
print("sample_name\tread_name\tfragment_start\tfragment_end\tfragment_length\tbase")
for alignment in bam_file.fetch(until_eof=True):
    base_at_position = "N"
    for (query_pos, reference_pos) in alignment.get_aligned_pairs():

        # pysam alignments are 0-based
        if reference_pos is not None and reference_pos + 1 == args.position and query_pos is not None:
            base_at_position = alignment.seq[query_pos]

    fragment_start = alignment.reference_start
    fragment_end = 0
    if not args.long_reads:
        tlen = alignment.template_length
        if tlen < 0:
            fragment_start = alignment.reference_end
        fragment_end = fragment_start + tlen
    else:
        fragment_end = alignment.reference_end
        if alignment.is_reverse:
            tmp = fragment_start
            fragment_start = fragment_end
            fragment_end = tmp
        tlen = 0
    print("%s\t%s\t%d\t%d\t%d\t%s" % (sample_name, alignment.query_name, fragment_start, fragment_end, tlen, base_at_position))

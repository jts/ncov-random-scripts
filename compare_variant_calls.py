import os
import csv
import sys
import pysam
import argparse

def load_variants(variants_filename):
    variant_map = dict()
    with open(variants_filename) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            variant_map[int(row['POS'])] = row
    return variant_map

# read in a .fofn file and map from sample -> variant call file
def read_fofn(fofn_filename):

    sample_to_file_map = dict()
    with open(fofn_filename) as fh:
        for path in fh:
            filename = os.path.basename(path)
            sample_name = filename.split(".")[0]
            sample_to_file_map[sample_name] = path.rstrip()
    return sample_to_file_map

parser = argparse.ArgumentParser()
parser.add_argument('--a-fofn', required=True, type=str)
parser.add_argument('--a-name', required=False, type=str, default="a")
parser.add_argument('--b-fofn', required=True, type=str)
parser.add_argument('--b-name', required=False, type=str, default="b")
args, extra = parser.parse_known_args()

a_map = read_fofn(args.a_fofn)
b_map = read_fofn(args.b_fofn)

print("\t".join(["sample", "contig", "position", "alt", "%s_depth" % (args.a_name), "%s_af" % (args.a_name), "%s_depth" % (args.b_name), "%s_af" % (args.b_name), "type"]))
for sample in a_map:
    if sample not in b_map:
        continue
    variants_a = load_variants(a_map[sample])
    variants_b = load_variants(b_map[sample])

    union_pos = list(set(variants_a.keys()) | set(variants_b.keys()))

    for pos in sorted(union_pos):
        a_depth = "na"
        a_af = 0
        b_depth = "na"
        b_af = 0
        qc_type = ""

        if pos in variants_a:
            a_depth = variants_a[pos]['TOTAL_DP']
            a_af = float(variants_a[pos]['ALT_FREQ'])
        if pos in variants_b:
            b_depth = variants_b[pos]['TOTAL_DP']
            b_af = float(variants_b[pos]['ALT_FREQ'])
        
        contig = ""
        var = None
        if a_depth == "na":
            qc_type = "%s_ONLY" % (args.b_name)
            var = variants_b[pos]
        elif b_depth == "na":
            qc_type = "%s_ONLY" % (args.a_name)
            var = variants_a[pos]
        elif abs(a_af - b_af) > 0.2:
            var = variants_a[pos]
            qc_type = "AF_CHANGE"
        
        if qc_type != "":
            contig = var["REGION"]
            alt = var["ALT"]
            print("%s\t%s\t%d\t%s\t%s\t%.2f\t%s\t%.2f\t%s" % (sample, contig, pos, alt, a_depth, a_af, b_depth, b_af, qc_type))


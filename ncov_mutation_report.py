#!/usr/bin/env python

import argparse
import pysam
import sys
import csv
import os
import re
from collections import defaultdict

class Sample:
    def __init__(self):
        self.lineage = "not_assigned"
        self.mutations = list()

def clean_sample_name(sample_name):
    sample_name = re.sub('^Consensus_', '', sample_name) # added by ivar
    sample_name = re.sub('.primertrimmed.consensus_threshold_0.75_quality_20', '', sample_name) # added by ivar
    sample_name = re.sub('_MN908947.3', '', sample_name) # added by pangolin
    sample_name = re.sub('/ARTIC/nanopolish', '', sample_name) # added by ARTIC
    sample_name = re.sub('/ARTIC/medaka', '', sample_name) # added by ARTIC
    sample_name = re.sub('.variants.tsv', '', sample_name) # added by ncov-watch (ivar)
    sample_name = re.sub('.variants.norm.vcf', '', sample_name) # added by ncov-watch (ivar)
    return sample_name

def load_watch_tsv(filename, data):
    with open(filename, 'r') as ifh:
        reader = csv.DictReader(ifh, delimiter='\t')
        for record in reader:
            sample = clean_sample_name(record['sample'])
            data[sample].mutations.append(record['mutation'])

def get_alt_for_type_variant(variant):
    fields = variant.split(":")
    if fields[0] == "del":
        return "del"
    else:
        assert(fields[0] == "aa" or fields[0] == "snp")
        match = re.search('(\D+)(\d+)(\D+)', fields[2])
        return match[3]

def load_type_variants(filename, data):
    with open(filename, 'r') as ifh:
        reader = csv.DictReader(ifh, delimiter=',')
        for record in reader:
            sample = clean_sample_name(record['query'])

            # this skips the non-genotype fields
            for variant, genotype in list(record.items())[5:]:
                variant_alt = get_alt_for_type_variant(variant)
                if genotype == variant_alt:
                    data[sample].mutations.append(variant)

def load_lineages(filename, data):
    with open(filename, 'r') as ifh:
        reader = csv.DictReader(ifh, delimiter=',')
        for record in reader:
            sample = clean_sample_name(record['taxon'])
            data[sample].lineage = record['lineage']

def main():

    description = 'Merge mutation and lineage results to produce a readable Report'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-l', '--lineage-data', required=True, help='the filename with pangolin output')
    parser.add_argument('-w', '--watch-data', default="", help='the filename with ncov-watch output')
    parser.add_argument('-t', '--type-variants-data', default="", help='the filename with type_variants.py output')
    parser.add_argument('-p', '--print-sample-names', action="store_true", help='print sample names for each record')
    args = parser.parse_args()

    if args.watch_data == "" and args.type_variants_data == "":
        sys.stderr.write("At least one of --watch-data or --type-variants-data must be provided\n")
        sys.exit(1)
    
    if args.watch_data != "" and args.type_variants_data != "":
        sys.stderr.write("At most one of --watch-data or --type-variants-data must be provided\n")
        sys.exit(1)

    data = defaultdict(Sample)

    if args.watch_data != "":
        load_watch_tsv(args.watch_data, data)
    else:
        load_type_variants(args.type_variants_data, data)
    
    load_lineages(args.lineage_data, data)

    lineage_count_by_mutation = defaultdict( lambda: defaultdict(int) )
    lineage_mutation_samples = defaultdict( lambda: defaultdict(list) )
    for sample_name, sample_data in data.items():
        for m in sample_data.mutations:
            lineage_count_by_mutation[m][sample_data.lineage] += 1
            lineage_mutation_samples[m][sample_data.lineage].append(sample_name)

    for m in lineage_count_by_mutation:
        print(m)
        for lineage, count in sorted(lineage_count_by_mutation[m].items(), key=lambda item: item[1], reverse=True):
            
            sample_name_str = ""
            if args.print_sample_names:
                sample_name_str = " (%s)" % (",".join(lineage_mutation_samples[m][lineage]))
            print("\t%s: %d%s" % (lineage, count, sample_name_str))
if __name__ == "__main__":
    main()

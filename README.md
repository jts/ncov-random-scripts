# ncov random scripts

Collection of random utility scripts for processing ARTIC sequencing results

## Scripts

### Extract evidence reads

This script will process each alignment in a BAM file and output which base each read supports at a particular reference position. Reads that do not cover the reference position are reported as `N`. I use this script to debug primer trimming and amplification artifacts. The coordinates output (`fragment_start`, `fragment_end`) are for the entire paired-end fragment, not the individual reads. This is to help visualize reads that cross amplicon boundaries.

Example usage:

`python extract_evidence_reads.py --bam input.bam --position 17747`

Output:

```
sample_name     read_name       fragment_start  fragment_end    fragment_length base
sample_abcd     read_1          17744           18103           359             T
```

### Compare variant calls between pipeline versions

This script compares collections of ivar `.variants.tsv` files to detect variants that are unique to one collection, or significantly changed variant allele frequency. I use this to evaluate the effects of changing pipeline versions.

Example usage:

```
find pipeline_v1.5_results -name "*.variants.tsv" > v1.5.fofn
find pipeline_v1.6_results -name "*.variants.tsv" > v1.6.fofn
python compare_variant_calls.py -a-fofn v1.5.fofn -b-fofn v1.6.fofn -a-name v1.5 -b-name v1.6
```

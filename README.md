# ncov random scripts

Collection of random utility scripts for processing ARTIC sequencing results

## Scripts

`extract_evidence_reads.py` - this script will process each alignment in a BAM file and output which base each read supports at a particular reference position. Reads that do not cover the reference position are reported as `N`. I use this script to debug primer trimming and amplification artifacts. The coordinates output (`fragment_start`, `fragment_end`) are for the entire paired-end fragment, not the individual reads. This is to help visualize reads that cross amplicon boundaries.

Example usage:

`python extract_evidence_reads.py --bam input.bam --position 17747`

Output:

```
sample_name     read_name       fragment_start  fragment_end    fragment_length base
sample_abcd     read_1          17744           18103           359             T
```

#! /bin/bash

IN_BAM=$1
PREFIX=`basename $IN_BAM .bam`
REF=nCoV-2019.reference.fasta

# the sed is to fix the header until a release is made with this fix
# https://github.com/freebayes/freebayes/pull/549
freebayes -p 1 -f $REF -F 0.2 -C 1 --pooled-continuous --min-coverage 10 \
    --gvcf --gvcf-dont-use-chunk true $IN_BAM | sed s/QR,Number=1,Type=Integer/QR,Number=1,Type=Float/ > $PREFIX.gvcf

# make depth mask, split variants into ambiguous/consensus
# NB: this has to happen before bcftools norm or else the depth mask misses any bases exposed during normalization
python process_gvcf.py -d 10 -m $PREFIX.mask.txt -a $PREFIX.ambiguous.vcf -c $PREFIX.consensus.vcf $PREFIX.gvcf

# normalize then gzip to make bcftools happy
for vt in "ambiguous" "consensus"; do
    bcftools norm -f $REF $PREFIX.$vt.vcf > $PREFIX.$vt.norm.vcf
    bgzip -f $PREFIX.$vt.norm.vcf
    tabix -f -p vcf $PREFIX.$vt.norm.vcf.gz
done

# apply ambiguous variants first using IUPAC codes. this variant set cannot contain indels or the subsequent step will break
bcftools consensus -f $REF -I $PREFIX.ambiguous.norm.vcf.gz > $PREFIX.ambiguous.fasta

# apply remaninng variants, including indels
bcftools consensus -f $PREFIX.ambiguous.fasta -m $PREFIX.mask.txt $PREFIX.consensus.norm.vcf.gz > $PREFIX.consensus.fasta

#TODO: rewrite header

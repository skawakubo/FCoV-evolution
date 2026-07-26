#!/bin/bash
set -euo pipefail

sample_ID="fecv_NARC3_S3"
ref="TRINITY_FCoV_FECV3_assembled_2.fas"
bam="FECV_NARC3_S3_mapped_on_assembled_v2.sort.bam"

samtools_cmd="/usr/local/package/samtools/1.9/bin/samtools"
bcftools_cmd="/usr/local/package/samtools/1.9/bin/bcftools"
vcfutils_pl_cmd="/usr/local/package/samtools/1.9/bin/vcfutils.pl"

mkdir -p bcf vcf

# 0. 参照のインデックス（未作成なら）
${samtools_cmd} faidx ${ref}

# 1. mpileup（深度上限を解除・BAQ無効化して取りこぼしを減らす）
${bcftools_cmd} mpileup -f ${ref} \
    -a FORMAT/AD,FORMAT/DP,INFO/AD \
    -d 1000000 \
    -q 0 -Q 0 -B \
    -Ob -o bcf/${sample_ID}_raw.bcf \
    ${bam}

# 2. call
${bcftools_cmd} call --ploidy 1 -m -v \
    -o vcf/${sample_ID}.pre.vcf \
    bcf/${sample_ID}_raw.bcf

# 3. filtering（最低深度1のみ、他の基準はすべて無効化）
${vcfutils_pl_cmd} varFilter \
    -d 1 \
    -D 1000000 \
    -a 1 \
    -Q 0 \
    -w 0 -W 0 \
    -1 0 -2 0 -3 0 -4 0 -e 0 \
    vcf/${sample_ID}.pre.vcf \
    > vcf/${sample_ID}.filtered.vcf

# 4. convert to easy-to-read table
python /home/shusuke/Izumi_FIP/vcf_to_table.py vcf/${sample_ID}.filtered.vcf > vcf/${sample_ID}.filtered_variant_table.txt

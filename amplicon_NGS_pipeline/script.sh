#!/bin/sh
#$ -S /bin/bash
#$ -pe def_slot 1
#$ -l s_vmem=300G
#$ #-o /dev/null
#$ -cwd
#$ -l lmem



READ_DIR=/home/skawa/Izumi_FIP/FCoV_seq_data  
REF=/home/skawa/Izumi_FIP/FCoV_seq_data_analysis_results/amplicon_base_mapping/FCoV1_KY566209.fasta         # 全長リファレンス
IDX=/home/skawa/Izumi_FIP/FCoV_seq_data_analysis_results/amplicon_base_mapping/FCoV1_ref               # bowtie2 index prefix
THREADS=8
OUT=all_amplicons_bt2_trim  
BED=/home/skawa/Izumi_FIP/FCoV_seq_data_analysis_results/amplicon_base_mapping/FCoV_amplicon_primer_v2.bed             
# TRIM_Q=20                   
# TRIM_MINLEN=30              
TRIM_DIR=trimmed            

bcftools_cmd="/usr/local/package/samtools/1.9/bin/bcftools"


bowtie2-build "${REF}" "${IDX}"
samtools faidx "${REF}"

mkdir -p bam_per_amplicon "${TRIM_DIR}"

# trim_galore → bowtie2 → BAM ====
for R1 in "${READ_DIR}"/v5-*_R1_*.fq.gz; do
  R2=${R1/_R1_/_R2_}
  SAMPLE=$(basename "${R1}" | sed 's/_R1_.*//')

  trim_galore --paired "${R1}" "${R2}" -o "${TRIM_DIR}"

  TR1="${TRIM_DIR}/$(basename "${R1}" .fq.gz)_val_1.fq.gz"
  TR2="${TRIM_DIR}/$(basename "${R2}" .fq.gz)_val_2.fq.gz"

  bowtie2 -p ${THREADS} --very-sensitive-local \
    -x "${IDX}" -1 "${TR1}" -2 "${TR2}" \
  | samtools view -@ ${THREADS} -b - \
  | samtools sort -@ ${THREADS} -o "bam_per_amplicon/${SAMPLE}.sorted.bam"

  samtools index "bam_per_amplicon/${SAMPLE}.sorted.bam"
done

# ==== merge, sort, index ====
samtools merge -@ ${THREADS} "${OUT}.merged.bam" bam_per_amplicon/*.sorted.bam
samtools sort  -@ ${THREADS} -o "${OUT}.merged.sorted.bam" "${OUT}.merged.bam"
samtools index "${OUT}.merged.sorted.bam"
rm -f "${OUT}.merged.bam"

# ==== trimming by iVar ====
ivar trim -i "${OUT}.merged.sorted.bam" -b "${BED}" -p "${OUT}.trimmed"
samtools sort -@ ${THREADS} -o "${OUT}.trimmed.sorted.bam" "${OUT}.trimmed.bam"
samtools index "${OUT}.trimmed.sorted.bam"

# ==== variant call ====
samtools depth -aa "${OUT}.trimmed.sorted.bam" > "${OUT}.depth.txt"

# ${bcftools_cmd} mpileup -f "${REF}" -Ou -a AD,DP -d 500000 "${OUT}.trimmed.sorted.bam" \
#   | ${bcftools_cmd} call -mv -Ob -o "bcf_tem"
#   | ${bcftools_cmd} filter -e 'QUAL<20 || DP<5' -Oz -o "${OUT}.variants.vcf.gz"
# ${bcftools_cmd} index "${OUT}.variants.vcf.gz"



${bcftools_cmd} mpileup -f "${REF}" -Ou -a AD,DP -d 500000 "${OUT}.trimmed.sorted.bam" \
  | ${bcftools_cmd} call -mv -Ob -o "bcf_tem"
${bcftools_cmd} index "bcf_tem"

${bcftools_cmd} view bcf_tem -Ou \
  | ${bcftools_cmd} filter -e 'QUAL<20 || FMT/DP<10' -Oz -o "${OUT}.variants.vcf.gz"



cat "${REF}" | ${bcftools_cmd} consensus -f - "${OUT}.variants.vcf.gz" > "${OUT}.consensus.fasta"

echo "[DONE]"
echo "  Final BAM:      ${OUT}.trimmed.sorted.bam"
echo "  Depth:          ${OUT}.depth.txt"
echo "  Variants:       ${OUT}.variants.vcf.gz"
echo "  Consensus:      ${OUT}.consensus.fasta"
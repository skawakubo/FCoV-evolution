#!/bin/sh
#$ -S /bin/bash
#$ -pe def_slot 1
#$ -l s_vmem=300G
#$ #-o /dev/null
#$ -cwd
#$ -l lmem



# set path for anaconda and activate conda environment
eval "$(/home/skawa/anaconda3/bin/conda shell.bash hook)"
conda activate FIP

# # make index in advance
# bwa_cmd="/usr/local/package/bwa/0.7.17/bin/bwa"
# ${bwa_cmd} index -p index_FCoV/FCoV_full_for_mappingidx index_FCoV/FCoV_full_for_mappingidx.fasta


run_id="v8"
OUT=${run_id}

READ_DIR=/home/skawa/Izumi_FIP/FCoV_seq_data  
# REF=/home/skawa/Izumi_FIP/FCoV_seq_data_analysis_results/amplicon_base_mapping/FCoV1_KY566209.fasta         # reference genome for bcftools consensus
IDX=/home/skawa/Izumi_FIP/amplicon_mapping/mapping/index_FCoV/FCoV_full_for_mappingidx              # bwa index prefix
MULTI_FASTA=/home/skawa/Izumi_FIP/amplicon_mapping/mapping/index_FCoV/FCoV_full_for_mappingidx.fasta
# OUT=all_amplicons_bt2_trim  
PRIMER_FA=/home/skawa/Izumi_FIP/amplicon_mapping/mapping/amplicon_primers.fa
# TRIM_Q=20                   
# TRIM_MINLEN=30              
TRIM_DIR=trimmed            
THREADS=8

bcftools_cmd="/usr/local/package/samtools/1.9/bin/bcftools"
bwa_cmd="/usr/local/package/bwa/0.7.17/bin/bwa"
samtools_cmd="/usr/local/package/samtools/1.9/bin/samtools"

mkdir -p ${run_id}
cd ${run_id}
mkdir -p bam_per_amplicon "${TRIM_DIR}"

# trim_galore → bwa → BAM ====
for R1 in "${READ_DIR}"/${run_id}-*_R1_*.fq.gz; do
  R2=${R1/_R1_/_R2_}
  # amplicon name
  RG_ID=$(basename "${R1}" | sed 's/_R1_.*//')

  trim_galore --paired "${R1}" "${R2}" -o "${TRIM_DIR}"

  CUT1="${TRIM_DIR}/$(basename "${R1}" .fq.gz)_val_1.fq.gz"
  CUT2="${TRIM_DIR}/$(basename "${R2}" .fq.gz)_val_2.fq.gz"
  CTO1="${TRIM_DIR}/$(basename "${R1}" .fq.gz)_val_1.prtrim.fq.gz"
  CTO2="${TRIM_DIR}/$(basename "${R2}" .fq.gz)_val_2.prtrim.fq.gz"

  echo "[INFO] Trimming primers (cutadapt, sequence-based): ${RG_ID}"
  cutadapt \
    -g file:"${PRIMER_FA}" -G file:"${PRIMER_FA}" \
    -a file:"${PRIMER_FA}" -A file:"${PRIMER_FA}" \
    --pair-filter=any --minimum-length 30 \
    -o "${CTO1}" -p "${CTO2}" "${CUT1}" "${CUT2}"

  echo "[INFO] Mapping to multi-ref index (bwa mem): ${RG_ID}"
  ${bwa_cmd} mem -t ${THREADS} \
    -R "@RG\tID:${RG_ID}\tSM:${run_id}\tPL:ILLUMINA" \
    ${IDX} "${CTO1}" "${CTO2}" \
  | ${samtools_cmd} view -@ ${THREADS} -b - \
  | ${samtools_cmd} sort -@ ${THREADS} -o "bam_per_amplicon/${RG_ID}.sorted.bam"

  ${samtools_cmd} index "bam_per_amplicon/${RG_ID}.sorted.bam"
done

# ==== merge, sort, index ====
${samtools_cmd} merge -@ ${THREADS} "${OUT}.merged.bam" bam_per_amplicon/*.sorted.bam
${samtools_cmd} sort  -@ ${THREADS} -o "${OUT}.merged.sorted.bam" "${OUT}.merged.bam"
${samtools_cmd} index "${OUT}.merged.sorted.bam"
rm -f "${OUT}.merged.bam"

# ==== choose best reference by coverage ====
# samtools coverage: column field 6  "coverage(%)" and 7 "mean coverage depth"
${samtools_cmd} coverage "${OUT}.merged.sorted.bam" > "${OUT}.coverage.txt"

best_ref=$(
  awk 'NR>1 {print $1, $6, $7}' "${OUT}.coverage.txt" \
  | sort -k2,2nr -k3,3nr \
  | head -1 | cut -d" " -f1
)
echo "[INFO] Best reference by coverage: ${best_ref}"

# extract best ref sequence from multi-FASTA
${samtools_cmd} faidx "${MULTI_FASTA}" "${best_ref}" > "${OUT}.bestref.fa"
REF="${OUT}.bestref.fa"


# extract reads mapped to the best reference
${samtools_cmd} view -@ ${THREADS} -b "${OUT}.merged.sorted.bam" "${best_ref}" > "${OUT}.bestref.bam"
${samtools_cmd} sort -@ ${THREADS} -o "${OUT}.bestref.sorted.bam" "${OUT}.bestref.bam"
${samtools_cmd} index "${OUT}.bestref.sorted.bam"
rm -f "${OUT}.bestref.bam"


# # ==== trimming by iVar ====
# ivar trim -i "${OUT}.bestref.sorted.bam" -b "${BED_USED}" -p "${OUT}.trimmed"
# ${samtools_cmd} sort -@ ${THREADS} -o "${OUT}.trimmed.sorted.bam" "${OUT}.trimmed.bam"
# ${samtools_cmd} index "${OUT}.trimmed.sorted.bam"
# rm -f "${OUT}.trimmed.bam"

# ==== depth ====
${samtools_cmd} depth -aa "${OUT}.trimmed.sorted.bam" > "${OUT}.depth.txt"

# ${bcftools_cmd} mpileup -f "${REF}" -Ou -a AD,DP -d 500000 "${OUT}.trimmed.sorted.bam" \
#   | ${bcftools_cmd} call -mv -Ob -o "bcf_tem"
#   | ${bcftools_cmd} filter -e 'QUAL<20 || DP<5' -Oz -o "${OUT}.variants.vcf.gz"
# ${bcftools_cmd} index "${OUT}.variants.vcf.gz"

# ==== variant call ====
${bcftools_cmd} mpileup -f "${REF}" -Ou \
  -a AD,DP,ADF,ADR \
  -d 500000 -q 20 -Q 20 \
  "${OUT}.trimmed.sorted.bam" \
| ${bcftools_cmd} call -mv --ploidy 1 -Ob -o "${OUT}.raw.bcf"

${bcftools_cmd} index "${OUT}.raw.bcf"

${bcftools_cmd} view "${OUT}.raw.bcf" -Ou \
| ${bcftools_cmd} filter -e 'QUAL<20 || FMT/DP<10' -Oz -o "${OUT}.variants.vcf.gz"

${bcftools_cmd} index "${OUT}.variants.vcf.gz"

# mask low coverage region (<10x) as N in consensus
awk '$3<10 {print $1"\t"$2-1"\t"$2}' "${OUT}.depth.txt" > "${OUT}.lowcov.bed"

# ==== write out consensus ====
${bcftools_cmd} consensus -f "${REF}" \
  --mask "${OUT}.lowcov.bed" \
  "${OUT}.variants.vcf.gz" > "${OUT}.consensus.fasta"

echo "[DONE]"
echo "  BestRef:        ${best_ref}"
echo "  BestRef FASTA:  ${REF}"
echo "  Coverage table: ${OUT}.coverage.txt"
echo "  Final BAM:      ${OUT}.trimmed.sorted.bam"
echo "  Depth:          ${OUT}.depth.txt"
echo "  Variants:       ${OUT}.variants.vcf.gz"
echo "  Lowcov BED:     ${OUT}.lowcov.bed"
echo "  Consensus:      ${OUT}.consensus.fasta"


#!/bin/sh
#$ -S /bin/bash
#$ -pe def_slot 1
#$ -l s_vmem=300G
#$ #-o /dev/null
#$ -cwd
#$ -l lmem

fasta_name="Alphacov1_ntblast_seq_cleaned_wg.fasta"
output_name="Alphacov1_ntblast_seq_cleaned_wg_mafftauto.fasta"

mafft_cmd="/usr/local/package/mafft/7.522/bin/mafft"

${mafft_cmd} --thread 8 --auto ${fasta_name} > ${output_name}

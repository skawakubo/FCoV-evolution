#!/bin/sh
#$ -S /bin/bash
#$ -pe def_slot 4
#$ -l s_vmem=480G
#$ #-o /dev/null
#$ -cwd
#$ -l lmem

#specify maximum run time to 2 month; 62 days; 1488 hours
#$ #-l d_rt=1488:00:00
#$ #-l s_rt=1488:00:00

# set path for anaconda
eval "$(/home/skawa/anaconda3/bin/conda shell.bash hook)"

# activate conda environment
conda activate blast

# run blast
blastn -db /usr/local/db/blast/ncbi/nt -query FCoV_full_NCBI_Izumi_250726.fasta -max_target_seqs 1000 -evalue 0.000001 -outfmt "6 qseqid qlen sseqid sacc slen pident length mismatch gapopen qstart qend sstart send qcovs evalue bitscore staxids sscinames scomnames" -out FCoV_full_NCBI_Izumi_250726_ntblast.out -num_threads 20

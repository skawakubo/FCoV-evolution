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
blastn -db /usr/local/db/blast/ncbi/nt -query FCoV_CCoV_FCoV23_TGEV_FRCoV_MCoV_full_nt_250818_cleaned_clu90.fasta -max_target_seqs 5000 -evalue 0.000001 -outfmt "6 qseqid qlen sseqid sacc slen pident length mismatch gapopen qstart qend sstart send qcovs evalue bitscore staxids sscinames scomnames" -out Alphacov1_ntblast.out -num_threads 40

#!/bin/bash
# 各サイトの塩基組成を「数えるだけ」で出力する（変異コールはしない）
set -euo pipefail

sample_ID="fecv_NARC3_S3"
ref="TRINITY_FCoV_FECV3_assembled_2.fas"
bam="FECV_NARC3_S3_mapped_on_assembled_v2.sort.bam"

samtools_cmd="/usr/local/package/samtools/1.9/bin/samtools"

mkdir -p tsv

# -B : BAQ無効化（IGVはBAQをかけない。ここが最大の差の原因）
# -Q 0: 塩基品質による除外なし（デフォルトは13で切られている）
# -q 0: マッピング品質による除外なし
# -d 1000000: 深度上限を実質解除（デフォルト8000）
# -a : カバレッジ0の位置も出力
# --ff UNMAP,SECONDARY,QCFAIL,DUP : IGVのデフォルト（重複除去ON）に合わせる
#      本当に全リードを見たいなら --ff UNMAP のみにする

${samtools_cmd} mpileup -f ${ref} \
    -B -Q 0 -q 0 -d 1000000 -a \
    --ff UNMAP,SECONDARY,QCFAIL,DUP \
    ${bam} \
  | python3 pileup_base_counts.py \
  > tsv/${sample_ID}.base_counts.tsv

# 2種類以上の塩基が観測されたサイトだけ抜き出す
python3 - <<'EOF'
import csv, sys
src = "tsv/fecv_NARC3_S3.base_counts.tsv"
dst = "tsv/fecv_NARC3_S3.polymorphic.tsv"
with open(src) as fi, open(dst, "w") as fo:
    r = csv.DictReader(fi, delimiter="\t")
    w = csv.DictWriter(fo, fieldnames=r.fieldnames, delimiter="\t")
    w.writeheader()
    n = 0
    for row in r:
        obs = sum(1 for b in "ACGT" if int(row[b]) > 0)
        if obs >= 2:
            w.writerow(row); n += 1
print(f"{n} polymorphic sites -> {dst}", file=sys.stderr)
EOF

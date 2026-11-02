#!/usr/bin/env python3
"""
samtools mpileup のテキスト出力を読んで、各サイトの塩基組成をTSVで出力する。
統計モデルは一切使わず、マップされたリードの塩基を数えるだけ。

usage:
  samtools mpileup -f ref.fas -B -q 0 -Q 0 -d 1000000 -a in.bam \
      | python3 pileup_base_counts.py > base_counts.tsv

  # 2種類以上の塩基が観測されたサイトだけ
  ... | python3 pileup_base_counts.py --min-alt-count 1 --min-alt-freq 0.0
"""
import argparse
import re
import sys

INDEL_RE = re.compile(r'^([+-])(\d+)')
BASES = ['A', 'C', 'G', 'T', 'N']


def parse_pileup_bases(s, ref_base):
    """mpileup の read base 文字列をパースして塩基カウントを返す。"""
    counts = {b: 0 for b in BASES}
    counts_fwd = {b: 0 for b in BASES}
    n_del = 0        # このリードでは欠失している（* / #）
    n_refskip = 0    # スプライス等のスキップ（> / <）
    ins_events = {}  # 挿入配列 -> 回数
    del_events = {}  # 欠失配列 -> 回数

    ref_base = ref_base.upper()
    i = 0
    L = len(s)
    while i < L:
        c = s[i]
        if c == '^':          # リード開始。次の1文字はMAPQ
            i += 2
            continue
        if c == '$':          # リード終端
            i += 1
            continue
        if c in '+-':         # インデル: +2AG / -3TTA
            m = INDEL_RE.match(s[i:])
            n = int(m.group(2))
            start = i + 1 + len(m.group(2))
            seq = s[start:start + n].upper()
            if c == '+':
                ins_events[seq] = ins_events.get(seq, 0) + 1
            else:
                del_events[seq] = del_events.get(seq, 0) + 1
            i = start + n
            continue
        if c in '.,':         # リファレンスと一致
            if ref_base in counts:
                counts[ref_base] += 1
                if c == '.':
                    counts_fwd[ref_base] += 1
            i += 1
            continue
        if c in '*#':         # 上流の欠失に覆われている
            n_del += 1
            i += 1
            continue
        if c in '><':         # reference skip
            n_refskip += 1
            i += 1
            continue
        u = c.upper()
        if u in counts:       # ミスマッチ塩基
            counts[u] += 1
            if c.isupper():
                counts_fwd[u] += 1
            i += 1
            continue
        i += 1                # 想定外の文字は無視
    return counts, counts_fwd, n_del, n_refskip, ins_events, del_events


def fmt_events(d):
    if not d:
        return '.'
    return ';'.join('%s:%d' % (k, v) for k, v in
                    sorted(d.items(), key=lambda x: -x[1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-alt-count', type=int, default=0,
                    help='非リファレンス塩基の最低リード数 (default: 0 = 全サイト出力)')
    ap.add_argument('--min-alt-freq', type=float, default=0.0,
                    help='非リファレンス塩基の最低頻度 (default: 0.0)')
    ap.add_argument('--min-depth', type=int, default=1)
    args = ap.parse_args()

    out = sys.stdout
    out.write('\t'.join([
        'chrom', 'pos', 'ref', 'depth_raw', 'depth_base',
        'A', 'C', 'G', 'T', 'N', 'del', 'refskip',
        'alt', 'alt_count', 'alt_freq',
        'A_fwd', 'C_fwd', 'G_fwd', 'T_fwd',
        'ins_events', 'del_events',
    ]) + '\n')

    for line in sys.stdin:
        f = line.rstrip('\n').split('\t')
        if len(f) < 5:
            continue
        chrom, pos, ref_base, depth_raw = f[0], f[1], f[2].upper(), int(f[3])
        bases = f[4] if len(f) > 4 else ''

        counts, fwd, n_del, n_skip, ins_ev, del_ev = parse_pileup_bases(bases, ref_base)
        depth_base = sum(counts.values())

        if depth_base < args.min_depth:
            continue

        # 最頻の非リファレンス塩基
        alts = [(b, c) for b, c in counts.items() if b != ref_base and c > 0]
        alts.sort(key=lambda x: -x[1])
        if alts:
            alt, alt_count = alts[0]
            alt = ','.join('%s' % b for b, _ in alts)
        else:
            alt, alt_count = '.', 0
        alt_freq = alt_count / depth_base if depth_base else 0.0

        if alt_count < args.min_alt_count:
            continue
        if alt_freq < args.min_alt_freq:
            continue

        out.write('\t'.join([
            chrom, pos, ref_base, str(depth_raw), str(depth_base),
            str(counts['A']), str(counts['C']), str(counts['G']),
            str(counts['T']), str(counts['N']), str(n_del), str(n_skip),
            alt, str(alt_count), '%.4f' % alt_freq,
            str(fwd['A']), str(fwd['C']), str(fwd['G']), str(fwd['T']),
            fmt_events(ins_ev), fmt_events(del_ev),
        ]) + '\n')


if __name__ == '__main__':
    main()

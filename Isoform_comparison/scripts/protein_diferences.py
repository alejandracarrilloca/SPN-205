"""
Author:
Dulce Alejandra Carrillo Carlos
Project: SPN-205 — Nonsense-mediated decay isoform expression analysis

Compares aligned protein sequences pairwise to identify and record regions of differences
between sequences, outputting the start and end positions along with the differing residues.
"""

from Bio import AlignIO
from itertools import combinations

aln_file = "results/alignments/spn.prot.mafft.fa"
aln = AlignIO.read(aln_file, "fasta")

def diff_ranges(seq1, seq2):
    ranges = []
    start = None
    current_diffs = []
    for i, (a, b) in enumerate(zip(seq1.seq, seq2.seq), start=1):
        if a != b:
            if start is None:
                start = i
            current_diffs.append((i, a, b))
        else:
            if start is not None:
                end = current_diffs[-1][0]
                ranges.append((start, end, ''.join([d[1] for d in current_diffs]),
                               ''.join([d[2] for d in current_diffs])))
                start = None
                current_diffs = []
    if start is not None:
        end = current_diffs[-1][0]
        ranges.append((start, end, ''.join([d[1] for d in current_diffs]),
                       ''.join([d[2] for d in current_diffs])))
    return ranges

out_file = "results/alignments/protein_comparison.tsv"
with open(out_file, "w") as out:
    out.write("transcript1\ttranscript2\tstart\tend\tseq1_diff\tseq2_diff\n")
    for a, b in combinations(aln, 2):
        ranges = diff_ranges(a, b)
        for start, end, s1, s2 in ranges:
            out.write(f"{a.id}\t{b.id}\t{start}\t{end}\t{s1}\t{s2}\n")
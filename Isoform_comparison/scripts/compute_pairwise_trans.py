"""
Author:
Dulce Alejandra Carrillo Carlos
Project: SPN-205 — Nonsense-mediated decay isoform expression analysis

Calculates pairwise percent identity between aligned cDNA sequences from a multiple sequence alignment.
"""

from Bio import AlignIO
from itertools import combinations

aln_file = "results/alignments/spn.cdna.mafft.fa"
aln = AlignIO.read(aln_file, "fasta")

def pid_aligned(seq1, seq2):
    s1 = str(seq1.seq)
    s2 = str(seq2.seq)
    assert len(s1) == len(s2), "Sequences must be aligned to same length"
    matches = sum(1 for a, b in zip(s1, s2) if a == b)
    return matches / len(s1) * 100

with open("results/alignments/cdna_pairwise.tsv", "w") as out:
    out.write("transcript1\ttranscript2\tpercent_identity\n")
    for a, b in combinations(aln, 2):
        out.write(f"{a.id}\t{b.id}\t{pid_aligned(a,b):.1f}\n")

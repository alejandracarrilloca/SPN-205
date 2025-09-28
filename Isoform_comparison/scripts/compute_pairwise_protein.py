from Bio import AlignIO
from itertools import combinations

# Load aligned protein sequences
aln_file = "results/alignments/spn.prot.mafft.fa"  # This must be an aligned FASTA
aln = AlignIO.read(aln_file, "fasta")

def pid_aligned(seq1, seq2):
    """Percent identity for aligned sequences including gaps"""
    s1 = str(seq1.seq)
    s2 = str(seq2.seq)
    assert len(s1) == len(s2), "Sequences must be aligned to same length"
    matches = sum(1 for a, b in zip(s1, s2) if a == b)
    return matches / len(s1) * 100

# Compute pairwise PID
with open("results/alignments/protein_pairwise.tsv", "w") as out:
    out.write("transcript1\ttranscript2\tpercent_identity\n")
    for a, b in combinations(aln, 2):
        out.write(f"{a.id}\t{b.id}\t{pid_aligned(a,b):.1f}\n")

print("Pairwise percent identity matrix saved to results/alignments/protein_pairwise.tsv")

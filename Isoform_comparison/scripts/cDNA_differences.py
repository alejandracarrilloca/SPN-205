from Bio import AlignIO
from itertools import combinations

# Load protein alignment
aln_file = "results/alignments/spn.cdna.mafft.fa"
aln = AlignIO.read(aln_file, "fasta")

# Function to find ranges of differences between two aligned sequences
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
    # catch trailing range
    if start is not None:
        end = current_diffs[-1][0]
        ranges.append((start, end, ''.join([d[1] for d in current_diffs]),
                       ''.join([d[2] for d in current_diffs])))
    return ranges

# Output file
out_file = "results/alignments/cDNA_comparison.tsv"
with open(out_file, "w") as out:
    out.write("transcript1\ttranscript2\tstart\tend\tseq1_diff\tseq2_diff\n")
    for a, b in combinations(aln, 2):
        ranges = diff_ranges(a, b)
        for start, end, s1, s2 in ranges:
            out.write(f"{a.id}\t{b.id}\t{start}\t{end}\t{s1}\t{s2}\n")

print(f"All pairwise protein difference ranges saved to {out_file}")

#Translate the sequence 

from Bio import SeqIO
from Bio.Seq import Seq
import sys

# Input and output files
input_fasta = "results/sequences/spn.cds.fa"
output_fasta = "results/sequences/spn.prot.fa"

# Open output file
with open(output_fasta, "w") as out_f:
    # Parse input CDS FASTA
    for rec in SeqIO.parse(input_fasta, "fasta"):
        seq = str(rec.seq).replace("\n", "")
        if not seq:
            continue  # skip empty sequences
        try:
            # Translate to protein (stop codons become '*')
            prot = str(Seq(seq).translate(to_stop=False))
            out_f.write(f">{rec.id}\n{prot}\n")
        except Exception as e:
            print(f"Error translating {rec.id}: {e}", file=sys.stderr)

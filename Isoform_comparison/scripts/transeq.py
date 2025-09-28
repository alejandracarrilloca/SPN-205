"""
Author:
Dulce Alejandra Carrillo Carlos
Project: SPN-205 — Nonsense-mediated decay isoform expression analysis

Translates coding DNA sequences (CDS) from a FASTA file into protein sequences,
writing the resulting amino acid sequences to a new FASTA file.
"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys

input_fasta = "results/sequences/spn.cds.fa"
output_fasta = "results/sequences/spn.prot.fa"

with open(output_fasta, "w") as out_f:
    for rec in SeqIO.parse(input_fasta, "fasta"):
        seq = str(rec.seq).replace("\n", "")
        if not seq:
            continue
        try:
            prot = str(Seq(seq).translate(to_stop=False))
            out_f.write(f">{rec.id}\n{prot}\n")
        except Exception as e:
            print(f"Error translating {rec.id}: {e}", file=sys.stderr)

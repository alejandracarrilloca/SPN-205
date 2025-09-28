
"""
Compute lengths of cDNA, CDS, protein sequences, and exon counts for SPN transcripts.
Output: results/sequences/spn_summary.tsv
Columns: transcript_id, cdna_len, cds_len, prot_len, exon_count
"""

from Bio import SeqIO
import pandas as pd
import subprocess

# Paths
cdna_fasta = "results/sequences/spn.cdna.fa"
cds_fasta  = "results/sequences/spn.cds.fa"
prot_fasta = "results/sequences/spn.prot.fa"
gtf_file   = "data/spn.gtf"
output_file = "results/sequences/spn_summary.tsv"

# --- Function to compute lengths from fasta ---
def fa_len(fasta_file):
    return {r.id: len(r.seq) for r in SeqIO.parse(fasta_file, "fasta")}

cdna_len = fa_len(cdna_fasta)
cds_len  = fa_len(cds_fasta)
prot_len = fa_len(prot_fasta)

# --- Extract transcript IDs from cDNA fasta ---
transcript_ids = list(cdna_len.keys())

# --- Compute exon counts from GTF ---
# Using bash commands via subprocess
cmd = f"""awk '$3=="exon"' {gtf_file} | sed 's/;/\\n/g' | grep transcript_id | awk '{{print $2}}' | sed 's/"//g' | sort | uniq -c | awk '{{print $2"\\t"$1}}'"""
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
exon_counts = {}
for line in result.stdout.strip().split("\n"):
    tid, count = line.split("\t")
    exon_counts[tid] = int(count)

# --- Create summary DataFrame ---
df = pd.DataFrame({
    "transcript_id": transcript_ids,
    "cdna_len": [cdna_len.get(tid,0) for tid in transcript_ids],
    "cds_len": [cds_len.get(tid,0) for tid in transcript_ids],
    "prot_len": [prot_len.get(tid,0) for tid in transcript_ids],
    "exon_count": [exon_counts.get(tid,0) for tid in transcript_ids]
})

# --- Save output ---
df.to_csv(output_file, sep="\t", index=False)
print(f"Written: {output_file}")
print(df.head())

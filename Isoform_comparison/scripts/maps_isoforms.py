
from Bio import SeqIO
import pandas as pd
from collections import defaultdict
import os

# Input files
gtf_file = "data/spn.gtf"
genome_file = "data/GRCh38.primary_assembly.genome.fa"
out_dir = "results/transcript_maps"

os.makedirs(out_dir, exist_ok=True)

# Load genome
genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

# Parse GTF
trans_features = defaultdict(list)
with open(gtf_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        cols = line.strip().split("\t")
        chrom, source, feature, start, end, score, strand, frame, attr = cols
        start, end = int(start), int(end)
        tid = None
        for part in attr.split(";"):
            if "transcript_id" in part:
                tid = part.split('"')[1]
        if tid:
            trans_features[tid].append({
                "feature": feature,
                "start": start,
                "end": end,
                "strand": strand
            })

# Process each transcript
for tid, feats in trans_features.items():
    # Order features by genomic start
    feats_sorted = sorted(feats, key=lambda x: x["start"])
    
    rows = []
    for f in feats_sorted:
        chrom_seq = genome[f["chrom"]].seq
        seq = chrom_seq[f["start"]-1:f["end"]]  # 1-based to 0-based
        if f["strand"] == "-":
            seq = seq.reverse_complement()
        rows.append({
            "feature": f["feature"],
            "genomic_start": f["start"],
            "genomic_end": f["end"],
            "sequence": str(seq)
        })
    
    df = pd.DataFrame(rows)
    df.to_csv(f"{out_dir}/{tid}_map.tsv", sep="\t", index=False)

print(f"Transcript maps saved to {out_dir}/")

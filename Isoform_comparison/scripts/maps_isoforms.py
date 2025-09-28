"""
Author:
Dulce Alejandra Carrillo Carlos
Project: SPN-205 — Nonsense-mediated decay isoform expression analysis

Extracts genomic sequences for each transcript from a GTF annotation and reference genome,
organizing features (exons, CDS, etc.) per transcript and saving them as TSV files.
"""

from Bio import SeqIO
import pandas as pd
from collections import defaultdict
import os

gtf_file = "data/spn.gtf"
genome_file = "data/GRCh38.primary_assembly.genome.fa"
out_dir = "results/transcript_maps"

os.makedirs(out_dir, exist_ok=True)

genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

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

for tid, feats in trans_features.items():
    feats_sorted = sorted(feats, key=lambda x: x["start"])
    
    rows = []
    for f in feats_sorted:
        chrom_seq = genome[f["chrom"]].seq
        seq = chrom_seq[f["start"]-1:f["end"]]
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
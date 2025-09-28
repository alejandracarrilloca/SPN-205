"""
Author:
Dulce Alejandra Carrillo Carlos
Project: SPN-205 — Nonsense-mediated decay isoform expression analysis

Description:
Filter transcript-level expression matrices (SPN-201 to SPN-205) 
based on TCGA sample ID files and create a single 
vertically formatted file per transcript.

Each sample entry includes:
- Sample ID
- Expression value (TPM)
- Cancer type (from file name)
- Tissue type: assumed "Tumor"
- Transcript: defined explicitly in the list

Inputs:
1. Expression matrix (TSV format):
    - Path: results/Transcripts/<transcript>/<transcript>_expression.tsv
    - Format: rows = transcripts, columns = sample IDs (TCGA-style)
    - Assumes only one transcript (one row)

2. Metadata files with sample IDs to extract:
    - Folder defined in `ids_folder`
    - Each file should:
        - End with `_ids.txt`
        - Begin with 4-letter cancer code (e.g., PAAD_ids.txt)
        - Contain sample IDs (first column)

Outputs:
- One combined file per transcript: <transcript>_TCGA_expression.tsv
- Format:
    sample_id | TPM | cancer_type | tissue_type | transcript
"""

import os
import pandas as pd

ids_folder = "results/IDS/TCGA_IDS"
transcripts_folder = "results/Transcripts"

transcripts = ["SPN-201", "SPN-202", "SPN-203", "SPN-204", "SPN-205"]

for transcript_name in transcripts:
    transcript_path = os.path.join(transcripts_folder, transcript_name)
    expression_file = os.path.join(transcript_path, f"{transcript_name}_expression.tsv")
    output_file = os.path.join(transcript_path, f"{transcript_name}_TCGA_expression.tsv")

    if not os.path.exists(expression_file):
        print(f"Expression file not found for {transcript_name}")
        continue

    expression_df = pd.read_csv(expression_file, sep="\t", index_col=0)

    combined_data = []

    for fname in os.listdir(ids_folder):
        if fname.endswith("_ids.txt") and fname[:4].isupper():
            cancer_code = fname[:4]
            id_file_path = os.path.join(ids_folder, fname)

            with open(id_file_path) as f:
                sample_ids = [line.strip().split('\t')[0] for line in f if line.strip()]

            filtered_df = expression_df.loc[:, expression_df.columns.intersection(sample_ids)]

            if filtered_df.empty:
                print(f"⚠ No matching samples for {cancer_code} in {transcript_name}")
                continue

            for sample in filtered_df.columns:
                tpm = filtered_df.at[expression_df.index[0], sample]
                combined_data.append({
                    "sample_id": sample,
                    "TPM": tpm,
                    "cancer_type": cancer_code,
                    "tissue_type": "Tumor",
                    "transcript": transcript_name
                })

            print(f"Processed {cancer_code} ({len(filtered_df.columns)} samples) for {transcript_name}")

    combined_df = pd.DataFrame(combined_data)
    if not combined_df.empty:
        combined_df.to_csv(output_file, sep="\t", index=False)
        print(f"\nFile saved: {output_file}\n")
"""
Author:
Dulce Alejandra Carrillo Carlos
Project: SPN-205 — Nonsense-mediated decay isoform expression analysis

Filter transcript-level expression matrices (SPN-201 to SPN-205) 
based on GTEx sample ID file and create a single 
vertically formatted file per transcript.
"""

import os
import pandas as pd

ids_folder = "results/IDS/GTEx_IDS"
transcripts_folder = "results/Transcripts"

transcripts = ["SPN-201", "SPN-202", "SPN-203", "SPN-204", "SPN-205"]

for transcript_name in transcripts:
    transcript_path = os.path.join(transcripts_folder, transcript_name)

    expr_files = [f for f in os.listdir(transcript_path) if f.endswith("_expression.tsv")]

    expression_file = os.path.join(transcript_path, expr_files[0])
    output_file = os.path.join(transcript_path, f"{transcript_name}_GTEx_expression.tsv")

    expression_df = pd.read_csv(expression_file, sep="\t", index_col=0)

    combined_data = []

    for fname in os.listdir(ids_folder):
        if fname.endswith("_ids.txt") and fname[:4].isupper():
            cancer_code = fname[:4]
            id_file_path = os.path.join(ids_folder, fname)

            sample_info = []
            with open(id_file_path) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) > 3:
                        sample_info.append((parts[0], parts[3]))
                    else:
                        sample_info.append((parts[0], "Unknown"))

            sample_to_tissue = dict(sample_info)

            valid_samples = set(expression_df.columns).intersection(sample_to_tissue.keys())
            filtered_df = expression_df.loc[:, list(valid_samples)]

            if filtered_df.empty:
                print(f"No matching samples for {cancer_code} in {transcript_name}")
                continue

            for sample in filtered_df.columns:
                tpm = filtered_df.at[expression_df.index[0], sample]
                tissue = sample_to_tissue.get(sample, "Unknown")
                combined_data.append({
                    "sample_id": sample,
                    "TPM": tpm,
                    "cancer_type": cancer_code,
                    "tissue_type": tissue,
                    "transcript": transcript_name
                })

            print(f"Processed {cancer_code} ({len(filtered_df.columns)} samples) for {transcript_name}")

    combined_df = pd.DataFrame(combined_data)
    if not combined_df.empty:
        combined_df.to_csv(output_file, sep="\t", index=False)
        print(f"\nGTEx file with tissue saved to: {output_file}\n")
"""
Author:
Dulce Alejandra Carrillo Carlos
Project: SPN-205 — Nonsense-mediated decay isoform expression analysis

Merges isoform-level expression data of SPN-201
and SPN-202 (principal SPN isoforms) into a single "canonical SPN" expression file for TCGA data.

"""

import pandas as pd

file1 = "results/Transcripts/SPN-201/SPN-201_TCGA_expression.tsv"
file2 = "results/Transcripts/SPN-202/SPN-202_TCGA_expression.tsv"
output_file = "results/Transcripts/SPN/SPN_TCGA_expression.tsv"

df1 = pd.read_csv(file1, sep="\t")
df2 = pd.read_csv(file2, sep="\t")
assert set(df1["sample_id"]) == set(df2["sample_id"]), "Las muestras no coinciden"
df_merged = pd.merge(df1, df2, on=["sample_id", "cancer_type", "tissue_type"], suffixes=('_201', '_202'))
df_merged["TPM"] = df_merged["TPM_201"] + df_merged["TPM_202"]
df_merged["Transcript"] = "SPN"
cols_to_save = ["sample_id", "TPM", "cancer_type", "tissue_type", "Transcript"]
df_merged[cols_to_save].to_csv(output_file, sep="\t", index=False)

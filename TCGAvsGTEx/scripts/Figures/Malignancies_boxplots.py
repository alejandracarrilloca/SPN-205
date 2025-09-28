"""
Author:
Dulce Alejandra Carrillo Carlos
Project: SPN-205 — Nonsense-mediated decay isoform expression analysis

Generate boxplots comparing expression levels of SPN transcripts (SPN and SPN-205)
between tumor (TCGA) and normal (GTEx) samples across multiple cancer types.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.patches import Patch

transcripts = ["SPN", "SPN-205"]

colors = {"Tumor": "#FC4E92", "Normal": "#7DEE75"}

base_dir = "results/Transcripts"
figures_dir = f"results/figures/Malignancies_boxplots"
os.makedirs(figures_dir, exist_ok=True)

df_all = []

for tr in transcripts:
    tumor_file = f"{base_dir}/{tr}/{tr}_TCGA_expression.tsv"
    normal_file = f"{base_dir}/{tr}/{tr}_GTEx_expression.tsv"

    tumor_df = pd.read_csv(tumor_file, sep="\t")
    normal_df = pd.read_csv(normal_file, sep="\t")

    tumor_df["tissue_type"] = "Tumor"
    normal_df["tissue_type"] = "Normal"

    df = pd.concat([tumor_df, normal_df], ignore_index=True)
    df["transcript"] = tr
    df_all.append(df)

df_all = pd.concat(df_all, ignore_index=True)
df_all["tissue_type"] = df_all["tissue_type"].str.strip().str.capitalize()
df_all["cancer_type"] = df_all["cancer_type"].str.strip()

df_all["cancer_type"] = df_all["cancer_type"].str[:4]

df_all["cancer_type"] = df_all["cancer_type"].str.rstrip("_")

for cancer in sorted(df_all["cancer_type"].unique()):
    sub_df = df_all[df_all["cancer_type"] == cancer]
    if sub_df.empty:
        continue

    positions = []
    data = []
    labels = []
    colors_list = []
    offset = 0

    for tr in transcripts:
        for tissue in ["Tumor", "Normal"]:
            values = sub_df[(sub_df["transcript"] == tr) & 
                            (sub_df["tissue_type"] == tissue)]["TPM"].values
            if len(values) == 0:
                continue
            data.append(values)
            positions.append(offset)
            labels.append(tr)
            colors_list.append(colors[tissue])
            offset += 1

    if len(data) < 2:
        print(f"{cancer}: insufficient data to plot.")
        continue

    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(5,5))

    # Boxplot
    box = ax.boxplot(data,
                     positions=positions,
                     patch_artist=True,
                     widths=0.6,
                     showfliers=False,
                     medianprops=dict(color='black', linewidth=1.2),
                     boxprops=dict(linewidth=1.5),
                     whiskerprops=dict(linewidth=1.5),
                     capprops=dict(linewidth=1.5))

    for patch, color in zip(box["boxes"], colors_list):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    for pos, values, color in zip(positions, data, colors_list):
        x_jitter = np.random.normal(pos, 0.05, len(values))
        ax.scatter(x_jitter, values, alpha=0.6, s=25, color=color,
                   edgecolors='black', linewidth=0.5)

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=12)
    ax.tick_params(axis='y', labelsize=12)
    ax.set_ylabel("log2(TPM+0.001)", fontsize=10, fontweight='bold')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)

    legend_elements = [Patch(facecolor=colors["Tumor"], alpha=0.7, label="Tumor"),
                       Patch(facecolor=colors["Normal"], alpha=0.7, label="Normal")]
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1.15),
              title='Tissue Type', fontsize=6, title_fontsize=8)

    ax.set_title(f"{cancer}", fontsize=12, fontweight="bold", pad=25)

    plt.tight_layout()
    save_path = os.path.join(figures_dir, f"{cancer}_comparison_boxplot.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
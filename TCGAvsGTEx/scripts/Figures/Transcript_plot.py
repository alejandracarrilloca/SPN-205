"""
Author:
Dulce Alejandra Carrillo Carlos
Project: SPN-205 — Nonsense-mediated decay isoform expression analysis

Compare expression levels of SPN gene transcripts between tumor and normal tissues
across multiple cancer types. Using the Mann-Whitney U test (Wilcoxon rank-sum) for each
transcript-cancer pair and visualize the results in a dot plot with significance and log2 fold 
change.

"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import mannwhitneyu
import seaborn as sns

transcripts = ["SPN-201", "SPN-202", "SPN-203", "SPN-204", "SPN-205"]

base_dir = "results/Transcripts"
figures_dir = "results/figures"
output_path = f"{base_dir}/wilcoxon_test.tsv"
os.makedirs(figures_dir, exist_ok=True)

results = []

for tr in transcripts:
    tumor_file = f"{base_dir}/{tr}/{tr}_TCGA_expression.tsv"
    normal_file = f"{base_dir}/{tr}/{tr}_GTEx_expression.tsv"

    tumor_df = pd.read_csv(tumor_file, sep="\t")
    normal_df = pd.read_csv(normal_file, sep="\t")
    normal_df["tissue_type"] = "Normal"

    combined_df = pd.concat([tumor_df, normal_df], ignore_index=True)
    combined_df["tissue_type"] = combined_df["tissue_type"].str.strip().str.capitalize()
    combined_df["cancer_type"] = combined_df["cancer_type"].str.strip()
    combined_df["transcript"] = tr

    for cancer in combined_df["cancer_type"].unique():
        sub = combined_df[combined_df["cancer_type"] == cancer]
        tumor_expr = sub[sub["tissue_type"] == "Tumor"]["TPM"]
        normal_expr = sub[sub["tissue_type"] == "Normal"]["TPM"]

        if len(tumor_expr) < 3 or len(normal_expr) < 3:
            print(f"{tr} - {cancer}: insufficient samples")
            continue

        stat, p_value = mannwhitneyu(tumor_expr, normal_expr, alternative="two-sided")

        results.append({
            "transcript": tr,
            "cancer_type": cancer,
            "n_tumor": len(tumor_expr),
            "mean_tumor": tumor_expr.mean(),
            "n_normal": len(normal_expr),
            "mean_normal": normal_expr.mean(),
            "p_value": p_value
        })

results_df = pd.DataFrame(results)
results_df["log2FC"] = results_df["mean_tumor"] - results_df["mean_normal"]
results_df["significance"] = 1 - results_df["p_value"]

threshold = results_df["log2FC"].abs().quantile(0.8)
results_df["extreme"] = results_df["log2FC"].abs() >= threshold

cancer_order = [
    "CHOL", "LIHC", "COAD", "READ", "STAD", "PAAD", "ESCA",
    "BLCA", "PRAD", "TGCT", "UCEC", "UCS", "CESC",
    "ACC", "PCPG", "THCA",
    "BRCA",
    "LUAD", "LUSC",
    "KICH", "KIRC", "KIRP",
    "GBM", "LGG",
    "LAML", "DLBC", "THYM",
    "HNSC", "SKCM"
]

results_df["cancer_code"] = results_df["cancer_type"].str[:4]

results_df["cancer_code"] = results_df["cancer_code"].str.rstrip("_")

results_df_plot = results_df[results_df["cancer_code"].isin(cancer_order)].copy()

cancer_to_x = {c: i for i, c in enumerate(cancer_order)}
transcripts_order = results_df_plot["transcript"].unique()
y_pos = np.arange(len(transcripts_order))
transcript_to_y = {tr: i for i, tr in enumerate(transcripts_order)}

plt.figure(figsize=(9,2))
ax = plt.gca()

row_colors = {
    0: "#75AAFF",
    1: "#75AAFF",
    2: "#FC4E92",
    3: "#B94EFC",
    4: "#FDD535"
}

ax.set_xlim(-0.5, len(cancer_order)-0.5)
ax.set_ylim(-0.5, len(transcripts_order)-0.5)

ax.add_patch(plt.Rectangle(
    (-0.5, -0.5),
    len(cancer_order),
    len(transcripts_order),
    fill=True,
    color='white',
    edgecolor='black',
    linewidth=1,
    zorder=-1
))

for i, tr in enumerate(transcripts_order):
    color = row_colors.get(i, "white")
    ax.add_patch(plt.Rectangle(
        (-0.5, i-0.5),
        len(cancer_order),
        1,
        facecolor=color,
        edgecolor='none',
        alpha=1,
        zorder=-1
    ))

for _, row in results_df_plot.iterrows():
    y = transcript_to_y[row["transcript"]]
    x = cancer_to_x[row["cancer_code"]]
    plt.scatter(
        x=x,
        y=y,
        s=row["significance"]*90,
        c=row["log2FC"],
        cmap="bwr",
        vmin=-results_df_plot["log2FC"].abs().max(),
        vmax=results_df_plot["log2FC"].abs().max(),
        edgecolor="black" if row["extreme"] else "none",
        linewidths=1.2 if row["extreme"] else 0
    )

plt.xticks(ticks=np.arange(len(cancer_order)), labels=cancer_order, rotation=90, ha="right", fontsize=9)
plt.yticks(y_pos, transcripts_order, fontsize=9)

cbar = plt.colorbar(label="log2FC", orientation='vertical', shrink=0.4)

p_thresholds = [1e-4, 0.05, 0.2, 0.5]
sizes_legend = [(1 - p) * 70 for p in p_thresholds]

ax_leg = plt.gcf().add_axes([0.70, 0.07, 0.08, 0.24])
ax_leg.set_xlim(-1.2, 1.2)
ax_leg.set_ylim(-0.5, len(p_thresholds) - 0.5)
ax_leg.axis('off')

for i, (s, p) in enumerate(zip(sizes_legend, p_thresholds)):
    ax_leg.scatter(0, i, s=s, color='black', edgecolors='black', linewidths=0.6, zorder=3)
    label = f"{p:.0e}" if p < 0.01 else f"{p:g}"
    ax_leg.text(0.6, i, label, ha='left', va='center', fontsize=7)

ax_leg.text(0.5, len(p_thresholds), "P-values", ha="center", fontsize=9, fontweight='bold')

ax.set_facecolor("white")
ax.patch.set_edgecolor("black")
ax.patch.set_linewidth(1)
plt.grid(False)
plt.tight_layout()
plt.subplots_adjust(right=0.85)

plt.savefig(os.path.join(figures_dir, "transcripts_dotplot.png"), dpi=300)
plt.show()

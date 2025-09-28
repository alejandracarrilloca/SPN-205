#!/bin/bash

#Author:
#Dulce Alejandra Carrillo Carlos
#Project: SPN-205 — Nonsense-mediated decay isoform expression analysis

# Extract normal sample IDs from the combined
# TCGA-GTEx phenotype file (`data/TcgaGTEX_phenotype.txt`) for
# each cancer type.


INPUT=data/TcgaGTEX_phenotype.txt
OUTDIR=results/IDS/GTEx_IDS
mkdir -p $OUTDIR

RESULT_FILE=results/IDS/Compare_GEPIA/GTEx_counts.txt
> $RESULT_FILE  

declare -A TISSUE_MAP=(
  [ACC]="Adrenal Gland"
  [BLCA]="Bladder"
  [BRCA]="Breast"
  [CESC]="Cervix Uteri"
  [CHOL]=""
  [COAD]="Colon"
  [DLBC]="Whole Blood"
  [ESCA]="Esophagus - Mucosa"
  [GBM]="Brain - Cortex"
  [HNSC]=""
  [KICH]="Kidney"
  [KIRC]="Kidney"
  [KIRP]="Kidney"
  [LAML]="Whole Blood"
  [LGG]="Brain - Cortex"
  [LIHC]="Liver"
  [LUAD]="Lung"
  [LUSC]="Lung"
  [MESO]=""
  [OV]="Ovary"
  [PAAD]="Pancreas"
  [PCPG]=""
  [PRAD]="Prostate"
  [READ]="Colon"
  [SARC]=""
  [SKCM]="Skin"
  [STAD]="Stomach"
  [TGCT]="Testis"
  [THCA]="Thyroid"
  [THYM]="Whole Blood"
  [UCEC]="Uterus"
  [UCS]="Uterus"
  [UVM]=""
)

declare -A CANCER_NAME=(
  [ACC]="Adrenocortical Carcinoma"
  [BLCA]="Bladder Urothelial Carcinoma"
  [BRCA]="Breast Invasive Carcinoma"
  [CESC]="Cervical & Endocervical Cancer"
  [CHOL]="Cholangiocarcinoma"
  [COAD]="Colon Adenocarcinoma"
  [DLBC]="Lymphoid Neoplasm Diffuse Large B-cell Lymphoma"
  [ESCA]="Esophageal Carcinoma"
  [GBM]="Glioblastoma Multiforme"
  [HNSC]="Head & Neck Squamous Cell Carcinoma"
  [KICH]="Kidney Chromophobe"
  [KIRC]="Kidney Clear Cell Carcinoma"
  [KIRP]="Kidney Papillary Cell Carcinoma"
  [LAML]="Acute Myeloid Leukemia"
  [LGG]="Brain Lower Grade Glioma"
  [LIHC]="Liver Hepatocellular Carcinoma"
  [LUAD]="Lung Adenocarcinoma"
  [LUSC]="Lung Squamous Cell Carcinoma"
  [MESO]="Mesothelioma"
  [OV]="Ovarian Serous Cystadenocarcinoma"
  [PAAD]="Pancreatic Adenocarcinoma"
  [PCPG]="Pheochromocytoma & Paraganglioma"
  [PRAD]="Prostate Adenocarcinoma"
  [READ]="Rectum Adenocarcinoma"
  [SARC]="Sarcoma"
  [SKCM]="Skin Cutaneous Melanoma"
  [STAD]="Stomach Adenocarcinoma"
  [TGCT]="Testicular Germ Cell Tumors"
  [THCA]="Thyroid Carcinoma"
  [THYM]="Thymoma"
  [UCEC]="Uterine Corpus Endometrioid Carcinoma"
  [UCS]="Uterine Carcinosarcoma"
  [UVM]="Uveal Melanoma"
)

for CANCER in "${!CANCER_NAME[@]}"; do
    TISSUE=${TISSUE_MAP[$CANCER]}
    NAME=${CANCER_NAME[$CANCER]}
    {
      head -n 1 $INPUT
      if [[ -n "$TISSUE" ]]; then
        grep -Ei "^GTEX.*${TISSUE}" $INPUT | grep -v "Cell Line"
      fi
      grep -Ei "^TCGA.*${NAME}.*Normal" $INPUT
    } > $OUTDIR/${CANCER}_normal_ids.txt
done

for FILE in "$OUTDIR"/*; do
    LINE_COUNT=$(wc -l < "$FILE")
    CANCER=$(basename "$FILE" | cut -d'_' -f1)
    echo "$CANCER : $LINE_COUNT" >> $RESULT_FILE
done

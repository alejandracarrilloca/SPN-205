#!/bin/bash

#Author:
#Dulce Alejandra Carrillo Carlos
#Project: SPN-205 — Nonsense-mediated decay isoform expression analysis

# Description:
# Extract TCGA sample IDs for each cancer type
# from the combined phenotype file `data/TcgaGTEX_phenotype.txt`.

(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Adrenocortical Cancer") > results/IDS/TCGA_IDS/ACC_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Bladder Urothelial Carcinoma") > results/IDS/TCGA_IDS/BLCA_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Breast Invasive Carcinoma") > results/IDS/TCGA_IDS/BRCA_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Cervical & Endocervical Cancer") > results/IDS/TCGA_IDS/CESC_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Cholangiocarcinoma") > results/IDS/TCGA_IDS/CHOL_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Colon Adenocarcinoma") > results/IDS/TCGA_IDS/COAD_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Diffuse Large B-Cell Lymphoma") > results/IDS/TCGA_IDS/DLBC_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Esophageal Carcinoma") > results/IDS/TCGA_IDS/ESCA_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Glioblastoma Multiforme") > results/IDS/TCGA_IDS/GBM_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Head & Neck Squamous Cell Carcinoma") > results/IDS/TCGA_IDS/HNSC_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Kidney Chromophobe") > results/IDS/TCGA_IDS/KICH_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Kidney Clear Cell Carcinoma") > results/IDS/TCGA_IDS/KIRC_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Kidney Papillary Cell Carcinoma") > results/IDS/TCGA_IDS/KIRP_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Acute Myeloid Leukemia") > results/IDS/TCGA_IDS/LAML_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Brain Lower Grade Glioma") > results/IDS/TCGA_IDS/LGG_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Liver Hepatocellular Carcinoma") > results/IDS/TCGA_IDS/LIHC_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Lung Adenocarcinoma") > results/IDS/TCGA_IDS/LUAD_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Lung Squamous Cell Carcinoma") > results/IDS/TCGA_IDS/LUSC_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Mesothelioma") > results/IDS/TCGA_IDS/MESO_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Ovarian Serous Cystadenocarcinoma") > results/IDS/TCGA_IDS/OV_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Pancreatic Adenocarcinoma") > results/IDS/TCGA_IDS/PAAD_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Pheochromocytoma & Paraganglioma") > results/IDS/TCGA_IDS/PCPG_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Prostate Adenocarcinoma") > results/IDS/TCGA_IDS/PRAD_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Rectum Adenocarcinoma") > results/IDS/TCGA_IDS/READ_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Sarcoma") > results/IDS/TCGA_IDS/SARC_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Skin Cutaneous Melanoma") > results/IDS/TCGA_IDS/SKCM_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Stomach Adenocarcinoma") > results/IDS/TCGA_IDS/STAD_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Testicular Germ Cell Tumor") > results/IDS/TCGA_IDS/TGCT_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Thyroid Carcinoma") > results/IDS/TCGA_IDS/THCA_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Thymoma") > results/IDS/TCGA_IDS/THYM_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Uterine Corpus Endometrioid Carcinoma") > results/IDS/TCGA_IDS/UCEC_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Uterine Carcinosarcoma") > results/IDS/TCGA_IDS/UCS_ids.txt
(head -n 1 data/TcgaGTEX_phenotype.txt; grep "^TCGA" data/TcgaGTEX_phenotype.txt | grep "Uveal Melanoma") > results/IDS/TCGA_IDS/UVM_ids.txt
#!/bin/bash

#Author:
#Dulce Alejandra Carrillo Carlos
#Project: SPN-205 — Nonsense-mediated decay isoform expression analysis

# Description:
# Extracts isoform-level expression data
# for SPN-201, SPN-202, SPN-203,
# SPN-204, SPN-205 from the GTEx/TCGA isoform dataset
# (TcgaGtex_rsem_isoform_tpm).


INPUT="data/TcgaGtex_rsem_isoform_tpm"

declare -A transcripts=(
    [SPN-201]="ENST00000360121.4"
    [SPN-202]="ENST00000395389.2"
    [SPN-203]="ENST00000436527.5"
    [SPN-204]="ENST00000561857.1"
    [SPN-205]="ENST00000563039.2"
)

header=$(head -n 1 "$INPUT")

for iso in "${!transcripts[@]}"; do
    enst="${transcripts[$iso]}"
    outdir="results/Transcripts/$iso"
    outfile="$outdir/${iso}_expression.tsv"

    mkdir -p "$outdir"

    {
        echo "$header"
        grep "^$enst" "$INPUT"
    } > "$outfile"

done

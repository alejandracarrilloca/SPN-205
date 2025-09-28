#Comparación de secuencias de transcritos 

>ENST00000360121.4 SPN-201 cdna:protein_coding - 
>ENST00000395389.2 SPN-202 cdna:protein_coding - 
>ENST00000436527.5 SPN-203 cdna:protein_coding - 
>ENST00000561857.1 SPN-204 cdna:retained_intron - Non protein coding 
>ENST00000563039.2 SPN-205 cdna:nonsense_mediated_decay - 
>ENST00000652691.1 SPN-206 cdna:protein_coding - 

Transcripts in my Gft comprehensive data 
ENST00000395389.2
ENST00000436527.5
ENST00000360121.4
ENST00000652691.1
ENST00000563039.2
ENST00000561857.1

1. Get files: https://www.gencodegenes.org/human/

Genome sequence, primary assembly (GRCh38)	PRI	
Nucleotide sequence of the GRCh38 primary genome assembly (chromosomes and scaffolds)
The sequence region names are the same as in the GTF/GFF3 files
Fasta - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz
 
Basic gene annotation	CHR	
It contains the basic gene annotation on the reference chromosomes only
This is a subset of the corresponding comprehensive annotation, including only those transcripts tagged as 'basic' in every gene
This is the main annotation file for most users
GTF
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.basic.annotation.gtf.gz

Extract metadata from gtf for SPN
# extract SPN lines from GTF
grep -w 'gene_name "SPN"' data/gencode.v48.annotation.gtf   > data/spn.gtf

# list unique transcript IDs + biotype + tags
awk '$3=="transcript" {print $0}' data/spn.gtf | sed 's/;/\n/g' | grep -E 'transcript_id|transcript_biotype|tag|gene_name' | paste - - - | sed 's/ gene_name "SPN"//'

#Extract sequences 

# cDNA sequences (full transcript including UTRs)
gffread data/spn.gtf -g data/GRCh38.primary_assembly.genome.fa -w results/sequences/spn.cdna.fa

# CDS sequences (coding regions only)
gffread data/spn.gtf -g data/GRCh38.primary_assembly.genome.fa -x results/sequences/spn.cds.fa

If the cluster has transeq

transeq -sequence results/sequences/spn.cds.fa -outseq results/sequences/spn.prot.fa

If not -> transeq.py

# Should each show 6 
grep -c '^>' results/sequences/spn.cdna.fa
grep -c '^>' results/sequences/spn.cds.fa
grep -c '^>' results/sequences/spn.prot.fa

#OUtput 
6
5
5
 -> One sequence is not portein coding 
 
Which one? 

grep '^>' results/sequences/spn.cds.fa | sed 's/>//' > results/sequences/cds_ids.txt
comm -23 <(sort results/sequences/transcript_ids.txt) <(sort cds_ids.txt)

ENST00000561857.1

# join counts into summary -> Isoform_comparison/scripts/spn_summary.py

       transcript_id  cdna_len  cds_len  prot_len  exon_count
0  ENST00000395389.2      1935     1203       401           2
1  ENST00000436527.5      1197     1065       355           2
2  ENST00000360121.4      6894     1203       401           2
3  ENST00000652691.1      7252     1203       401           2
4  ENST00000563039.2      1853     1203       401           3
5  ENST00000561857.1      2077        0         0           1


Protein aligment - Determinar que tan iguales son sus proteinas 

mafft --auto results/sequences/spn.prot.fa > results/alignments/spn.prot.mafft.fa

# produce a percent-identity pairwise matrix 
-> compute_pairwise_protein.py


# Diferencias en secuencia de proteina usando protein_diferences.py 

out_file = "results/alignments/protein_comparison.tsv"

# cDNA alignment - Determinar que tan diferentes son las secuencias de cDNA usando 

mafft --auto results/sequences/spn.prot.fa > results/alignments/spn.prot.mafft.fa


-> compute_pairwise_trans.py 


# Ahora armamos los mapas de cada transcrito 
-> Maps isoforms 

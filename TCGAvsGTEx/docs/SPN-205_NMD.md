## Antes de: No se encontraron datos usando BioLinks en R de TCGA, en ninguno de los proyectos, no directamente

1. Descargar datos de USC Xena: 
https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_rsem_isoform_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=http%3A%2F%2F127.0.0.1%3A7222

Nota: son datos de este proyecto https://www.nature.com/articles/nbt.3772

SPN-205: ENST00000563039.2
SPN-201: ENST00000360121.4
SPN-202: ENST00000395389.2
SPN-203: ENST00000436527.5
SPN-204
SPN-206


2. Verificar que el transcrito se encuentre en la lista de ids de el archivo antes de descargar: grep -q "ENSG00000197471" TcgaTargetGtex_rsem_isopct && echo "Sí se encontró" || echo "No se encontró"

3. Descargar 
4. Filtrar para SPN-205: 
(cat data/TcgaGtex_rsem_isoform_tpm | head -n 1 && grep "^ENST00000360121.4" data/TcgaGtex_rsem_isoform_tpm) > /home/carrillo/Documentos/SPN-205_NMD/results/SPN-201/SPN-201_expression.tsv

5. Encontrar proyectos que me son relevantes: 
acute myeloid leukemia (LAML), thymoma (THYM), Testicular germ cell tumor (TGCT), diffuse large B-cell lymphoma (DLBC), pancreatic adenocarcinoma (PAAD), and stomach adenocarcinoma (STAD). 

6. Para encotrarlos hay que descargar el archivo de fenpotipos para saber que IDs pertenecen a que tipo de cancer: https://xenabrowser.net/datapages/?dataset=TcgaTargetGTEX_phenotype.txt&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=http%3A%2F%2F127.0.0.1%3A7222

7. Debo filtrar los canceres que me interesan -> en el filtro se debe especificar que pertenezcan a TCGA y que no sean de tejidos normales 

(head -n 1 "$INPUT_FILE" && \
 grep -Ei '^TCGA.*Stomach Adenocarcinoma' "$INPUT_FILE" | grep -vi 'Normal') \
> "$OUTPUT_DIR/STAD_ids.txt"

8. Ahora tengo un archivo con los ids que me interesan, vamos a separarlos por tipo de cancer con sus respectivas expresiones 

Comandos individuales por tipo de cancer: 

(head -n 1 docs/TcgaTargetGTEX_phenotype.txt && grep -Ei 'Acute Myeloid Leukemia' docs/TcgaTargetGTEX_phenotype.txt) > results/LAML_ids.txt

(head -n 1 docs/TcgaTargetGTEX_phenotype.txt && grep -Ei 'Thymoma' docs/TcgaTargetGTEX_phenotype.txt) > results/THYM_ids.txt

(head -n 1 docs/TcgaTargetGTEX_phenotype.txt & grep -Ei 'Testicular Germ Cell Tumor' docs/TcgaTargetGTEX_phenotype.txt) > results/TGCT_ids.txt

(head -n 1 docs/TcgaTargetGTEX_phenotype.txt && grep -Ei 'Diffuse Large B-Cell Lymphoma' docs/TcgaTargetGTEX_phenotype.txt) > results/DLBC_ids.txt

(head -n 1 docs/TcgaTargetGTEX_phenotype.txt && grep -Ei 'Pancreatic Adenocarcinoma' docs/TcgaTargetGTEX_phenotype.txt) > results/PAAD_ids.txt

(head -n 1 docs/TcgaTargetGTEX_phenotype.txt && grep -Ei 'Stomach Adenocarcinoma' docs/TcgaTargetGTEX_phenotype.txt) > results/STAD_ids.txt

10. Ahora vamos a generar archivos a partir de (SPN 205 expression data) que sean exclusivos de cada cancer, filtrando con los ids que ya habiamos encontrado -> Filter_exp_data_by_TCGA-id.py 

9. Ahora necesito los tejidos normales: Se verifica en GEPIA para ver los tejidos nromales que se usan para la comparacion contra tejidos tumorales, tener cuidado de añadir los tejidos normales incluidos en TCGA 

10. Ahora si filtramos en el archivo de SPN-205 tejidos normales y hacemos un archivo por cancer

11. Conseguir los datos de SPN normal para la comparacion, mismo archivo de isoformas de UCSC Xena 

Se va a llevar a cabo el análisis usando dos de los transcritos principales, la expresion del tercer transcrito identificado como principal no fue encontrada en las bases de datos, no se encontro en los proyectos de TCGA ni GTEx, tampoco se encuentra en GEPIA2. 

Se van a usar los otros dos transcritos que muestran  una expresion elevada (SPN-202 y SPN-201), la cual se pudo observar en GEPIA2, siendo el transcrito SPN-201 el de la expresion mas significativa.  

12. Repetir el analisis anterior ahora para SPN-202 y SPN-201

Se invesitgó que transcriitos de SPN eran los principales en APPRIS: https://appris.bioinfo.cnio.es/#/database/id/homo_sapiens/ENSG00000197471?as=hg38&sc=ensembl&ds=e114v50

SPN-201 
SPN-202
SPN-203
SPN-206 -> No hay datos de RNA-seq 


Antes de llevar a cabo el analisis vamos a definir la expresion canonica de SPN, sumando la expresion de SPN-201 y SPN-202 y SPN-203, tenemos que verificar que estamos sumando muestras iguales 

Ahora tenenmos archivos con la expresion canonica de SPN 

13. Ya cuando se tienen los datos de los otros tres transcritos, se debe llevar a cabo el analisis de expresion diferencial

-> Ahora si vamos a comparar SPN con SPN-205 

Contamos con TPM, por lo tanto no se puede usar Deseq2 ya que necesitaríamos counts, por lo que vamos a usar un analisis no parametrico 

Se utilizó la prueba de rango con signo de Wilcoxon (específicamente la versión para muestras independientes, conocida como prueba de Mann–Whitney U) debido a que los datos de expresión génica, particularmente en estudios transcriptómicos como este, no suelen seguir una distribución normal. Este tipo de datos puede estar sesgado, contener valores extremos (outliers) y presentar varianza heterogénea entre grupos, lo que invalida los supuestos de normalidad requeridos por pruebas paramétricas como el *t*-test. La prueba de Wilcoxon es una alternativa no paramétrica robusta que permite comparar las distribuciones de expresión entre muestras tumorales y normales sin asumir una forma específica de la distribución, haciendo de ella una opción adecuada y confiable para detectar diferencias significativas en este contexto.

Se calculo el Fold change para ver la expresion de cada transcrito individual e idenitificar los canceres con sobre expresiones positivas o negativas. A partir de eso se tomo los canceres con FC significativo en SPN-205. Se realizaron box plots para poder vizualizar la expresion de los transcritos marcados como principales contra SPN-205


Se hizo un bubble plot con el
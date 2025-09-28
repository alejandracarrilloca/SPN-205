Exploring the role of SPN-205: an nonsense-mediated decay targeted CD43 isoform
Canva: https://www.canva.com/design/DAGvO9HsPp8/B84BcWVRRvncNn-DMlo9Qw/edit?utm_content=DAGvO9HsPp8&utm_campaign=designshare&utm_medium=link2&utm_source=sharebutton
No se encontraron datos directamente relacionados con SPN-205 en los proyectos de TCGA. Esto llevó a explorar fuentes alternativas como los datasets procesados por TOIL en UCSC Xena, los cuales integran datos de TCGA, TARGET y GTEx.

Recursos utilizados:
(https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_rsem_isoform_tpm&host=https%3A%2F%2Ftoil.xenahubs.net)

Paso 1: Descarga de datos

Se descargaron los datos del siguiente enlace:
(https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_rsem_isoform_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=http%3A%2F%2F127.0.0.1%3A7222)

Nota: estos datos provienen del proyecto publicado en *Nature Biotechnology* (2017): (https://www.nature.com/articles/nbt.3772)

Transcrito de interés: SPN-205 ID: `ENST00000563039.2`

Paso 3: Descarga y filtrado

Se filtraron los datos para SPN-205 con este comando:

```bash
(cat data/TcgaGtex_rsem_isoform_tpm | head -n 1 && grep "^ENST00000360121.4" data/TcgaGtex_rsem_isoform_tpm) > /home/carrillo/Documentos/SPN-205_NMD/results/SPN-205/SPN-205_expression.tsv
```

Paso 4: Selección de cánceres relevantes

Se seleccionaron los siguientes proyectos para análisis, basados en observaciones preliminares:

* Acute Myeloid Leukemia (LAML)
* Thymoma (THYM)
* Testicular Germ Cell Tumor (TGCT)
* Diffuse Large B-Cell Lymphoma (DLBC)
* Pancreatic Adenocarcinoma (PAAD)
* Stomach Adenocarcinoma (STAD)

Como nuestro archivo contiene experimentos con IDs no descriptivos debemos filtrar los experimentos que son de nuestro interes, por lo tanto se debe obtener la anotacion de cada experimento, descargando el archivo de fenotipos: 

Paso 5: Obtener IDs de muestra por tipo de cáncer

Se descargó el archivo de fenotipos desde:
[https://xenabrowser.net/datapages/?dataset=TcgaTargetGTEX\_phenotype.txt\&host=https%3A%2F%2Ftoil.xenahubs.net](https://xenabrowser.net/datapages/?dataset=TcgaTargetGTEX_phenotype.txt&host=https%3A%2F%2Ftoil.xenahubs.net)

Después, se filtraron las muestras por tipo de cáncer con los siguientes comandos:

```bash
(head -n 1 docs/TcgaTargetGTEX_phenotype.txt && grep -Ei 'Acute Myeloid Leukemia' docs/TcgaTargetGTEX_phenotype.txt) > results/LAML_ids.txt
(head -n 1 docs/TcgaTargetGTEX_phenotype.txt && grep -Ei 'Thymoma' docs/TcgaTargetGTEX_phenotype.txt) > results/THYM_ids.txt
(head -n 1 docs/TcgaTargetGTEX_phenotype.txt & grep -Ei 'Testicular Germ Cell Tumor' docs/TcgaTargetGTEX_phenotype.txt) > results/TGCT_ids.txt
(head -n 1 docs/TcgaTargetGTEX_phenotype.txt && grep -Ei 'Diffuse Large B-Cell Lymphoma' docs/TcgaTargetGTEX_phenotype.txt) > results/DLBC_ids.txt
(head -n 1 docs/TcgaTargetGTEX_phenotype.txt && grep -Ei 'Pancreatic Adenocarcinoma' docs/TcgaTargetGTEX_phenotype.txt) > results/PAAD_ids.txt
(head -n 1 docs/TcgaTargetGTEX_phenotype.txt && grep -Ei 'Stomach Adenocarcinoma' docs/TcgaTargetGTEX_phenotype.txt) > results/STAD_ids.txt
```

Paso 6: Separar datos por tipo de cáncer

Se usarán los archivos de expresión filtrados y los IDs para generar archivos por transcrito utilizando el script TCGA_filter_exp_data_IDS, para los datos tumorales 

## Paso 7: Identificar tejidos normales equivalentes

Se identificaron muestras normales por tejido para comparación tumoral:

Gepia               	Yo 

DLBC        	Blood               	Blood
LAML        	Bone Marrow         	Blood
PAAD        	Pancreas            	Pancreas
STAD        	Stomach             	Stomach
TCGT        	Testis              	Testis
THYM        	Blood               	Blood


## Paso 8: Filtrar datos para muestras normales

Por ejemplo, para muestras normales de sangre se usó:

```bash
(head -n 1 data/TcgaGTEX_phenotype.txt && grep -Ei 'Whole Blood' data/TcgaGTEX_phenotype.txt) > /home/carrillo/Documentos/SPN-205_NMD/results/SPN-205/GTEx_IDS/LAML_GTEx_ids.txt
```
Comparación de las muestras de este análisis contra las disponibles en GEPIA 2 

###         Gepia                                   Yo
                                  
            T       N(TCGA)     N(GTEx)             T       N(TCGA)     N(GTEx)
DLBX        47      -           337                 47      -           337
LAML        173     -           70                  173     -           337
PAAD        179     4           167                 179     4           167
STAD        408     36          175                 414     36          175
TCGT        137     -           165                 154     -           165
THYM        118     2           337                 119     2           337


Paso 9: Comparación con transcrito SPN canónico

Para definir la expresión del gen canónico se debe identificar los transcritos principales y sumar su expresión 

Paso 10: Validación de transcritos principales

Se investigaron los transcritos principales de SPN usando APPRIS:
[https://appris.bioinfo.cnio.es/#/database/id/homo\_sapiens/ENSG00000197471?as=hg38\&sc=ensembl\&ds=e114v50](https://appris.bioinfo.cnio.es/#/database/id/homo_sapiens/ENSG00000197471?as=hg38&sc=ensembl&ds=e114v50)

Se encontraron tres transcritos marcados como principales. Uno de ellos no se encontró en los datasets de TCGA ni GTEx, ni tampoco en GEPIA2, por lo cual fue descartado. Los otros dos transcritos presentaron niveles de expresión elevados y fueron visibles tanto en GEPIA2 como en los archivos descargados. De estos, el transcrito SPN-201 fue el que mostró la expresión más significativa.

SPN-201 
SPN-202
SPN-203
SPN-206 -> No hay datos de RNA-seq 

Se lleva acabo el mismo analisis que se hizo para SPN-205

Antes de llevar a cabo el analisis vamos a definir la expresion canonica de SPN, sumando la expresion de SPN-201 y SPN-202 y SPN-203, es importante verificar que se estan sumando muestras iguales 

Ahora tenemos archivos con la expresión canónica de SPN 

13. Ya cuando se tienen los datos de los otros tres transcritos, se debe llevar a cabo el análisis de expresión diferencial

-> Ahora si vamos a comparar SPN con SPN-205 

Contamos con TPM, por lo tanto no se puede usar Deseq2 ya que necesitaríamos counts, por lo que vamos a usar un análisis no paramétrico para comparar muestras independientes 

Se utilizó la prueba de rango con signo de Wilcoxon (específicamente la versión para muestras independientes, conocida como prueba de Mann–Whitney U) debido a que los datos de expresión génica, particularmente en estudios transcriptómicos como este, no suelen seguir una distribución normal. Este tipo de datos puede estar sesgado, contener valores extremos (outliers) y presentar varianza heterogénea entre grupos, lo que invalida los supuestos de normalidad requeridos por pruebas paramétricas como el *t*-test. La prueba de Wilcoxon es una alternativa no paramétrica robusta que permite comparar las distribuciones de expresión entre muestras tumorales y normales sin asumir una forma específica de la distribución, haciendo de ella una opción adecuada y confiable para detectar diferencias significativas en este contexto.

A partir del archivo generado por el wilcoxon test, se va a generar boxplots de la expresión de SPN y SPN-205 y después una comparación de la expresión de ambos transcritos por cada cáncer. 



















DEGs detections between Control and Vent
----------------------------------------
working directory: ~/Documents/2021/White_island/reads_number_matrix 
# Detect the Outliers in each species
***
## Blenny
Get the infor and raw read nb of all of the Blenny inds   
```bash
less coldata_rename_trait.txt|perl -alne 'print if /^\s+/ || /^Blenny/' >coldata_blenny.txt
extract_reads_nb --matrix all_species_raw_rename_matrix.xls --samples coldata_blenny.txt > Blenny_raw_rename_matrix.xls
mv coldata_blenny.txt Blenny_raw_rename_matrix.xls Blenny/
```
### 1. pca based on all genes (cd Blenny)  
Result: Blenny_all_gene.pdf
```bash
mpca --matrix Blenny_raw_rename_matrix.xls \
--samples coldata_blenny.txt \
--column Site_2 \
--title Blenny \
--label \
--prefix Blenny_all_gene
```
***
### 2. pca based on top 1000 variance genes (DESeq2)   
Result: Blenny_top1000.pdf   
```bash
mpca_rna --matrix Blenny_raw_rename_matrix.xls \
--samples coldata_blenny.txt \
--column Site_2 \
--title Blenny \
--prefix Blenny_top1000
```
***
### 3. WCGNA (cd Blenny) 
```bash
extract_reads_nb --matrix ../White_island.TPM.TMM.sqrt.rename.matrix --samples coldata_blenny.txt >Blenny_normalized_matrix.xls
```
```R
setwd("~/Documents/2021/White_island/reads_number_matrix/Blenny")
library(WGCNA)
options(stringsAsFactors = FALSE)
# read expression data
####################################################
white_island = read.table("Blenny_raw_rename_matrix.xls", header = T, row.names = 1)
# transposed the expression data for the hclust
datExpr0=t(white_island)
# remove the bad genes and bad smaples
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr0), method = "average")
# plot the dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
```
Result: Blenny_wgcna_outlier.pdf (outlier: Blenny_Vn_3).
***
### 4. Hclust 
```bash
mhclust --matrix Blenny_raw_rename_matrix.xls --samples coldata_blenny.txt --column Site_2 --title Blenny --prefix Blenny_hclust
```
Result: Blenny_hclust.pdf (outlier: Blenny_Vn_3).
***
# DEGs detection
## Blenny
remove the outlier from the information and matrix of Blenny  
coldata_blenny_remove_outlier.txt (remove Blenny_Vn_3)  
```bash
extract_reads_nb --matrix Blenny_raw_rename_matrix.xls --samples coldata_blenny_remove_outlier.txt >Blenny_raw_rename_matrix_remove_outlier.xls
DESeq --matrix Blenny_raw_rename_matrix_remove_outlier.xls --samples coldata_blenny_remove_outlier.txt --column Site_2 --prefix Blenny
extract_anno --genes Blenny_Control_Vent.DEGs.txt --anno ../unprot_name_description_orthgroup.txt --col 1 >Blenny_Control_Vent.DEGs.ano.txt
```
Results files: Blenny_Control_Vent.DEGs.txt; Blenny_Control_Vent.DEGs.ano.txt; Blenny_Control_Vent.csv (19 DEGs)  
***
## Blue_eyed
Result: Blue_eyed_wgcna_outlier.pdf; Blue_eyed_hclust.pdf.   
Outliers: Blue_eyed_Cs_3; Blue_eyed_Vn_2.
### DEGs detection
remove Blue_eyed_Cs_3 and Blue_eyed_Vn_2 (information table: coldata_blue_eyed_remove_outlier.txt)  




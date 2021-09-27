Detect the gene expression correlation with trait data
------------------------------------------------------
matrix: ~/Documents/2021/White_island/reads_number_matrix/White_island.TPM.TMM.sqrt.rename.filtered.matrix    
## Based on all species using normalized reads nb matrix
```bash
# ~/Documents/2021/White_island/reads_number_matrix
cp coldata_rename_trait_recode_remove_outlier.txt ~/Documents/2021/White_island/WGCNA/coldata_rename_trait_recode_remove_outlier_wgcna.txt
cp White_island.TPM.TMM.sqrt.rename.filtered.matrix ~/Documents/2021/White_island/WGCNA/
```
### trait data: coldata_rename_trait_recode_remove_outlier_wgcna.txt  
Remove some traits in coldata_rename_trait_recode_remove_outlier_wgcna.txt  
keep: Site (control, vent); Species (Blenny, Common, Yaldwyn, Blue_eyed); pH (pH_Mett);	Salinity;	Length.   
remove the ind if it has "NA" in any one of these traits (Blenny: 17, Common: 9, Yaldwyn: 20, Blue_eyed: 9).   
### Get the normalized reads matrix
```bash
extract_reads_nb --matrix White_island.TPM.TMM.sqrt.rename.filtered.matrix \
--samples coldata_rename_trait_recode_remove_outlier_wgcna.txt \
>White_island.TPM.TMM.sqrt.rename.filtered.wgcna.matrix
```
### transport the matrix of expression and trait to working station
working dir: /media/HDD/white_island/wgcna  
```bash
scp coldata_rename_trait_recode_remove_outlier_wgcna.txt White_island.TPM.TMM.sqrt.rename.filtered.wgcna.matrix Kang@147.8.76.231:/media/HDD/white_island/wgcna
```

***
## Blenny  
working dir: /media/HDD/white_island/wgcna/Blenny  
The matrix of normalized reads nb and trait data: Blenny_matrix_wgcna.xls; coldata_Blenny_wgcna.txt.  
### Run wgcna_no_merge.R
**Result file**:  
|**File Name**|**File Description**|
|:---:|:---:|
|dynamicColors.csv|gene nb in each color modules|
|textMatrix_1.csv|the correlation bettwen trait and module|
### Check the sig. modules per trait
```bash
sig_module >sig_module.txt
```
### Save the correlation info. between trait and module per gene
out the code and run the code in R console
```bash
code_trait_gene_cor --traits pH Salinity Length
```
**Result**: geneInfo_Length.csv; geneInfo_pH.csv; geneInfo_Salinity.csv.  

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
Get the matrix of normalized reads nb and trait data  

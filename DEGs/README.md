DEGs detections between Control and Vent
----------------------------------------
working directory: ~/Documents/2021/White_island/reads_number_matrix 
### Detect the Outliers in each species
***
#### Blenny
Get the infor and raw read nb of all of the Blenny inds   
```bash
less coldata_rename_trait.txt|perl -alne 'print if /^\s+/ || /^Blenny/' >coldata_blenny.txt
extract_reads_nb --matrix all_species_raw_rename_matrix.xls --samples coldata_blenny.txt > Blenny_raw_rename_matrix.xls
mv coldata_blenny.txt Blenny_raw_rename_matrix.xls Blenny/
```
1. pca based on all genes (cd Blenny) 
Result: Blenny_all_gene.pdf
```bash
mpca --matrix Blenny_raw_rename_matrix.xls \
--samples coldata_blenny.txt \
--column Site_2 \
--title Blenny \
--label \
--prefix Blenny_all_gene
```

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
### Check the sig. modules and genes in each sig. modules per trait
```bash
sig_module >sig_module.txt
```
### Save the correlation info. between trait and module per gene
out the code and run the code in R console
```bash
code_trait_gene_cor --traits pH Salinity Length
```
**Result**: geneInfo_Length.csv; geneInfo_pH.csv; geneInfo_Salinity.csv.  
### Based on "sig_module.txt" to extract the sig. genes in sig. correlated modules（the correlation >=0.6 or <=-0.6）
working dir: /media/HDD/white_island/wgcna/Blenny    
**Notice**: the "sig_module.txt" should be the current working dir  
```bash
sig_module_trait --traits geneInfo_pH.csv geneInfo_Length.csv geneInfo_Salinity.csv
```
output files: geneInfo_Length_sig.txt, geneInfo_pH_sig.txt, geneInfo_Salinity_sig.txt;
|**Trait**|**module nb**|**gene nb**|
|:---:|:---:|:---:|
|pH|module_nb(15)|gene_nb(144)|
|Length|module_nb(43)|gene_nb(948)|
|Salinity|module_nb(2)|gene_nb(24)|
### transform the expression data to ratio, and plot these gene expression data according to positive or negative correlated
```bash
less coldata_Blenny_wgcna.txt|perl -alne 'if (/^\s+/){print}else{($F[1]==1)?($F[1]="vent"):($F[1]="control");print join("\t",@F)}' >coldata_Blenny.txt
less geneInfo_pH_sig.txt|perl -alne 'next if /^\s+/;print $F[0] if $F[2]>0' >geneInfo_pH_sig_pos.txt
less geneInfo_pH_sig.txt|perl -alne 'next if /^\s+/;print $F[0] if $F[2]<0' >geneInfo_pH_sig_neg.txt
extract_gene_expression_plot --matrix expression_ratio.txt --sample coldata_Blenny.txt --gene geneInfo_pH_sig_pos.txt --col1 Site --order1 control vent >pH_sig_pos.plot.txt
extract_gene_expression_plot --matrix expression_ratio.txt --sample coldata_Blenny.txt --gene geneInfo_pH_sig_neg.txt --col1 Site --order1 control vent >pH_sig_neg.plot.txt
```
**accomplish by connecting the above codes**  
Results: geneInfo_pH_sig_pos_plot.txt, geneInfo_pH_sig_neg_plot.txt
```bash
prep_exp_wgcna_plot --wgcna Blenny_matrix_wgcna.xls --trait geneInfo_pH_sig.txt --sample coldata_Blenny.txt
```
***
## Common

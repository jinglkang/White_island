# Restart the analysis of white island transcriptomes
**Combine all individuals sampled from CO2 vents (south or north) as Vent, and from Control sites (south or north) as Control**   
## PCA for all species
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 二 12 20 16:35:07 ~/Documents/2022/White_island_Restart
extract_reads_nb --matrix all_species_raw_rename_matrix.xls --samples coldata_rename_trait.txt > all_species_raw_matrix.xls
rm all_species_raw_rename_matrix.xls
# PCA plot with label in points: PCA_point_label.R
# PCA plot without label in points: PCA_point_label.R
# Remove: Blue_eye_Cs_3, Yaldwyn_Cn_4, Yaldwyn_Cs_2
```
## DESeq2 for DEGs detection  
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 三 12 21 12:08:25 ~/Documents/2022/White_island_Restart/reads_number_matrix/Blue_eyed
extract_reads_nb --matrix Blue_eyed_raw_rename_matrix.xls --samples coldata_Blue_eyed_remove_outlier.txt > Blue_eyed_remove_outlier_matrix.xls
# kangjingliang@kangjingliangdeMacBook-Pro 三 12 21 12:18:36 ~/Documents/2022/White_island_Restart/reads_number_matrix/Yaldwyn
extract_reads_nb --matrix Yaldwyn_raw_rename_matrix.xls --samples coldata_Yaldwyn_remove_outlier.txt > Yaldwyn_remove_outlier_matrix.xls

# Blenny: 1 DEG
# kangjingliang@kangjingliangdeMacBook-Pro 三 12 21 14:48:23 ~/Documents/2022/White_island_Restart/reads_number_matrix/Blenny
DESeq --matrix Blenny_raw_rename_matrix.xls --samples coldata_blenny.txt --column Site_2 --prefix Blenny
extract_anno --genes Blenny_Control_Vent.DEGs.txt --anno ../unprot_name_description_orthgroup.txt --col 1 >Blenny_Control_Vent.DEGs.ano.txt

# Blue_eyed: 112 DEGs
# kangjingliang@kangjingliangdeMacBook-Pro 三 12 21 14:53:15 ~/Documents/2022/White_island_Restart/reads_number_matrix/Blue_eyed
DESeq --matrix Blue_eyed_remove_outlier_matrix.xls --samples coldata_Blue_eyed_remove_outlier.txt --column Site_2 --prefix Blue_eyed
extract_anno --genes Blue_eyed_Control_Vent.DEGs.txt --anno ../unprot_name_description_orthgroup.txt --col 1 >Blue_eyed_Control_Vent.DEGs.ano.txt

# Common: 4 DEGs
# kangjingliang@kangjingliangdeMacBook-Pro 三 12 21 15:07:08 ~/Documents/2022/White_island_Restart/reads_number_matrix/Common
DESeq --matrix Common_raw_rename_matrix.xls --samples coldata_Common.txt --column Site_2 --prefix Common
extract_anno --genes Common_Control_Vent.DEGs.txt --anno ../unprot_name_description_orthgroup.txt --col 1 >Common_Control_Vent.DEGs.ano.txt

# Yaldwyn: 50 DEGs
# kangjingliang@kangjingliangdeMacBook-Pro 三 12 21 15:09:53 ~/Documents/2022/White_island_Restart/reads_number_matrix/Yaldwyn
DESeq --matrix Yaldwyn_remove_outlier_matrix.xls --samples coldata_Yaldwyn_remove_outlier.txt --column Site_2 --prefix Yaldwyn
extract_anno --genes Yaldwyn_Control_Vent.DEGs.txt --anno ../unprot_name_description_orthgroup.txt --col 1 >Yaldwyn_Control_Vent.DEGs.ano.txt
```
## WGCNA for DEGs detection  
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 三 12 21 16:19:43 ~/Documents/2021/White_island/reads_number_matrix
cp White_island.TPM.TMM.sqrt.rename.matrix ~/Documents/2022/White_island_Restart/wgcna/wgcna_tpm_matrix.txt

# kangjingliang@kangjingliangdeMacBook-Pro 三 12 21 16:30:15 ~/Documents/2022/White_island_Restart/wgcna
cp ../reads_number_matrix/Blenny/coldata_blenny.txt ./coldata_blenny_recode.txt
# coldata_blenny_recode.txt: "Site_2" (Vent => 1; Control => 2)
cp ../reads_number_matrix/Blue_eyed/coldata_Blue_eyed_remove_outlier.txt coldata_Blue_eyed_remove_outlier_recode.txt

mv coldata_Common_recode.txt Common_inds.txt
mv coldata_Blue_eyed_remove_outlier_recode.txt Blue_eyed_inds.txt
mv coldata_Yaldwyn_remove_outlier_recode.txt Yaldwyn_inds.txt
mv coldata_blenny_recode.txt Blenny_inds.txt

extract_reads_nb --matrix wgcna_tpm_matrix.txt --samples Blenny_inds.txt > Blenny_tpm_matrix.xls
extract_reads_nb --matrix wgcna_tpm_matrix.txt --samples Blue_eyed_inds.txt > Blue_eyed_tpm_matrix.xls
extract_reads_nb --matrix wgcna_tpm_matrix.txt --samples Yaldwyn_inds.txt > Yaldwyn_tpm_matrix.xls
extract_reads_nb --matrix wgcna_tpm_matrix.txt --samples Common_inds.txt > Common_tpm_matrix.xls

# Run WGCNA
# Run_Blenny_wgcna.R; Run_Blue_eyed_wgcna.R; Run_Common_wgcna.R; Run_Yaldwyn_wgcna.R
# Identify genes that have a high significance for pH

# Blenny: 27 DEG
# kangjingliang@kangjingliangdeMacBook-Pro 一  1 02 12:44:07 ~/Documents/2022/White_island_Restart/wgcna
less Blenny_wgcna_results.csv|perl -alne 's/\"//g;@a=split /\,/;print "$a[0]\t$a[1]\t$a[2]\t$a[3]" if abs($a[2])>=0.7 && $a[3]<=0.05' > Blenny_wgcna_DEGs.txt
cut -f 1 Blenny_wgcna_DEGs.txt >Blenny_wgcna_DEGs_id.txt
extract_anno --genes Blenny_wgcna_DEGs_id.txt --anno ../reads_number_matrix/unprot_name_description_orthgroup.txt --col 1 >Blenny_wgcna_DEGs_ano.txt


# Blue_eyed: 82 DEG
less Blue_eyed_wgcna_results.csv|perl -alne 's/\"//g;@a=split /\,/;print "$a[0]\t$a[1]\t$a[2]\t$a[3]" if abs($a[2])>=0.7 && $a[3]<=0.05' > Blue_eyed_wgcna_DEGs.txt
cut -f 1 Blue_eyed_wgcna_DEGs.txt >Blue_eyed_wgcna_DEGs_id.txt
extract_anno --genes Blue_eyed_wgcna_DEGs_id.txt --anno ../reads_number_matrix/unprot_name_description_orthgroup.txt --col 1 >Blue_eyed_wgcna_DEGs_ano.txt


# Yaldwyn: 65 DEG
less Yaldwyn_wgcna_results.csv|perl -alne 's/\"//g;@a=split /\,/;print "$a[0]\t$a[1]\t$a[2]\t$a[3]" if abs($a[2])>=0.7 && $a[3]<=0.05' > Yaldwyn_wgcna_DEGs.txt
cut -f 1 Yaldwyn_wgcna_DEGs.txt >Yaldwyn_wgcna_DEGs_id.txt
extract_anno --genes Yaldwyn_wgcna_DEGs_id.txt --anno ../reads_number_matrix/unprot_name_description_orthgroup.txt --col 1 >Yaldwyn_wgcna_DEGs_ano.txt

# Common: 14 DEG
less Common_wgcna_results.csv|perl -alne 's/\"//g;@a=split /\,/;print "$a[0]\t$a[1]\t$a[2]\t$a[3]" if abs($a[2])>=0.7 && $a[3]<=0.05' > Common_wgcna_DEGs.txt
cut -f 1 Common_wgcna_DEGs.txt >Common_wgcna_DEGs_id.txt
extract_anno --genes Common_wgcna_DEGs_id.txt --anno ../reads_number_matrix/unprot_name_description_orthgroup.txt --col 1 >Common_wgcna_DEGs_ano.txt
```

## Combine DEGs from DEseq2 and WGCNA
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 一  1 02 13:14:18 ~/Documents/2022/White_island_Restart
mkdir combine; cd combine

# kangjingliang@kangjingliangdeMacBook-Pro 一  1 02 13:16:36 ~/Documents/2022/White_island_Restart/combine
cat ../wgcna/Blenny_wgcna_DEGs_ano.txt ../reads_number_matrix/Blenny/Blenny_Control_Vent.DEGs.ano.txt|sort -u > Blenny_DEGs.txt # 27 DEGs
cat ../wgcna/Blue_eyed_wgcna_DEGs_ano.txt ../reads_number_matrix/Blue_eyed/Blue_eyed_Control_Vent.DEGs.ano.txt|sort -u > Blue_eyed_DEGs.txt # 162 DEGs
cat ../wgcna/Common_wgcna_DEGs_ano.txt ../reads_number_matrix/Common/Common_Control_Vent.DEGs.ano.txt|sort -u > Common_DEGs.txt # 16 DEGs
cat ../wgcna/Yaldwyn_wgcna_DEGs_ano.txt ../reads_number_matrix/Yaldwyn/Yaldwyn_Control_Vent.DEGs.ano.txt|sort -u > Yaldwyn_DEGs.txt # 75 DEGs

# extract genes underlying specific GOs
# kangjingliang@kangjingliangdeMacBook-Pro 一  1 02 21:24:58 ~/Documents/2022/White_island_Restart/combine
extract_gene_functions -i *_enrichment.txt -a ../reads_number_matrix/unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions pH_funcs.txt --output pH_DEGs
extract_gene_functions -i *_enrichment.txt -a ../reads_number_matrix/unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions ion_funcs.txt --output ion_DEGs
extract_gene_functions -i *_enrichment.txt -a ../reads_number_matrix/unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions immune_funcs.txt --output immune_DEGs
extract_gene_functions -i *_enrichment.txt -a ../reads_number_matrix/unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions GABA_funcs.txt --output GABA_DEGs
extract_gene_functions -i *_enrichment.txt -a ../reads_number_matrix/unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions CR_funcs.txt --output CR_DEGs
```

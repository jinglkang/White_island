## DESeq2 for DEGs detection (previous)  
### South
```bash
# Cs vs. Vs
# Blenny: 108
# kangjingliang@kangjinangdeMBP 二  1 03 18:25:56 ~/Documents/2021/White_island/reads_number_matrix/Blenny
extract_anno --genes Blenny_Cs_Vs.DEGs.txt --anno ../unprot_name_description_orthgroup.txt --col 1 >Blenny_Cs_Vs_DEGs_ano.txt

# Blueeyed: 31
# kangjingliang@kangjinangdeMBP 二  1 03 20:04:35 ~/Documents/2021/White_island/reads_number_matrix/Blue_eyed
extract_anno --genes Blueeyed_Cs_Vs.DEGs.txt --anno ../unprot_name_description_orthgroup.txt --col 1 >Blueeyed_Cs_Vs_DEGs_ano.txt

# Common: 429
# kangjingliang@kangjinangdeMBP 二  1 03 20:07:23 ~/Documents/2021/White_island/reads_number_matrix/Common
extract_anno --genes Common_Cs_Vs.DEGs.txt --anno ../unprot_name_description_orthgroup.txt --col 1 >Common_Cs_Vs_DEGs_ano.txt

# Yaldwyn: 143
# kangjingliang@kangjinangdeMBP 二  1 03 20:09:37 ~/Documents/2021/White_island/reads_number_matrix/Yaldwyn
extract_anno --genes Yaldwyn_Cs_Vs.DEGs.txt --anno ../unprot_name_description_orthgroup.txt --col 1 >Yaldwyn_Cs_Vs_DEGs_ano.txt

# kangjingliang@kangjinangdeMBP 二  1 03 20:20:43 ~/Documents/2021/White_island/reads_number_matrix
mkdir South_enrichment
# kangjingliang@kangjinangdeMBP 二  1 03 20:21:26 ~/Documents/2021/White_island/reads_number_matrix/South_enrichment
extract_gene_functions -i *_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions pH_funcs.txt --output pH_DEGs
extract_gene_functions -i *_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions ion_funcs.txt --output ion_DEGs
extract_gene_functions -i *_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions immune_funcs.txt --output immune_DEGs
extract_gene_functions -i *_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions GABA_funcs.txt --output GABA_DEGs
extract_gene_functions -i *_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions CR_funcs.txt --output CR_DEGs
```
### North
```bash
# Blenny: 36
# Blueeyed: 5
# Common: 93
# Yaldwyn: 24
# kangjingliang@kangjinangdeMBP 三  1 04 16:13:04 ~/Documents/2021/White_island/reads_number_matrix/North_enrichment
extract_gene_functions -i *_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions pH_funcs.txt --output pH_DEGs
extract_gene_functions -i *_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions ion_funcs.txt --output ion_DEGs
extract_gene_functions -i *_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions immune_funcs.txt --output immune_DEGs
extract_gene_functions -i *_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions GABA_funcs.txt --output GABA_DEGs
extract_gene_functions -i *_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions CR_funcs.txt --output CR_DEGs
```

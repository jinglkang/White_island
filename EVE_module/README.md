Detect the divergent and plastic gene
-------------------------------------
# Used the species tree from orthofinder
/media/HDD/white_island/orthologue/input_pep/orthofinder_input/OrthoFinder/Results_Mar29/Species_Tree/SpeciesTree_rooted.txt  
working dir: /media/HDD/white_island/EVE_module  
```bash
cp /media/HDD/white_island/orthologue/input_pep/orthofinder_input/OrthoFinder/Results_Mar29/Species_Tree/SpeciesTree_rooted.txt ./
cp SpeciesTree_rooted.txt 1_Tree.newick
```
**SpeciesTree_rooted.txt**:   
(Blenny-1:0.185616,(Common-1:0.0978023,(Yaldwyn-1:0.0499941,Blue_eyed-1:0.0536857)0.597529:0.0950789)1:0.185616);    
**1_Tree.newick**:    
4    
(1:0.185616,(2:0.0978023,(3:0.0499941,4:0.0536857):0.0950789):0.185616);
|**Species**|**Recode Nb**|
|:---:|:---:|
|Blenny|1|
|Common|2|
|Yaldwyn|3|
|Blue_eyed|4|
***
# Prepare the expression data
all_species_raw_rename_matrix.xls: raw reads number matrix  
## all ind should not be less than 1 mapped read, mean reads nb should be more than 10
working dir: ~/Documents/2021/White_island/reads_number_matrix  
```bash
less all_species_raw_rename_matrix.xls|perl -alne 'print if /^\s+/; next if  /^\s+/;my $sum;for (my $i=1;$i < @F; $i++){if ($F[$i]==0){$sum=0;last};$sum+=$F[$i]};$mean=$sum/(@F-1);print "$_" unless $mean<=10' > all_species_raw_rename_matrix_filtered.xls
```
Result file: all_species_raw_rename_matrix_filtered.xls (15,908 genes) 
***
## keep these genes (filtered_genes.txt)
```bash
less all_species_raw_rename_matrix_filtered.xls|perl -alne 'next if /^\s+/;print $F[0]' >filtered_genes.txt
```
***
## keep inds (remove outliers: Blenny_Vn_3, Blue_eyed_Cs_3, Blue_eyed_Vn_2)
individuals should be ranked acoording to species
(1: Blenny; 2: Common; 3: Yaldwyn; 4: Blue_eyed)  
Result file: coldata_rename_trait_recode_remove_outlier.txt  
***
## Prepare normalized expression data based on filtered_genes.txt and coldata_rename_trait_recode_remove_outlier.txt  
```bash
extract_reads_nb --matrix White_island.TPM.TMM.sqrt.rename.matrix \
--genes filtered_genes.txt \
--samples coldata_rename_trait_recode_remove_outlier.txt \
>White_island.TPM.TMM.sqrt.rename.filtered.matrix
```
Result: White_island.TPM.TMM.sqrt.rename.filtered.matrix (will used in the EVE and WCGNA analysis)  
***
## Prepare 2_Nindivs.indiv: ind nb per species  
|**Species**|**Ind nb**|
|:---:|:---:|
|Blenny|17|
|Common|20|
|Yaldwyn|20|
|Blue_eyed|18|
2_Nindivs.indiv: 17 20 20 18  
### head -n 2 3_sampleExpr.dat
```bash
15908
OG0038649	4.213	2.936	2.707	2.080	3.054	3.879	3.778	3.501	3.315
```
***
# RUN eve module to detect expression shif in each species
```bash
nohup ~/software/EVE_release/EVEmodel -O -o 1 -n 15908 -t 1_Tree.newick -i 2_Nindivs.indiv -d 3_sampleExpr.dat -f _Blenny -v 50 > Blenny_eve.process.txt 2>&1 & 
nohup ~/software/EVE_release/EVEmodel -O -o 2 -n 15908 -t 1_Tree.newick -i 2_Nindivs.indiv -d 3_sampleExpr.dat -f _Common -v 50 > Common_eve.process.txt 2>&1 &
nohup ~/software/EVE_release/EVEmodel -O -o 3 -n 15908 -t 1_Tree.newick -i 2_Nindivs.indiv -d 3_sampleExpr.dat -f _Yaldwyn -v 50 > Yaldwyn_eve.process.txt 2>&1 &
nohup ~/software/EVE_release/EVEmodel -O -o 4 -n 15908 -t 1_Tree.newick -i 2_Nindivs.indiv -d 3_sampleExpr.dat -f _Blue_eyed -v 50 > Blue_eve.process.txt 2>&1 &
```
# RUN eve module on all inds
```bash
nohup ~/software/EVE_release/EVEmodel -S -n 15908 -t 1_Tree.newick -i 2_Nindivs.indiv -d 3_sampleExpr.dat -f _EVE -v 50 > EVE.process.txt 2>&1 &
```
***
# Estimate diverge and diverse gene of each species
working dir: /media/HDD/white_island/EVE_module     
```bash
less White_island.TPM.TMM.sqrt.rename.filtered.matrix|perl -alne 'print $F[0] if /^OG/' >allgenes_names.txt
```
run EVEresults.R    

Detect the divergent and plastic gene
-------------------------------------
# Used the species tree from orthofinder
/media/HDD/white_island/orthologue/input_pep/orthofinder_input/OrthoFinder/Results_Mar29/Species_Tree/SpeciesTree_rooted.txt  
working directory: /media/HDD/white_island/EVE_module  
```bash
cp /media/HDD/white_island/orthologue/input_pep/orthofinder_input/OrthoFinder/Results_Mar29/Species_Tree/SpeciesTree_rooted.txt ./
cp SpeciesTree_rooted.txt 1_Tree.newick
```
**SpeciesTree_rooted.txt**:   
(Blenny-1:0.185616,(Common-1:0.0978023,(Yaldwyn-1:0.0499941,Blue_eyed-1:0.0536857)0.597529:0.0950789)1:0.185616);    
**1_Tree.newick**:    
(1:0.185616,(2:0.0978023,(3:0.0499941,4:0.0536857):0.0950789):0.185616);
|**Species**|**Recode Nb**|
|:---:|:---:|
|Blenny|1|
|Common|2|
|Yaldwyn|3|
|Blue_eyed|4|

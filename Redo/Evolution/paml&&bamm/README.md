# White Island
## Evolutionary rate comparisons between the six reference fish species and six coral fish species in gcb
### 1. OrthoFinder for orthologous genes detection
**The protein sequences of ORFs detected by OrthoFinder as input**   
```bash
# (base) kang1234@celia-PowerEdge-T640 Mon Mar 20 09:16:09 ~/white_island/Compevo/Input_pep
nohup orthofinder -f ./ -a 20 >orthofinder-process 2>&1 &
# [1] 26527

# (base) kang1234@celia-PowerEdge-T640 Wed Mar 22 04:25:35 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21/Orthogroups
wc -l Orthogroups_SingleCopyOrthologues.txt # 47 single copy ortholgous

# Change the header: "Blue_eyed_1" --> "Blueeyed_1"; "ENSXMAG00000021355" -> "Platyfish_ENSXMAG00000021355"
# (base) kang1234@celia-PowerEdge-T640 Wed Mar 22 04:48:55 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21
perl ChangeHeader.pl # the changed header files will be in "Orth_Seq_changeHeader/"

# concatenate the single copy ortholgous for phylogeny
# copy the single copy fasta file to a new dir "single_copy" and change the name
# (base) kang1234@celia-PowerEdge-T640 Wed Mar 22 04:57:16 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21
perl Change_SCGheader.pl # the 46 single copy ortholgous will be in "single_copy/"
# (base) kang1234@celia-PowerEdge-T640 Wed Mar 22 05:03:00 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21
perl Align_trim_conca.pl > single_copy.concatenated.fasta
fasta2phy.pl single_copy.concatenated.fasta >single_copy.concatenated.phy
# (base) kang1234@celia-PowerEdge-T640 Wed Mar 22 05:04:18 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21
nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s single_copy.concatenated.phy -o Spottedgar -n single_copy.concatenated -T 24 > raxml.process 2>&1 &
# [1] 17749
```
### 2. Divid the orthogroups into small orthogroups
```bash
# Kang@fishlab3 Wed Mar 22 05:40:34 /media/HDD/white_island/Compevo/Results_Mar21
conda activate possvm
export DISPLAY=:0.0
nohup perl build_sub_orth.pl > build_sub_orth.process 2>&1 &
# [1] 31664
# Kang@fishlab3 Thu Mar 23 06:15:17 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth
less sub_orth_genecount.txt|perl -alne 'my $j;for (my $i = 1; $i <= 16; $i++){$j++ if $F[$i]>=1};print if $j==16'|wc -l # 5216 sub_orth
```
### 3. Run PAML
```bash
# Anotate the proteins of all species
# (base) kang1234@celia-PowerEdge-T640 Thu Mar 23 06:29:17 ~/white_island/Compevo/Input_pep
cat *.fasta > all.fasta
nohup diamond blastp -q all.fasta -e 1e-5 --sensitive -k 1 -d ~/swiss-prot/uniprot-filtered-reviewed_yes.fasta --out all_swissprot_diamond_ano.txt  > diamond_blastp.process 2>&1 &
# [1] 29772
# Extract gene list for PAML
# (base) kang1234@celia-PowerEdge-T640 Thu Mar 23 09:42:39 ~/white_island/Compevo/Input_pep
rm all.fasta diamond_blastp.process

# only use one seq per species
# 1. select the subgroup if 16 species have sequences in the subgroup
# 2. the seq of these selected species should be annotated by the same gene name
# 3. select the one with higher blast score if a species has mutiple seqs with the same gene name

# Kang@fishlab3 Thu Mar 23 09:52:19 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth
scp kang1234@147.8.76.177:~/white_island/Compevo/Input_pep/all_swissprot_diamond_ano.txt ./
# Kang@fishlab3 Thu Mar 23 14:32:52 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth
perl sub_orth_same_name.pl # output: sub_orth_genecount_final.txt; sub_orth_id_final.txt # 4564 sub_orth
# only keep the one with highest blastp score if a species with mutiple sequences in the sub_orth
perl Generate_paml_list.pl >sub_orth_id_paml.txt # the genes in "sub_orth_id_paml.txt" will be inputs of paml

# prepare paml input
# Kang@fishlab3 Thu Mar 23 15:28:26 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth
mkdir paml_input
# Kang@fishlab3 Thu Mar 23 15:32:43 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
cp ../sub_orth_id_paml.txt ./ # remove the header

######################################################
######################################################
# the pep fasta files
# Kang@fishlab3 Mon Mar 20 23:05:34 ~/Desktop/PapueNewGuinea-new/longest_pep/input_pep
# these pep files were come from "~/Desktop/PapueNewGuinea-new/longest_pep/input_pep"
# longest transcripts
cp /media/HDD/white_island/Compevo/*.fasta ./
less sub_orth_id_paml.txt|perl -alne 's/Blue_eyed/Blueeyed/;print' >sub_orth_id_paml.txt.1; mv sub_orth_id_paml.txt.1 sub_orth_id_paml.txt
# (base) kang1234@celia-PowerEdge-T640 Thu Mar 23 15:35:44 ~/white_island/Compevo/Input_pep
scp Medaka.fasta Platyfish.fasta Spottedgar.fasta Stickleback.fasta Zebrafish.fasta Fugu.fasta Kang@147.8.76.229:/media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
# Kang@fishlab3 Thu Mar 23 15:47:50 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
less Medaka.fasta|perl -alne 'if (/>/){s/>//;$nm="Medaka_".$_;print ">$nm"}else{print}' >Medaka.fasta.1; mv Medaka.fasta.1 Medaka.fasta
less Platyfish.fasta|perl -alne 'if (/>/){s/>//;$nm="Platyfish_".$_;print ">$nm"}else{print}' >Platyfish.fasta.1; mv Platyfish.fasta.1 Platyfish.fasta
less Spottedgar.fasta|perl -alne 'if (/>/){s/>//;$nm="Spottedgar_".$_;print ">$nm"}else{print}' >Spottedgar.fasta.1; mv Spottedgar.fasta.1 Spottedgar.fasta
less Zebrafish.fasta|perl -alne 'if (/>/){s/>//;$nm="Zebrafish_".$_;print ">$nm"}else{print}' >Zebrafish.fasta.1; mv Zebrafish.fasta.1 Zebrafish.fasta
less Fugu.fasta|perl -alne 'if (/>/){s/>//;$nm="Fugu_".$_;print ">$nm"}else{print}' >Fugu.fasta.1; mv Fugu.fasta.1 Fugu.fasta
less Stickleback.fasta|perl -alne 'if (/>/){s/>//;$nm="Stickleback_".$_;print ">$nm"}else{print}' >Stickleback.fasta.1; mv Stickleback.fasta.1 Stickleback.fasta
ll *.fasta|perl -alne '($spe)=$F[-1]=~/(.*)\.fasta/;$new=$spe.".pep.fasta";print "mv $F[-1] $new"'
mv Acura.fasta Acura.pep.fasta
mv Apoly.fasta Apoly.pep.fasta
mv Blenny.fasta Blenny.pep.fasta
mv Blueeyed.fasta Blueeyed.pep.fasta
mv Common.fasta Common.pep.fasta
mv Daru.fasta Daru.pep.fasta
mv Fugu.fasta Fugu.pep.fasta
mv Medaka.fasta Medaka.pep.fasta
mv Ocomp.fasta Ocomp.pep.fasta
mv Padel.fasta Padel.pep.fasta
mv Platyfish.fasta Platyfish.pep.fasta
mv Pmol.fasta Pmol.pep.fasta
mv Spottedgar.fasta Spottedgar.pep.fasta
mv Stickleback.fasta Stickleback.pep.fasta
mv Yaldwyn.fasta Yaldwyn.pep.fasta
mv Zebrafish.fasta Zebrafish.pep.fasta
######################################################
######################################################

######################################################
######################################################
# the nuc fasta file
# longest transcripts
# the cds sequences of six reference fish species
# Kang@fishlab3 Thu Mar 23 16:30:00 /media/HDD/cleaner_fish/genome/gene_family_3/longest_cds
cp Longest_Medaka_cds.fasta Longest_Platyfish_cds.fasta Longest_Spottedgar_cds.fasta Longest_Stickleback_cds.fasta Longest_Zebrafish_cds.fasta Longest_Fugu_cds.fasta /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
# Kang@fishlab3 Thu Mar 23 16:31:54 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
ll Longest_*.fasta|perl -alne '($spe)=$F[-1]=~/Longest_(.*)_cds\.fasta/;$new=$spe.".cds.fasta";print "mv $F[-1] $new"'
mv Longest_Fugu_cds.fasta Fugu.cds.fasta
mv Longest_Medaka_cds.fasta Medaka.cds.fasta
mv Longest_Platyfish_cds.fasta Platyfish.cds.fasta
mv Longest_Spottedgar_cds.fasta Spottedgar.cds.fasta
mv Longest_Stickleback_cds.fasta Stickleback.cds.fasta
mv Longest_Zebrafish_cds.fasta Zebrafish.cds.fasta

# gcb species
# Kang@fishlab3 Thu Mar 23 16:44:24 ~/Desktop/PapueNewGuinea-new/orthologue/orthofinder_input_nuc
cp *.fa /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
# Kang@fishlab3 Thu Mar 23 16:49:05 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
ll *_nuc.fa|perl -alne '($spe)=$F[-1]=~/(.*)_nuc\.fa/;$new=$spe.".cds.fasta";print "mv $F[-1] $new"'
mv Acura_nuc.fa Acura.cds.fasta
mv Apoly_nuc.fa Apoly.cds.fasta
mv Daru_nuc.fa Daru.cds.fasta
mv Ocomp_nuc.fa Ocomp.cds.fasta
mv Padel_nuc.fa Padel.cds.fasta
mv Pmol_nuc.fa Pmol.cds.fasta

# four fish species this study
# Kang@fishlab3 Thu Mar 23 16:54:31 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
cp /media/HDD/white_island/orthologue/input_nuc/orf_nuc/rename/Blenny.fa Blenny.cds.fasta
cp /media/HDD/white_island/orthologue/input_nuc/orf_nuc/rename/Blue_eyed.fa Blueeyed.cds.fasta
cp /media/HDD/white_island/orthologue/input_nuc/orf_nuc/rename/Common.fa Common.cds.fasta
cp /media/HDD/white_island/orthologue/input_nuc/orf_nuc/rename/Yaldwyn.fa Yaldwyn.cds.fasta
# change the header of "Blue_eyed" to "Blueeyed"
# Kang@fishlab3 Thu Mar 23 17:03:15 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
less Blueeyed.pep.fasta|perl -alne 'if (/>/){s/>//;$_=~s/\_//;print ">$_"}else{print}' >Blueeyed.pep.fasta.1; mv Blueeyed.pep.fasta.1 Blueeyed.pep.fasta
less Blueeyed.cds.fasta|perl -alne 'if (/>/){s/>//;$_=~s/\_//;print ">$_"}else{print}' >Blueeyed.cds.fasta.1; mv Blueeyed.cds.fasta.1 Blueeyed.cds.fasta

######################################################
# Ready
######################################################
cp /media/HDD/cleaner_fish/genome/gene_family_2/paml_input/prepare_input_paml.pl ./
cp /media/HDD/cleaner_fish/genome/gene_family_2/paml_input/prepare_input_paml_parallel.pl ./
less sub_orth_id_paml.txt|head -n 1|perl -alne 'for (my $i = 1; $i < @F; $i++){($spe)=$F[$i]=~/(.*)\_/;$pep=$spe.".pep.fasta";$nuc=$spe.".cds.fasta";print "$spe\t$pep\t$nuc"}' >correlation.txt
# Kang@fishlab3 Thu Mar 23 17:27:23 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
mkdir paml_files
vi prepare_input_paml_parallel.pl # change "--output ./" to "--output paml_files/"
nohup perl prepare_input_paml_parallel.pl sub_orth_id_paml.txt >prepare_input_paml.process 2>&1 &
# [1] 2794
# Kang@fishlab3 Thu Mar 23 22:20:44 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
vi spe.tre
# (((Ocomp,((Stickleback,Fugu),(((Daru,((Acura,Apoly),(Pmol,Padel))),(Blenny,((Yaldwyn,Blueeyed),Common))),(Platyfish,Medaka)))),Zebrafish),Spottedgar);

# 3826 sub_orth will be used in the paml analysis
###############
# 1. free-ratio
###############
# "my $cmd="perl codeml.pl --input temp/$temp --model free-ratio --dir paml_files --tree spe.tre --icode 0 --omega 1.2";"
# Kang@fishlab3 Sat Mar 25 07:30:51 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
vi codeml_parallel.pl
nohup perl codeml_parallel.pl final_orth_input_paml.txt >free_ratio.process 2>&1 &
# [1] 24546
# plot the free-ratio results
# Kang@fishlab3 Sun Mar 26 22:35:35 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
cat free-ratio/*.txt >free_ratio_result.txt
perl make_bin_plot_total.pl >free_ratio_bin_plot_total.txt
# kangjingliang@kangjingliangdeMacBook-Pro 一  3 27 01:02:21 ~/Desktop
less free_ratio_result.txt|perl -alne 'if (/^Or/i){print}else{s/\s+/\t/g;print}' >free_ratio_result.txt.1
mv free_ratio_result.txt.1 free_ratio_result.txt
# kangjingliang@kangjingliangdeMacBook-Pro 一  3 27 01:07:15 ~/Desktop
less free_ratio_result.txt|perl -alne 'if (/^Or/i){print}elsif($F[1]=~/[a-z]+/){print}' >free_ratio_result.txt.1
mv free_ratio_result.txt.1 free_ratio_result.txt
# Plot
# kangjingliang@kangjingliangdeMacBook-Pro 一  4 17 10:03:54 ~/Documents/2023/WI/paml
scp Kang@147.8.76.229:/media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input/free_ratio_result.txt ./
less free_ratio_result.txt|perl -alne 'print if $F[1]=~/[A-Z]+/'|perl -alne '$j++;if ($j==1){print "Orth_id\tspecies\tlenth\tbranch\tt\tN\tS\tdN/dS\tdN\tdS\tN*dN\tS*dS";s/\s+/\t/g;print}else{s/\s+/\t/g;print}' > free_ratio_result_plot.txt
perl temp1.pl > free_ratio_result_plot.txt.1
mv free_ratio_result_plot.txt.1 free_ratio_result_plot.txt

# select the species highest dNdS in each ortholgous genes
# kangjingliang@kangjingliangdeMacBook-Pro 一  4 17 14:50:13 ~/Documents/2023/WI/paml
perl select_max_dNdS.pl |grep -i 'common'|less # CLOCK;

# some branches have very few free ratio curated results left
##################################
# From ziheng: free ratio model
# That model is parameter rich and the estimates may involve large sampling errors. A very large omega for a very short branch, for example, does not mean much
# if you concatenate the 42 genes and then use the model, the estimates will much smaller sampling errors.  obviously they only mean averages over all codons and over all genes
##################################
# concatenate the sequences of all the orthologous genes as the input of free ratio model
# for 3826 orthologous genes in "final_orth_input_paml.txt"
# Kang@fishlab3 Mon Mar 27 09:35:58 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
perl conca_orth_paml.pl >Orth_conca_paml.fasta
fasta2phy.pl Orth_conca_paml.fasta >Orth_conca_paml.phy
# Run free ratio in SNORLAX
# (base) kang1234@celia-PowerEdge-T640 Mon Mar 27 09:38:52 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21
mkdir Orth_conca; cd Orth_conca
# Kang@fishlab3 Mon Mar 27 09:41:21 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input/paml_files/OG0000030_OG13
scp free-ratio.ctr kang1234@147.8.76.177:~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21/Orth_conca
scp spe.tre kang1234@147.8.76.177:~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21/Orth_conca
# (base) kang1234@celia-PowerEdge-T640 Mon Mar 27 09:42:14 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21/Orth_conca
ls # free-ratio.ctr free-ratio-result.txt Orth_conca_paml.phy
codeml free-ratio.ctr
# w ratios as labels for TreeView:
# (((Ocomp #0.0740261 , ((Stickleback #0.156886 , Fugu #0.121413 ) #0.139846 , (((Daru #0.0969731 , ((Acura #0.0903161 , Apoly #0.176893 ) #0.102446 , (Pmol #0.123169 , Padel #0.117265 ) #0.0980881 ) #0.0986756 ) #0.0760218 , (Blenny #0.113977 , ((Yaldwyn #0.165079 , Blueeyed #0.209281 ) #0.122188 , Common #0.125025 ) #0.134158 ) #0.107869 ) #0.214283 , (Platyfish #0.112533 , Medaka #0.118092 ) #0.141452 ) #0.096879 ) #0.143607 ) #0.126492 , Zebrafish #0.0865591 ) #0.128264 , Spottedgar #0.0001 );
# Common and Blenny didn't show a more rapid evolutionary rate than other species

####################################
# Postive selection: transfer to HPC
####################################
# Kang@fishlab3 Sat Apr 15 18:06:46 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth
nohup scp -r paml_input/ jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/Orth10_reference/ > nohup.out 2>&1
# control + z
bg
########
# Common
########
# jlkang@hpc2021 Sat Apr 15 23:14:17 /lustre1/g/sbs_schunter/Kang/Orth10_reference/paml_input
vi spe_Common.tre
# (((Ocomp,((Stickleback,Fugu),(((Daru,((Acura,Apoly),(Pmol,Padel))),(Blenny,((Yaldwyn,Blueeyed),Common #1))),(Platyfish,Medaka)))),Zebrafish),Spottedgar);
# perl codeml.pl --input temp/$temp --model branch-site --dir paml_files --output_suf Common --tree spe_Common.tre --icode 0 --omega 1.2
vi codeml_parallel.pl

# Kang@fishlab3 Sat Apr 15 23:35:34 ~/software
scp -r paml4.9j/ jlkang@hpc2021-io1.hku.hk:~/software
wc -l final_orth_input_paml.txt # 3,826 lines => split to 5 files
split -l 766 final_orth_input_paml.txt final_orth_input_paml_split_

# final_orth_input_paml_split_a[a-e]
# perl codeml_parallel.pl final_orth_input_paml_split_aa >codeml_aa.process
# perl codeml_parallel.pl final_orth_input_paml_split_ab >codeml_ab.process
# perl codeml_parallel.pl final_orth_input_paml_split_ac >codeml_ac.process
# perl codeml_parallel.pl final_orth_input_paml_split_ad >codeml_ad.process
# perl codeml_parallel.pl final_orth_input_paml_split_ae >codeml_ae.process

# install perl module
# jlkang@hpc2021 Sun Apr 16 01:06:11 /lustre1/g/sbs_schunter/Kang/Orth10_reference/paml_input
mkdir -p /home/jlkang/perl/lib/5.34.0
module load perl/5.34.0 perl-lib/5.34.0
cpanm install Parallel::ForkManager

# jlkang@hpc2021 Sat Apr 15 23:58:38 /lustre1/g/sbs_schunter/Kang/Orth10_reference/paml_input
cp ../Blueeyed_RSEM_output/script1.cmd ./ ; vi script1.cmd
vi script2.cmd; vi script3.cmd; vi script4.cmd; vi script5.cmd
sbatch script1.cmd # Submitted batch job 1108160
sbatch script2.cmd # Submitted batch job 1108161
sbatch script3.cmd # Submitted batch job 1108162
sbatch script4.cmd # Submitted batch job 1108163
sbatch script5.cmd # Submitted batch job 1108164

########
# Blenny
########
# jlkang@hpc2021 Sun Apr 16 07:46:55 /lustre1/g/sbs_schunter/Kang/Orth10_reference/paml_input
module load perl/5.34.0 perl-lib/5.34.0
vi spe_Blenny.tre
# (((Ocomp,((Stickleback,Fugu),(((Daru,((Acura,Apoly),(Pmol,Padel))),(Blenny #1,((Yaldwyn,Blueeyed),Common))),(Platyfish,Medaka)))),Zebrafish),Spottedgar);
# jlkang@hpc2021 Sun Apr 16 07:51:50 /lustre1/g/sbs_schunter/Kang/Orth10_reference/paml_input
cp codeml_parallel.pl codeml_parallel_Blenny.pl; vi codeml_parallel_Blenny.pl
vi script1_Bleeny.cmd; vi script2_Bleeny.cmd; vi script3_Bleeny.cmd; vi script4_Bleeny.cmd; vi script5_Bleeny.cmd
sbatch script1_Bleeny.cmd # Submitted batch job 1108175
sbatch script2_Bleeny.cmd # Submitted batch job 1108176
sbatch script3_Bleeny.cmd # Submitted batch job 1108177
sbatch script4_Bleeny.cmd # Submitted batch job 1108178
sbatch script5_Bleeny.cmd # Submitted batch job 1108179

##########
# Blueeyed
##########
# jlkang@hpc2021 Sun Apr 16 13:13:11 /lustre1/g/sbs_schunter/Kang/Orth10_reference/paml_input
vi spe_Blueeyed.tre
# (((Ocomp,((Stickleback,Fugu),(((Daru,((Acura,Apoly),(Pmol,Padel))),(Blenny,((Yaldwyn,Blueeyed #1),Common))),(Platyfish,Medaka)))),Zebrafish),Spottedgar);
cp codeml_parallel_Blenny.pl codeml_parallel_Blueeyed.pl; vi codeml_parallel_Blueeyed.pl
cp script1_Bleeny.cmd script1_Blueeyed.cmd ; vi script1_Blueeyed.cmd
cp script1_Blueeyed.cmd script2_Blueeyed.cmd
cp script1_Blueeyed.cmd script3_Blueeyed.cmd
cp script1_Blueeyed.cmd script4_Blueeyed.cmd
cp script1_Blueeyed.cmd script5_Blueeyed.cmd

sbatch script1_Blueeyed.cmd # Submitted batch job 1108660
sbatch script2_Blueeyed.cmd # Submitted batch job 1108661
sbatch script3_Blueeyed.cmd # Submitted batch job 1108662
sbatch script4_Blueeyed.cmd # Submitted batch job 1108663
sbatch script5_Blueeyed.cmd # Submitted batch job 1108664

##########
# Yaldwyn
##########
# jlkang@hpc2021 Sun Apr 16 17:22:17 /lustre1/g/sbs_schunter/Kang/Orth10_reference/paml_input
vi spe_Yaldwyn.tre
# (((Ocomp,((Stickleback,Fugu),(((Daru,((Acura,Apoly),(Pmol,Padel))),(Blenny,((Yaldwyn #1,Blueeyed),Common))),(Platyfish,Medaka)))),Zebrafish),Spottedgar);
cp codeml_parallel_Blenny.pl codeml_parallel_Yaldwyn.pl; vi codeml_parallel_Yaldwyn.pl
cp script1_Bleeny.cmd script1_Yaldwyn.cmd ; vi script1_Yaldwyn.cmd
cp script1_Yaldwyn.cmd script2_Yaldwyn.cmd
cp script1_Yaldwyn.cmd script3_Yaldwyn.cmd
cp script1_Yaldwyn.cmd script4_Yaldwyn.cmd
cp script1_Yaldwyn.cmd script5_Yaldwyn.cmd

sbatch script1_Yaldwyn.cmd # Submitted batch job 1108780
sbatch script2_Yaldwyn.cmd # Submitted batch job 1108782
sbatch script3_Yaldwyn.cmd # Submitted batch job 1108783
sbatch script4_Yaldwyn.cmd # Submitted batch job 1108785
sbatch script5_Yaldwyn.cmd # Submitted batch job 1108786
```

```bash
# Kang@fishlab3 Tue Apr 18 13:33:18 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
scp -r jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/Orth10_reference/paml_input/postively_selected_genes ./
scp kang1234@147.8.76.177:~/genome/paml_new/paml_input/Extract_PSGs.pl ./
scp kang1234@147.8.76.177:~/swiss-prot/unprot_name_description.txt ./
```

```extract_ano.pl
#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

my $TANO="../all_swissprot_diamond_ano.txt";
my %tano;
open TANO, $TANO or die "can not open $TANO\n";
while (<TANO>) {
	chomp;
	my @a=split;
	$tano{$a[0]}=$a[1]; # Acura_1 => sp|Q498W5|T198B_DANRE
}

my $UNIP="unprot_name_description.txt";
my %unip;
open UNIP, $UNIP or die "can not open $UNIP\n";
while (<UNIP>) {
	chomp;
	my @a=split /\t/;
	$unip{$a[0]}=$a[1]; # sp|P30604|CHS2_XYLBA => Chitin synthase 2 (Fragment)
}

my $SUB="../sub_orth_id_final.txt";
open SUB, $SUB or die "can not open $SUB\n";
while (<SUB>) {
	chomp;
	my @a=split /\t/;
	if (/^Suborth/) {
		next;
	} else {
		my @b=split /\,/, $a[1];
		my $orth=$a[0];
		my $name=$tano{$b[0]};
		my $desc=$unip{$name};
#		print "$orth\t$b[0]\t$name\t$a[-1]\t$desc\n";
		print "$orth\t$name\t$a[-1]\t$desc\n";
	}
}
```

```bash
# Kang@fishlab3 Tue Apr 18 15:09:13 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
perl extract_ano.pl > ortho16_ano.txt
perl Extract_PSGs.pl ortho16_ano.txt >orth16_PSGs.txt
```


## BAMM estimate macroevolution
### 1. Estimate the divergence time (MCMCTree)
```bash
# the input is the coding nucleotide sequences alignments of the 47 single copy genes across 16 fish species
# extract the cds seq of 47 single copy genes
# Kang@fishlab3 Fri Mar 24 19:06:12 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
# perl Get_single_copy_cds.pl > sing_copy_nuc.fasta

# use "/media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input/Orth_conca_paml.phy" as the input of MCMCTree
# use the nucleotide sequences of single copy genes as the input of MCMCTree "/media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input/sing_copy_nuc.phy"
# Kang@fishlab3 Mon Mar 27 09:45:45 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
mkdir mcmctree; cd mcmctree
# Kang@fishlab3 Mon Mar 27 10:19:28 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input/mcmctree
vi baseml.tree
# 16 1
# (((Ocomp,((Stickleback,Fugu),(((Daru,((Acura,Apoly),(Pmol,Padel))),(Blenny,((Yaldwyn,Blueeyed),Common))),(Platyfish,Medaka)))'@(.969,1.509)'),Zebrafish)'@(1.4985, 1.652)',Spottedgar);
mkdir baseml; cd baseml
# Kang@fishlab3 Mon Mar 27 10:22:42 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input/mcmctree/baseml
cp /media/HDD/cleaner_fish/genome/gene_family_2/single_copy_id/mcmctree/baseml/baseml.ctl ./
# Kang@fishlab3 Mon Mar 27 10:28:07 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input/mcmctree
cd baseml; baseml # run baseml
```
**Substitution rate is per time unit**(in ~/Desktop/PapueNewGuinea-new/longest_pep_sec/input_pep/mcmctree/baseml/mlb)   
```bash
# in "mlb"
# Substitution rate is per time unit
#    0.019895 +-     -nan  # (use this value to set the prior for the mean substitution rate in the Bayesian analysis)
# a = (m/s)^2 and b = m/s^2 (m=s=0.019895)
# rgene_gamma = 1 50.3
```
**Estimation of the gradient and Hessian**   
```bash
# Kang@fishlab3 Mon Mar 27 11:07:47 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input/mcmctree
vi mcmc.tree
# 16 2
# (((Ocomp,((Stickleback,Fugu),(((Daru,((Acura,Apoly),(Pmol,Padel))),(Blenny,((Yaldwyn,Blueeyed),Common))),(Platyfish,Medaka)))'B(.969, 1.509)'),Zebrafish)'B(1.4885, 1.652)',Spottedgar);
# Kang@fishlab3 Mon Mar 27 11:15:20 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input/mcmctree/baseml
cp ~/Desktop/PapueNewGuinea-new/longest_pep_sec/input_pep/mcmctree/baseml/mcmctree*.ctl ./
# in the control file (mcmctree.ctl): set "usedata = 3", "ndata = 1" (no partition)
# Kang@fishlab3 Mon Mar 27 11:22:38 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input/mcmctree/baseml
nohup ~/software/paml4.9j/src/mcmctree mcmctree1.ctl > mcmctree1.process 2>&1 &
# [1] 31395
mv out.BV in.BV
cat in.BV rst2 > 1.txt
mv 1.txt in.BV

nohup ~/software/paml4.9j/src/mcmctree mcmctree2.ctl >mcmctree2.process 2>&1 &
# [1] 31639
```
### 2. Orthofinder for the ten fish species sampled from CO2 seep
```bash
# Kang@fishlab3 Tue Apr 04 07:15:47 /media/HDD/white_island/Compevo
mkdir Orth_ten_CO2; cd Orth_ten_CO2
# Kang@fishlab3 Tue Apr 04 07:16:33 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
cp Acura.pep.fasta Apoly.pep.fasta Daru.pep.fasta Ocomp.pep.fasta Padel.pep.fasta Pmol.pep.fasta Blenny.pep.fasta Blueeyed.pep.fasta Common.pep.fasta Yaldwyn.pep.fasta /media/HDD/white_island/Compevo/Orth_ten_CO2

# move to SNORLAX for orthologous genes detection
# (base) kang1234@celia-PowerEdge-T640 Tue Apr 04 07:19:49 ~/white_island/Compevo
scp -r Kang@147.8.76.229:/media/HDD/white_island/Compevo/Orth_ten_CO2 ./
# (base) kang1234@celia-PowerEdge-T640 Tue Apr 04 07:22:21 ~/white_island/Compevo/Orth_ten_CO2
nohup orthofinder -f ./ -a 22 >orthofinder-process 2>&1 &
# [1] 9623

# 252 orthologous genes among the ten fish species: all orthologous genes should be longer than 400
# Kang@fishlab3 Tue Apr 04 14:41:06 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
scp Acura.cds.fasta Apoly.cds.fasta Daru.cds.fasta Ocomp.cds.fasta Padel.cds.fasta Pmol.cds.fasta Blenny.cds.fasta Blueeyed.cds.fasta Common.cds.fasta Yaldwyn.cds.fasta kang1234@147.8.76.177:~/white_island/Compevo/Orth_ten_CO2
# extract the nucleotide sequences of 50 orthologous genes with pep length > 400 as bite sequences
```

```extract_50_single_copy_400peps.pl
use strict;
use warnings;
use File::Basename;

my $i;
my %nuc;
my @cds =</home/kang1234/white_island/Compevo/Orth_ten_CO2/*.cds.fasta>;
foreach my $cds (@cds) {
        open CDS, $cds or die "can not open $cds\n";
        my $name;
        while (<CDS>) {
                chomp;
                if (/>/) {
                        s/\>//;
                        $name=$_;
                } else {
                        $nuc{$name}.=$_;
                }
        }
        close CDS;
}
my @spes=qw(Acura Apoly Blenny Blueeyed Common Daru Ocomp Padel Pmol Yaldwyn);
my @fas=<*.fa>;
foreach my $fa (@fas) {
        open FA, $fa or die "can not open $fa\n";
        my %hash;
        my ($len, $ort, $gene, $spe, $info);
        ($ort)=$fa=~/(.*)\.fa/;
        while (<FA>) {
                chomp;
                if (/^>/) {
                        s/>//;
                        ($spe)=$_=~/(.*)\_.*/;
                        $gene=$ort."-".$_;
                        $hash{$spe}={
                                'GENE'=> $gene,
                                'NAME'=> $_
                        };
                } else {
                        $len=length($_);
                        last if $len<=400;
                }
        }
        if ($len>400) {
                $i++;
                if ($i<=50) {
                        foreach my $spe (@spes) {
                                my $file=$spe."_50_single_copy.fa";
                                open FILE, ">>$file" or die "can not create $file\n";
                                my $gene=$hash{$spe}->{'GENE'};
                                my $name=$hash{$spe}->{'NAME'};
                                my $cds =$nuc{$name};
                                print FILE ">$gene\n$cds\n";
                        }
                        close FILE;
                }
        }
        close FA;
}
```

```bash
# (base) kang1234@celia-PowerEdge-T640 Tue Apr 04 15:06:19 ~/white_island/Compevo/Orth_ten_CO2/OrthoFinder/Results_Apr04/Single_Copy_Orthologue_Sequences
perl extract_50_single_copy_400peps.pl
ls *_50_single_copy.fa
# Acura_50_single_copy.fa  Blenny_50_single_copy.fa    Common_50_single_copy.fa  Ocomp_50_single_copy.fa  Pmol_50_single_copy.fa
# Apoly_50_single_copy.fa  Blueeyed_50_single_copy.fa  Daru_50_single_copy.fa    Padel_50_single_copy.fa  Yaldwyn_50_single_copy.fa
# (base) kang1234@celia-PowerEdge-T640 Tue Apr 04 15:09:39 ~/white_island/Compevo/Orth_ten_CO2/OrthoFinder/Results_Apr04
mkdir gene_capture_single_copy_nuc
mv Single_Copy_Orthologue_Sequences/*_50_single_copy.fa gene_capture_single_copy_nuc/
# (base) kang1234@celia-PowerEdge-T640 Tue Apr 04 17:12:31 ~/white_island/kraken
vi Ind_info.txt
# B1	Common	Vs
# B2	Common	Vs
# B3	Common	Vs
# B4	Common	Vs

# (base) kang1234@celia-PowerEdge-T640 Tue Apr 04 17:34:42 ~/white_island/kraken
cp ~/white_island/Compevo/Orth_ten_CO2/OrthoFinder/Results_Apr04/gene_capture_single_copy_nuc/*_50_single_copy.fa ./
makeblastdb -in Blenny_50_single_copy.fa -out Blenny -dbtype nucl
makeblastdb -in Blueeyed_50_single_copy.fa -out Blueeyed -dbtype nucl
makeblastdb -in Common_50_single_copy.fa -out Common -dbtype nucl
makeblastdb -in Yaldwyn_50_single_copy.fa -out Yaldwyn -dbtype nucl

# mkdir for gene capture
# change the fastq name, and copy to results/
mkdir gene_capture
```

```cp_change_fq_nm.pl
#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

my (%hash1, %hash2);
my $indi="Ind_info.txt";
open INDI, $indi or die "can not open $indi\n";
while (<INDI>) {
	chomp;
	my @a=split;
	$hash1{$a[1]}++;
	my ($ind, $spe, $id);
	$spe=$a[1];
	$ind=$spe.$hash1{$a[1]};
	$id =$a[0];
	$hash2{$id}={
		'IND' => $ind,
		'SPE' => $spe
	};
}

my @fq=<*.fastq.gz>;
my @cmds;
foreach my $fq (@fq) {
	my ($id, $info)=$fq=~/(.*)_(\d+\.fastq\.gz)/;
	if ($hash2{$id}) {
		my $spe=$hash2{$id}->{'SPE'};
		my $ind=$hash2{$id}->{'IND'};

		my $new=$ind."_".$info;
		my $cmd="cp $fq gene_capture/$new";
		push @cmds, $cmd;
	}
}

my $manager = new Parallel::ForkManager(5);
foreach my $cmd (@cmds) {
    $manager->start and next;
    system($cmd);
    $manager->finish;
}
$manager -> wait_all_children;
```

```bash
# (base) kang1234@celia-PowerEdge-T640 Tue Apr 04 18:21:24 ~/white_island/kraken
nohup perl cp_change_fq_nm.pl >cp_change.process 2>&1 &
# [1] 6608
# (base) kang1234@celia-PowerEdge-T640 Tue Apr 04 22:52:05 ~/white_island/kraken/gene_capture
mv ../Blenny_50_single_copy.fa ../Blueeyed_50_single_copy.fa ../Yaldwyn_50_single_copy.fa ../Common_50_single_copy.fa ./
cp ~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/*.pm ./

# subset reads (0.1)
# (base) kang1234@celia-PowerEdge-T640 Wed Apr 05 09:31:32 ~/white_island/kraken/gene_capture
nohup perl subset_fq.pl >subset_fq.process 2>&1 &
# [1] 18562

# 1. rmrep.pl
# (base) kang1234@celia-PowerEdge-T640 Wed Apr 05 15:27:58 ~/white_island/kraken/gene_capture
cp ~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/rmrep.pl ./
nohup perl gene_capture_parallel.pl >gene_capture.process 2>&1 &
# [1] 27548

# 2. bandp.pl
# (base) kang1234@celia-PowerEdge-T640 Wed Apr 05 16:49:26 ~/white_island/kraken/gene_capture
cp ~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/bandp.pl ./
nohup perl run_bandp_parallel.pl > run_bandp_parallel.process 2>&1 &
# [1] 6384

# 3. Trinity.pl
# (base) kang1234@celia-PowerEdge-T640 Thu Apr 06 06:40:51 ~/white_island/kraken/gene_capture
cp ~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/Trinity_new.pl ./
export SHARED_DIR=$PWD
nohup ./Trinity_new.pl -species=Blenny -output=trinity >Blenny_trinity.process 2>&1 &
# [2] 13030
nohup ./Trinity_new.pl -species=Blueeyed -output=trinity >Blueeyed_trinity.process 2>&1 &
# [3] 3415
nohup ./Trinity_new.pl -species=Common -output=trinity >Common_trinity.process 2>&1 &
# [4] 4660
nohup ./Trinity_new.pl -species=Yaldwyn -output=trinity >Yaldwyn_trinity.process 2>&1 &
# [5] 7148

# 4. run_getbest.pl
# (base) kang1234@celia-PowerEdge-T640 Thu Apr 06 10:59:47 ~/white_island/kraken/gene_capture
cp ~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/run_getbest.pl ./
cp ~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/getbest_new.pl ./
mkdir result

ll -d Blenny*_results|perl -alne '(my $ind)=$F[-1]=~/(.*)\_/;$info.="$ind ";END{print $info}'|less
ll -d Blueeyed*_results|perl -alne '(my $ind)=$F[-1]=~/(.*)\_/;$info.="$ind ";END{print $info}'|less
ll -d Common*_results|perl -alne '(my $ind)=$F[-1]=~/(.*)\_/;$info.="$ind ";END{print $info}'|less
ll -d Yaldwyn*_results|perl -alne '(my $ind)=$F[-1]=~/(.*)\_/;$info.="$ind ";END{print $info}'|less
```

```run_getbest.sh
perl run_getbest.pl -query="Blenny_50_single_copy" -subject="Blenny10 Blenny11 Blenny12 Blenny13 Blenny14 Blenny15 Blenny16 Blenny1 Blenny2 Blenny3 Blenny4 Blenny5 Blenny6 Blenny7 Blenny8 Blenny9"
perl run_getbest.pl -query="Blueeyed_50_single_copy" -subject="Blueeyed10 Blueeyed11 Blueeyed12 Blueeyed13 Blueeyed14 Blueeyed15 Blueeyed1 Blueeyed2 Blueeyed3 Blueeyed4 Blueeyed5 Blueeyed6 Blueeyed7 Blueeyed8 Blueeyed9"
perl run_getbest.pl -query="Common_50_single_copy" -subject="Common10 Common11 Common12 Common13 Common14 Common15 Common1 Common2 Common3 Common4 Common5 Common6 Common7 Common8 Common9"
perl run_getbest.pl -query="Yaldwyn_50_single_copy" -subject="Yaldwyn10 Yaldwyn11 Yaldwyn12 Yaldwyn13 Yaldwyn14 Yaldwyn15 Yaldwyn1 Yaldwyn2 Yaldwyn3 Yaldwyn4 Yaldwyn5 Yaldwyn6 Yaldwyn7 Yaldwyn8 Yaldwyn9"
```

```bash
# (base) kang1234@celia-PowerEdge-T640 Thu Apr 06 16:42:46 ~/white_island/kraken/gene_capture
nohup sh run_getbest.sh > getbest.process 2>&1 &
# [1] 18317

# 5. reblast.pl
# (base) kang1234@celia-PowerEdge-T640 Thu Apr 06 19:45:22 ~/white_island/kraken/gene_capture/result
cp ../*_50_single_copy.fa ./
cp ~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/result/reblast.pl ./
makeblastdb -in Blenny_50_single_copy.fa -out Blenny -dbtype nucl
makeblastdb -in Blueeyed_50_single_copy.fa -out Blueeyed -dbtype nucl
makeblastdb -in Common_50_single_copy.fa -out Common -dbtype nucl
makeblastdb -in Yaldwyn_50_single_copy.fa -out Yaldwyn -dbtype nucl
```

```run_reblast.sh
perl reblast.pl -query Blenny_50_single_copy.resultnf_new -database Blenny
perl reblast.pl -query Blueeyed_50_single_copy.resultnf_new -database Blueeyed
perl reblast.pl -query Common_50_single_copy.resultnf_new -database Common
perl reblast.pl -query Yaldwyn_50_single_copy.resultnf_new -database Yaldwyn
```

```bash
# (base) kang1234@celia-PowerEdge-T640 Thu Apr 06 21:56:19 ~/white_island/kraken/gene_capture/result
nohup sh run_reblast.sh > run_reblast.process 2>&1 &
# [1] 22169
```
### RNA-seq data of PNG individuals transfer to SNORALX (CO2 seep: 5 inviduals; Control: 5 inviduals)
```bash
# Kang@fishlab3 Thu Apr 06 11:28:51 ~/Desktop/PapueNewGuinea-new/merge_clean
vi Indviduals_transfer.txt

# (base) kang1234@celia-PowerEdge-T640 Thu Apr 06 18:04:49 ~/white_island/kraken
mkdir gene_capture_PNG
# Kang@fishlab3 Thu Apr 06 18:08:54 ~/Desktop/PapueNewGuinea-new/merge_clean
less Indviduals_transfer.txt|perl -alne '$ind.=$F[0]."*.fastq.gz ";END{$ind=~s/\s+$//;print "nohup scp $ind kang1234\@147.8.76.177:~/white_island/kraken/gene_capture_PNG > nohup.out 2>&1"}'|less
nohup scp Apoly21*.fastq.gz Apoly23*.fastq.gz Apoly24*.fastq.gz Apoly26*.fastq.gz Apoly27*.fastq.gz Apoly44*.fastq.gz Apoly45*.fastq.gz Apoly46*.fastq.gz Apoly47*.fastq.gz Apoly48*.fastq.gz Acura3*.fastq.gz Acura4*.fastq.gz Acura5*.fastq.gz Acura6*.fastq.gz Acura7*.fastq.gz Acura24*.fastq.gz Acura25*.fastq.gz Acura26*.fastq.gz Acura27*.fastq.gz Acura28*.fastq.gz Ocomp1*.fastq.gz Ocomp2*.fastq.gz Ocomp3*.fastq.gz Ocomp4*.fastq.gz Ocomp5*.fastq.gz Ocomp25*.fastq.gz Ocomp26*.fastq.gz Ocomp27*.fastq.gz Ocomp28*.fastq.gz Ocomp29*.fastq.gz Daru1*.fastq.gz Daru2*.fastq.gz Daru3*.fastq.gz Daru4*.fastq.gz Daru6*.fastq.gz Daru25*.fastq.gz Daru26*.fastq.gz Daru27*.fastq.gz Daru29*.fastq.gz Daru30*.fastq.gz Padel6*.fastq.gz Padel7*.fastq.gz Padel8*.fastq.gz Padel9*.fastq.gz Padel10*.fastq.gz Padel23*.fastq.gz Padel24*.fastq.gz Padel25*.fastq.gz Padel26*.fastq.gz Padel28*.fastq.gz Pmol7*.fastq.gz Pmol8*.fastq.gz Pmol9*.fastq.gz Pmol10*.fastq.gz Pmol11*.fastq.gz Pmol21*.fastq.gz Pmol22*.fastq.gz Pmol23*.fastq.gz Pmol24*.fastq.gz Pmol25*.fastq.gz kang1234@147.8.76.177:~/white_island/kraken/gene_capture_PNG > nohup.out 2>&1
# Ctrl+Z
bg
# (base) kang1234@celia-PowerEdge-T640 Fri Apr 07 06:42:45 ~/white_island/kraken/gene_capture_PNG
cp ~/white_island/kraken/gene_capture/*.pl ./
cp ~/white_island/kraken/gene_capture/*.pm ./
```

```bash
# (base) kang1234@celia-PowerEdge-T640 Fri Apr 07 06:44:22 ~/white_island/kraken/gene_capture_PNG
cp ~/white_island/Compevo/Orth_ten_CO2/OrthoFinder/Results_Apr04/gene_capture_single_copy_nuc/*_50_single_copy.fa ./
rm Blenny_50_single_copy.fa Blueeyed_50_single_copy.fa Common_50_single_copy.fa Yaldwyn_50_single_copy.fa
makeblastdb -in Acura_50_single_copy.fa -out Acura -dbtype nucl
makeblastdb -in Apoly_50_single_copy.fa -out Apoly -dbtype nucl
makeblastdb -in Daru_50_single_copy.fa -out Daru -dbtype nucl
makeblastdb -in Ocomp_50_single_copy.fa -out Ocomp -dbtype nucl
makeblastdb -in Padel_50_single_copy.fa -out Padel -dbtype nucl
makeblastdb -in Pmol_50_single_copy.fa -out Pmol -dbtype nucl

# 1. rmrep.pl
# (base) kang1234@celia-PowerEdge-T640 Fri Apr 07 06:51:05 ~/white_island/kraken/gene_capture_PNG
nohup perl gene_capture_parallel.pl >gene_capture.process 2>&1 &
# [1] 4075

# 2. bandp.pl
# (base) kang1234@celia-PowerEdge-T640 Fri Apr 07 21:35:52 ~/white_island/kraken/gene_capture_PNG
nohup perl run_bandp_parallel.pl > run_bandp_parallel.process 2>&1 &
# [1] 12987

# 3. Trinity.pl
# (base) kang1234@celia-PowerEdge-T640 Sat Apr 08 21:17:04 ~/white_island/kraken/gene_capture_PNG
export SHARED_DIR=$PWD
nohup ./Trinity_new.pl -species=Acura -output=trinity >Acura_trinity.process 2>&1 &
# [1] 10133
nohup ./Trinity_new.pl -species=Apoly -output=trinity >Blueeyed_trinity.process 2>&1 &
# [2] 19002
nohup ./Trinity_new.pl -species=Daru -output=trinity >Daru_trinity.process 2>&1 &
# [3] 22622
nohup ./Trinity_new.pl -species=Ocomp -output=trinity >Ocomp_trinity.process 2>&1 &
# [4] 7487
nohup ./Trinity_new.pl -species=Padel -output=trinity >Padel_trinity.process 2>&1 &
# [5] 8385
nohup ./Trinity_new.pl -species=Pmol -output=trinity >Pmol_trinity.process 2>&1 &
# [6] 17781

# 4. run_getbest.pl
# (base) kang1234@celia-PowerEdge-T640 Sat Apr 08 21:42:27 ~/white_island/kraken/gene_capture_PNG
mkdir result

ll -d Acura*_results|perl -alne '(my $ind)=$F[-1]=~/(.*)\_/;$info.="$ind ";END{print $info}'|less
ll -d Apoly*_results|perl -alne '(my $ind)=$F[-1]=~/(.*)\_/;$info.="$ind ";END{print $info}'|less
ll -d Daru*_results|perl -alne '(my $ind)=$F[-1]=~/(.*)\_/;$info.="$ind ";END{print $info}'|less
ll -d Ocomp*_results|perl -alne '(my $ind)=$F[-1]=~/(.*)\_/;$info.="$ind ";END{print $info}'|less
ll -d Padel*_results|perl -alne '(my $ind)=$F[-1]=~/(.*)\_/;$info.="$ind ";END{print $info}'|less
ll -d Pmol*_results|perl -alne '(my $ind)=$F[-1]=~/(.*)\_/;$info.="$ind ";END{print $info}'|less
```

```run_getbest.sh
perl run_getbest.pl -query="Acura_50_single_copy" -subject="Acura24 Acura25 Acura26 Acura27 Acura28 Acura3 Acura4 Acura5 Acura6 Acura7"
perl run_getbest.pl -query="Apoly_50_single_copy" -subject="Apoly21 Apoly23 Apoly24 Apoly26 Apoly27 Apoly44 Apoly45 Apoly46 Apoly47 Apoly48"
perl run_getbest.pl -query="Daru_50_single_copy" -subject="Daru1 Daru25 Daru26 Daru27 Daru29 Daru2 Daru30 Daru3 Daru4 Daru6"
perl run_getbest.pl -query="Ocomp_50_single_copy" -subject="Ocomp1 Ocomp25 Ocomp26 Ocomp27 Ocomp28 Ocomp29 Ocomp2 Ocomp3 Ocomp4 Ocomp5"
perl run_getbest.pl -query="Padel_50_single_copy" -subject="Padel10 Padel23 Padel24 Padel25 Padel26 Padel28 Padel6 Padel7 Padel8 Padel9"
perl run_getbest.pl -query="Pmol_50_single_copy" -subject="Pmol10 Pmol11 Pmol21 Pmol22 Pmol23 Pmol24 Pmol25 Pmol7 Pmol8 Pmol9"
```

```bash
# (base) kang1234@celia-PowerEdge-T640 Sat Apr 08 21:49:57 ~/white_island/kraken/gene_capture_PNG
nohup sh run_getbest.sh > getbest.process 2>&1 &
# [1] 4554

# 4. reblast.pl
# (base) kang1234@celia-PowerEdge-T640 Sun Apr 09 19:09:10 ~/white_island/kraken/gene_capture_PNG/result
cp ../*_50_single_copy.fa ./
cp ~/CO2-seeps/high_index/paired/kraken/merge/gene_capture/result/reblast.pl ./
makeblastdb -in Acura_50_single_copy.fa -out Acura -dbtype nucl
makeblastdb -in Apoly_50_single_copy.fa -out Apoly -dbtype nucl
makeblastdb -in Daru_50_single_copy.fa -out Daru -dbtype nucl
makeblastdb -in Ocomp_50_single_copy.fa -out Ocomp -dbtype nucl
makeblastdb -in Padel_50_single_copy.fa -out Padel -dbtype nucl
makeblastdb -in Pmol_50_single_copy.fa -out Pmol -dbtype nucl
```

```run_reblast.sh
perl reblast.pl -query Acura_50_single_copy.resultnf_new -database Acura
perl reblast.pl -query Apoly_50_single_copy.resultnf_new -database Apoly
perl reblast.pl -query Daru_50_single_copy.resultnf_new -database Daru
perl reblast.pl -query Ocomp_50_single_copy.resultnf_new -database Ocomp
perl reblast.pl -query Padel_50_single_copy.resultnf_new -database Padel
perl reblast.pl -query Pmol_50_single_copy.resultnf_new -database Pmol
```

```bash
# (base) kang1234@celia-PowerEdge-T640 Sun Apr 09 19:15:31 ~/white_island/kraken/gene_capture_PNG/result
nohup sh run_reblast.sh > run_reblast.process 2>&1 &
# [1] 9136

# copy white island individuals to here
# (base) kang1234@celia-PowerEdge-T640 Sun Apr 09 19:40:35 ~/white_island/kraken/gene_capture_PNG/result/genebin
cp ~/white_island/kraken/gene_capture/result/genebin/*.fas ./
ll *.fas|perl -alne 'if ($F[4]==0){`rm $F[-1]`}' # delete the empty file
# cat all sequences with the same orthogroup id together 
ll *.fas|perl -alne '(my $na)=$F[-1]=~/(.*)-/;$hash{$na}++;push @na, $na if $hash{$na}==1;END{foreach my $na (@na){print"$na\t$hash{$na}"}}'|perl -alne '`cat $F[0]*.fas > $F[0].fas`'
rm *-*_*.fas
```

**align, trim, concatenate**   
```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;

# muscle -in OG0008797.fas -out OG0008797-align.fas
# trimal -in OG0008797-align.fas -out OG0008797-align-trim.fas -gt 0.9 -st 0.001 -cons 60

my @fasta=<*.fas>;
foreach my $fasta (@fasta) {
	(my $name)=$fasta=~/(.*)\.fas/;
	my $align=$name."-align.fas";
	my $trim=$name."-trim.fas";
	my $cmd1="muscle -in $fasta -out $align";
	my $cmd2="trimal -in $align -out $trim -gt 0.9 -st 0.001 -cons 60";
	system($cmd1);
	system($cmd2);
}

my %hash; 
my $gene; my @genes;
my $i;
my @trims=<*-trim.fas>;
foreach my $trim (@trims) {
	$i++;
	open TRIM, "$trim";
	while (<TRIM>) {
		chomp;
		if (/>/) {
			s/>//;
			$gene=$_;
			push @genes, $gene if $i==1;
		} else {
			$hash{$gene}.=$_;
		}
	}
}

open FASTA, ">All_gene_concatenated.fasta";
foreach my $gen (@genes) {
	my $seq=$hash{$gen};
	print FASTA ">$gen\n$seq\n";
}
```

```bash
# (base) kang1234@celia-PowerEdge-T640 Sun Apr 09 20:57:03 ~/white_island/kraken/gene_capture_PNG/result/genebin
perl Check_capture_info.pl > Capture_info.txt
# Check which individual didn't capture the orthologous genes
less Capture_info.txt|perl -alne '$j++;@head=@F if $j==1;my $a;for (my $i = 1; $i < @F; $i++){$a.=$head[$i].";" if $F[$i]==0};print "$F[0]\t$a" if $j>1'|less
# remove: Daru3; Blenny6; Blueeyed8; Blueeyed12; Common5; Common7
# (base) kang1234@celia-PowerEdge-T640 Sun Apr 09 22:02:00 ~/white_island/kraken/gene_capture_PNG/result/genebin
vi selected_orth.txt # 43 orthlogous genes; "remmove OG0061615 OG0061823 OG0061956 OG0062083 OG0062498 OG0062526"
mkdir select
```

```select.pl
#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

my (%qind, %qorth);
# don't use the sequences of these individuals as they didn't captured genes in enough orthologous genes
my @qinds=qw(Daru3 Blenny6 Blueeyed8 Blueeyed12 Common5 Common7);
foreach my $ind (@qinds) {
	$qind{$ind}++;
}

# don't use the sequences of these orthologous genes as they are not captured in enough indviduals
my @qorth=qw(OG0061615 OG0061823 OG0061956 OG0062083 OG0062498 OG0062526);
foreach my $orth (@qorth) {
	$qorth{$orth}++;
}

my @fas=<*.fas>;
foreach my $fa (@fas) {
	(my $orth)=$fa=~/(.*)\.fas/;
	if ($qorth{$orth}) {
		next;
	} else {
		my $ind;
		my %hash;
		open FAS, $fa or die "can not open $fa\n";
		while (<FAS>) {
			chomp;
			if (/\>/) {
				s/\>//;
				$ind=$_;
			} else {
				$hash{$ind}.=$_;
			}
		}
		close FAS;

		my $newfas="select/$fa";
		open NEWFAS, ">$newfas" or die "can not create $newfas\n";
		foreach my $ind (sort keys %hash) {
			if ($qind{$ind}) {
				next;
			} else {
				my $seq=$hash{$ind};
				print NEWFAS ">$ind\n$seq\n";
			}
		}
		close NEWFAS;
	}
}
```

```bash
# (base) kang1234@celia-PowerEdge-T640 Tue Apr 11 11:42:16 ~/white_island/kraken/gene_capture_PNG/result/genebin
perl select.pl # resulted fasta files in select/

# (base) kang1234@celia-PowerEdge-T640 Tue Apr 11 11:44:08 ~/white_island/kraken/gene_capture_PNG/result/genebin/select
cp ../temp1.pl ./
perl temp1.pl
fasta2phy.pl All_gene_concatenated.fasta >All_gene_concatenated.phy
```

## Detect ORF and then align/trim
**Transfer to my own work station for next step**    
```bash
# Kang@fishlab3 Tue Apr 11 13:51:07 /media/HDD/white_island/Compevo
mkdir genecapture; cd genecapture/
scp kang1234@147.8.76.177:~/white_island/kraken/gene_capture_PNG/result/ORFextract ./
scp kang1234@147.8.76.177:~/white_island/kraken/gene_capture_PNG/result/ORFextract.sh ./

# Kang@fishlab3 Tue Apr 11 13:53:18 /media/HDD/white_island/Compevo/genecapture
# get the longest ORF of all transcripts of each gene for all individuals
# ORFextract.sh
perl ORFextract --dir Acura_50_single_copy.resultnf_new --out Acura_orf
perl ORFextract --dir Apoly_50_single_copy.resultnf_new --out Apoly_orf
perl ORFextract --dir Daru_50_single_copy.resultnf_new --out Daru_orf
perl ORFextract --dir Ocomp_50_single_copy.resultnf_new --out Ocomp_orf
perl ORFextract --dir Padel_50_single_copy.resultnf_new --out Padel_orf
perl ORFextract --dir Pmol_50_single_copy.resultnf_new --out Pmol_orf
perl ORFextract --dir Blenny_50_single_copy.resultnf_new --out Blenny_orf
perl ORFextract --dir Blueeyed_50_single_copy.resultnf_new --out Blueeyed_orf
perl ORFextract --dir Common_50_single_copy.resultnf_new --out Common_orf
perl ORFextract --dir Yaldwyn_50_single_copy.resultnf_new --out Yaldwyn_orf

# Kang@fishlab3 Tue Apr 11 13:53:18 /media/HDD/white_island/Compevo/genecapture
nohup sh ORFextract.sh > run_ORFextract.process 2>&1 &
# [1] 13853
ll Ocomp_orf/*.fas|perl -alne '(my $n)=$F[-1]=~/\/(.*)-/;print $n' >1.txt
```

```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

my @spes=qw(Daru Ocomp Pmol Padel Acura Apoly Blenny Blueeyed Common Yaldwyn);
open TXT, "1.txt" or die "can not open 1.txt";
while (<TXT>) {
	chomp;
	my $cmd2="cat ";
	my $gene=$_;
	foreach my $spe (@spes) {
		my $dir=$spe."_orf";
		my $gene1=$dir."/".$gene."\*".".fas";
		$cmd2.=$gene1." ";
	}
	$cmd2.=">$gene.fas";
	`mkdir Total_orf` unless -d "Total_orf";
	system($cmd2);
	system("mv $gene.fas Total_orf/");
}
```

```bash
# Kang@fishlab3 Tue Apr 11 14:06:02 /media/HDD/white_island/Compevo/genecapture
perl temp1.pl # the results in "Total_orf/"
# cat: 'Apoly_orf/OG0061823*.fas': No such file or directory
# cat: 'Common_orf/OG0061823*.fas': No such file or directory

# Kang@fishlab3 Tue Apr 11 14:10:56 /media/HDD/white_island/Compevo/genecapture/Total_orf
scp kang1234@147.8.76.177:~/white_island/kraken/gene_capture_PNG/result/genebin/Check_capture_info.pl ./
perl Check_capture_info.pl > Capture_info.txt
less Capture_info.txt|perl -alne '$j++;@head=@F if $j==1;my $a;for (my $i = 1; $i < @F; $i++){$a.=$head[$i].";" if $F[$i]==0};print "$F[0]\t$a" if $j>1'|less
# if the orthologous gene didn't capture more than two individuals, remove this ortholgous gene
less Capture_info.txt|perl -alne '$j++;@head=@F if $j==1;my $a;for (my $i = 1; $i < @F; $i++){$a.=$head[$i].";" if $F[$i]==0};print "$F[0]\t$a" if $j>1'|perl -alne '@b=split /\;/, $F[1];print $F[0] if @b>1'|less
# remove orthologous genes (21 ortholgous genes will be removed)
# OG0061742 OG0061823 OG0061909 OG0061954 OG0061956 OG0062169 OG0062324 
# OG0062344 OG0062473 OG0062487 OG0062488 OG0062498 OG0062526 OG0062611 
# OG0062625 OG0062672 OG0062769 OG0062782 OG0062880 OG0062896 OG0062897

# remove individuals (three individuals will be removed)
less Capture_info.txt|perl -alne '$j++;@head=@F if $j==1;my $a;for (my $i = 1; $i < @F; $i++){$a.=$head[$i].";" if $F[$i]==0};print "$F[0]\t$a" if $j>1'|perl -alne '@b=split /\;/, $F[1];print if @b<=1'|perl -alne '@a=split /\;/, $F[1];foreach my $in (@a){print $in}'|sort -u
# Apoly48 Daru3 Ocomp4

# Kang@fishlab3 Tue Apr 11 14:34:10 /media/HDD/white_island/Compevo/genecapture/Total_orf
mkdir select
perl select.pl
# align and trim
# Kang@fishlab3 Tue Apr 11 14:36:44 /media/HDD/white_island/Compevo/genecapture/Total_orf/select
scp kang1234@147.8.76.177:~/white_island/kraken/gene_capture_PNG/result/genebin/select/temp1.pl ./
perl temp1.pl 
fasta2phy.pl All_gene_concatenated.fasta >All_gene_concatenated.phy
# Transfer to SNORALX for the phylogenetic tree construction
# Kang@fishlab3 Tue Apr 11 15:02:35 /media/HDD/white_island/Compevo/genecapture/Total_orf/select
scp All_gene_concatenated.phy kang1234@147.8.76.177:~/white_island/kraken/gene_capture_PNG/
# set all Ocomp individuals as root
nohup raxmlHPC -f a -x 12345 -p 12345 -# 100 -m GTRCAT -T 24 -s All_gene_concatenated.phy -o Ocomp1,Ocomp2,Ocomp25,Ocomp26,Ocomp27,Ocomp28,Ocomp29,Ocomp3,Ocomp5 -n All_gene_concatenated > Raxml.process 2>&1 &
# error
# because OG0061786 remove Acura3 after trim
# remove OG0061786

# remove orthologous genes (22 ortholgous genes will be removed)
# OG0061742 OG0061823 OG0061909 OG0061954 OG0061956 OG0062169 OG0062324 
# OG0062344 OG0062473 OG0062487 OG0062488 OG0062498 OG0062526 OG0062611 
# OG0062625 OG0062672 OG0062769 OG0062782 OG0062880 OG0062896 OG0062897
# OG0061786

# Kang@fishlab3 Tue Apr 11 16:00:05 /media/HDD/white_island/Compevo/genecapture/Total_orf/select
rm OG0061786-align.fas OG0061786.fas OG0061786-trim.fas
perl temp1.pl
scp All_gene_concatenated.phy kang1234@147.8.76.177:~/white_island/kraken/gene_capture_PNG/
# set all Ocomp individuals as root
# (base) kang1234@celia-PowerEdge-T640 Tue Apr 11 15:51:49 ~/white_island/kraken/gene_capture_PNG
nohup raxmlHPC -f a -x 12345 -p 12345 -# 100 -m GTRCAT -T 24 -s All_gene_concatenated.phy -o Ocomp1,Ocomp2,Ocomp25,Ocomp26,Ocomp27,Ocomp28,Ocomp29,Ocomp3,Ocomp5 -n All_gene_concatenated > Raxml.process 2>&1 &
# [1] 18512

# iTol to colour the label and node
# Kang@fishlab3 Tue Apr 11 23:03:07 /media/HDD/white_island/Compevo/genecapture/Total_orf/select
perl temp2.pl|less
```

```temp2.pl
#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

my $fa="OG0061615.fas";
my %col=(
	'Blenny'=>'#00468B',
	'Blueeyed'=>'#ED0000',
	'Common'=>'#0099B4',
	'Yaldwyn'=>'#925E9F',
	'Acura'=>'#FDAF91',
	'Apoly'=>'#ADB6B6',
	'Daru'=>'#1B19FF',
	'Ocomp'=>'#EFC000',
	'Padel'=>'#FF7F0E',
	'Pmol'=>'#E377C2',
	);
open FA, $fa or die "can not open $fa\n";
while (<FA>) {
	chomp;
	if (/>/) {
		s/>//;
		(my $spe)=$_=~/(\D+)\d+/;
		my $color=$col{$spe};
		print "$_,1,10,$color,1,1\n";
	}
}
```

```bash
# Kang@fishlab3 Wed Apr 12 00:01:51 /media/HDD/white_island/Compevo/genecapture/Total_orf/select
scp kang1234@147.8.76.177:~/white_island/kraken/gene_capture_PNG/*All_gene_concatenated ./
cp /media/HDD/cleaner_fish/genome/gene_family/cafetutorial_prep_r8s.py ./
# Kang@fishlab3 Wed Apr 12 00:27:01 /media/HDD/white_island/Compevo/genecapture/Total_orf/select
python cafetutorial_prep_r8s.py -i RAxML_bestTree.All_gene_concatenated -o r8s_ctl_file_1.txt -s 638853 -p 'Ocomp1,Ocomp2,Ocomp25,Ocomp26,Ocomp27,Ocomp28,Ocomp29,Ocomp3,Ocomp5,Daru1,Daru2,Daru25,Daru26,Daru27,Daru29,Daru30,Daru4,Daru6' -c '108.78'
nohup r8s -b -f r8s_ctl_file_1.txt > r8s_tmp_1.txt 2>&1 &
# [1] 9185
tail -n 1 r8s_tmp_1.txt | cut -c 16- > r8s_ultrametric.txt
# "r8s_ultrametric.txt" will be the input phylogeny of bamm
```

## Orth_ten_CO2: Divid the orthogroups into small orthogroups
```bash
# Kang@fishlab3 Thu Apr 06 09:44:45 /media/HDD/white_island/Compevo/Orth_ten_CO2
mkdir OrthoFinder; cd OrthoFinder
# Need: 1. "Gene_Trees/"; 2. "Orthogroup_Sequences/"; 3. "Orthogroups/"
# Kang@fishlab3 Thu Apr 06 09:50:31 /media/HDD/white_island/Compevo/Orth_ten_CO2/OrthoFinder
scp -r kang1234@147.8.76.177:~/white_island/Compevo/Orth_ten_CO2/OrthoFinder/Results_Apr04/Gene_Trees/ ./
scp -r kang1234@147.8.76.177:~/white_island/Compevo/Orth_ten_CO2/OrthoFinder/Results_Apr04/Orthogroup_Sequences/ ./
scp -r kang1234@147.8.76.177:~/white_island/Compevo/Orth_ten_CO2/OrthoFinder/Results_Apr04/Orthogroups/ ./

# Kang@fishlab3 Thu Apr 06 10:06:15 /media/HDD/white_island/Compevo/Orth_ten_CO2/OrthoFinder
vi build_sub_orth.pl
conda activate possvm
export DISPLAY=:0.0
# (possvm) Kang@fishlab3 Thu Apr 06 13:47:01 /media/HDD/white_island/Compevo/Orth_ten_CO2/OrthoFinder
nohup perl build_sub_orth.pl > build_sub_orth.process 2>&1 &
# [1] 22919

# Construct the phylogenetic tree of the ten fish species
# (base) kang1234@celia-PowerEdge-T640 Thu Apr 06 23:05:21 ~/white_island/Compevo/Orth_ten_CO2/OrthoFinder/Results_Apr04
mkdir single_copy; cd single_copy
# (base) kang1234@celia-PowerEdge-T640 Thu Apr 06 23:05:33 ~/white_island/Compevo/Orth_ten_CO2/OrthoFinder/Results_Apr04/single_copy
cp ../Single_Copy_Orthologue_Sequences/*.fa ./
# (base) kang1234@celia-PowerEdge-T640 Thu Apr 06 23:11:49 ~/white_island/Compevo/Orth_ten_CO2/OrthoFinder/Results_Apr04
cp ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21/Align_trim_conca.pl ./
perl Align_trim_conca.pl > single_copy.concatenated.fasta
fasta2phy.pl single_copy.concatenated.fasta >single_copy.concatenated.phy
nohup raxmlHPC -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s single_copy.concatenated.phy -o Ocomp -n single_copy.concatenated -T 24 > raxml.process 2>&1 &
# [1] 1161

# Kang@fishlab3 Tue Apr 18 16:04:36 /media/HDD/white_island/Compevo/Orth_ten_CO2/OrthoFinder/sub_orth
cp /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input/extract_ano.pl ./
cp /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input/unprot_name_description.txt ./
perl extract_ano.pl >ortho10_ano.txt
```

## Redo the analysis paml and bamm
```bash
####################################################################################################
# 1. make sure all species have at least one transcript and only kept the genes with same gene name;
# 2. and select the longest one if there are more than one transcript for a species
####################################################################################################

# Step 1: 
# (base) kang1234@celia-PowerEdge-T640 Tue May 09 14:34:46 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21/Orthogroups
perl select_same_gene_nm.pl|perl -alne '$nb=@F;print if $nb==18' >Orth_same_gene_nm_16spe.txt

# Step 2
# (base) kang1234@celia-PowerEdge-T640 Tue May 09 16:03:09 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21/Orthogroups
perl select_longest.pl > Orth_same_gene_nm_16spe_longest.txt

# Step 3: recode the orthologous id and build the ref sequences to get the number of mapped reads
# (base) kang1234@celia-PowerEdge-T640 Tue May 09 16:54:54 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21/Orthogroups
perl record_build_ref.pl

# Build the matrix of mapped reads number and do PAML again
# (base) kang1234@celia-PowerEdge-T640 Wed May 10 10:44:21 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21/Orthogroups
perl change_name.pl > orth_id_paml.txt
# (base) kang1234@celia-PowerEdge-T640 Mon May 15 17:11:17 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21/Orthogroups
perl temp1.pl > orth_id_paml_ano.txt # annotation to each orthologous

# Kang@fishlab3 Wed May 10 10:48:38 /media/HDD/white_island/Compevo
mkdir orth16_new
# Kang@fishlab3 Wed May 10 10:49:52 /media/HDD/white_island/Compevo/orth16_new
scp kang1234@147.8.76.177:~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21/Orthogroups/orth_id_paml.txt ./
cp /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input/*.fasta ./
cp /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input/prepare_input_paml* ./
cp /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input/correlation.txt ./

nohup perl prepare_input_paml_parallel.pl orth_id_paml.txt >prepare_input_paml.process 2>&1 &
# [1] 11368
```

**Transfer to HPC to prepare paml input files**   
```bash
# jlkang@hpc2021 Wed May 10 21:53:46 /lustre1/g/sbs_schunter/Kang
mkdir orth16_new

# Kang@fishlab3 Wed May 10 21:56:34 /media/HDD/white_island/Compevo/orth16_new
scp *.fasta jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new
scp *.txt jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new
scp *.pl jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new

module load perl/5.34.0 perl-lib/5.34.0
# "orth_id_paml.txt" split into six files
split -l 1473 orth_id_paml.txt orth_id_paml_split_
sbatch script_aa.cmd # Submitted batch job 1163359: amd
sbatch script_ab.cmd # Submitted batch job 1163376: amd
sbatch script_ac.cmd # Submitted batch job 1163375: intel
sbatch script_ad.cmd # Submitted batch job 1163380: intel
sbatch script_ae.cmd # Submitted batch job 1163378: intel
sbatch script_af.cmd # Submitted batch job 1163381: amd

# check some orthogroups not in the final list
split -l 603 left_paml.txt orth_id_paml_split_

########################################
# free ratio estimate evolutionary rate
########################################
# 1. concatenate the sequences to estimate the evolutionary rate
# jlkang@hpc2021 Thu May 11 09:32:12 /lustre1/g/sbs_schunter/Kang/orth16_new
perl conca_orth_paml.pl >Orth_conca_paml.fasta
# Kang@fishlab3 Thu May 11 19:15:02 /media/HDD/white_island/Compevo/orth16_new
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/Orth_conca_paml.fasta ./
scp -r kang1234@147.8.76.177:~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21/Orth_conca ./
rm 2NG.dN 2NG.dS 2NG.t 4fold.nuc free-ratio-result.txt lnf Orth_conca_paml.phy rst rst1 rub
mv ../Orth_conca_paml.phy ./
nohup codeml free-ratio.ctr > free-ratio.process 2>&1 &
# [1] 20354

# Kang@fishlab3 Wed Jun 21 13:42:34 /media/HDD/white_island/Compevo/orth16_new/Orth_conca
less free-ratio-result.txt|grep '^#'  # species number
# copy the estimation results to "1.txt"
# (base) kang1234@celia-PowerEdge-T640 Wed Jun 21 13:37:55 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21/Orth_conca
less 1.txt|perl -alne 'if (/\d+/ || /\D+/){print}'|perl -alne 's/^\s+//;s/\s+/\t/g;print' # copy to excel and edit
###################################################################################################
# branch  t       N       S       dN/dS   dN      dS      N*dN    S*dS
# 17..18  0.739   5124968.1       1915572.9       0.1263  0.0855  0.6770  438266.4        1296777.0
# 18..19  0.428   5124968.1       1915572.9       0.1249  0.0491  0.3932  251685.8        753247.6
# 19..9   0.332   5124968.1       1915572.9       0.0768  0.0259  0.3376  132912.5        646709.6
###################################################################################################

# 2. each gene
# Kang@fishlab3 Sat May 13 18:48:48 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
scp codeml.pl codeml_parallel.pl spe.tre jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new
split -l 1955 final_orth_input_paml.txt final_orth_input_paml_split_
# perl codeml_parallel.pl final_orth_input_paml_split_aa
# jlkang@hpc2021 Sat May 13 21:43:29 /lustre1/g/sbs_schunter/Kang/orth16_new
module load perl/5.34.0 perl-lib/5.34.0
sbatch script_aa.cmd # Submitted batch job 1173142
sbatch script_ab.cmd # Submitted batch job 1173143
sbatch script_ac.cmd # Submitted batch job 1173144
sbatch script_ad.cmd # Submitted batch job 1173145

split -l 800 left_free_ratio.txt final_orth_input_paml_split_

# jlkang@hpc2021 Mon May 15 16:14:36 /lustre1/g/sbs_schunter/Kang/orth16_new
cat free-ratio/*.txt >free_ratio_result.txt
# kangjingliang@kangjingliangdeMacBook-Pro 一  5 15 16:29:42 ~/Documents/2023/WI/paml
mkdir free-ratio_new
# kangjingliang@kangjingliangdeMacBook-Pro 一  5 15 16:30:52 ~/Documents/2023/WI/paml/free-ratio_new
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/free_ratio_result.txt ./
less free_ratio_result.txt|perl -alne 'if (/^Or/i){print}else{s/\s+/\t/g;print}' >free_ratio_result.txt.1
mv free_ratio_result.txt.1 free_ratio_result.txt
less free_ratio_result.txt|perl -alne 'if (/^Or/i){print}elsif($F[1]=~/[a-z]+/){print}' >free_ratio_result.txt.1
mv free_ratio_result.txt.1 free_ratio_result.txt
perl temp1.pl > free_ratio_result_plot.txt

# maximum dN/dS per gene
# kangjingliang@kangjingliangdeMacBook-Pro 一  5 15 17:15:06 ~/Documents/2023/WI/paml/free-ratio_new
perl select_max_dNdS.pl > max_dNdS_per_orth.txt
########################################
# Positive selected genes
########################################
# Common
# jlkang@hpc2021 Mon May 15 15:01:13 /lustre1/g/sbs_schunter/Kang/orth16_new
vi spe_Common.tre
# (((Ocomp,((Stickleback,Fugu),(((Daru,((Acura,Apoly),(Pmol,Padel))),(Blenny,((Yaldwyn,Blueeyed),Common #1))),(Platyfish,Medaka)))),Zebrafish),Spottedgar);
vi codeml_parallel.pl
# perl codeml.pl --input temp/$temp --model branch-site --dir paml_files --output_suf Common --tree spe_Common.tre --icode 0 --omega 1.2
split -l 1183 final_orth_input_paml.txt final_orth_input_paml_split_
# jlkang@hpc2021 Mon May 15 16:10:50 /lustre1/g/sbs_schunter/Kang/orth16_new
sbatch script_aa.cmd # Submitted batch job 1175377
sbatch script_ab.cmd # Submitted batch job 1175392
sbatch script_ac.cmd # Submitted batch job 1175394
sbatch script_ad.cmd # Submitted batch job 1175395
sbatch script_ae.cmd # Submitted batch job 1175396
sbatch script_af.cmd # Submitted batch job 1175397

########
# Blenny
########
# jlkang@hpc2021 Tue May 16 14:53:52 /lustre1/g/sbs_schunter/Kang/orth16_new
vi spe_Blenny.tre
# (((Ocomp,((Stickleback,Fugu),(((Daru,((Acura,Apoly),(Pmol,Padel))),(Blenny #1,((Yaldwyn,Blueeyed),Common))),(Platyfish,Medaka)))),Zebrafish),Spottedgar);
vi codeml_parallel.pl
# perl codeml.pl --input temp/$temp --model branch-site --dir paml_files --output_suf Blenny --tree spe_Blenny.tre --icode 0 --omega 1.2
sbatch script_aa.cmd
sbatch script_ab.cmd
sbatch script_ac.cmd
sbatch script_ad.cmd
sbatch script_ae.cmd
sbatch script_af.cmd

###########
# Blue-eyed
###########
# jlkang@hpc2021 Wed May 17 16:11:03 /lustre1/g/sbs_schunter/Kang/orth16_new
vi spe_Blueeyed.tre
# (((Ocomp,((Stickleback,Fugu),(((Daru,((Acura,Apoly),(Pmol,Padel))),(Blenny,((Yaldwyn,Blueeyed #1),Common))),(Platyfish,Medaka)))),Zebrafish),Spottedgar);
vi codeml_parallel.pl
# perl codeml.pl --input temp/$temp --model branch-site --dir paml_files --output_suf Blueeyed --tree spe_Blueeyed.tre --icode 0 --omega 1.2
sbatch script_aa.cmd
sbatch script_ab.cmd
sbatch script_ac.cmd
sbatch script_ad.cmd
sbatch script_ae.cmd
sbatch script_af.cmd

###########
# Yaldwyn
###########
# jlkang@hpc2021 Thu May 18 15:01:39 /lustre1/g/sbs_schunter/Kang/orth16_new
vi spe_Yaldwyn.tre
# (((Ocomp,((Stickleback,Fugu),(((Daru,((Acura,Apoly),(Pmol,Padel))),(Blenny,((Yaldwyn #1,Blueeyed),Common))),(Platyfish,Medaka)))),Zebrafish),Spottedgar);
vi codeml_parallel.pl
# perl codeml.pl --input temp/$temp --model branch-site --dir paml_files --output_suf Yaldwyn --tree spe_Yaldwyn.tre --icode 0 --omega 1.2
sbatch script_aa.cmd
sbatch script_ab.cmd
sbatch script_ac.cmd
sbatch script_ad.cmd
sbatch script_ae.cmd
sbatch script_af.cmd

###########################################
# HYPHY to detect positively selected genes
###########################################
# Blenny
# perl run_hyphy.pl final_orth_input_paml_split_aa spe_Blenny_hyphy.tre Blenny
# jlkang@hpc2021 Tue Jun 20 15:31:03 /lustre1/g/sbs_schunter/Kang/orth16_new
module load hyphy; module load perl/5.34.0 perl-lib/5.34.0
sbatch script_aa.cmd; # ... ... Blenny

# keep in mind: change the name in the result
# jlkang@hpc2021 Mon Jun 26 09:57:01 /lustre1/g/sbs_schunter/Kang/orth16_new
perl temp5.pl Blenny  # final_alignment.fa.BUSTED.json => Blenny_BUSTED.json

# mv final_alignment.fa.BUSTED.json final_alignment.fa.Blenny.BUSTED.json
# jlkang@hpc2021 Mon Jun 26 09:27:00 /lustre1/g/sbs_schunter/Kang/orth16_new/Hyphy
cp /lustre1/g/sbs_schunter/Kang/Ldim_revision/Hyphy/temp1.pl ./
vi temp1.pl # 
perl temp1.pl Blenny > Blenny_hyphy_PSGs.txt # my ($orth, $spe)=$txt=~/(.*)\_(.*)\.txt/;
module load R

R
#############
# fdr correct
p_apoly<-read.table(file="Blenny_hyphy_PSGs.txt")
p_apoly$fdr<- p.adjust(p_apoly$V3,method="fdr",length(p_apoly$V3))
write.table(p_apoly, file="Blenny_hyphy_PSGs_fdr.txt",row.names=F,col.names=F,quote=F,sep="\t")
#############
# jlkang@hpc2021 Mon Jun 26 09:39:16 /lustre1/g/sbs_schunter/Kang/orth16_new/Hyphy
mv Blenny_hyphy_PSGs_fdr.txt ../Hyphy_PSGs/
# jlkang@hpc2021 Mon Jun 26 09:41:07 /lustre1/g/sbs_schunter/Kang/orth16_new/Hyphy_PSGs
less Blenny_hyphy_PSGs_fdr.txt|perl -alne 'print if $F[-1]<=0.05'|wc -l # 103 PSGs in Blenny; 99 are consistent with paml

# Common
# perl run_hyphy.pl final_orth_input_paml_split_aa spe_Common_hyphy.tre Common
# jlkang@hpc2021 Mon Jun 26 10:06:59 /lustre1/g/sbs_schunter/Kang/orth16_new
sbatch script_aa.cmd # ... ... Common
# jlkang@hpc2021 Thu Jun 29 10:10:42 /lustre1/g/sbs_schunter/Kang/orth16_new/Hyphy
perl temp1.pl Common > Common_hyphy_PSGs.txt
module load R
R
#############
# fdr correct
p_apoly<-read.table(file="Common_hyphy_PSGs.txt")
p_apoly$fdr<- p.adjust(p_apoly$V3,method="fdr",length(p_apoly$V3))
write.table(p_apoly, file="Common_hyphy_PSGs_fdr.txt",row.names=F,col.names=F,quote=F,sep="\t")
#############
mv Common_hyphy_PSGs_fdr.txt ../Hyphy_PSGs/
# jlkang@hpc2021 Thu Jun 29 10:16:38 /lustre1/g/sbs_schunter/Kang/orth16_new/Hyphy_PSGs
less Common_hyphy_PSGs_fdr.txt|perl -alne 'print if $F[-1]<=0.05'|wc -l # 236 PSGs in Common; 217 are consistent with paml

# change the resulted file name
# jlkang@hpc2021 Mon Jun 26 09:57:01 /lustre1/g/sbs_schunter/Kang/orth16_new
perl temp5.pl Common  # final_alignment.fa.BUSTED.json => Common_BUSTED.json

# Blueeyed
# perl run_hyphy.pl final_orth_input_paml_split_aa spe_Blueeyed_hyphy.tre Blueeyed
# divide to 10 files to run
split -l 710 final_orth_input_paml.txt final_orth_input_paml_split_
sbatch script_aa.cmd
sbatch script_ab.cmd
sbatch script_ac.cmd
sbatch script_ad.cmd
sbatch script_ae.cmd
sbatch script_af.cmd
sbatch script_ag.cmd
sbatch script_ah.cmd
sbatch script_ai.cmd
sbatch script_aj.cmd

# the process ag is wrong
split -l 142 final_orth_input_paml_split_ag final_orth_input_paml_split_
sbatch script_aa.cmd
sbatch script_ab.cmd
sbatch script_ac.cmd
sbatch script_ad.cmd
sbatch script_ae.cmd

# jlkang@hpc2021 Sat Jul 01 10:14:55 /lustre1/g/sbs_schunter/Kang/orth16_new/Hyphy
perl temp1.pl Blueeyed > Blueeyed_hyphy_PSGs.txt
module load R
R
#############
# fdr correct
p_apoly<-read.table(file="Blueeyed_hyphy_PSGs.txt")
p_apoly$fdr<- p.adjust(p_apoly$V3,method="fdr",length(p_apoly$V3))
write.table(p_apoly, file="Blueeyed_hyphy_PSGs_fdr.txt",row.names=F,col.names=F,quote=F,sep="\t")
#############
mv Blueeyed_hyphy_PSGs_fdr.txt ../Hyphy_PSGs/
less Blueeyed_hyphy_PSGs_fdr.txt|perl -alne 'print if $F[-1]<=0.05'|wc -l # 414 PSGs by hyhpy; 349 are consistent with paml
# change the resulted file name
# jlkang@hpc2021 Mon Jun 26 09:57:01 /lustre1/g/sbs_schunter/Kang/orth16_new
perl temp5.pl Blueeyed  # final_alignment.fa.BUSTED.json => Blueeyed_BUSTED.json

# Yaldwyn
# perl run_hyphy.pl final_orth_input_paml_split_aa spe_Yaldwyn_hyphy.tre Yaldwyn
sbatch script_aa.cmd ... ...

# jlkang@hpc2021 Tue Jul 04 02:41:06 /lustre1/g/sbs_schunter/Kang/orth16_new/Hyphy
perl temp1.pl Yaldwyn > Yaldwyn_hyphy_PSGs.txt
#############
# fdr correct
p_apoly<-read.table(file="Yaldwyn_hyphy_PSGs.txt")
p_apoly$fdr<- p.adjust(p_apoly$V3,method="fdr",length(p_apoly$V3))
write.table(p_apoly, file="Yaldwyn_hyphy_PSGs_fdr.txt",row.names=F,col.names=F,quote=F,sep="\t")
#############
mv Yaldwyn_hyphy_PSGs_fdr.txt ../Hyphy_PSGs/
# jlkang@hpc2021 Tue Jul 04 02:43:27 /lustre1/g/sbs_schunter/Kang/orth16_new/Hyphy_PSGs
less Yaldwyn_hyphy_PSGs_fdr.txt|perl -alne 'print if $F[-1]<=0.05'|wc -l # 324 PSGs by hyphy; 269 are consistent with paml
# change the resulted file name
# jlkang@hpc2021 Mon Jun 26 09:57:01 /lustre1/g/sbs_schunter/Kang/orth16_new
perl temp5.pl Yaldwyn  # final_alignment.fa.BUSTED.json => Yaldwyn_BUSTED.json

# only the PSGs detected by both hyphy and paml were considered as the final PSGs
# Kang@fishlab3 Tue Jul 04 02:51:50 /media/HDD/white_island/Compevo/orth16_new
mkdir PSGs_paml_hyphy; mv *_distinc_info.txt PSGs_paml_hyphy/
# Kang@fishlab3 Tue Jul 04 03:29:05 /media/HDD/white_island/Compevo/orth16_new/PSGs_paml_hyphy
perl temp1.pl Blenny_hyphy_PSGs_fdr.txt Blenny_PSGs_distinc_info.txt >Blenny_PSGs_final.txt
perl temp1.pl Blueeyed_hyphy_PSGs_fdr.txt Blueeyed_PSGs_distinc_info.txt >Blueeyed_PSGs_final.txt
perl temp1.pl Common_hyphy_PSGs_fdr.txt Common_PSGs_distinc_info.txt >Common_PSGs_final.txt
perl temp1.pl Yaldwyn_hyphy_PSGs_fdr.txt Yaldwyn_PSGs_distinc_info.txt >Yaldwyn_PSGs_final.txt
```


**Map to get the reads number matrix**   
```temp2.pl
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
        unless ($a[2] eq "Cs") {
                $hash2{$id}={
                        'IND' => $ind,
                        'SPE' => $spe
                };
        }
}

my @fq=<*.fastq.gz>;
my @cmds;
foreach my $fq (@fq) {
        my ($id, $info)=$fq=~/(.*)_(\d+\.fastq\.gz)/;
        if ($hash2{$id}) {
                my $spe=$hash2{$id}->{'SPE'};
                my $ind=$hash2{$id}->{'IND'};

                my $new=$ind."_".$info;
                my $cmd="mv $fq $new";
                system($cmd);
        }
}
```

```bash
# Use the whole nucleotide per species as the reference to get the reads number matrix
# Do not select the individuals from Cs and change the name of fastq

# (base) kang1234@celia-PowerEdge-T640 Thu May 11 11:31:21 ~/white_island/kraken
perl temp2.pl
# jlkang@hpc2021 Thu May 11 11:46:21 /lustre1/g/sbs_schunter/Kang/orth16_new
mkdir fastqs

# Transfer the fastq data to HPC from SNORLAX
# WI
# (base) kang1234@celia-PowerEdge-T640 Thu May 11 14:12:31 ~/white_island/kraken
nohup scp Blenny*_1.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Blenny_1.out 2>&1
nohup scp Blenny*_2.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Blenny_2.out 2>&1
nohup scp Blueeyed*_1.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Blueeyed_1.out 2>&1
nohup scp Blueeyed*_2.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Blueeyed_2.out 2>&1
nohup scp Common*_1.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Common_1.out 2>&1
nohup scp Common*_2.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Common_2.out 2>&1
nohup scp Yaldwyn*_1.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Yaldwyn_1.out 2>&1
nohup scp Yaldwyn*_2.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Yaldwyn_2.out 2>&1

# Transfer the fastq data to HPC from my own workstation
# PNG
# Kang@fishlab3 Thu May 11 14:18:16 ~/Desktop/PapueNewGuinea-new/merge_clean
nohup scp Apoly21*.fastq.gz Apoly23*.fastq.gz Apoly24*.fastq.gz Apoly26*.fastq.gz Apoly27*.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Apoly_1.out 2>&1
nohup scp Apoly44*.fastq.gz Apoly45*.fastq.gz Apoly46*.fastq.gz Apoly47*.fastq.gz Apoly48*.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Apoly_2.out 2>&1
nohup scp Acura3*.fastq.gz Acura4*.fastq.gz Acura5*.fastq.gz Acura6*.fastq.gz Acura7*.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Acura_1.out 2>&1
nohup scp Acura24*.fastq.gz Acura25*.fastq.gz Acura26*.fastq.gz Acura27*.fastq.gz Acura28*.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Acura_2.out 2>&1
nohup scp Ocomp1*.fastq.gz Ocomp2*.fastq.gz Ocomp3*.fastq.gz Ocomp4*.fastq.gz Ocomp5*.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Ocomp_1.out 2>&1
nohup scp Ocomp25*.fastq.gz Ocomp26*.fastq.gz Ocomp27*.fastq.gz Ocomp28*.fastq.gz Ocomp29*.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Ocomp_2.out 2>&1
nohup scp Daru1*.fastq.gz Daru2*.fastq.gz Daru3*.fastq.gz Daru4*.fastq.gz Daru6*.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Daru_1.out 2>&1
nohup scp Daru25*.fastq.gz Daru26*.fastq.gz Daru27*.fastq.gz Daru29*.fastq.gz Daru30*.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Daru_2.out 2>&1
nohup scp Padel6*.fastq.gz Padel7*.fastq.gz Padel8*.fastq.gz Padel9*.fastq.gz Padel10*.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Padel_1.out 2>&1
nohup scp Padel23*.fastq.gz Padel24*.fastq.gz Padel25*.fastq.gz Padel26*.fastq.gz Padel28*.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Padel_2.out 2>&1
nohup scp Pmol7*.fastq.gz Pmol8*.fastq.gz Pmol9*.fastq.gz Pmol10*.fastq.gz Pmol11*.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Pmol_1.out 2>&1
nohup scp Pmol21*.fastq.gz Pmol22*.fastq.gz Pmol23*.fastq.gz Pmol24*.fastq.gz Pmol25*.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs > Pmol_2.out 2>&1

# Create the index
module load RSEM/1.3.3; module load bowtie2; module load bowtie

rsem-prepare-reference --bowtie2 Blenny.cds.fasta Blenny --bowtie2
rsem-prepare-reference --bowtie2 Blueeyed.cds.fasta Blueeyed --bowtie2
rsem-prepare-reference --bowtie2 Common.cds.fasta Common --bowtie2
rsem-prepare-reference --bowtie2 Yaldwyn.cds.fasta Yaldwyn --bowtie2
rsem-prepare-reference --bowtie2 Acura.cds.fasta Acura --bowtie2
rsem-prepare-reference --bowtie2 Apoly.cds.fasta Apoly --bowtie2
rsem-prepare-reference --bowtie2 Daru.cds.fasta Daru --bowtie2
rsem-prepare-reference --bowtie2 Ocomp.cds.fasta Ocomp --bowtie2
rsem-prepare-reference --bowtie2 Padel.cds.fasta Padel --bowtie2
rsem-prepare-reference --bowtie2 Pmol.cds.fasta Pmol --bowtie2

sbatch script.cmd # Submitted batch job 1164549

# Blenny
# jlkang@hpc2021 Thu May 11 22:36:29 /lustre1/g/sbs_schunter/Kang/orth16_new/fastqs
mkdir Blenny_RSEM_output
# jlkang@hpc2021 Thu May 11 23:09:40 /lustre1/g/sbs_schunter/Kang/orth16_new/fastqs/Blenny_RSEM_output
ll ../Blenny*_1.fastq.gz|perl -alne '($nm)=$F[-1]=~/\.\.\/(.*)/;$info.=$nm.",";END{$info=~s/\,$//;print "Blenny\t$info"}' > Blenny_info.txt
cp /lustre1/g/sbs_schunter/Kang/Orth10_reference/Blueeyed_RSEM_output/RSEM_align.pl ./
cp /lustre1/g/sbs_schunter/Kang/Orth10_reference/Blueeyed_RSEM_output/script*.cmd ./
perl RSEM_align.pl Blenny_info.txt # rsem_1.sh; rsem_2.sh
module load RSEM/1.3.3; module load bowtie2; module load bowtie
sbatch script1.cmd # Submitted batch job 1165739
sbatch script2.cmd # Submitted batch job 1165740

# Blueeyed
# jlkang@hpc2021 Thu May 11 23:40:55 /lustre1/g/sbs_schunter/Kang/orth16_new/fastqs
mkdir Blueeyed_RSEM_output
# jlkang@hpc2021 Thu May 11 23:41:22 /lustre1/g/sbs_schunter/Kang/orth16_new/fastqs/Blueeyed_RSEM_output
cp ../Blenny_RSEM_output/RSEM_align.pl ./
cp ../Blenny_RSEM_output/script*.cmd ./
ll ../Blueeyed*_1.fastq.gz|perl -alne '($nm)=$F[-1]=~/\.\.\/(.*)/;$info.=$nm.",";END{$info=~s/\,$//;print "Blueeyed\t$info"}' > Blueeyed_info.txt
perl RSEM_align.pl Blueeyed_info.txt
sbatch script1.cmd # Submitted batch job 1165770
sbatch script2.cmd # Submitted batch job 1165774

# Common
# jlkang@hpc2021 Fri May 12 00:06:25 /lustre1/g/sbs_schunter/Kang/orth16_new/fastqs
mkdir Common_RSEM_output
# jlkang@hpc2021 Thu May 11 23:56:12 /lustre1/g/sbs_schunter/Kang/orth16_new/fastqs/Common_RSEM_output
cp ../Blenny_RSEM_output/RSEM_align.pl ./
cp ../Blenny_RSEM_output/script*.cmd ./
ll ../Common*_1.fastq.gz|perl -alne '($nm)=$F[-1]=~/\.\.\/(.*)/;$info.=$nm.",";END{$info=~s/\,$//;print "Common\t$info"}' > Common_info.txt
perl RSEM_align.pl Common_info.txt
sbatch script1.cmd # Submitted batch job 1165782
sbatch script2.cmd # Submitted batch job 1165786

# Yaldwyn
# jlkang@hpc2021 Fri May 12 00:06:25 /lustre1/g/sbs_schunter/Kang/orth16_new/fastqs
mkdir Yaldwyn_RSEM_output
cp ../Blenny_RSEM_output/RSEM_align.pl ./
cp ../Blenny_RSEM_output/script*.cmd ./
ll ../Yaldwyn*_1.fastq.gz|perl -alne '($nm)=$F[-1]=~/\.\.\/(.*)/;$info.=$nm.",";END{$info=~s/\,$//;print "Yaldwyn\t$info"}' > Yaldwyn_info.txt
perl RSEM_align.pl Yaldwyn_info.txt
sbatch script1.cmd # Submitted batch job 1165790
sbatch script2.cmd # Submitted batch job 1165791

# Acura: in my own workstation
# Kang@fishlab3 Fri May 12 09:12:14 ~/Desktop/PapueNewGuinea-new/merge_clean
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/Acura.* ./
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/Apoly.* ./
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/Ocomp.* ./
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/Daru.* ./
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/Pmol.* ./
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/Padel.* ./
mkdir RSEM_output
# jlkang@hpc2021 Fri May 12 09:21:28 /lustre1/g/sbs_schunter/Kang/orth16_new/fastqs/Blenny_RSEM_output
ll ../Acura*_1.fastq.gz|perl -alne '($nm)=$F[-1]=~/\.\.\/(.*)/;$info.=$nm.",";END{$info=~s/\,$//;print "Acura\t$info"}' > Acura_info.txt
ll ../Apoly*_1.fastq.gz|perl -alne '($nm)=$F[-1]=~/\.\.\/(.*)/;$info.=$nm.",";END{$info=~s/\,$//;print "Apoly\t$info"}' > Apoly_info.txt
ll ../Ocomp*_1.fastq.gz|perl -alne '($nm)=$F[-1]=~/\.\.\/(.*)/;$info.=$nm.",";END{$info=~s/\,$//;print "Ocomp\t$info"}' > Ocomp_info.txt
ll ../Daru*_1.fastq.gz|perl -alne '($nm)=$F[-1]=~/\.\.\/(.*)/;$info.=$nm.",";END{$info=~s/\,$//;print "Daru\t$info"}' > Daru_info.txt
ll ../Pmol*_1.fastq.gz|perl -alne '($nm)=$F[-1]=~/\.\.\/(.*)/;$info.=$nm.",";END{$info=~s/\,$//;print "Pmol\t$info"}' > Pmol_info.txt
ll ../Padel*_1.fastq.gz|perl -alne '($nm)=$F[-1]=~/\.\.\/(.*)/;$info.=$nm.",";END{$info=~s/\,$//;print "Padel\t$info"}' > Padel_info.txt
# # Kang@fishlab3 Fri May 12 09:28:34 ~/Desktop/PapueNewGuinea-new/merge_clean/RSEM_output
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs/Blenny_RSEM_output/*_info.txt ./
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs/Blenny_RSEM_output/RSEM_align.pl ./
perl RSEM_align.pl Acura_info.txt > rsem.sh
nohup sh rsem.sh > rsem.process 2>&1 &
# [1] 2231

perl RSEM_align.pl Apoly_info.txt > rsem.sh
nohup sh rsem.sh > rsem.process 2>&1 &
# [1] 30918

perl RSEM_align.pl Ocomp_info.txt > rsem.sh
perl RSEM_align.pl Daru_info.txt >> rsem.sh
perl RSEM_align.pl Pmol_info.txt >> rsem.sh
perl RSEM_align.pl Padel_info.txt >> rsem.sh
# [1] 26870

# Creat the reads number matrix
# Kang@fishlab3 Tue May 23 10:07:17 ~/Desktop/PapueNewGuinea-new/merge_clean/RSEM_output
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs/Blenny_RSEM_output/*.genes.results ./
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs/Blueeyed_RSEM_output/*.genes.results ./
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs/Common_RSEM_output/*.genes.results ./
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/fastqs/Yaldwyn_RSEM_output/*.genes.results ./
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/orth_id_paml.txt ./
mkdir orth_genes_paml
# Kang@fishlab3 Tue May 23 10:46:00 ~/Desktop/PapueNewGuinea-new/merge_clean/RSEM_output
perl temp1.pl # the results in orth_genes_paml/
# Kang@fishlab3 Tue May 23 11:24:44 ~/Desktop/PapueNewGuinea-new/merge_clean/RSEM_output/orth_genes_paml
RNAnorm --method TPM --norm --prefix paml_Orth16 *genes.results
mkdir bamm
cp ../paml_Orth16.TPM.TMM.sqrt.matrix ./
cp /media/HDD/white_island/Compevo/genecapture/Total_orf/select/r8s_ultrametric.txt ./
cp ~/Desktop/PapueNewGuinea-new/merge_clean/gene_capture/result/Total_orf/Final_ind_2/bamm/temp1.pl ./
# CLOCK: OG0000954
perl temp1.pl OG0000954 > OG0000954.txt
mkdir OG0000954_trait; mv OG0000954.txt OG0000954_trait/
cd OG0000954_trait/
cp ../r8s_ultrametric.txt ./
bamm -c traitcontrol.txt
# kangjingliang@kangjingliangdeMacBook-Pro 二  5 23 11:45:38 ~/Documents/2023/WI/paml/free-ratio_new/bamm
scp -r Kang@147.8.76.229:~/Desktop/PapueNewGuinea-new/merge_clean/RSEM_output/orth_genes_paml/bamm/OG0000954_trait ./

# Select the core CR, vision, and pH genes to test
# Kang@fishlab3 Tue May 23 11:45:19 ~/Desktop/PapueNewGuinea-new/merge_clean/RSEM_output/orth_genes_paml/bamm
cp /media/HDD/white_island/Compevo/Orth_ten_CO2/bamm/temp2.pl ./
# BMAL1: OG0002046
perl temp2.pl OG0002046
cd OG0002046_trait
# Kang@fishlab3 Tue May 23 11:58:56 ~/Desktop/PapueNewGuinea-new/merge_clean/RSEM_output/orth_genes_paml/bamm/OG0002046_trait
bamm -c traitcontrol.txt

# Select the target genes
# kangjingliang@kangjingliangdeMacBook-Pro 二  5 23 12:34:23 ~/Documents/2023/WI/paml/free-ratio_new
scp orth_id_paml_ano.txt Kang@147.8.76.229:~/Desktop/PapueNewGuinea-new/merge_clean/RSEM_output/orth_genes_paml/bamm
# revise temp2.pl
# Kang@fishlab3 Tue May 23 12:47:32 ~/Desktop/PapueNewGuinea-new/merge_clean/RSEM_output/orth_genes_paml/bamm
perl temp3.pl # create "orth_id_paml_ano.txt" && "run_bamm.sh"
sh run_bamm.sh

# use the previous result? just white island individuals
# Kang@fishlab3 Tue May 23 16:02:29 /media/HDD/white_island/Compevo/Orth_ten_CO2/bamm_WI
cp ~/Desktop/PapueNewGuinea-new/merge_clean/RSEM_output/orth_genes_paml/bamm/bamm_genes_plot.txt ./
cp ~/Desktop/PapueNewGuinea-new/merge_clean/RSEM_output/orth_genes_paml/bamm/temp3.pl ./
# kangjingliang@kangjingliangdeMacBook-Pro 二  5 23 16:05:44 ~/Documents/2023/WI/DEGs/Enrichment
scp unprot_name_description_orthgroup.txt Kang@147.8.76.229:/media/HDD/white_island/Compevo/Orth_ten_CO2/bamm_WI
# Kang@fishlab3 Tue May 23 16:09:48 /media/HDD/white_island/Compevo/Orth_ten_CO2/bamm_WI
perl temp3.pl
# Kang@fishlab3 Tue May 23 16:23:07 /media/HDD/white_island/Compevo/Orth_ten_CO2/bamm_WI
nohup sh run_bamm.sh > run_bamm.process 2>&1 &
# [1] 18125

# Test the stress response genes
# OG0045765: HS90A; OG0085246: FOS; OG0068475: JUN; OG0083480: JUN
# Kang@fishlab3 Fri Jun 23 11:54:31 /media/HDD/white_island/Compevo/Orth_ten_CO2/bamm_WI
vi bamm_genes_plot.txt # put the gene name inside
perl temp3.pl
cat bamm_genes_plot_ano.txt
# OG0007248       sp|P53450|FOS_TAKRU     Proto-oncogene c-Fos    FOS
# OG0044371       sp|P12981|JUN_COTJA     Transcription factor AP-1       JUN
# OG0045765       sp|Q4R4P1|HS90A_MACFA   Heat shock protein HSP 90-alpha HS90A # DE
# OG0050100       sp|P54864|JUN_SERCA     Transcription factor AP-1       JUN
# OG0068475       sp|P18870|JUN_CHICK     Transcription factor AP-1       JUN # DE
# OG0083480       sp|P54864|JUN_SERCA     Transcription factor AP-1       JUN # DE
# OG0085246       sp|P79702|FOS_CYPCA     Proto-oncogene c-Fos    FOS # DE
sh run_bamm.sh
# kangjingliang@kangjingliangdeMacBook-Pro 五  6 23 12:01:08 ~/Documents/2023/WI/bamm_WI
scp -r Kang@147.8.76.229:/media/HDD/white_island/Compevo/Orth_ten_CO2/bamm_WI/OG0045765_trait ./
scp -r Kang@147.8.76.229:/media/HDD/white_island/Compevo/Orth_ten_CO2/bamm_WI/OG0068475_trait ./
scp -r Kang@147.8.76.229:/media/HDD/white_island/Compevo/Orth_ten_CO2/bamm_WI/OG0083480_trait ./
scp -r Kang@147.8.76.229:/media/HDD/white_island/Compevo/Orth_ten_CO2/bamm_WI/OG0085246_trait ./
```
## Annotation
```
# Annotate the new orthologous genes
# (base) kang1234@celia-PowerEdge-T640 Mon May 22 12:38:18 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21/Orthogroups
perl temp3.pl > final_orth_input_paml.fa
# kangjingliang@kangjingliangdeMacBook-Pro 一  5 22 12:41:29 ~/Documents/2023/WI/paml/free-ratio_new/annotate
scp kang1234@147.8.76.177:~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21/Orthogroups/final_orth_input_paml.fa ./

# 2.txt is the functional table
# kangjingliang@kangjingliangdeMacBook-Pro 一  5 22 15:49:11 ~/Documents/2023/WI/paml/free-ratio_new
less 2.txt|perl -alne 'my @F=split /\t/;my $nb=@F;if ($nb==12){my $info;for (my $i = 0; $i < @F-2; $i++){$info.=$F[$i]."\t";};$info=~s/\s+$//;my $ge=$F[-2].";".$F[-1];print "$info\t$ge"}else{print}' >3.txt
# Total_enrichment.txt is all of the input paml genes in each function
vi 3.txt; mv 3.txt Total_enrichment.txt; rm 1.txt 2.txt
less Total_enrichment.txt|perl -alne 'my @F=split /\t/;my $info;for (my $i = 0; $i < @F; $i++){$F[$i]=~s/^;//;$info.=$F[$i]."\t";};$info=~s/\s+$//;print $info' > Total_enrichment.txt.1
mv Total_enrichment.txt.1 Total_enrichment.txt
# kangjingliang@kangjinangdeMBP 一  5 22 20:34:04 ~/Documents/2023/WI/DEGs/Enrichment
cp CR_funcs.txt Vision_funcs.txt pH_funcs.txt ~/Documents/2023/WI/paml/free-ratio_new

# Check the gene of paml underlying functions
# kangjingliang@kangjingliangdeMacBook-Pro 一  5 22 22:57:43 ~/Documents/2023/WI/paml/free-ratio_new
perl temp2.pl Max_dNdS_common_enrichment.txt pH_funcs.txt > paml_pH.txt
perl temp2.pl Max_dNdS_common_enrichment.txt Vision_funcs.txt > paml_vision.txt
perl temp2.pl Max_dNdS_common_enrichment.txt CR_funcs.txt

# Positivie selected genes
# kangjingliang@kangjingliangdeMacBook-Pro 一  5 22 16:52:41 ~/Documents/2023/WI/paml/free-ratio_new
scp orth_id_paml_ano.txt jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new
# Kang@fishlab3 Mon May 22 16:57:48 /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/paml_input
scp Extract_PSGs.pl jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new
# jlkang@hpc2021 Mon May 22 17:02:29 /lustre1/g/sbs_schunter/Kang/orth16_new
perl Extract_PSGs.pl orth_id_paml_ano.txt > ortho16_new_PSGs.txt
```
## Give up evolution, change to cytoscape
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 四  5 25 09:49:02 ~/Documents/2023/WI/DEGs/Enrichment
less unprot_name_description_orthgroup.txt|perl -alne 'next if /^Orth_id/;(my $nm)=$F[1]=~/sp\|.*\|(.*)\_.*/;my $in=$nm."_HUMAN";print $in' > 1.txt
wc -l 1.txt # 22570
less 1.txt |sort -u|wc -l # 12574
# https://www.uniprot.org/uploadlists/ transvert to string id
# 12574 => 11416 string id (all_orth_stringid.txt)

# Transfer the paml results to my own station
# Kang@fishlab3 Mon May 29 13:49:20 /media/HDD/white_island/Compevo/orth16_new
scp -r jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files ./

###############################################
# Check the positive selection in pep sqeuences
###############################################
# Blenny
# OG0000120_2: sp|Q14831|GRM7_HUMAN	Metabotropic; glutamate receptor 7
# jlkang@hpc2021 Mon May 29 14:29:14 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000120_2
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
# jlkang@hpc2021 Mon May 29 14:30:29 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000120_2
Show_pos_in_seq.pl final_alignment_pep.fa 226
# distinct position: 829, 830, 837, 838, 841
# auto detect the positive selection one time
# jlkang@hpc2021 Mon May 29 15:20:06 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000120_2
Show_pos_in_seq_auto.pl Blenny-branch-site-alt-result.txt Blenny
# 829	830	837	838	841	Species
# A	A	N	G	K	Acura
# A	A	N	G	K	Apoly
# L	R	C	H	G	Blenny
# A	A	N	G	K	Blueeyed
# A	A	N	G	K	Common
# A	A	N	G	K	Daru
# A	A	N	G	K	Fugu
# A	A	N	G	K	Medaka
# A	A	N	G	K	Ocomp
# A	A	N	G	K	Padel
# A	A	N	G	K	Platyfish
# A	A	N	G	K	Pmol
# A	A	N	G	K	Spottedgar
# A	A	N	G	K	Stickleback
# A	A	N	G	K	Yaldwyn
# A	A	N	G	K	Zebrafish
perl temp1.pl > OG0000120_2_GRM7.fa
# kangjingliang@kangjingliangdeMacBook-Pro 一  5 29 16:42:55 ~/Documents/2023/WI/paml
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000120_2/OG0000120_2_GRM7.fa ./

# OG0000134_3: sp|P19493|GRIA4_RAT	Glutamate receptor 4
# jlkang@hpc2021 Mon May 29 14:33:27 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000134_3
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Blenny-branch-site-alt-result.txt Blenny
# pos: 407	409	412	413	415	416	420	421	429	430	433	436
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0000134_3_GRIA4.fa
# kangjingliang@kangjingliangdeMacBook-Pro 一  5 29 16:52:51 ~/Documents/2023/WI/paml
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000134_3/OG0000134_3_GRIA4.fa ./

# OG0000252_1:	sp|Q6U841|S4A10_HUMAN	Sodium-driven chloride bicarbonate exchanger
# jlkang@hpc2021 Mon May 29 15:22:17 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000134_3
cd ../OG0000252_1
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Blenny-branch-site-alt-result.txt Blenny
# pos: 711	712	713	716	726	732
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0000252_1_S4A10.fa
# kangjingliang@kangjingliangdeMacBook-Pro 一  5 29 16:58:39 ~/Documents/2023/WI/paml
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000252_1/OG0000252_1_S4A10.fa ./

# OG0000289_2	sp|Q38PU3|GRIK2_MACFA	Glutamate receptor ionotropic, kainate 2
# jlkang@hpc2021 Mon May 29 15:25:32 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000289_2
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Blenny-branch-site-alt-result.txt Blenny
# pos: 8	10	13	14	15	16	19	20	26	28	30	31	32	35	38	40	42	43	45
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0000289_2_GRIK2.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000289_2/OG0000289_2_GRIK2.fa ./


# OG0000493	sp|P51791|CLCN3_MOUSE	H(+)/Cl(-) exchange transporter 3
# jlkang@hpc2021 Mon May 29 15:26:58 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000493
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Blenny-branch-site-alt-result.txt Blenny
# pos: 16	40	81	86	89	93	96	103	105	108	112	168	170	175
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0000493_CLCN3.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000493/OG0000493_CLCN3.fa ./


# OG0004840	sp|O88602|CCG2_MOUSE	Voltage-dependent calcium channel gamma-2 subunit
# jlkang@hpc2021 Mon May 29 15:28:55 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0004840
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Blenny-branch-site-alt-result.txt Blenny
# pos: 59
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0004840_CCG2.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0004840/OG0004840_CCG2.fa ./


# OG0008000	sp|Q03064|MOT1_CRILO	Monocarboxylate transporter 1
# jlkang@hpc2021 Mon May 29 15:29:50 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0008000
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Blenny-branch-site-alt-result.txt Blenny
# pos: 270
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0008000_MOT1.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0008000/OG0008000_MOT1.fa ./


# OG0013997	sp|Q5ZJ75|SL9A8_CHICK	Sodium/hydrogen exchanger 8
# jlkang@hpc2021 Mon May 29 15:30:54 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0013997
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Blenny-branch-site-alt-result.txt Blenny
# pos: 318	323	324	325	328	329	330	331	339	341	344	345	347	348	349	353
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0013997_SL9A8.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0013997/OG0013997_SL9A8.fa ./


# OG0014857	sp|Q80T41|GABR2_MOUSE	Gamma-aminobutyric acid type B receptor subunit 2
# jlkang@hpc2021 Mon May 29 15:32:42 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0014857
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Blenny-branch-site-alt-result.txt Blenny
# pos: 1	3	8	9	11
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0014857_GABR2.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0014857/OG0014857_GABR2.fa ./


# Common
# OG0000133_1	sp|Q9H2X9|S12A5_HUMAN	Solute carrier family 12 member 5
# jlkang@hpc2021 Mon May 29 15:36:03 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000133_1
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Common-branch-site-alt-result.txt Common
# pos: 154	173	205	207	322	409	435	482
cd ../OG0000133_1
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0000133_1_S12A5.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000133_1/OG0000133_1_S12A5.fa ./


# OG0000493	sp|P51791|CLCN3_MOUSE	H(+)/Cl(-) exchange transporter 3
# jlkang@hpc2021 Mon May 29 15:38:21 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000493
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Common-branch-site-alt-result.txt Common
# pos: 332	395	416	422	423	434
cd ../OG0000493
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0000493_CLCN3.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000493/OG0000493_CLCN3.fa ./


# OG0000658	sp|Q61626|GRIK5_MOUSE	Glutamate receptor ionotropic, kainate 5
# jlkang@hpc2021 Mon May 29 15:39:19 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000658
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Common-branch-site-alt-result.txt Common
# pos: 1	7	9	10	12	13	14	15	16	19	20	22	23	32	33	52	84	166	167	168	170	175	177	187	188
cd ../OG0000658
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0000658_GRIK5.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000658/OG0000658_GRIK5.fa ./

# OG0000853	sp|P47869|GBRA2_HUMAN	Gamma-aminobutyric acid receptor subunit alpha-2
# jlkang@hpc2021 Mon May 29 15:40:34 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000853
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Common-branch-site-alt-result.txt Common
# pos: 99	105	112	113	114	115	116	117	118	120	121	122	124	130	131	132	134	135	138	141	142	144	145	146	148	152
cd ../OG0000853
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0000853_GBRA2.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000853/OG0000853_GBRA2.fa ./


# OG0002710	sp|Q00961|NMDE3_RAT	Glutamate receptor ionotropic, NMDA 2C
# jlkang@hpc2021 Mon May 29 15:42:03 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0002710
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Common-branch-site-alt-result.txt Common
# pos: 1	6	11	13	19	20	21	22	24	26
cd ../OG0002710
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0002710_NMDE3.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0002710/OG0002710_NMDE3.fa ./


# OG0007113	sp|P62956|CCG7_MOUSE	Voltage-dependent calcium channel gamma-7 subunit
# jlkang@hpc2021 Mon May 29 15:43:43 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0007113
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Common-branch-site-alt-result.txt Common
# pos: 174
cd ../OG0007113
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0007113_CCG7.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0007113/OG0007113_CCG7.fa ./

# OG0014857	sp|Q80T41|GABR2_MOUSE	Gamma-aminobutyric acid type B receptor subunit 2
# jlkang@hpc2021 Mon May 29 15:44:46 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0014857
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Common-branch-site-alt-result.txt Common
# pos: 240
cd ../OG0014857
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0014857_GABR2.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0014857/OG0014857_GABR2.fa ./


# OG0017327	sp|Q9XT54|PDE6D_CANLF	Retinal rod rhodopsin-sensitive cGMP 3',5'-cyclic phosphodiesterase subunit delta
# jlkang@hpc2021 Mon May 29 15:45:46 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0017327
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Common-branch-site-alt-result.txt Common
# pos: 1	2	6	7	9	11	13	14	17
cd ../OG0017327
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0017327_PDE6D.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0017327/OG0017327_PDE6D.fa ./


# OG0026786	sp|Q16864|VATF_HUMAN	V-type proton ATPase subunit F
# jlkang@hpc2021 Mon May 29 15:46:32 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0026786
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
Show_pos_in_seq_auto.pl Common-branch-site-alt-result.txt Common
# pos: 103	104	105	106	107
cd ../OG0026786
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0026786_VATF.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0026786/OG0026786_VATF.fa ./
```

## detect the position of all PSGs at one time
```detect_positive_site.pl
#!/usr/bin/perl
use strict;
use warnings;

my $File="final_alignment_pep.fa";
# my $numb=$ARGV[1];
open FILE, $File or die "can not open $File\n"; # pep fasta file
my (%hash, %hash1);
my $spe;
my @spes;

while (<FILE>) {
    chomp;
    if (/^\>/) {
        s/\>//;
        $spe=$_;
        push @spes, $spe;
    } else {
        my $seq=$_;
        my $len=length($seq);
        for (my $i = 0; $i < $len; $i++) {
                my $spepos=substr($seq,$i,1);
                $hash{$spe}->{$i}=$spepos;
        }
    }
}

my @poss;
my $resu=$ARGV[0];
open RESU, $resu or die "can not open $resu\n";
while (<RESU>) {
        chomp;
        if (/^Bayes Empirical Bayes \(BEB\)/) {
                while (<RESU>) {
                        chomp;
                        if (/\*/) {
                                s/^\s+//g;
                                my @a=split;
                                push @poss, $a[0];
                        }
                }
        } elsif (/The grid \(see ternary graph/) {
                last;
        }
}

my $info; my @poss2;
my $speci=$ARGV[1];
foreach my $pos (@poss) {
        my $pos1=$pos-1;
        my %hash2;
        foreach my $spe1 (@spes) {
                next if $spe1 eq $speci;
                my $spepos=$hash{$spe1}->{$pos1};
                $hash2{$spepos}++;
        }
        my @b=keys %hash2;
        push @poss2, $pos if @b==1 && ($hash{$speci}->{$pos1} ne $hash{"Zebrafish"}->{$pos1});
}

my $orth1=$ARGV[2];

unless (@poss2) {
	die "no distinct position in $orth1\n";
}

foreach my $pos (@poss2) {
        $info.=$pos.";";
}
$info=~s/\;$//;
print "$speci\t$orth1\t$info\n" if @poss2;
```

```detect_positive_site_all.pl
#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

# cds2pep.pl final_alignment.fa > final_alignment_pep.fa
# Show_pos_in_seq_auto.pl Blenny-branch-site-alt-result.txt Blenny

my $genes=$ARGV[0]; # positively selected genes
my $resul=$ARGV[1]; # positive selection results
my $spe  =$ARGV[2]; # species
my $pwd = getcwd;

open GENES, $genes or die "can not open $genes\n";
while (<GENES>) {
	chomp;
	my $orth=$_;
	my $dir="paml_files/$orth";
	system("cp detect_positive_site.pl $dir/");
	chdir $dir;
	system("cds2pep.pl final_alignment.fa > final_alignment_pep.fa");
	system("perl detect_positive_site.pl $resul $spe $orth");
	chdir $pwd;
}
```

```bash
# Postitively selected genes detected by PAML
# jlkang@hpc2021 Thu Jun 29 11:11:34 /lustre1/g/sbs_schunter/Kang/orth16_new
cat postively_selected_genes/Blenny-*.txt > Blenny_paml_PSGs.txt # 789 PSGs in blenny
# fdr
R
#############
# fdr correct
p_apoly<-read.table(file="Blenny_paml_PSGs.txt")
p_apoly$fdr<- p.adjust(p_apoly$V5,method="fdr",length(p_apoly$V5))
write.table(p_apoly, file="Blenny_paml_PSGs_fdr.txt",row.names=F,col.names=F,quote=F,sep="\t")
#############
# 789 PSGs in blenny after fdr

cat postively_selected_genes/Blueeyed-*.txt > Blueeyed_paml_PSGs.txt # 680 PSGs in Blueeyed
# fdr
R
#############
# fdr correct
p_apoly<-read.table(file="Blueeyed_paml_PSGs.txt")
p_apoly$fdr<- p.adjust(p_apoly$V5,method="fdr",length(p_apoly$V5))
write.table(p_apoly, file="Blueeyed_paml_PSGs_fdr.txt",row.names=F,col.names=F,quote=F,sep="\t")
#############
less Blueeyed_paml_PSGs_fdr.txt|cut -f 5|perl -alne 'print if $_<=0.05'|wc -l # 680 PSGs in Blueeyed
less Blueeyed_paml_PSGs_fdr.txt|cut -f 6|perl -alne 'print if $_<=0.05'|wc -l # 680 PSGs in Blueeyed after fdr

cat postively_selected_genes/Common-*.txt > Common_paml_PSGs.txt # 582 PSGs in Common
# fdr
R
#############
# fdr correct
p_apoly<-read.table(file="Common_paml_PSGs.txt")
p_apoly$fdr<- p.adjust(p_apoly$V5,method="fdr",length(p_apoly$V5))
write.table(p_apoly, file="Common_paml_PSGs_fdr.txt",row.names=F,col.names=F,quote=F,sep="\t")
#############
less Common_paml_PSGs_fdr.txt|cut -f 5|perl -alne 'print if $_<=0.05'|wc -l # 582 PSGs in Common
less Common_paml_PSGs_fdr.txt|cut -f 6|perl -alne 'print if $_<=0.05'|wc -l # 582 PSGs in Common after fdr

cat postively_selected_genes/Yaldwyn-*.txt > Yaldwyn_paml_PSGs.txt # 502 PSGs in Yaldwyn
# fdr
R
#############
# fdr correct
p_apoly<-read.table(file="Yaldwyn_paml_PSGs.txt")
p_apoly$fdr<- p.adjust(p_apoly$V5,method="fdr",length(p_apoly$V5))
write.table(p_apoly, file="Yaldwyn_paml_PSGs_fdr.txt",row.names=F,col.names=F,quote=F,sep="\t")
#############
less Yaldwyn_paml_PSGs_fdr.txt|cut -f 5|perl -alne 'print if $_<=0.05'|wc -l # 502 PSGs in Yaldwyn
less Yaldwyn_paml_PSGs_fdr.txt|cut -f 6|perl -alne 'print if $_<=0.05'|wc -l # 502 PSGs in Yaldwyn after fdr

# jlkang@hpc2021 Thu Jun 29 11:30:25 /lustre1/g/sbs_schunter/Kang/orth16_new
mkdir PAML_PSGs
mv *_paml_PSGs.txt *_paml_PSGs_fdr.txt PAML_PSGs/
# Kang@fishlab3 Thu Jun 29 11:12:31 /media/HDD/white_island/Compevo/orth16_new
scp -r jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/PAML_PSGs/ ./
less PAML_PSGs/Blueeyed_paml_PSGs_fdr.txt|cut -f 1 > Blueeyed_PSGs.txt
less PAML_PSGs/Yaldwyn_paml_PSGs_fdr.txt|cut -f 1 > Yaldwyn_PSGs.txt

# annotations to orthologous genes
# Kang@fishlab3 Thu Jun 29 11:54:03 /media/HDD/white_island/Compevo/orth16_new
scp kang1234@147.8.76.177:~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar21/Orthogroups/orth_id_paml_ano.txt ./

vi detect_positive_site.pl
# Kang@fishlab3 Tue May 30 23:50:52 /media/HDD/white_island/Compevo/orth16_new
perl detect_positive_site_all.pl Blenny_PSGs.txt Blenny-branch-site-alt-result.txt Blenny > Blenny_PSGs_distinc.txt # 789 => 540
perl detect_positive_site_all.pl Common_PSGs.txt Common-branch-site-alt-result.txt Common > Common_PSGs_distinc.txt # 582 => 389
perl detect_positive_site_all.pl Blueeyed_PSGs.txt Blueeyed-branch-site-alt-result.txt Blueeyed > Blueeyed_PSGs_distinc.txt # 689 => 479
perl detect_positive_site_all.pl Yaldwyn_PSGs.txt Yaldwyn-branch-site-alt-result.txt Yaldwyn > Yaldwyn_PSGs_distinc.txt # 502 => 361
# add the annotation information to *PSGs_distinc.txt
perl temp4.pl orth_id_paml_ano.txt Blenny_PSGs_distinc.txt Blenny > Blenny_PSGs_distinc_info.txt
perl temp4.pl orth_id_paml_ano.txt Common_PSGs_distinc.txt Common > Common_PSGs_distinc_info.txt
perl temp4.pl orth_id_paml_ano.txt Blueeyed_PSGs_distinc.txt Blueeyed > Blueeyed_PSGs_distinc_info.txt
perl temp4.pl orth_id_paml_ano.txt Yaldwyn_PSGs_distinc.txt Yaldwyn > Yaldwyn_PSGs_distinc_info.txt
```

## Detect convergent
```detect_convergent_site.pl
#!/usr/bin/perl
use strict;
use warnings;

my $File="final_alignment_pep.fa";
# my $numb=$ARGV[1];
open FILE, $File or die "can not open $File\n"; # pep fasta file
my (%hash, %hash1);
my $spe;
my @spes;

while (<FILE>) {
    chomp;
    if (/^\>/) {
        s/\>//;
        $spe=$_;
        push @spes, $spe;
    } else {
        my $seq=$_;
        my $len=length($seq);
        for (my $i = 0; $i < $len; $i++) {
                my $spepos=substr($seq,$i,1);
                $hash{$spe}->{$i}=$spepos;
        }
    }
}

my @poss;
my $resu=$ARGV[0];
open RESU, $resu or die "can not open $resu\n";
while (<RESU>) {
        chomp;
        if (/^Bayes Empirical Bayes \(BEB\)/) {
                while (<RESU>) {
                        chomp;
                        if (/\*/) {
                                s/^\s+//g;
                                my @a=split;
                                push @poss, $a[0];
                        }
                }
        } elsif (/The grid \(see ternary graph/) {
                last;
        }
}

my $info; my @poss2;
my $speci=$ARGV[1];
my $spec2=$ARGV[2];
foreach my $pos (@poss) {
        my $pos1=$pos-1;
        my %hash2;
        foreach my $spe1 (@spes) {
                next if (($spe1 eq $speci) || ($spe1 eq $spec2));
                my $spepos=$hash{$spe1}->{$pos1};
                $hash2{$spepos}++;
        }
        my @b=keys %hash2;
        push @poss2, $pos if @b==1 && ($hash{$speci}->{$pos1} eq $hash{$spec2}->{$pos1}) && ($hash{$speci}->{$pos1} ne $hash{"Zebrafish"}->{$pos1});
}

my $orth1=$ARGV[3];

unless (@poss2) {
        die "no distinct position in $orth1\n";
}

foreach my $pos (@poss2) {
        $info.=$pos.";";
}
$info=~s/\;$//;
print "convergent\t$speci\t$orth1\t$info\n" if @poss2;
```

```detect_convergent_site_all.pl
#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

# cds2pep.pl final_alignment.fa > final_alignment_pep.fa
# Show_pos_in_seq_auto.pl Blenny-branch-site-alt-result.txt Blenny
# example: perl detect_convergent_site_all.pl Blenny_PSGs.txt Blenny-branch-site-alt-result.txt Blenny Common
# example: perl detect_convergent_site.pl Blenny-branch-site-alt-result.txt Blenny Common $orth

my $genes=$ARGV[0]; # positively selected genes
my $resul=$ARGV[1]; # positive selection results
my $spe1  =$ARGV[2]; # species
my $spe2  =$ARGV[3]; # species
my $pwd = getcwd;

open GENES, $genes or die "can not open $genes\n";
while (<GENES>) {
        chomp;
        my $orth=$_;
        my $dir="paml_files/$orth";
        system("cp detect_convergent_site.pl $dir/");
        chdir $dir;
        system("cds2pep.pl final_alignment.fa > final_alignment_pep.fa");
        system("perl detect_convergent_site.pl $resul $spe1 $spe2 $orth");
        chdir $pwd;
}
```

```bash
# Kang@fishlab3 Wed May 31 00:34:00 /media/HDD/white_island/Compevo/orth16_new
perl detect_convergent_site_all.pl Blenny_PSGs.txt Blenny-branch-site-alt-result.txt Blenny Common > Convergent_PSGs.txt
perl detect_convergent_site_all.pl Common_PSGs.txt Common-branch-site-alt-result.txt Common Blenny >> Convergent_PSGs.txt

# sort the convergent genes
```

```temp2.pl
#!/usr/bin/perl
use strict;
use warnings;

my $conv="Convergent_PSGs.txt";
my (%hash1, %hash2, %hash3);
my @orths;
open CONV, $conv or die "can not open $conv\n";
while (<CONV>) {
	chomp;
	my @a=split /\t/;
	my ($orth, $pos)=($a[2], $a[3]);
	$hash3{$orth}++;
	push @orths, $orth if $hash3{$orth}==1;
	my @b=split /\;/, $pos;
	foreach my $pos1 (@b) {
		$hash1{$orth}->{$pos1}++;
		if ($hash1{$orth}->{$pos1} == 1) { # "==1" combine all site; "==2" same site by both of two species
			$hash2{$orth}=[] unless $hash2{$orth};
			push @{$hash2{$orth}}, $pos1;
		}
	}
}

foreach my $orth (@orths) {
	next unless $hash2{$orth};
	my @pos=@{$hash2{$orth}};
	my $info;
	foreach my $pos (sort {$a<=>$b} @pos) {
		$info.=$pos.";";
	}
	$info=~s/\;$//;
	print "$orth\t$info\n";
}
```

```bash
# Kang@fishlab3 Wed May 31 01:22:13 /media/HDD/white_island/Compevo/orth16_new
perl temp2.pl # "==1"
# OG0001301	93;255;453;489;509;532;635;648;697;700;743
# OG0002007	458;460;461;462	Transient receptor potential cation channel subfamily M member 7
# OG0007782	351
# OG0016105	164
# OG0022156	402;407;409;412
# OG0023496	7;39;60;68;88;92;102;128
# OG0023786	588;590;594;595;596;599;600;601;604;606;608;617;628
# OG0014462	273
# OG0015063	122

perl temp2.pl # "==2"
# OG0001301	93;255;453;509;635;697
# OG0002007	458;460;461
# OG0007782	351
# OG0016105	164
# OG0022156	402;407;412
# OG0023496	7;88;92;102
# OG0023786	599;604;628
python ~/software/Busco/scripts/generate_plot.py -wd ./

# Find the id of the single copy genes in orth4 for bamm based on orth10 single cope genes
# orth10: Kang@fishlab3 Thu Jul 20 12:07:59 /media/HDD/white_island/Compevo/Orth_ten_CO2/OrthoFinder/Orthogroups
# The selected capture orthogroups: Kang@fishlab3 Thu Jul 20 11:59:24 /media/HDD/white_island/Compevo/genecapture/Total_orf/select_WI
# orth4: Kang@fishlab3 Tue Jul 04 14:45:07 /media/HDD/white_island/orthologue/input_pep/orthofinder_input/OrthoFinder/Results_Mar29/final_reference
vi orth10_capture_Yaldwyn_id.txt
# two orth in orth10 were not in orth4
# Yaldwyn_230592 (OG0061954); Yaldwyn_39008 (OG0062118)
# Kang@fishlab3 Thu Jul 20 13:23:14 /media/HDD/white_island/Compevo/Orth_ten_CO2/OrthoFinder/Orthogroups
grep -i 'Yaldwyn_230592' Orthogroups.txt; grep -i 'Yaldwyn_39008' Orthogroups.txt
# orth10 annotation: Kang@fishlab3 Thu Jul 20 13:39:59 /media/HDD/white_island/Compevo/Orth_ten_CO2/OrthoFinder/sub_orth/ortho10_ano.txt
# Kang@fishlab3 Thu Jul 20 13:40:00 /media/HDD/white_island/Compevo/Orth_ten_CO2/OrthoFinder/sub_orth
vi gene_capture_id.txt
perl temp1.pl

# WI PSGs plot
# kangjingliang@kangjingliangdeMacBook-Pro 二  7 04 20:34:35 ~/Desktop
mkdir WI_PSGs_plot; cd WI_PSGs_plot

# Blenny
# OG0013997: Sodium/hydrogen exchanger 8; SL9A8
# kangjingliang@kangjingliangdeMacBook-Pro 二  7 04 20:39:30 ~/Desktop/WI_PSGs_plot
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0013997/Blenny_BUSTED.json Blenny_BUSTED_OG0013997_SL9A8.json
scp Kang@147.8.76.229:/media/HDD/white_island/Compevo/orth16_new/paml_files/OG0013997/final_alignment_pep.fa Blenny_OG0013997_SL9A8_pep.fa
less Blenny_BUSTED_OG0013997_SL9A8.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my ($i, $j);foreach my $a (@a){$i++;$j.=$i.";" if $a >= 10};$j=~s/\;$//;print $j'
# 316;318;324;326;328;329;330;339;343;345;353
# select the distinct sites under positive selection by both PAML and Hyphy
# PAML
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:00:00 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Blenny_PSGs_distinc.txt|grep 'OG0013997'|perl -alne '@a=split /\;/, $F[2]; foreach my $a(@a){print $a}'
# Hyphy
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:03:57 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Blenny_BUSTED_OG0013997_SL9A8.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my $i;foreach my $a (@a){$i++; print $i if $a >= 10}'
# venn
# 318;324;328;329;330;339;345;353

# Blue-eyed
# Blue-eyed: the results were convered because i forgot to rename the result json; need rerun
# OG0000770_2: Sodium/hydrogen exchanger 7; SL9A7
# jlkang@hpc2021 Sun Jul 09 16:48:01 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000770_2
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
# kangjingliang@kangjingliangdeMacBook-Pro 日  7 09 17:14:56 ~/Desktop/WI_PSGs_plot
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000770_2/Blueeyed_BUSTED.json Blueeyed_BUSTED_OG0000770_2_SL9A7.json
scp Kang@147.8.76.229:/media/HDD/white_island/Compevo/orth16_new/paml_files/OG0000770_2/final_alignment_pep.fa Blueeyed_OG0000770_2_SL9A7_pep.fa
less Blueeyed_BUSTED_OG0000770_2_SL9A7.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my ($i, $j);foreach my $a (@a){$i++;$j.=$i.";" if $a >= 10};$j=~s/\;$//;print $j'
# 298;299;321;322;323;324;325;326;327;328;329;330;331;333;334;335;337;338;339;340
# select the distinct sites under positive selection by both PAML and Hyphy
# PAML
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:00:00 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Blueeyed_PSGs_distinc.txt|grep 'OG0000770_2'|perl -alne '@a=split /\;/, $F[2]; foreach my $a(@a){print $a}'
# Hyphy
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:03:57 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Blueeyed_BUSTED_OG0000770_2_SL9A7.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my $i;foreach my $a (@a){$i++; print $i if $a >= 10}'
# venn
# 298;299;321;322;323;324;325;326;327;328;329;330;331;333;337;339;340


# OG0004342       CAH14   Carbonic anhydrase 14
# kangjingliang@kangjingliangdeMacBook-Pro 日  7 09 17:14:56 ~/Desktop/WI_PSGs_plot
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0004342/Blueeyed_BUSTED.json Blueeyed_BUSTED_OG0004342_CAH14.json
scp Kang@147.8.76.229:/media/HDD/white_island/Compevo/orth16_new/paml_files/OG0004342/final_alignment_pep.fa Blueeyed_OG0004342_CAH14_pep.fa
less Blueeyed_BUSTED_OG0004342_CAH14.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my ($i, $j);foreach my $a (@a){$i++;$j.=$i.";" if $a >= 10};$j=~s/\;$//;print $j'
# 186;187;188;189;190;194;196;198;199;200;201;202;203;204;205;206;207;208;210;211;212;213;214
# select the distinct sites under positive selection by both PAML and Hyphy
# PAML
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:00:00 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Blueeyed_PSGs_distinc.txt|grep 'OG0004342'|perl -alne '@a=split /\;/, $F[2]; foreach my $a(@a){print $a}'
# Hyphy
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:03:57 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Blueeyed_BUSTED_OG0004342_CAH14.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my $i;foreach my $a (@a){$i++; print $i if $a >= 10}'
# venn
# 188;196;206;211;214

# OG0005553     VATE1   V-type proton ATPase subunit E 1
# jlkang@hpc2021 Sun Jul 09 17:15:49 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0005553
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0005553/Blueeyed_BUSTED.json Blueeyed_BUSTED_OG0005553_VATE1.json
scp Kang@147.8.76.229:/media/HDD/white_island/Compevo/orth16_new/paml_files/OG0005553/final_alignment_pep.fa Blueeyed_OG0005553_VATE1_pep.fa
less Blueeyed_BUSTED_OG0005553_VATE1.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my ($i, $j);foreach my $a (@a){$i++;$j.=$i.";" if $a >= 10};$j=~s/\;$//;print $j'
# 206;209;210;211;212;215;216;217;218;219;220;221;222
# select the distinct sites under positive selection by both PAML and Hyphy
# PAML
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:00:00 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Blueeyed_PSGs_distinc.txt|grep 'OG0005553'|perl -alne '@a=split /\;/, $F[2]; foreach my $a(@a){print $a}'
# Hyphy
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:03:57 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Blueeyed_BUSTED_OG0005553_VATE1.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my $i;foreach my $a (@a){$i++; print $i if $a >= 10}'
# venn
# 206;211;215;216;218;219;220;221;222


# OG0010541     CLCN7   H(+)/Cl(-) exchange transporter 7
# jlkang@hpc2021 Sun Jul 09 17:19:31 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0010541
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0010541/Blueeyed_BUSTED.json Blueeyed_BUSTED_OG0010541_CLCN7.json
scp Kang@147.8.76.229:/media/HDD/white_island/Compevo/orth16_new/paml_files/OG0010541/final_alignment_pep.fa Blueeyed_OG0010541_CLCN7_pep.fa
less Blueeyed_BUSTED_OG0010541_CLCN7.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my ($i, $j);foreach my $a (@a){$i++;$j.=$i.";" if $a >= 10};$j=~s/\;$//;print $j'
# 304;305;307;308;309;310
# select the distinct sites under positive selection by both PAML and Hyphy
# PAML
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:00:00 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Blueeyed_PSGs_distinc.txt|grep 'OG0010541'|perl -alne '@a=split /\;/, $F[2]; foreach my $a(@a){print $a}'
# Hyphy
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:03:57 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Blueeyed_BUSTED_OG0010541_CLCN7.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my $i;foreach my $a (@a){$i++; print $i if $a >= 10}'
# venn
# 305;309;310

# Common
# OG0026786     VATF    V-type proton ATPase subunit F
# jlkang@hpc2021 Sun Jul 09 17:33:55 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0026786
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0026786/Common_BUSTED.json Common_BUSTED_OG0026786_VATF.json
scp Kang@147.8.76.229:/media/HDD/white_island/Compevo/orth16_new/paml_files/OG0026786/final_alignment_pep.fa Common_OG0026786_VATF_pep.fa
less Common_BUSTED_OG0026786_VATF.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my ($i, $j);foreach my $a (@a){$i++;$j.=$i.";" if $a >= 10};$j=~s/\;$//;print $j'
# 103;104;105;106;107
# select the distinct sites under positive selection by both PAML and Hyphy
# PAML
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:00:00 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Common_PSGs_distinc.txt|grep 'OG0026786'|perl -alne '@a=split /\;/, $F[2]; foreach my $a(@a){print $a}'
# Hyphy
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:03:57 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Common_BUSTED_OG0026786_VATF.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my $i;foreach my $a (@a){$i++; print $i if $a >= 10}'
# venn
# 103;104;105;106;107

# Yaldwin
# OG0000466: Chloride channel protein 2; CLCN2
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000466/Yaldwyn_BUSTED.json Yaldwyn_BUSTED_OG0000466_CLCN2.json
scp Kang@147.8.76.229:/media/HDD/white_island/Compevo/orth16_new/paml_files/OG0000466/final_alignment_pep.fa Yaldwyn_OG0000466_CLCN2_pep.fa
less Yaldwyn_BUSTED_OG0000466_CLCN2.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my ($i, $j);foreach my $a (@a){$i++;$j.=$i.";" if $a >= 10};$j=~s/\;$//;print $j'
# 366;382;394;486
# select the distinct sites under positive selection by both PAML and Hyphy
# PAML
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:00:00 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Yaldwyn_PSGs_distinc.txt|grep 'OG0000466'|perl -alne '@a=split /\;/, $F[2]; foreach my $a(@a){print $a}'
# Hyphy
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:03:57 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Yaldwyn_BUSTED_OG0000466_CLCN2.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my $i;foreach my $a (@a){$i++; print $i if $a >= 10}'
# venn
# 382;486

# OG0001941_1     CLIC4   Chloride intracellular channel protein 4
# jlkang@hpc2021 Sun Jul 09 17:25:20 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0001941_1
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
# kangjingliang@kangjingliangdeMacBook-Pro 日  7 09 17:21:50 ~/Desktop/WI_PSGs_plot
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0001941_1/Yaldwyn_BUSTED.json Yaldwyn_BUSTED_OG0001941_1_CLIC4.json
scp Kang@147.8.76.229:/media/HDD/white_island/Compevo/orth16_new/paml_files/OG0001941_1/final_alignment_pep.fa Yaldwyn_OG0001941_1_CLIC4_pep.fa
less Yaldwyn_BUSTED_OG0001941_1_CLIC4.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my ($i, $j);foreach my $a (@a){$i++;$j.=$i.";" if $a >= 10};$j=~s/\;$//;print $j'
# 1;2;3;4;5;6;8;11;15;17;20;23;24;26;27;28;29;30;31;32;34;35;36;37;38;39;40;41;43;44;45;46;48;51;52;53;54;58;59
# select the distinct sites under positive selection by both PAML and Hyphy
# PAML
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:00:00 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Yaldwyn_PSGs_distinc.txt|grep 'OG0001941_1'|perl -alne '@a=split /\;/, $F[2]; foreach my $a(@a){print $a}'
# Hyphy
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:03:57 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Yaldwyn_BUSTED_OG0001941_1_CLIC4.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my $i;foreach my $a (@a){$i++; print $i if $a >= 10}'
# venn
# 1;2;3;4;5;6;8;11;15;17;20;23;24;26;27;28;29;30;31;32;34;35;36;38;39;40;41;43;44;45;46;48;51;52;53;54;58;59


# OG0019184     CLIC2   Chloride intracellular channel protein 2
# jlkang@hpc2021 Sun Jul 09 17:29:34 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0019184
cds2pep.pl final_alignment.fa > final_alignment_pep.fa
# kangjingliang@kangjingliangdeMacBook-Pro 日  7 09 17:21:50 ~/Desktop/WI_PSGs_plot
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0019184/Yaldwyn_BUSTED.json Yaldwyn_BUSTED_OG0019184_CLIC2.json
scp Kang@147.8.76.229:/media/HDD/white_island/Compevo/orth16_new/paml_files/OG0019184/final_alignment_pep.fa Yaldwyn_OG0019184_CLIC2_pep.fa
less Yaldwyn_BUSTED_OG0019184_CLIC2.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my ($i, $j);foreach my $a (@a){$i++;$j.=$i.";" if $a >= 10};$j=~s/\;$//;print $j'
# 133;134;135;136;137
# select the distinct sites under positive selection by both PAML and Hyphy
# PAML
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:00:00 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Yaldwyn_PSGs_distinc.txt|grep 'OG0019184'|perl -alne '@a=split /\;/, $F[2]; foreach my $a(@a){print $a}'
# Hyphy
# kangjingliang@kangjingliangdeMacBook-Pro 一  7 10 17:03:57 ~/Documents/2023/WI/paml/WI_PSGs_plot
less Yaldwyn_BUSTED_OG0019184_CLIC2.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my $i;foreach my $a (@a){$i++; print $i if $a >= 10}'
# venn
# 134

################
# GABA receptors
################
# Crested blenny
# OG0014857       GABR2   Gamma-aminobutyric acid type B receptor subunit 2
# kangjingliang@kangjingliangdeMacBook-Pro 三  7 19 10:09:28 ~/Documents/2023/WI/paml/WI_PSGs_plot
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0014857/Blenny_BUSTED.json Blenny_BUSTED_OG0014857_GABR2.json
scp Kang@147.8.76.229:/media/HDD/white_island/Compevo/orth16_new/paml_files/OG0014857/final_alignment_pep.fa Blenny_OG0014857_GABR2_pep.fa
perl temp1.pl Blenny_OG0014857_GABR2_pep.fa > Blenny_OG0014857_GABR2_pep.fa.1
mv Blenny_OG0014857_GABR2_pep.fa.1 Blenny_OG0014857_GABR2_pep.fa
# Hyphy
less Blenny_BUSTED_OG0014857_GABR2.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my ($i, $j);foreach my $a (@a){$i++;$j.=$i.";" if $a >= 10};$j=~s/\;$//;print $j'
# 3;6;8
less blenny_PSGs_distinc.txt|grep 'OG0014857'|perl -alne '@a=split /\;/, $F[2]; foreach my $a(@a){print $a}'
# 1;3;8;9;11
# common position under positive selection
# 3;8

# Blue-eyed
# OG0000523_1   GBRB2   Gamma-aminobutyric acid receptor subunit beta-2
# kangjingliang@kangjingliangdeMacBook-Pro 三  7 19 10:09:28 ~/Documents/2023/WI/paml/WI_PSGs_plot
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000523_1/Blueeyed_BUSTED.json Blueeyed_BUSTED_OG0000523_1_GBRB2.json
scp Kang@147.8.76.229:/media/HDD/white_island/Compevo/orth16_new/paml_files/OG0000523_1/final_alignment_pep.fa Blueeyed_OG0000523_1_GBRB2_pep.fa
perl temp1.pl Blueeyed_OG0000523_1_GBRB2_pep.fa > Blueeyed_OG0000523_1_GBRB2_pep.fa.1
mv Blueeyed_OG0000523_1_GBRB2_pep.fa.1 Blueeyed_OG0000523_1_GBRB2_pep.fa
#  Hyphy
less Blueeyed_BUSTED_OG0000523_1_GBRB2.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my ($i, $j);foreach my $a (@a){$i++;print $i if $a >= 10}'
# 1;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23
# paml
less Blueeyed_PSGs_distinc.txt|grep 'OG0000523_1'|perl -alne '@a=split /\;/, $F[2]; foreach my $a(@a){print $a}'
# common position under positive selection
# 1;3;4;5;7;8;10;11;12;14;15;16;17;18;19;20;21;22;23

# Common
# OG0000853       GBRA2   Gamma-aminobutyric acid receptor subunit alpha-2
# kangjingliang@kangjingliangdeMacBook-Pro 三  7 19 10:09:28 ~/Documents/2023/WI/paml/WI_PSGs_plot
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000853/Common_BUSTED.json Common_BUSTED_OG0000853_GBRA2.json
scp Kang@147.8.76.229:/media/HDD/white_island/Compevo/orth16_new/paml_files/OG0000853/final_alignment_pep.fa Common_OG0000853_GBRA2_pep.fa
perl temp1.pl Common_OG0000853_GBRA2_pep.fa > Common_OG0000853_GBRA2_pep.fa.1
mv Common_OG0000853_GBRA2_pep.fa.1 Common_OG0000853_GBRA2_pep.fa
# Hyphy
less Common_BUSTED_OG0000853_GBRA2.json|head -n 4 |tail -n 1|perl -alne 's/\[//g;s/\]//g;s/\,//g;my @a=split;my ($i, $j);foreach my $a (@a){$i++;print $i if $a >= 10}'
less Common_PSGs_distinc.txt|grep 'OG0000853'|perl -alne '@a=split /\;/, $F[2]; foreach my $a(@a){print $a}'
# common position under positive selection
# 112;113;115;116;118;121;122;124;130;131;132;141;144;145;146;148;152

# Confirm the change of fastq name
# Kang@fishlab3 Fri Jul 21 10:15:32 /media/HDD/white_island/paired
perl temp1.pl > Fq_nm_change.txt
perl temp2.pl > final_nm_fq.txt
perl temp3.pl >filtering_info_final.txt # filtering information needed
```

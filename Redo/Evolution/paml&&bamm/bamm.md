```bash
# Kang@fishlab3 Wed Apr 12 09:36:12 /media/HDD/white_island/Compevo/Orth_ten_CO2/OrthoFinder/sub_orth
less sub_orth_genecount.txt|perl -alne 'my $j;for (my $i = 1; $i <= 10; $i++){$j++ if $F[$i]>=1};print if $j==10'|wc -l
# 9068 sub_orth

# Kang@fishlab3 Wed Apr 12 09:39:16 /media/HDD/white_island/Compevo/Orth_ten_CO2/OrthoFinder/sub_orth
cp /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/all_swissprot_diamond_ano.txt ./
cp /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/sub_orth_same_name.pl ./
cp /media/HDD/white_island/Compevo/Results_Mar21/sub_orth/Generate_paml_list.pl ./

less all_swissprot_diamond_ano.txt|perl -alne 's/Blue_eyed/Blueeyed/g;print' >all_swissprot_diamond_ano_1.txt
mv all_swissprot_diamond_ano_1.txt all_swissprot_diamond_ano.txt
perl sub_orth_same_name.pl # output: sub_orth_genecount_final.txt; sub_orth_id_final.txt # 7508 sub_orth
perl Generate_paml_list.pl >sub_orth_id_paml.txt

# Create the reference nucleotide sequences to map against
scp kang1234@147.8.76.177:~/white_island/Compevo/Orth_ten_CO2/*.cds.fasta ./
```

```Create_ref.pl
use strict;
use warnings;

my $allcds="all.cds.fasta";
system("cat *.cds.fasta > $allcds");
my %nuc;
my $nm;
open ALL, $allcds or die "can not open $allcds\n";
while (<ALL>) {
	chomp;
	if (/\>/) {
		s/\>//;
		$nm=$_;
	} else {
		$nuc{$nm}.=$_;
	}
}
system("rm $allcds");

my $orthid="sub_orth_id_paml.txt";
open ORTH, $orthid or die "can not open $orthid\n";
while (<ORTH>) {
	chomp;
	next if /^Suborth/;
	my @a=split;
	for (my $i = 1; $i < @a; $i++) {
		my $id=$a[0]."||".$a[$i];
		my $nm=$a[$i];
		(my $spe)=$a[$i]=~/(.*)\_\d+/;
		my $outp=$spe."_refnuc.fasta";
		open OUTP, ">>$outp" or die "can not create $outp\n";
		my $seq=$nuc{$nm};
		print OUTP ">$id\n$seq\n";
	}
}
```

```bash
# Kang@fishlab3 Wed Apr 12 10:50:41 /media/HDD/white_island/Compevo/Orth_ten_CO2/OrthoFinder/sub_orth
perl Create_ref.pl # resulted *_refnuc.fasta as the reference per species to map

# (base) kang1234@celia-PowerEdge-T640 Wed Apr 12 10:58:04 ~/white_island/kraken
mkdir Orth10_reference; cd Orth10_reference
# Kang@fishlab3 Wed Apr 12 11:01:59 /media/HDD/white_island/Compevo/Orth_ten_CO2/OrthoFinder/sub_orth
scp *_refnuc.fasta kang1234@147.8.76.177:~/white_island/kraken/Orth10_reference
```

```Create_ref.sh
rsem-prepare-reference --bowtie2 Blenny_refnuc.fasta Blenny --bowtie2
rsem-prepare-reference --bowtie2 Blueeyed_refnuc.fasta Blueeyed --bowtie2
rsem-prepare-reference --bowtie2 Common_refnuc.fasta Common --bowtie2
rsem-prepare-reference --bowtie2 Yaldwyn_refnuc.fasta Yaldwyn --bowtie2
rsem-prepare-reference --bowtie2 Acura_refnuc.fasta Acura --bowtie2
rsem-prepare-reference --bowtie2 Apoly_refnuc.fasta Apoly --bowtie2
rsem-prepare-reference --bowtie2 Daru_refnuc.fasta Daru --bowtie2
rsem-prepare-reference --bowtie2 Ocomp_refnuc.fasta Ocomp --bowtie2
rsem-prepare-reference --bowtie2 Padel_refnuc.fasta Padel --bowtie2
rsem-prepare-reference --bowtie2 Pmol_refnuc.fasta Pmol --bowtie2
```

```bash
# download rsem_perl_utils.pm
# (base) kang1234@celia-PowerEdge-T640 Wed Apr 12 11:22:29 ~/white_island/kraken/Orth10_reference
sudo mv rsem_perl_utils.pm /usr/bin
# (base) kang1234@celia-PowerEdge-T640 Wed Apr 12 11:14:30 ~/white_island/kraken/Orth10_reference
nohup sh Create_ref.sh > Create_ref.process 2>&1 &
cp ~/white_island/kraken/final_reference/*_info.txt ./
```

```RSEM_align.pl
use strict;
use warnings;
open "fil", "$ARGV[0]" or die "can not open $ARGV[0]";
while (<fil>) {
        chomp;
        my @a=split;
        my $spe=$a[0];
        my $result_dir=$spe."_RSEM_output";
        mkdir $result_dir unless (-d $result_dir);
        chdir "./$result_dir";
#       system 'pwd';
        my @R1=split /\,/, $a[1];
        foreach my $R1 (@R1) {
                my ($sample)=$R1=~/(.*)\_1\.fastq\.gz/;
                my $R2=$sample."_2.fastq.gz";
                my $cmd="rsem-calculate-expression -p 24 --bowtie2 --paired-end ../../$R1 ../../$R2 ../$spe $sample";
                system "$cmd";
        }
        chdir "../";
#       system 'pwd';
}
```

```bash
# (base) kang1234@celia-PowerEdge-T640 Wed Apr 12 11:29:50 ~/white_island/kraken/Orth10_reference
nohup perl RSEM_align.pl Blenny_info.txt >Blenny_RSEM_process.txt 2>&1 &
# [1] 27792
# Too slow, transfer some data to HPC do RSEM
# (base) kang1234@celia-PowerEdge-T640 Thu Apr 13 15:09:00 ~/white_island/kraken
less Orth10_reference/Blue_eyed_info.txt|perl -alne '@a=split /\,/, $F[1];foreach my $a (@a){(my $ind)=$a=~/(.*)_1\.fastq\.gz/;$info.=$a." ".$ind."_2.fastq.gz "};END{$info=~s/\s+$//;print $info}'
scp B21_1.fastq.gz B21_2.fastq.gz B22_1.fastq.gz B22_2.fastq.gz B23_1.fastq.gz B23_2.fastq.gz B24_1.fastq.gz B24_2.fastq.gz B25_1.fastq.gz B25_2.fastq.gz B26_1.fastq.gz B26_2.fastq.gz B27_1.fastq.gz B27_2.fastq.gz B28_1.fastq.gz B28_2.fastq.gz B29_1.fastq.gz B29_2.fastq.gz B30_1.fastq.gz B30_2.fastq.gz B51_1.fastq.gz B51_2.fastq.gz B52_1.fastq.gz B52_2.fastq.gz B53_1.fastq.gz B53_2.fastq.gz B54_1.fastq.gz B54_2.fastq.gz B55_1.fastq.gz B55_2.fastq.gz B56_1.fastq.gz B56_2.fastq.gz B57_1.fastq.gz B57_2.fastq.gz B58_1.fastq.gz B58_2.fastq.gz B59_1.fastq.gz B59_2.fastq.gz B60_1.fastq.gz B60_2.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang
# (base) kang1234@celia-PowerEdge-T640 Thu Apr 13 15:28:43 ~/white_island/kraken/Orth10_reference
scp Blueeyed* Blue_eyed_info.txt RSEM_align.pl jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/Orth10_reference

# jlkang@hpc2021 Thu Apr 13 15:19:57 /lustre1/g/sbs_schunter/Kang
module load RSEM/1.3.3; module load bowtie2; module load bowtie
mkdir Orth10_reference; cd Orth10_reference
mkdir Blueeyed_RSEM_output; cd Blueeyed_RSEM_output
# jlkang@hpc2021 Thu Apr 13 15:51:29 /lustre1/g/sbs_schunter/Kang/Orth10_reference/Blueeyed_RSEM_output
mv ../RSEM_align.pl ./; mv ../Blue_eyed_info.txt ./; mv Blue_eyed_info.txt Blueeyed_info.txt
perl RSEM_align.pl Blueeyed_info.txt
```

```rsem_1.sh
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B21_1.fastq.gz ../../B21_2.fastq.gz ../Blueeyed B21
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B22_1.fastq.gz ../../B22_2.fastq.gz ../Blueeyed B22
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B23_1.fastq.gz ../../B23_2.fastq.gz ../Blueeyed B23
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B24_1.fastq.gz ../../B24_2.fastq.gz ../Blueeyed B24
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B25_1.fastq.gz ../../B25_2.fastq.gz ../Blueeyed B25
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B26_1.fastq.gz ../../B26_2.fastq.gz ../Blueeyed B26
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B27_1.fastq.gz ../../B27_2.fastq.gz ../Blueeyed B27
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B28_1.fastq.gz ../../B28_2.fastq.gz ../Blueeyed B28
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B29_1.fastq.gz ../../B29_2.fastq.gz ../Blueeyed B29
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B30_1.fastq.gz ../../B30_2.fastq.gz ../Blueeyed B30
```

```rsem_2.sh
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B51_1.fastq.gz ../../B51_2.fastq.gz ../Blueeyed B51
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B52_1.fastq.gz ../../B52_2.fastq.gz ../Blueeyed B52
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B53_1.fastq.gz ../../B53_2.fastq.gz ../Blueeyed B53
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B54_1.fastq.gz ../../B54_2.fastq.gz ../Blueeyed B54
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B55_1.fastq.gz ../../B55_2.fastq.gz ../Blueeyed B55
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B56_1.fastq.gz ../../B56_2.fastq.gz ../Blueeyed B56
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B57_1.fastq.gz ../../B57_2.fastq.gz ../Blueeyed B57
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B58_1.fastq.gz ../../B58_2.fastq.gz ../Blueeyed B58
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B59_1.fastq.gz ../../B59_2.fastq.gz ../Blueeyed B59
rsem-calculate-expression -p 64 --bowtie2 --paired-end ../../B60_1.fastq.gz ../../B60_2.fastq.gz ../Blueeyed B60
```

```script1.cmd
#!/bin/bash
#SBATCH --job-name=rsem-calculate-expression        # 1. Job name
#SBATCH --mail-type=BEGIN,END,FAIL    # 2. Send email upon events (Options: NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlkang@hku.hk     #    Email address to receive notification
#SBATCH --partition=amd               # 3. Request a partition
#SBATCH --qos=normal                  # 4. Request a QoS
#SBATCH --nodes=1                     #    Request number of node(s)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem-per-cpu=1G
#SBATCH --time=7-00:00:00             # 7. Job execution duration limit day-hour:min:sec
#SBATCH --output=%x_%j.out            # 8. Standard output log as $job_name_$job_id.out
#SBATCH --error=%x_%j.err             #    Standard error log as $job_name_$job_id.err

# print the start time
date
sh rsem_1.sh
# print the end time
date
```

```bash
# jlkang@hpc2021 Thu Apr 13 16:08:01 /lustre1/g/sbs_schunter/Kang/Orth10_reference/Blueeyed_RSEM_output
vi script1.cmd
sbatch script1.cmd # Submitted batch job 1106084
sbatch script2.cmd # Submitted batch job 1106085

# jlkang@hpc2021 Thu Apr 13 20:08:50 /lustre1/g/sbs_schunter/Kang/Orth10_reference/
mkdir Common_RSEM_output; cd Common_RSEM_output
# jlkang@hpc2021 Thu Apr 13 20:10:37 /lustre1/g/sbs_schunter/Kang/Orth10_reference/Common_RSEM_output
cp ../Blueeyed_RSEM_output/RSEM_align.pl ./ ; mv ../common_info.txt ./ # divide into four nodes to run
perl RSEM_align.pl common_info.txt
vi rsem_1.sh; vi rsem_2.sh; vi rsem_3.sh; vi rsem_4.sh
cp ../Blueeyed_RSEM_output/script1.cmd ./
cp script1.cmd script2.cmd; cp script1.cmd script3.cmd; cp script1.cmd script4.cmd
vi script1.cmd; vi script2.cmd; vi script3.cmd; vi script4.cmd

sbatch script1.cmd # Submitted batch job 1106166
sbatch script2.cmd # Submitted batch job 1106167
sbatch script3.cmd # Submitted batch job 1106168
sbatch script4.cmd # Submitted batch job 1106170

# jlkang@hpc2021 Thu Apr 13 20:23:17 /lustre1/g/sbs_schunter/Kang/Orth10_reference
mkdir Yaldwyn_RSEM_output; cd Yaldwyn_RSEM_output
# jlkang@hpc2021 Thu Apr 13 20:23:52 /lustre1/g/sbs_schunter/Kang/Orth10_reference/Yaldwyn_RSEM_output
cp ../Blueeyed_RSEM_output/RSEM_align.pl ./ ; mv ../Yaldwyn_info.txt ./
sbatch script1.cmd # Submitted batch job 1106171
sbatch script2.cmd # Submitted batch job 1106172
sbatch script3.cmd # Submitted batch job 1106173
sbatch script4.cmd # Submitted batch job 1106174
```

```bash
# align in my own workstation
# Kang@fishlab3 Wed Apr 12 15:06:36 ~/Desktop/PapueNewGuinea-new/merge_clean
mkdir Orth10_reference; cd Orth10_reference

# Kang@fishlab3 Wed Apr 12 15:08:01 ~/Desktop/PapueNewGuinea-new/merge_clean/Orth10_reference
scp kang1234@147.8.76.177:~/white_island/kraken/Orth10_reference/Acura* ./
scp kang1234@147.8.76.177:~/white_island/kraken/Orth10_reference/Apoly* ./
scp kang1234@147.8.76.177:~/white_island/kraken/Orth10_reference/Ocomp* ./
scp kang1234@147.8.76.177:~/white_island/kraken/Orth10_reference/Daru* ./
scp kang1234@147.8.76.177:~/white_island/kraken/Orth10_reference/Pmol* ./
scp kang1234@147.8.76.177:~/white_island/kraken/Orth10_reference/Padel* ./

scp kang1234@147.8.76.177:~/white_island/kraken/Orth10_reference/RSEM_align.pl ./
# Kang@fishlab3 Wed Apr 12 15:18:56 /media/HDD/white_island/Compevo/genecapture/Total_orf/select
less OG0061615.fas|perl -alne 'if (/\>/){s/\>//;(my $spe)=$_=~/(\D+)\d+/;$hash{$spe}.=$_."_1.fastq.gz," if /Acura/ || /Apoly/ || /Daru/ || /Pmol/ || /Padel/ || /Ocomp/};END{foreach my $key (sort keys %hash){$info=$hash{$key};$info=~s/\,$//;print "$key\t$info"}}' >PNG_ind.info
# Kang@fishlab3 Wed Apr 12 15:19:43 ~/Desktop/PapueNewGuinea-new/merge_clean/Orth10_reference
cp /media/HDD/white_island/Compevo/genecapture/Total_orf/select/PNG_ind.info ./
nohup perl RSEM_align.pl PNG_ind.info > RSEM_process.txt 2>&1 &
# [1] 25470

##################################################################################################
######################################    FINISH    ##############################################
# Kang@fishlab3 Sat Apr 15 10:59:08 ~/Desktop/PapueNewGuinea-new/merge_clean/Orth10_reference
mkdir Total_RSEM_output
cp Acura_RSEM_output/*.genes.results ./Total_RSEM_output
cp Apoly_RSEM_output/*.genes.results ./Total_RSEM_output
cp Daru_RSEM_output/*.genes.results ./Total_RSEM_output
cp Ocomp_RSEM_output/*.genes.results ./Total_RSEM_output
cp Pmol_RSEM_output/*.genes.results ./Total_RSEM_output
cp Padel_RSEM_output/*.genes.results ./Total_RSEM_output

# Kang@fishlab3 Sat Apr 15 11:13:03 ~/Desktop/PapueNewGuinea-new/merge_clean/Orth10_reference/Total_RSEM_output
scp kang1234@147.8.76.177:~/white_island/kraken/Orth10_reference/Blenny_RSEM_output/*.genes.results ./
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/Orth10_reference/Blueeyed_RSEM_output/*.genes.results ./
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/Orth10_reference/Common_RSEM_output/*.genes.results ./
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/Orth10_reference/Yaldwyn_RSEM_output/*.genes.results ./

scp kang1234@147.8.76.177:~/white_island/kraken/Ind_info.txt ./
```

```change_nm.pl
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

my @fq=<*.genes.results>;
my @cmds;
foreach my $fq (@fq) {
        my ($id)=$fq=~/(.*)\.genes\.results/;
        if ($hash2{$id}) {
                my $spe=$hash2{$id}->{'SPE'};
                my $ind=$hash2{$id}->{'IND'};

                my $new=$ind.".genes.results";
                system("mv $fq $new");
        }
}
```

```bash
# Kang@fishlab3 Sat Apr 15 11:44:59 ~/Desktop/PapueNewGuinea-new/merge_clean/Orth10_reference/Total_RSEM_output
perl change_nm.pl
scp kang1234@147.8.76.177:~/CO2-seeps/high_index/paired/kraken/merge/merge_RSEM_frag_counts_single_table.pl ./
perl -e '@files=<*.genes.results>;print "merge_RSEM_frag_counts_single_table.pl";foreach $file(@files){print " $file "};END{print "> total_gene_matrix.txt"}'
mkdir new
```

```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

my @fq=<*.genes.results>;
foreach my $fq (@fq) {
        open FQ, $fq or die "can not open $fq\n";
        my $new="new/$fq";
        open NEW, ">$new" or die "can not create $new\n";
        while (<FQ>) {
                chomp;
                my @a=split;
                if (/^gene_id/) {
                        print NEW "$_\n";
                } else {
                        my $info;
                        for (my $i = 0; $i < @a; $i++) {
                                $a[$i]=~s/\|\|.*$//;
                                $info.=$a[$i]."\t";
                        }
                        $info=~s/\s+$//;
                        print NEW "$info\n";
                }
        }
        close FQ; close NEW;
}
```

```bash
# Kang@fishlab3 Sat Apr 15 12:21:59 ~/Desktop/PapueNewGuinea-new/merge_clean/Orth10_reference/Total_RSEM_output
mkdir new; cd new/
# Kang@fishlab3 Sat Apr 15 17:55:29 ~/Desktop/PapueNewGuinea-new/merge_clean/Orth10_reference/Total_RSEM_output/new
RNAnorm --method TPM --norm --prefix Orth10 *genes.results
# output: Orth10.TPM.raw.matrix; Orth10.TPM.TMM.matrix; Orth10.TPM.TMM.sqrt.matrix (will use Orth10.TPM.TMM.matrix)
```

## run bamm
```bash
# Kang@fishlab3 Sun Apr 16 18:00:48 /media/HDD/white_island/Compevo/Orth_ten_CO2
mkdir bamm; cd bamm
# Kang@fishlab3 Sun Apr 16 18:01:48 /media/HDD/white_island/Compevo/Orth_ten_CO2/bamm
cp ~/Desktop/PapueNewGuinea-new/merge_clean/Orth10_reference/Total_RSEM_output/new/Orth10.TPM.TMM.sqrt.matrix ./
cp /media/HDD/white_island/Compevo/genecapture/Total_orf/select/r8s_ultrametric.txt ./
cp ~/Desktop/PapueNewGuinea-new/merge_clean/gene_capture/result/Total_orf/Final_ind_2/bamm/temp1.pl ./

perl temp1.pl OG0015225_OG1 > OG0015225_OG1.txt # CAH10
mkdir OG0015225_OG1_trait; mv OG0015225_OG1.txt OG0015225_OG1_trait/
cd OG0015225_OG1_trait/
cp ../r8s_ultrametric.txt ./
cp ~/Desktop/PapueNewGuinea-new/merge_clean/gene_capture/result/Total_orf/Final_ind_2/bamm/OG0002874_trait/traitcontrol.txt ./
vi traitcontrol.txt
bamm -c traitcontrol.txt

# Plot
# kangjingliang@kangjingliangdeMacBook-Pro 日  4 16 18:53:08 ~/Documents/2023/WI/bamm
less label.txt|perl -alne 'my $info; for ($i=0; $i<@F-1;$i++){$info.="which(tree\$tip.label == \"$F[$i]\")\, "};$info=~s/\,\s+$//;print $info'|less

# Kang@fishlab3 Sun Apr 16 19:37:11 /media/HDD/white_island/Compevo/Orth_ten_CO2/OrthoFinder/sub_orth
less sub_orth_genecount_final.txt|perl -alne 'print if $F[-1] eq S12A5'
# OG0007634_OG0	2	3	2	1	1	2	1	1	2	1	S12A5
# Kang@fishlab3 Sun Apr 16 19:32:50 /media/HDD/white_island/Compevo/Orth_ten_CO2/bamm
perl temp2.pl OG0007634_OG0
cd OG0007634_OG0_trait/
bamm -c traitcontrol.txt

# OG0001269_OG6: BMAL1
perl temp2.pl OG0001269_OG6; cd OG0001269_OG6_trait/; bamm -c traitcontrol.txt
# OG0005083_OG0: BMAL2
perl temp2.pl OG0005083_OG0; cd OG0005083_OG0_trait/; bamm -c traitcontrol.txt
# OG0014481_OG1: OPN5
perl temp2.pl OG0014481_OG1; cd OG0014481_OG1_trait/; bamm -c traitcontrol.txt
# OG0015590_OG0: OPN4A
perl temp2.pl OG0015590_OG0; cd OG0015590_OG0_trait/; bamm -c traitcontrol.txt
# OG0008973_OG3: CAH1
perl temp2.pl OG0008973_OG3; cd OG0008973_OG3_trait/; bamm -c traitcontrol.txt
# OG0038714_OG0: VATF
perl temp2.pl OG0038714_OG0; cd OG0038714_OG0_trait/; bamm -c traitcontrol.txt
```

## Only select individuals from white island
```bash
# Kang@fishlab3 Mon Apr 17 16:06:15 /media/HDD/white_island/Compevo/genecapture/Total_orf
vi WI_inds.txt
extract_reads_nb --matrix Capture_info.txt --samples WI_inds.txt > Capture_info_WI.txt

# Select the genes captured all individuals: 37 individuals
less Capture_info_WI.txt|perl -alne '$j++;@head=@F if $j==1;my $a;for (my $i = 1; $i < @F; $i++){$a.=$head[$i].";" if $F[$i]==0};print "$F[0]\t$a" if $j>1'|perl -alne 'print "$F[0]" if @F==1'
mkdir select_WI
perl select_WI.pl; cd select_WI
cp ../select/temp1.pl ./ ; perl temp1.pl
# Kang@fishlab3 Mon Apr 17 17:22:16 /media/HDD/white_island/Compevo/genecapture/Total_orf/select_WI
fasta2phy.pl All_gene_concatenated.fasta >All_gene_concatenated_orth10WI.phy
scp All_gene_concatenated_orth10WI.phy kang1234@147.8.76.177:~/white_island/kraken/gene_capture_PNG/
nohup raxmlHPC -f a -x 12345 -p 12345 -# 100 -m GTRCAT -T 24 -s All_gene_concatenated_orth10WI.phy -o Blenny2,Blenny5,Blenny3,Blenny12,Blenny16,Blenny6,Blenny14,Blenny10,Blenny7,Blenny9,Blenny13,Blenny4,Blenny15,Blenny1,Blenny11,Blenny8 -n All_gene_concatenated_orth10WI > Raxml.process 2>&1 &
# [1] 23051
scp kang1234@147.8.76.177:~/white_island/kraken/gene_capture_PNG/RAxML_bestTree.All_gene_concatenated_orth10WI ./
cp /media/HDD/cleaner_fish/genome/gene_family/cafetutorial_prep_r8s.py ./
less All_gene_concatenated.fasta |perl -alne 'if (/>/){s/\>//;print if /Blenny/ || /Common/}'|sort|perl -alne '$info.=$_.",";END{$info=~s/\,$//;print $info}'
python cafetutorial_prep_r8s.py -i RAxML_bestTree.All_gene_concatenated_orth10WI -o r8s_ctl_file_1.txt -s 638853 -p 'Blenny1,Blenny10,Blenny11,Blenny12,Blenny13,Blenny14,Blenny15,Blenny16,Blenny2,Blenny3,Blenny4,Blenny5,Blenny6,Blenny7,Blenny8,Blenny9,Common1,Common10,Common11,Common12,Common13,Common14,Common15,Common2,Common3,Common4,Common5,Common6,Common7,Common8,Common9' -c '67.17'
nohup r8s -b -f r8s_ctl_file_1.txt > r8s_tmp_1.txt 2>&1 &
# [1] 10279
tail -n 1 r8s_tmp_1.txt | cut -c 16- > r8s_ultrametric.txt
# "r8s_ultrametric.txt" will be the input phylogeny of bamm for the WI individuals

# Kang@fishlab3 Mon Apr 17 20:07:24 /media/HDD/white_island/Compevo/Orth_ten_CO2
mkdir bamm_WI; cd bamm_WI

# (base) kang1234@celia-PowerEdge-T640 Mon Apr 17 20:15:59 ~/white_island/kraken/final_reference
mkdir all_RSEM_output
cp Blenny_RSEM_output/*.genes.results all_RSEM_output/
cp Blue_eyed_RSEM_output/*.genes.results all_RSEM_output/
cp Common_RSEM_output/*.genes.results all_RSEM_output/
cp Yaldwyn_RSEM_output/*.genes.results all_RSEM_output/

# Kang@fishlab3 Mon Apr 17 20:31:49 /media/HDD/white_island/Compevo/Orth_ten_CO2/bamm_WI
scp change_nm.pl  Ind_info.txt kang1234@147.8.76.177:~/white_island/kraken/final_reference/all_RSEM_output
# (base) kang1234@celia-PowerEdge-T640 Mon Apr 17 20:26:49 ~/white_island/kraken/final_reference/all_RSEM_output
RNAnorm --method TPM --norm --prefix Orth4 *.genes.results
# Kang@fishlab3 Mon Apr 17 20:39:22 /media/HDD/white_island/Compevo/Orth_ten_CO2/bamm_WI
scp kang1234@147.8.76.177:~/white_island/kraken/final_reference/all_RSEM_output/Orth4.TPM.TMM.sqrt.matrix ./
cp /media/HDD/white_island/Compevo/genecapture/Total_orf/select_WI/r8s_ultrametric.txt ./
cp /media/HDD/white_island/Compevo/Orth_ten_CO2/bamm/*.pl ./

# kangjingliang@kangjingliangdeMacBook-Pro 一  4 17 20:43:06 ~/Documents/2023/WI/DEGs/Enrichment
less CR_DEGs_regulation_coreCR.txt|grep -i 'clock'|sort -u
# Common	circadian rhythm	OG0008988	sp|Q6YGZ4|CLOCK_TYTAL	Circadian locomoter output cycles protein kaput	Common_Cs_Vn_up;Common_Vn_Vs_dw;	CLOCK_2
# Common	circadian rhythm	OG0012851	sp|Q6YGZ4|CLOCK_TYTAL	Circadian locomoter output cycles protein kaput	Common_Cs_Vn_up;	CLOCK_1
# Common	circadian rhythm	OG0013172	sp|Q8QGQ6|CLOCK_CHICK	Circadian locomoter output cycles protein kaput	Common_Cs_Vn_up;Common_Vn_Vs_dw;	CLOCK_4
# Common	circadian rhythm	OG0024759	sp|Q8QGQ6|CLOCK_CHICK	Circadian locomoter output cycles protein kaput	Common_Cs_Vn_up;Common_Vn_Vs_dw;	CLOCK_3

# OG0012851: CLOCK_1
# Kang@fishlab3 Mon Apr 17 20:46:17 /media/HDD/white_island/Compevo/Orth_ten_CO2/bamm_WI
perl temp1.pl OG0012851 > OG0012851.txt
mkdir OG0012851_trait; mv OG0012851.txt OG0012851_trait/
cd OG0012851_trait/
cp ../r8s_ultrametric.txt ./
cp ~/Desktop/PapueNewGuinea-new/merge_clean/gene_capture/result/Total_orf/Final_ind_2/bamm/OG0002874_trait/traitcontrol.txt ./
vi traitcontrol.txt
bamm -c traitcontrol.txt

less label.txt|perl -alne 'my $info; for ($i=0; $i<@F-1;$i++){$info.="which(tree\$tip.label == \"$F[$i]\")\, "};$info=~s/\,\s+$//;print $info'|less

# OG0008988: CLOCK_2
perl temp2.pl OG0008988; cd OG0008988_trait/; bamm -c traitcontrol.txt

# OG0024759: CLOCK_3
perl temp2.pl OG0024759; cd OG0024759_trait/; bamm -c traitcontrol.txt

# OG0013172: CLOCK_4
perl temp2.pl OG0013172; cd OG0013172_trait/; bamm -c traitcontrol.txt

# OG0017066: S4A8
perl temp2.pl OG0017066; cd OG0017066_trait/; bamm -c traitcontrol.txt

# OG0085219: Rhodopsin
perl temp2.pl OG0085219; cd OG0085219_trait/; bamm -c traitcontrol.txt

# OG0005900: Melatonin receptor type 1B-B (MR1BB)
perl temp2.pl OG0005900; cd OG0005900_trait/; bamm -c traitcontrol.txt

# OG0074362: CAH12
perl temp2.pl OG0074362; cd OG0074362_trait/; bamm -c traitcontrol.txt

# OG0072207: VATF
perl temp2.pl OG0072207; cd OG0072207_trait/; bamm -c traitcontrol.txt
```

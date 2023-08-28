## WI: Mutation rate calculation
```bash
# Kang@fishlab3 Mon Aug 07 15:55:54 /media/HDD/white_island/Compevo/orth16_new/Orth_conca
less Orth_conca_paml.phy|perl -alne 'if (/^\d+/){print}else{my $a=sprintf("%-10s",$F[0]);print "$a\t$F[1]"}' > Orth_conca_paml_2.phy
phylip dnadist # calculate based on Jukes-Cantor
# kangjingliang@kangjingliangdeMacBook-Pro 一  8 07 16:35:26 ~/Desktop
scp Kang@147.8.76.229:/media/HDD/white_island/Compevo/orth16_new/Orth_conca/outfile ./Distance.txt
less Distance.txt|perl -alne 's/\s+/\t/g;print' >Distance.txt.1
mv Distance.txt.1 Distance.txt

# Use the average mutation rate of fish
# fishes (average of all species 5.97e−9, 95% CI of the mean 4.39e10−9 to 7.55e10−9)
# transfer the fastq file to hpc: control+z; bg
# (base) kang1234@celia-PowerEdge-T640 Mon Aug 07 18:34:02 ~/white_island/kraken
nohup scp Blenny*_1.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/WI > Blenny_1.out 2>&1
nohup scp Blenny*_2.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/WI > Blenny_2.out 2>&1
nohup scp Blueeyed*_1.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/WI > Blueeyed_1.out 2>&1
nohup scp Blueeyed*_2.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/WI > Blueeyed_2.out 2>&1
nohup scp Common*_1.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/WI > Common_1.out 2>&1
nohup scp Common*_2.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/WI > Common_2.out 2>&1
nohup scp Yaldwyn*_1.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/WI > Yaldwyn_1.out 2>&1
nohup scp Yaldwyn*_2.fastq.gz jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/WI > Yaldwyn_2.out 2>&1

# use the de novo assembly as reference for mapping
# transfer the de novo assembly from DRAP to HPC
# (base) kang1234@celia-PowerEdge-T640 Mon Aug 07 18:50:37 ~/white_island/kraken/Blenny_runDrap/e-rmbt_editing
scp all_contigs.second_pass.fa jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/WI/Blenny_denovo.fa
# (base) kang1234@celia-PowerEdge-T640 Mon Aug 07 18:53:29 ~/white_island/kraken/Blue_eyed_runDrap/e-rmbt_editing
scp all_contigs.second_pass.fa jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/WI/Blueeyed_denovo.fa
# (base) kang1234@celia-PowerEdge-T640 Mon Aug 07 18:55:47 ~/white_island/kraken/Common_runDrap/e-rmbt_editing
scp all_contigs.second_pass.fa jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/WI/Common_denovo.fa
# (base) kang1234@celia-PowerEdge-T640 Mon Aug 07 18:59:17 ~/white_island/kraken/Yaldwyn_runDrap/e-rmbt_editing
scp all_contigs.second_pass.fa jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/WI/Yaldwyn_denovo.fa

module load perl/5.34.0 perl-lib/5.34.0
module load bowtie2/2.4.4
```

```script_index_bowtie2.cmd
# build the index for each species
#!/bin/bash
#SBATCH --job-name=Index_bowtie2        # 1. Job name
#SBATCH --mail-type=BEGIN,END,FAIL    # 2. Send email upon events (Options: NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlkang@hku.hk     #    Email address to receive notification
#SBATCH --partition=amd               # 3. Request a partition
#SBATCH --qos=normal                  # 4. Request a QoS
#SBATCH --nodes=1                     #    Request number of node(s)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --mem-per-cpu=1G
#SBATCH --time=7-00:00:00             # 7. Job execution duration limit day-hour:min:sec
#SBATCH --output=%x_%j.out            # 8. Standard output log as $job_name_$job_id.out
#SBATCH --error=%x_%j.err             #    Standard error log as $job_name_$job_id.err

# print the start time
date
bowtie2-build Blenny_denovo.fa Blenny
bowtie2-build Blueeyed_denovo.fa Blueeyed
bowtie2-build Common_denovo.fa Common
bowtie2-build Yaldwyn_denovo.fa Yaldwyn
# print the end time
date
```

```bash
# jlkang@hpc2021 Mon Aug 07 22:25:38 /lustre1/g/sbs_schunter/Kang/WI
sbatch script_index_bowtie2.cmd
# Submitted batch job 1368203s
module load samtools/1.14
# bowtie2 -x Blenny -1 Blenny10_1.fastq.gz -2 Blenny10_2.fastq.gz -S Blenny10.sam -p 32 # 1_align
# samtools view -bS -@ 32 Blenny10.sam -o Blenny10.bam # sam to bam
# samtools sort -@ 32 Blenny10.bam -o Blenny10_sorted.bam # sorted bam
# samtools index -@ 32 Blenny10_sorted.bam # index bam
# rm -rf Blenny10.sam Blenny10.bam
```

```bowtie2_align.pl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Parallel::ForkManager;

# bowtie2 -x Blenny -1 Blenny10_1.fastq.gz -2 Blenny10_2.fastq.gz -S Blenny10.sam -p 32 # align
# samtools view -bS -@ 32 Blenny10.sam -o Blenny10.bam # sam to bam
# samtools sort -@ 32 Blenny10.bam -o Blenny10_sorted.bam # sorted bam
# samtools index -@ 32 Blenny10_sorted.bam # index bam
# rm -rf Blenny10.sam Blenny10.bam

my $spe =$ARGV[0];
my @fqs=<$spe*_1.fastq.gz>;
foreach my $fq (@fqs) {
	(my $nm)=$fq=~/(.*)_1\.fastq\.gz/;
	my ($fqr1, $fqr2);
	$fqr1=$fq;
	$fqr2=$nm."_2.fastq.gz";
	my $sam=$nm.".sam";
	my $bam=$nm.".bam";
	my $sob=$nm."_sorted.bam";
	# algin
	my $cmd1="bowtie2 -x $spe -1 $fqr1 -2 $fqr2 -S $sam -p 32";
	# sam => bam
	my $cmd2="samtools view -bS -@ 32 $sam -o $bam";
	# sort
	my $cmd3="samtools sort -@ 32 $bam -o $sob";
	# index
	my $cmd4="samtools index -@ 32 $sob";
	# delete
	my $cmd5="rm $sam $bam";
#	print "$cmd1\n$cmd2\n$cmd3\n$cmd4\n$cmd5\n";
	system($cmd1);
	system($cmd2);
	system($cmd3);
	system($cmd4);
	system($cmd5);
}
```

```script_align_bowtie2_Blenny.cmd
#!/bin/bash
#SBATCH --job-name=Align_bowtie2_Blenny        # 1. Job name
#SBATCH --mail-type=BEGIN,END,FAIL    # 2. Send email upon events (Options: NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jlkang@hku.hk     #    Email address to receive notification
#SBATCH --partition=amd               # 3. Request a partition
#SBATCH --qos=normal                  # 4. Request a QoS
#SBATCH --nodes=1                     #    Request number of node(s)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1G
#SBATCH --time=7-00:00:00             # 7. Job execution duration limit day-hour:min:sec
#SBATCH --output=%x_%j.out            # 8. Standard output log as $job_name_$job_id.out
#SBATCH --error=%x_%j.err             #    Standard error log as $job_name_$job_id.err

# print the start time
date
perl bowtie2_align.pl Blenny
# print the end time
date
```

```bash
sbatch script_align_bowtie2_Blenny.cmd
sbatch script_align_bowtie2_Blue.cmd
sbatch script_align_bowtie2_Common.cmd
sbatch script_align_bowtie2_Yaldwyn.cmd

# psmc
# jlkang@hpc2021 Wed Aug 09 10:55:34 /lustre1/g/sbs_schunter/Kang/WI
module load psmc
module load bamtools/2.5.2
module load perl/5.34.0 perl-lib/5.34.0
module load bowtie2/2.4.4
module load samtools/1.14
module load gnuplot
module load bcftools
# jlkang@hpc2021 Wed Aug 09 14:20:40 /lustre1/g/sbs_schunter/Kang/WI
sbatch script_psmc_Yaldwyn.cmd # perl To_psmc.pl Yaldwyn
```

```To_psmc.pl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Parallel::ForkManager;

# get "*.psmc"

# samtools mpileup -C50 -uf ./Apoly.fa Apoly10.split.bam|bcftools call -c - |vcfutils.pl vcf2fq -d 10 -D 100 -Q 30|gzip > Apoly10.split.fq.gz
# fq2psmcfa -q20 Apoly10.split.fq.gz >Apoly10.split.psmc.fa 
# psmc -N50 -t15 -r5 -p "4+25*2+4+6" -o Apoly10.split.psmc Apoly10.split.psmcfa

my $spe=$ARGV[0];
my @cmds;
my $ref=$spe."_denovo.fa";
my @bams=<$spe*_sorted.bam>;
foreach my $bam (@bams) {
	(my $ind)=$bam=~/(.*)\_sorted\.bam/;
	my $fq=$ind.".fq.gz";
	my $cmd1="bcftools mpileup -C50 -f $ref $bam \| bcftools call -c - \| vcfutils.pl vcf2fq -d 10 -D 100 -Q 30 \| gzip > $fq";
	my $psmcfa=$ind.".psmc.fa";
	my $cmd2="fq2psmcfa -q20 $fq > $psmcfa";

	my $psmc=$ind.".psmc";
	my $cmd3="psmc -N50 -t15 -r5 -p \"4+25*2+4+6\" -o $psmc $psmcfa";
	my $cmd =$cmd1.";".$cmd2.";".$cmd3;
	push @cmds, $cmd;
#	print "$cmd\n";
#	print "$cmd1\n$cmd2\n$cmd3\n";
#	system($cmd1);
#	system($cmd2);
#	system($cmd3);
}

my $manager = new Parallel::ForkManager(4);
foreach my $cmd (@cmds) {
    $manager->start and next;
    system($cmd);
    $manager->finish;
}
$manager -> wait_all_children;
```

```bash
# my $cmd2="fq2psmcfa -q20 $fq > $psmcfa";
# my $cmd3="psmc -N50 -t15 -r5 -p \"4+25*2+4+6\" -o $psmc $psmcfa";
# Kang@fishlab3 Thu Aug 10 19:52:46 /media/HDD/white_island
mkdir psmc; cd psmc
# Kang@fishlab3 Thu Aug 10 20:09:50 /media/HDD/white_island/psmc
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/WI/Blenny*.fq.gz ./
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/WI/Blueeyed*.fq.gz ./
scp jlkang@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/WI/Yaldwyn*.fq.gz ./
```

```psmc.pl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Parallel::ForkManager;

# get "*.psmc"

# samtools mpileup -C50 -uf ./Apoly.fa Apoly10.split.bam|bcftools call -c - |vcfutils.pl vcf2fq -d 10 -D 100 -Q 30|gzip > Apoly10.split.fq.gz
# fq2psmcfa -q20 Apoly10.split.fq.gz >Apoly10.split.psmc.fa 
# psmc -N50 -t15 -r5 -p "4+25*2+4+6" -o Apoly10.split.psmc Apoly10.split.psmcfa

my @cmds;
my @fqs=<*.fq.gz>;
foreach my $fq (@fqs) {
	(my $ind)=$fq=~/(.*)\.fq\.gz/;
	my $psmcfa=$ind.".psmc.fa";
	next if (-z $psmcfa);
	my $psmc=$ind.".psmc";
	#my $cmd1="bcftools mpileup -C50 -f $ref $bam \| bcftools call -c - \| vcfutils.pl vcf2fq -d 10 -D 100 -Q 30 \| gzip > $fq";
	my $cmd2="~/software/psmc/utils/fq2psmcfa -q20 $fq > $psmcfa";
	my $cmd3="~/software/psmc/psmc -N50 -t15 -r5 -p \"4+25*2+4+6\" -o $psmc $psmcfa";
	system($cmd2);
	system($cmd3);
}
```

```bash
# Kang@fishlab3 Thu Aug 10 20:35:48 /media/HDD/white_island/psmc
nohup perl psmc.pl > psmc.process 2>&1 &
# [1] 13118

# plot
# 1. Concatenate the individual .psmc files you wish to plot
# Kang@fishlab3 Thu Aug 10 21:07:55 /media/HDD/white_island/psmc
cat Blueeyed*.psmc > Blueeyed_combined.psmc
# 2. Plot using multiline mode
ll Blueeyed*.fq.gz|perl -alne '(my $ind)=$F[-1]=~/(.*)\.fq\.gz/;$info.=$ind.",";END{$info=~s/\,$//;print $info}'
# Blueeyed10,Blueeyed11,Blueeyed12,Blueeyed13,Blueeyed14,Blueeyed15,Blueeyed1,Blueeyed2,Blueeyed3,Blueeyed4,Blueeyed5,Blueeyed6,Blueeyed7,Blueeyed8
~/software/psmc/utils/psmc_plot.pl -p -u 7.37e-09 -g 1 -M 'Blueeyed10,Blueeyed11,Blueeyed12,Blueeyed13,Blueeyed14,Blueeyed15,Blueeyed1,Blueeyed2,Blueeyed3,Blueeyed4,Blueeyed5,Blueeyed6,Blueeyed7,Blueeyed8' Blueeyed_combined Blueeyed_combined.psmc
~/software/psmc/utils/psmc_plot.pl -p -u 7.37e-09 -g 1 Blueeyed_combined Blueeyed_combined.psmc

cat Blenny*.psmc > Blenny_combined.psmc
psmc_plot.pl -p -u 7.37e-09 -g 1 -Y 100 Blenny_combined Blenny_combined.psmc
cat Blenny1.psmc Blenny3.psmc Blenny6.psmc Blenny8.psmc Blenny9.psmc Blenny15.psmc > Blenny_combined.psmc
~/software/psmc/utils/psmc_plot.pl -p -u 7.37e-09 -g 1 -Y 100 Blenny_combined Blenny_combined.psmc

cat Yaldwyn*.psmc > Yaldwyn_combined.psmc
~/software/psmc/utils/psmc_plot.pl -p -u 7.37e-09 -g 1 -Y 100 Yaldwyn_combined Yaldwyn_combined.psmc
# Kang@fishlab3 Fri Aug 11 10:43:09 /media/HDD/white_island/psmc
ll Yaldwyn*.fq.gz|perl -alne '(my $ind)=$F[-1]=~/(.*)\.fq\.gz/;$info.=$ind.".psmc ";END{$info=~s/\,$//;print $info}'
# Yaldwyn10.psmc Yaldwyn11.psmc Yaldwyn12.psmc Yaldwyn13.psmc Yaldwyn14.psmc Yaldwyn15.psmc Yaldwyn1.psmc Yaldwyn2.psmc Yaldwyn3.psmc Yaldwyn4.psmc Yaldwyn5.psmc Yaldwyn6.psmc Yaldwyn7.psmc Yaldwyn8.psmc
# don't include Yaldwyn6.psmc
cat Yaldwyn10.psmc Yaldwyn11.psmc Yaldwyn12.psmc Yaldwyn13.psmc Yaldwyn14.psmc Yaldwyn15.psmc Yaldwyn1.psmc Yaldwyn2.psmc Yaldwyn3.psmc Yaldwyn4.psmc Yaldwyn5.psmc Yaldwyn7.psmc Yaldwyn8.psmc > Yaldwyn_combined.psmc
~/software/psmc/utils/psmc_plot.pl -p -u 7.37e-09 -g 1 Yaldwyn_combined Yaldwyn_combined.psmc

cat Common*.psmc > Common_combined.psmc
~/software/psmc/utils/psmc_plot.pl -p -u 7.37e-09 -g 0.5 Common_combined Common_combined.psmc

# plot for each individual to make sure which one is outlier
perl plot.pl
```

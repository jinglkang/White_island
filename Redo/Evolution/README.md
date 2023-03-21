# White Island
## Evolutionary rate comparisons between fish species exposed to elevated pCO2 using RNA-seq
### 1. Download the data
```bash
# (base) kang1234@celia-PowerEdge-T640 Thu Mar 16 16:13:08 ~/white_island
mkdir Compevo; cd Compevo

# Download the genome sequences
# 1. Dicentrarchus labrax (European seabass); PRJEB43042; GCA_905237075.1 (latest); genome assembly; RefSeq
# 2. Tripterygion delaisi (Black-faced blenny): PRJNA186408; GAJK00000000.1; RNA-seq; Genebank;
# 3. Salaria pavo (Peacock blenny): PRJNA329073; GenBank: GEVH01.1.gbff.gz; RNA-seq;
# 4. Lates calcarifer (Asian Seabass): genome; 
# 5. Gadus morhua (Atlantic cod): genome;
# 6. Amphiprion percula (orange clownfish): genebank; genome;
# 7. Trematomus bernacchii (Emerald rockcod): genebank; genome;
# 8. Larimichthys crocea (Large yellow croaker): genebank; genome;
# 9. Oncorhynchus kisutch (coho salmon): PRJNA352719; GenBank assembly accession (GCA_002021735.2; latest); genome

# Download the gtf
# 1. Dicentrarchus_labrax => GCF_905237075.1-RS_2023_02
# 2. 
# 3. 
# 4. Lates_calcarifer => GCF_001640805.2-RS_2023_02
# 5. Gadus_morhua => GCF_902167405.1_gadMor3.0
# 6. 
# 7. Trematomus_bernacchii => GCF_902827165.1_fTreBer1.1
# 8. Larimichthys_crocea => GCF_000972845.2_L_crocea_2.0
# 9. Oncorhynchus_kisutch => GCF_002021735.2_Okis_V2
```
### 2. Run TransDecoder for ORFs per species
```Run_TransDecoder.pl
#!/usr/bin/perl
use strict;
use warnings;

# nohup TransDecoder.LongOrfs -t Acura_tra_nuc.fa -O Acura_orf >Acura_transdecoder.process 2>&1 &
my $pref="pre_TransDecoder.txt";
open PREF, $pref or die "can not open $pref\n";
while (<PREF>) {
	chomp;
	my @a=split;
	my ($nm1, $nm2)=($a[0], $a[1]);
	#my $proc=$nm2."_transdecoder.process";
	#$nm2.="_orf";
	my $cmd="TransDecoder.LongOrfs -t $nm1";
	#print "$cmd\n";
	system($cmd);
}
```

```bash
# nohup TransDecoder.LongOrfs -t Acura_tra_nuc.fa -O Acura_orf >Acura_transdecoder.process 2>&1 &
# (base) kang1234@celia-PowerEdge-T640 Fri Mar 17 15:07:29 ~/white_island/Compevo
nohup perl Run_TransDecoder.pl >TransDecoder.LongOrfs.process 2>&1 &
# [1] 23387
```
### 3. OrthoFinder for orthologous genes detection
**The protein sequences of ORFs detected by OrthoFinder as input**   
```pre_OrthoFinder.pl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

my $file="pre_TransDecoder.txt";
my @dirs; my %hash;
open FILE, $file or die "can not open $file\n";
while (<FILE>) {
	chomp;
	my @a=split;
	my $file_old=$a[0].".transdecoder_dir/longest_orfs.pep";
	my $file_new="Input_pep/".$a[1].".fasta";
	open OLD, $file_old or die "can not open $file_old\n";
	open NEW, ">$file_new" or die "can not create $file_new\n";
	my $i;
	while (<OLD>) {
		chomp;
		if (/>/) {
			$i++;
			my $gene=$a[1]."_".$i;
			print NEW ">$gene\n";
		} else {
			print NEW "$_\n";
		}
	}
}
```

```bash
# (base) kang1234@celia-PowerEdge-T640 Fri Mar 17 17:11:59 ~/white_island/Compevo
mkdir Input_pep
perl pre_OrthoFinder.pl

# Time for the six species in gcb
# Kang@fishlab3 Fri Mar 17 18:17:05 ~/Desktop/PapueNewGuinea-new/orthologue/orthofinder_input_pep/all_input
scp *.fasta kang1234@147.8.76.177:~/white_island/Compevo/Input_pep

# Time for the four species in this study
# Kang@fishlab3 Fri Mar 17 18:21:06 /media/HDD/white_island/orthologue/input_pep/orthofinder_input
scp *.fa kang1234@147.8.76.177:~/white_island/Compevo/Input_pep

# Run OrthoFinder
# Transfer to my own workstation to run OrthoFinder
# Kang@fishlab3 Fri Mar 17 18:28:06 /media/HDD/white_island
mkdir Compevo
scp kang1234@147.8.76.177:~/white_island/Compevo/Input_pep/*fasta ./Compevo
# Kang@fishlab3 Fri Mar 17 18:30:54 /media/HDD/white_island/Compevo
nohup orthofinder -f ./ -a 20 >orthofinder-process 2>&1 &
# [1] 26906

# Crush
# (base) kang1234@celia-PowerEdge-T640 Mon Mar 20 09:16:09 ~/white_island/Compevo/Input_pep
nohup orthofinder -f ./ -a 20 >orthofinder-process 2>&1 &
# [1] 10342

# add another six ref-species into the orthologous detection
# (base) kang1234@celia-PowerEdge-T640 Mon Mar 20 23:04:23 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar20
mkdir Extra_species
# Kang@fishlab3 Mon Mar 20 23:05:34 ~/Desktop/PapueNewGuinea-new/longest_pep/input_pep
scp Fugu.fasta Medaka.fasta Platyfish.fasta Spottedgar.fasta Stickleback.fasta Zebrafish.fasta kang1234@147.8.76.177:~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar20/Extra_species
# (base) kang1234@celia-PowerEdge-T640 Mon Mar 20 23:06:59 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar20
nohup orthofinder -b WorkingDirectory/ -f Extra_species >Extra_orthofinder.process 2>&1 &
# [1] 20768

# Each species only keep one sequence per orthogroup; for phylogenetic tree construction
# (base) kang1234@celia-PowerEdge-T640 Tue Mar 21 11:17:38 ~/white_island/Compevo/Input_pep/OrthoFinder/Results_Mar20/WorkingDirectory/OrthoFinder/Results_Mar20
perl Orth_1spe_1seq.pl
```

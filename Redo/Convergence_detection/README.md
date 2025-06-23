# Detect the convergence selection in the common triplefin and crested blenny
## 1. translate the nucleotide sequences of all orthologous genes "in final_orth_input_paml.txt"
```perl
# translate_OG_pep.pl
#!/usr/bin/perl
use strict;
use warnings;

my $orth="final_orth_input_paml.txt";
open ORTH, $orth or die "can not open $orth\n";
while (<ORTH>) {
	chomp;
	my $orthdir="paml_files/".$_;
	chdir "$orthdir";
	my $cmd1="translateDna.pl -i final_alignment.fa > final_alignment_pep.fa";
	system($cmd1);
	chdir "/data2/jlkang/White_island/paml_input";
}
```

```bash
# (base) jlkang@hnu2024 Mon May 26 2025 15:35:14 /data2/jlkang/White_island/paml_input
perl translate_OG_pep.pl
```

```perl
# HKU HPC
# translate_OG_pep.pl
#!/usr/bin/perl
use strict;
use warnings;

my $orth="final_orth_input_paml.txt";
open ORTH, $orth or die "can not open $orth\n";
while (<ORTH>) {
	chomp;
	my $orthdir="paml_files/".$_;
	chdir "$orthdir";
	my $cmd1="/lustre1/g/sbs_schunter/Kang/orth16_new/translateDna.pl -i final_alignment.fa > final_alignment_pep.fa";
	system($cmd1);
	chdir "/lustre1/g/sbs_schunter/Kang/orth16_new/";
}
```

```bash
# (base) jlkang@hnu2024 Mon May 26 2025 15:35:14 /data2/jlkang/White_island/paml_input
perl translate_OG_pep.pl

# HKU HPC
# (base) romeodan@hpc2021 Mon May 26 19:49:26 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files
perl translate_OG_pep.pl
```

## 2. detect the putative convergent sites
```perl
#!/usr/bin/perl
# Detect_Nons.pl
use strict;
use warnings;
use Array::Utils qw(:all);

my %seq; my $name;
my $orth=$ARGV[0];
my $fas="paml_files/$orth/final_alignment_pep.fa";
open FAS, $fas or die "can not open $fas\n";
while (<FAS>) {
        chomp;
        if (/^>/) {
                s/\>//;
                $name=$_;
        } else {
                $seq{$name}.=$_;
        }
}

my %cleaner=(
        'Blenny'=> 1,
        'Common'=> 1,
        );

# compare the nonsynonymous position pep sequences one by one
my %hash1;
my @poss;
my @nocls=qw(Acura Apoly Blueeyed Daru Fugu Medaka Ocomp Padel Platyfish Pmol Spottedgar Stickleback Yaldwyn Zebrafish);
my @cleas=qw(Blenny Common);
my @aspes=qw(Acura Apoly Blenny Blueeyed Common Daru Fugu Medaka Ocomp Padel Platyfish Pmol Spottedgar Stickleback Yaldwyn Zebrafish);

&Build_pos_hash(\@aspes);

sub Build_pos_hash {
        my ($grp)=@_;
        my @grp=@{$grp};
        my $len;
        foreach my $spe (@grp) {
                my $seq=$seq{$spe};
                $len=length($seq);
                for (my $i = 0; $i < $len; $i++) {
                        my $spepos=substr($seq,$i,1);
                        $hash1{$spe}->{$i}=$spepos;
                }
        }

        for (my $i = 0; $i < $len; $i++) {
                my (%hash2,%hash3);
                my $pos=$i;
                my $newp=$pos+1;
                my $info=$newp.":";
                foreach my $spe (@aspes) {
                        my $spepos=$hash1{$spe}->{$pos};
                        $info.=$spe."($spepos);";
                }

                my (@cleas_pos, @nocls_pos);
                foreach my $spe (@cleas) {
                        my $spepos=$hash1{$spe}->{$pos};
                        $hash2{$spepos}++;
                        push @cleas_pos, $spepos;
                }
                foreach my $spe (@nocls) {
                        my $spepos=$hash1{$spe}->{$pos};
                        $hash3{$spepos}++;
                        push @nocls_pos, $spepos;
                }
                #my %array1 = map { $_ => 1 } @cleas_pos;
				#my @isect = grep { $array1{$_} } @nocls_pos;
                my @isect = intersect(@cleas_pos, @nocls_pos);
                my $numb2=keys %hash2;
                my $numb3=keys %hash3;
                unless (@isect) {
                        print "$orth\t$numb2\t$numb3\t$info\n" if $numb2==1 && $numb3==1;
                }
        }
}
```

```perl
#!/usr/bin/perl
# Detect_Nons_all.pl
use strict;
use warnings;

my $orth="final_orth_input_paml.txt";
open ORTH, $orth or die "can not open $orth\n";
while (<ORTH>) {
        chomp;
        system("perl Detect_Nons.pl $_");
}
```

```bash
# crested blenny and common triplefin are the same, and all the other species are the same but not same with the two species
# (base) jlkang@hnu2024 Mon May 26 2025 16:08:19 /data2/jlkang/White_island/paml_input
perl Detect_Nons_all.pl > convergent_evo_genes.txt
# 125 convergent sites
```

```perl
#!/usr/bin/perl
# HPC: Detect_Nons.pl
use strict;
use warnings;
#use Array::Utils qw(:all);

my %seq; my $name;
my $orth=$ARGV[0];
my $fas="paml_files/$orth/final_alignment_pep.fa";
open FAS, $fas or die "can not open $fas\n";
while (<FAS>) {
        chomp;
        if (/^>/) {
                s/\>//;
                $name=$_;
        } else {
                $seq{$name}.=$_;
        }
}

my %cleaner=(
        'Blenny'=> 1,
        'Common'=> 1,
        );

# compare the nonsynonymous position pep sequences one by one
my %hash1;
my @poss;
my @nocls=qw(Acura Apoly Blueeyed Daru Fugu Medaka Ocomp Padel Platyfish Pmol Spottedgar Stickleback Yaldwyn Zebrafish);
my @cleas=qw(Blenny Common);
my @aspes=qw(Acura Apoly Blenny Blueeyed Common Daru Fugu Medaka Ocomp Padel Platyfish Pmol Spottedgar Stickleback Yaldwyn Zebrafish);

&Build_pos_hash(\@aspes);

sub Build_pos_hash {
        my ($grp)=@_;
        my @grp=@{$grp};
        my $len;
        foreach my $spe (@grp) {
                my $seq=$seq{$spe};
                $len=length($seq);
                for (my $i = 0; $i < $len; $i++) {
                        my $spepos=substr($seq,$i,1);
                        $hash1{$spe}->{$i}=$spepos;
                }
        }

        for (my $i = 0; $i < $len; $i++) {
                my (%hash2,%hash3);
                my $pos=$i;
                my $newp=$pos+1;
                my $info=$newp.":";
                foreach my $spe (@aspes) {
                        my $spepos=$hash1{$spe}->{$pos};
                        $info.=$spe."($spepos);";
                }

                my (@cleas_pos, @nocls_pos);
                foreach my $spe (@cleas) {
                        my $spepos=$hash1{$spe}->{$pos};
                        $hash2{$spepos}++;
                        push @cleas_pos, $spepos;
                }
                foreach my $spe (@nocls) {
                        my $spepos=$hash1{$spe}->{$pos};
                        $hash3{$spepos}++;
                        push @nocls_pos, $spepos;
                }
                my %array1 = map { $_ => 1 } @cleas_pos;
				my @isect = grep { $array1{$_} } @nocls_pos;
                #my @isect = intersect(@cleas_pos, @nocls_pos);
                my $numb2=keys %hash2;
                my $numb3=keys %hash3;
                unless (@isect) {
                        print "$orth\t$numb2\t$numb3\t$info\n" if $numb2==1 && $numb3==1;
                }
        }
}
```

```perl
#!/usr/bin/perl
# Detect_Nons_all.pl
use strict;
use warnings;

my $orth="final_orth_input_paml.txt";
open ORTH, $orth or die "can not open $orth\n";
while (<ORTH>) {
        chomp;
        system("perl Detect_Nons.pl $_");
}
```

```bash
# (base) romeodan@hpc2021 Mon May 26 20:28:16 /lustre1/g/sbs_schunter/Kang/orth16_new
perl Detect_Nons_all.pl > convergent_evo_genes.txt
# annotation
perl temp7.pl > convergent_evo_genes_ano.txt

# kangjingliang@KangdeMacBook-Pro-2 一  5 26 2025 20:58:56 ~/Documents/2025/WI/Convergence
less convergent_evo_genes_ano.txt|perl -alne 'my @a=split /\t/;$hash1{$a[0]}++;$hash2{$a[0]}=$a[-2]."\t".$a[-1];END{foreach my $key(keys %hash1){my $nb=$hash1{$key};my $ano=$hash2{$key};print "$key\t$nb\t$ano"}}' > convergent_evo_genes_ano_times.txt
```

## 3. Detect the false and random convergence sites
### 3.1. codeml estimation branch length etc.
```bash
# use OG0000493: CLCN3 as an example, which have been detected with eight convergent sites
# 1. use codeml to estimate the branch lengths, amino acid frequencies and the best shape parameter for variable rates among sites (alpha) based on the amino acid sequences
# (base) jlkang@hnu2024 Tue May 27 2025 17:35:35 /data2/jlkang/WI_convergence/OG0000030_2
vi spe.tre
# (((Ocomp,((Stickleback,Fugu),(((Daru,((Acura,Apoly),(Pmol,Padel))),(Blenny,((Yaldwyn,Blueeyed),Common))),(Platyfish,Medaka)))),Zebrafish),Spottedgar);
fasta2phy.pl final_alignment_pep.fa # phylip format: final_alignment_pep.fa.phy
cp /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input/OG0009073/estimate_before_stimulation.ctr ./
codeml estimate_before_stimulation.ctr # it doesn't matter if occurs "error: end of tree file."
```

### 3.2. evolver for amino acid sequence simulation to detect random convergence
#### Based on "mlc" (result file)
```branch_length
# branch_length
(((Ocomp: 0.004178, ((Stickleback: 0.037526, Fugu: 0.015085): 0.001538, (((Daru: 0.000004, ((Acura: 0.004078, Apoly: 0.051951): 0.000004, (Pmol: 0.000004, Padel: 0.00
0004): 0.001356): 0.000004): 0.001413, (Blenny: 0.160510, ((Yaldwyn: 0.001546, Blueeyed: 0.160354): 0.000004, Common: 0.203980): 0.006755): 0.000004): 0.000004, (Plat
yfish: 0.009600, Medaka: 0.012345): 0.001408): 0.002567): 0.000004): 0.016736, Zebrafish: 0.029813): 0.014303, Spottedgar: 0.019210);
```

```alpha_gamma
# alpha_gamma
0.34650 3
```

```amino_acid_freq
# amino_acid_freq
0.08651 0.04727 0.02702 0.03832 0.02385 0.02075 0.05915 0.07823 0.01975 0.06459 0.10534 0.04192 0.02652 0.05848 0.04577 0.06384 0.05848
 0.02527 0.03506 0.07388
A R N D C Q E G H I L K M F P S T W Y V
```

```bash
# create "MCaa.dat" in the current directory and run "evolver" which will use MCaa automatically
# (base) jlkang@hnu2024 Wed Apr 30 2025 20:24:53 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input/OG0009073
evolver 7
ls -lt # list the files according the date (the lastest to earliest)
# the out put is: ancestral.txt; mc.txt; siterates.txt
# And then estimate the ratio of convergence and non-convergence in the selected site in the 1000 replicates
# Check whether the putative convergent sites are convergent in all 1000 replicates (in mc.txt)

# 1. create directory for each replicate
# (base) jlkang@hnu2024 Wed May 28 2025 01:56:12 /data2/jlkang/WI_convergence/OG0000030_2
perl create_rep.pl
```

```perl
#!/usr/bin/perl
use strict;
use warnings;

my $rep=">rep.txt"; my $i;
open REP, $rep or die "can not open $rep\n";
for (my $i = 1; $i <= 1000; $i++) {
	my $nm="rep".$i;
	print REP "$nm\n";
}
```

```bash
# 2. make sure all the putative convengent sites that are consistent with all replicates
# 2.1 create "rep.txt"
# (base) jlkang@hnu2024 Wed May 28 2025 10:59:23 /data2/jlkang/WI_convergence/OG0000030_2
perl temp1.pl
# (base) jlkang@hnu2024 Wed May 28 2025 11:05:44 /data2/jlkang/WI_convergence/OG0000030_2
vi Detect_Nons.pl; vi Detect_Nons_all.pl
perl Detect_Nons_all.pl > convergent_evo_genes.txt
# no site is same with the putative convergent site
```

```bash
# OG0000493: sp|P51791|CLCN3_MOUSE   H(+)/Cl(-) exchange transporter 3
# (base) jlkang@hnu2024 Wed May 28 2025 11:17:46 /data2/jlkang/WI_convergence/OG0000493
cp /data2/jlkang/WI_convergence/OG0000030_2/estimate_before_stimulation.ctr ./
cp /data2/jlkang/WI_convergence/OG0000030_2/MCaa.dat ./
cp /data2/jlkang/WI_convergence/OG0000030_2/create_rep.pl ./
cp /data2/jlkang/WI_convergence/OG0000030_2/Detect_Nons.pl ./
cp /data2/jlkang/WI_convergence/OG0000030_2/Detect_Nons_all.pl ./
```

```branch_length
# branch_length
(((Ocomp: 0.387232, ((Stickleback: 0.178151, Fugu: 0.302299): 0.012656, (((Daru: 0.286117, ((Acura: 0.185769, Apoly: 0.530571): 0.000004, (Pmol: 0.004223, Padel: 0.009603): 0.249941): 0.014376): 0.000004, (Blenny: 0.119831, ((Yaldwyn: 0.207710, Blueeyed: 0.546447): 0.034464, Common: 0.236088): 0.007761): 0.047420): 0.000004, (Platyfish: 0.115291, Medaka: 0.139867): 0.443068): 0.000004): 0.158934): 0.000004, Zebrafish: 0.353532): 0.147484, Spottedgar: 0.131182);
```

```alpha_gamma
# alpha_gamma
0.48235 3
```

```amino_acid_freq
# amino_acid_freq
0.05940 0.05592 0.04209 0.07285 0.02196 0.03589 0.07672 0.07624 0.02051 0.05012 0.08969 0.05128 0.01964 0.06153 0.03038 0.07837 0.04992
 0.02138 0.02486 0.06124
A R N D C Q E G H I L K M F P S T W Y V
```


### 3.3. rebuild the ancestral sequences of each node to detect false convergence
```bash
# (base) jlkang@hnu2024 Tue May 27 2025 21:15:25 /data2/jlkang/WI_convergence/OG0000030_2
codeml Ancestral.ctr
# the sequences of each node in the rst (output file)
# compare the site with the same site in their ancsetral node sequence  
```

### 4. plot the convergent genes
```bash
# kangjingliang@KangdeMacBook-Pro-2 五  5 30 2025 17:51:22 ~/Documents/2025/WI/Convergence
scp -r jlkang@10.33.247.14:/data2/jlkang/Convergence/OG0004578/final_alignment_pep.fa OG0004578_GNRR2_pep.fa
perl temp2.pl OG0004578_GNRR2_pep.fa > OG0004578_GNRR2_pep.fa.1;mv OG0004578_GNRR2_pep.fa.1 OG0004578_GNRR2_pep.fa

# vision
scp -r jlkang@10.33.247.14:/data2/jlkang/Convergence/OG0011926/final_alignment_pep.fa OG0011926_AR2BP_pep.fa
perl temp2.pl OG0011926_AR2BP_pep.fa > OG0011926_AR2BP_pep.fa.1;mv OG0011926_AR2BP_pep.fa.1 OG0011926_AR2BP_pep.fa
scp -r jlkang@10.33.247.14:/data2/jlkang/Convergence/OG0023960/final_alignment_pep.fa OG0023960_XRP2_pep.fa
perl temp2.pl OG0023960_XRP2_pep.fa > OG0023960_XRP2_pep.fa.1;mv OG0023960_XRP2_pep.fa.1 OG0023960_XRP2_pep.fa

# CR
scp -r jlkang@10.33.247.14:/data2/jlkang/Convergence/OG0020719/final_alignment_pep.fa OG0020719_BHE40_pep.fa
perl temp2.pl OG0020719_BHE40_pep.fa > OG0020719_BHE40_pep.fa.1;mv OG0020719_BHE40_pep.fa.1 OG0020719_BHE40_pep.fa

scp -r jlkang@10.33.247.14:/data2/jlkang/Convergence/OG0001325/final_alignment_pep.fa OG0001325_TEF_pep.fa
perl temp2.pl OG0001325_TEF_pep.fa > OG0001325_TEF_pep.fa.1;mv OG0001325_TEF_pep.fa.1 OG0001325_TEF_pep.fa

# Ca2+
scp -r jlkang@10.33.247.14:/data2/jlkang/Convergence/OG0002007/final_alignment_pep.fa OG0002007_TRPM7_pep.fa
perl temp2.pl OG0002007_TRPM7_pep.fa > OG0002007_TRPM7_pep.fa.1;mv OG0002007_TRPM7_pep.fa.1 OG0002007_TRPM7_pep.fa
# OG0000169_2
scp -r jlkang@10.33.247.14:/data2/jlkang/Convergence/OG0000169_2/final_alignment_pep.fa OG0000169_2_KCC1D_pep.fa
perl temp2.pl OG0000169_2_KCC1D_pep.fa > OG0000169_2_KCC1D_pep.fa.1;mv OG0000169_2_KCC1D_pep.fa.1 OG0000169_2_KCC1D_pep.fa
# OG0024824, TM38A
scp -r jlkang@10.33.247.14:/data2/jlkang/Convergence/OG0024824/final_alignment_pep.fa OG0024824_TM38A_pep.fa
perl temp2.pl OG0024824_TM38A_pep.fa > OG0024824_TM38A_pep.fa.1;mv OG0024824_TM38A_pep.fa.1 OG0024824_TM38A_pep.fa


# Energy
# OG0011056, FAKD4
scp -r jlkang@10.33.247.14:/data2/jlkang/Convergence/OG0011056/final_alignment_pep.fa OG0011056_FAKD4_pep.fa
perl temp2.pl OG0011056_FAKD4_pep.fa > OG0011056_FAKD4_pep.fa.1;mv OG0011056_FAKD4_pep.fa.1 OG0011056_FAKD4_pep.fa
# OG0000435_2, KCRM
scp -r jlkang@10.33.247.14:/data2/jlkang/Convergence/OG0000435_2/final_alignment_pep.fa OG0000435_2_KCRM_pep.fa
perl temp2.pl OG0000435_2_KCRM_pep.fa > OG0000435_2_KCRM_pep.fa.1;mv OG0000435_2_KCRM_pep.fa.1 OG0000435_2_KCRM_pep.fa
# OG0000842, ODO1
scp -r jlkang@10.33.247.14:/data2/jlkang/Convergence/OG0000842/final_alignment_pep.fa OG0000842_ODO1_pep.fa
perl temp2.pl OG0000842_ODO1_pep.fa > OG0000842_ODO1_pep.fa.1;mv OG0000842_ODO1_pep.fa.1 OG0000842_ODO1_pep.fa
# OG0023454, PPR3B
scp -r jlkang@10.33.247.14:/data2/jlkang/Convergence/OG0023454/final_alignment_pep.fa OG0023454_PPR3B_pep.fa
perl temp2.pl OG0023454_PPR3B_pep.fa > OG0023454_PPR3B_pep.fa.1;mv OG0023454_PPR3B_pep.fa.1 OG0023454_PPR3B_pep.fa

# Heat protein
# OG0026576, DNJC3
scp -r jlkang@10.33.247.14:/data2/jlkang/Convergence/OG0026576/final_alignment_pep.fa OG0026576_DNJC3_pep.fa
perl temp2.pl OG0026576_DNJC3_pep.fa > OG0026576_DNJC3_pep.fa.1;mv OG0026576_DNJC3_pep.fa.1 OG0026576_DNJC3_pep.fa
# OG0020983, DJC10
scp -r jlkang@10.33.247.14:/data2/jlkang/Convergence/OG0020983/final_alignment_pep.fa OG0020983_DJC10_pep.fa
perl temp2.pl OG0020983_DJC10_pep.fa > OG0020983_DJC10_pep.fa.1;mv OG0020983_DJC10_pep.fa.1 OG0020983_DJC10_pep.fa
# OG0026587, DNJB1
scp -r jlkang@10.33.247.14:/data2/jlkang/Convergence/OG0026587/final_alignment_pep.fa OG0026587_DNJB1_pep.fa
perl temp2.pl OG0026587_DNJB1_pep.fa > OG0026587_DNJB1_pep.fa.1;mv OG0026587_DNJB1_pep.fa.1 OG0026587_DNJB1_pep.fa
```

### 5. Detect the specific response in the common triplefin and crested blenny
#### A. V1 vs. C
##### 5.1 Detect the specific responses only in the common triplefin
```specific_sig_funcs_commontriplefin.pl
#!/usr/bin/perl
use strict;
use warnings;

my (%yal, %blu, %ble);
my $Yaldwyn="YaldwynNofilter_enrichment.txt";
open YAL, $Yaldwyn or die "can not open $Yaldwyn\n";
while (<YAL>) {
	chomp;
	next if /^Tags/;
	my @a=split /\t/;
	$yal{$a[2]}++ if $a[6]==0;
}
my $Blueyed="BlueeyedNofilter_enrichment.txt";
open BLU, $Blueyed or die "can not open $Blueyed\n";
while (<BLU>) {
	chomp;
	next if /^Tags/;
	my @a=split /\t/;
	$blu{$a[2]}++ if $a[6]==0;
}

my $blenny="BlennyNofilter_enrichment.txt";
open BLE, $blenny or die "can not open $blenny\n";
while (<BLE>) {
	chomp;
	next if /^Tags/;
	my @a=split /\t/;
	$ble{$a[2]}++ if $a[6]==0;
}

my $common="CommonNofilter_enrichment.txt";
open COM, $common or die "can not open $common\n";
while (<COM>) {
	chomp;
	next if /^Tags/;
	my @a=split /\t/;
	print "$a[2]\n" if $yal{$a[2]} && $blu{$a[2]} && $ble{$a[2]} && $a[4]<=0.05;
}
```

```bash
# select the significant enriched functions in the common triplefin
# and make sure the other three species without any genes enriched in these functions
# kangjingliang@KangdeMacBook-Pro-2 三  6 04 2025 11:02:34 ~/Documents/2023/WI/Enrichment_VsCs
perl specific_sig_funcs_commontriplefin.pl > specific_sig_funcs_commontriplefin_V1C.txt # 278 functions
```
##### 5.2 Detect the specific responses in both the common triplefin and crested blenny
```both_specific_sig_funcs_commonblenny.pl
#!/usr/bin/perl
use strict;
use warnings;

my (%yal, %blu, %ble, %com);
my $Yaldwyn="YaldwynNofilter_enrichment.txt";
open YAL, $Yaldwyn or die "can not open $Yaldwyn\n";
while (<YAL>) {
	chomp;
	next if /^Tags/;
	my @a=split /\t/;
	$yal{$a[2]}++ if $a[6]==0;
}
my $Blueyed="BlueeyedNofilter_enrichment.txt";
open BLU, $Blueyed or die "can not open $Blueyed\n";
while (<BLU>) {
	chomp;
	next if /^Tags/;
	my @a=split /\t/;
	$blu{$a[2]}++ if $a[6]==0;
}

my $blenny="BlennyNofilter_enrichment.txt";
open BLE, $blenny or die "can not open $blenny\n";
while (<BLE>) {
	chomp;
	next if /^Tags/;
	my @a=split /\t/;
	$ble{$a[2]}++ if $a[6]>0 && $a[4]>0.05;
}

my $common="CommonNofilter_enrichment.txt";
open COM, $common or die "can not open $common\n";
while (<COM>) {
	chomp;
	next if /^Tags/;
	my @a=split /\t/;
	$com{$a[2]}++ if $a[6]>0 && $a[4]>0.05;
	print "common\t$a[2]\n" if $yal{$a[2]} && $blu{$a[2]} && $ble{$a[2]} && $a[4]<=0.05;
}

open BLE, $blenny or die "can not open $blenny\n";
while (<BLE>) {
	chomp;
	next if /^Tags/;
	my @a=split /\t/;
	print "blenny\t$a[2]\n" if $yal{$a[2]} && $blu{$a[2]} && $com{$a[2]} && $a[4]<=0.05;
}
```

```bash
# select the significant enriched functions in the common triplefin and crested blenny
# and make sure the other two species without any genes enriched in these functions
# kangjingliang@KangdeMacBook-Pro-2 三  6 04 2025 11:02:34 ~/Documents/2023/WI/Enrichment_VsCs
perl both_specific_sig_funcs_commonblenny.pl > both_specific_sig_funcs_commonblenny_V1C.txt #  132 functions
```
##### 5.3 Detect the specific responses only in the crested blenny
```specific_sig_funcs_blenny.pl
#!/usr/bin/perl
use strict;
use warnings;

my (%yal, %blu, %com);
my $Yaldwyn="YaldwynNofilter_enrichment.txt";
open YAL, $Yaldwyn or die "can not open $Yaldwyn\n";
while (<YAL>) {
	chomp;
	next if /^Tags/;
	my @a=split /\t/;
	$yal{$a[2]}++ if $a[6]==0;
}
my $Blueyed="BlueeyedNofilter_enrichment.txt";
open BLU, $Blueyed or die "can not open $Blueyed\n";
while (<BLU>) {
	chomp;
	next if /^Tags/;
	my @a=split /\t/;
	$blu{$a[2]}++ if $a[6]==0;
}

my $common="CommonNofilter_enrichment.txt";
open COM, $common or die "can not open $common\n";
while (<COM>) {
	chomp;
	next if /^Tags/;
	my @a=split /\t/;
	$com{$a[2]}++ if $a[6]==0;
}

my $blenny="BlennyNofilter_enrichment.txt";
open BLE, $blenny or die "can not open $blenny\n";
while (<BLE>) {
	chomp;
	next if /^Tags/;
	my @a=split /\t/;
	print "$a[2]\n" if $yal{$a[2]} && $blu{$a[2]} && $com{$a[2]} && $a[4]<=0.05;
}
```

```bash
# select the significant enriched functions in the crested blenny
# and make sure the other three species without any genes enriched in these functions
# kangjingliang@KangdeMacBook-Pro-2 三  6 04 2025 11:02:34 ~/Documents/2023/WI/Enrichment_VsCs
perl specific_sig_funcs_blenny.pl|wc -l # 0 gene
```
#### A. V2 vs. C
##### 5.1 Detect the specific responses only in the common triplefin
```bash
# kangjingliang@KangdeMacBook-Pro-2 三  6 04 2025 11:07:14 ~/Documents/2023/WI/Enrichment_VnCs
perl specific_sig_funcs_commontriplefin.pl > specific_sig_funcs_commontriplefin_V2C.txt # 46 functions
```
##### 5.2 Detect the specific responses in both the common triplefin and crested blenny
```bash
# kangjingliang@KangdeMacBook-Pro-2 三  6 04 2025 11:07:14 ~/Documents/2023/WI/Enrichment_VnCs
perl both_specific_sig_funcs_commonblenny.pl > both_specific_sig_funcs_commonblenny_V2C.txt # 2 functions
```
##### 5.3 Detect the specific responses only in the crested blenny
```bash
# kangjingliang@KangdeMacBook-Pro-2 三  6 04 2025 11:07:14 ~/Documents/2023/WI/Enrichment_VnCs
perl specific_sig_funcs_blenny.pl > specific_sig_funcs_blenny_V2C.txt # 4 functions
```

### Plot these functions
```generate_funcs_plot.pl
#!/usr/bin/perl
use strict;
use warnings;

my $head="Tags\tGO_ID\tGO_Name\tGO_Category\tFDR\tP-Value\tNr_Test\tNr_Reference\tNon_Annot_Test\tNon_Annot_Reference\tType\tSpecies";
print "$head\n";
my (%yal, %blu, %ble, %com);

my $Yaldwyn="YaldwynNofilter_enrichment.txt";
my $Blueyed="BlueeyedNofilter_enrichment.txt";
my $blenny="BlennyNofilter_enrichment.txt";
my $common="CommonNofilter_enrichment.txt";

%yal=&Build_hash($Yaldwyn);
%blu=&Build_hash($Blueyed);
%ble=&Build_hash($blenny);
%com=&Build_hash($common);

my $funcs=$ARGV[0];
open FUNCS, $funcs or die "can not open $funcs\n";
while (<FUNCS>) {
	chomp;
	s/\r//g;
	next if /^Function/;
	my @a=split /\t/;
	my $yalInfo=$yal{$a[0]};
	my $bluInfo=$blu{$a[0]};
	my $bleInfo=$ble{$a[0]};
	my $comInfo=$com{$a[0]};
	print "$yalInfo\t$a[1]\tYaldwin\n";
	print "$bluInfo\t$a[1]\tblue-eyed\n";
	print "$bleInfo\t$a[1]\tblenny\n";
	print "$comInfo\t$a[1]\tcommon\n";
}


sub Build_hash {
	my ($file)=@_; my %hash;
	open FILE, $file or die "can not open $file\n";
	while (<FILE>) {
		chomp;
		my @a=split /\t/;
		if (/^Tags/) {
			next;			
		} else {
			my $line;
			for (my $i = 0; $i < 10; $i++) {
				unless ($a[6]) {
					$a[6]="NA";
				}
				$line.=$a[$i]."\t";
			}
			$line=~s/\s+$//;
			$hash{$a[2]}=$line;
		}
	}
	return %hash;
}
```

```bash
# V1 vs. C
# specific functions in the common triplefin
# kangjingliang@KangdeMacBook-Pro-2 三  6 04 2025 16:17:20 ~/Documents/2023/WI/Enrichment_VsCs
perl generate_funcs_plot.pl specific_sig_funcs_commontriplefin_V1C_plot.txt > specific_sig_funcs_commontriplefin_V1C_plot_final.txt

# functions in both the common triplefin and crested blenny
# kangjingliang@KangdeMacBook-Pro-2 三  6 04 2025 16:52:12 ~/Documents/2023/WI/Enrichment_VsCs
perl generate_funcs_plot.pl both_specific_sig_funcs_commonblenny_V1C_plot.txt > both_specific_sig_funcs_commonblenny_V1C_plot_final.txt

# V2 vs. C
# specific functions in the common triplefin
# kangjingliang@KangdeMacBook-Pro-2 三  6 04 2025 19:53:58 ~/Documents/2023/WI/Enrichment_VnCs
perl generate_funcs_plot.pl specific_sig_funcs_commontriplefin_V2C_plot.txt > specific_sig_funcs_commontriplefin_V2C_plot_final.txt

# functions in both the common triplefin and crested blenny
# kangjingliang@KangdeMacBook-Pro-2 三  6 04 2025 20:13:53 ~/Documents/2023/WI/Enrichment_VnCs
perl generate_funcs_plot.pl both_specific_sig_funcs_commonblenny_V2C.txt > both_specific_sig_funcs_commonblenny_V2C_plot_final.txt

# functions only in the crested blenny
# kangjingliang@KangdeMacBook-Pro-2 三  6 04 2025 20:43:47 ~/Documents/2023/WI/Enrichment_VnCs
perl generate_funcs_plot.pl specific_sig_funcs_blenny_V2C.txt > specific_sig_funcs_blenny_V2C_plot_final.txt
```

### Extract the genes underlying the specific functions
```bash
# V1 vs. C
# Extract the genes underlying the specific functions in the common triplefin
# kangjingliang@KangdeMacBook-Pro-2 一  6 09 2025 20:03:20 ~/Documents/2023/WI/Enrichment_VsCs
less specific_sig_funcs_commontriplefin_V1C_plot.txt|perl -alne 'my @a=split /\t/;next if /Function/i;print $a[0]' > specific_sig_funcs_commontriplefin_V1C_extract.txt
extract_gene_functions -i CommonNofilter_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions specific_sig_funcs_commontriplefin_V1C_extract.txt --output specific_sig_funcs_commontriplefin_V1C_genes
# Extract the genes underlying functions shared by both the common triplefin and crested blenny
less both_specific_sig_funcs_commonblenny_V1C_plot.txt|perl -alne 'my @a=split /\t/;next if /Function/i;print $a[0]' > both_specific_sig_funcs_commonblenny_V1C_extract.txt
extract_gene_functions -i CommonNofilter_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions both_specific_sig_funcs_commonblenny_V1C_extract.txt --output specific_sig_funcs_sharedCommontriplefin_V1C_genes
extract_gene_functions -i BlennyNofilter_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions both_specific_sig_funcs_commonblenny_V1C_extract.txt --output specific_sig_funcs_sharedBlenny_V1C_genes
# add the categoried information
# specific in the common triplefin
# kangjingliang@KangdeMacBook-Pro-2 五  6 13 2025 17:36:32 ~/Documents/2023/WI/Enrichment_VsCs
perl temp2.pl specific_sig_funcs_commontriplefin_V1C_plot.txt specific_sig_funcs_commontriplefin_V1C_genes.txt > specific_sig_funcs_commontriplefin_V1C_genes_categoried.txt
# in the common triplefin shared with the crested blenny
perl temp2.pl both_specific_sig_funcs_commonblenny_V1C_plot.txt specific_sig_funcs_sharedCommontriplefin_V1C_genes.txt > specific_sig_funcs_sharedCommontriplefin_V1C_genes_categoried.txt
# in the crested blenny shared with the common triplefin
perl temp2.pl both_specific_sig_funcs_commonblenny_V1C_plot.txt specific_sig_funcs_sharedBlenny_V1C_genes.txt > specific_sig_funcs_sharedBlenny_V1C_genes_categoried.txt

# Search all the genes underlying these functions
# kangjingliang@KangdeMacBook-Pro-2 二  6 10 2025 10:36:42 ~/Documents/2023/WI/Enrichment_VsCs
less specific_sig_funcs_commontriplefin_V1C_genes.txt|perl -alne '@a=split /\t/;my ($nm)=$a[-2]=~/sp\|.*\|(.*)_.*/;;print "$nm\t$a[-1]"'|sort -u|less

# V2 vs. C
# Extract the genes underlying the specific functions in the common triplefin
# kangjingliang@KangdeMacBook-Pro-2 一  6 09 2025 17:26:02 ~/Documents/2023/WI/Enrichment_VnCs
less specific_sig_funcs_commontriplefin_V2C_plot.txt|perl -alne 'my @a=split /\t/;next if /Function/i;print $a[0]' > specific_sig_funcs_commontriplefin_V2C_extract.txt
extract_gene_functions -i CommonNofilter_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions specific_sig_funcs_commontriplefin_V2C_extract.txt --output specific_sig_funcs_commontriplefin_V2C_genes
# check the genes in specific_sig_funcs_commontriplefin_V2C_genes.txt
# kangjingliang@KangdeMacBook-Pro-2 一  6 16 2025 10:22:01 ~/Documents/2023/WI/Enrichment_VnCs
cp ../Enrichment_VsCs/temp2.pl ./
perl temp2.pl specific_sig_funcs_commontriplefin_V2C_plot.txt specific_sig_funcs_commontriplefin_V2C_genes.txt > specific_sig_funcs_commontriplefin_V2C_genes_categoried.txt
# extract the genes underlying shared functions
extract_gene_functions -i CommonNofilter_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions both_specific_sig_funcs_commonblenny_V2C.txt --output specific_sig_funcs_sharedCommontriplefin_V2C_genes
extract_gene_functions -i BlennyNofilter_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions both_specific_sig_funcs_commonblenny_V2C.txt --output specific_sig_funcs_sharedBlenny_V2C_genes

# Extract the genes underlying functions only in the crested blenny
extract_gene_functions -i BlennyNofilter_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions specific_sig_funcs_blenny_V2C.txt --output specific_sig_funcs_blenny_V2C_genes
```

### The functions that were enriched by DETs in all species and were significantly enriched by DETs in at least one species
```functions_allspes.pl
#!/usr/bin/perl
use strict;
use warnings;

my (%yal, %blu, %ble, %com);

my $Yaldwyn="YaldwynNofilter_enrichment.txt";
my $Blueyed="BlueeyedNofilter_enrichment.txt";
my $blenny="BlennyNofilter_enrichment.txt";
my $common="CommonNofilter_enrichment.txt";

%yal=&Build_hash($Yaldwyn);
%blu=&Build_hash($Blueyed);
%ble=&Build_hash($blenny);
%com=&Build_hash($common);

foreach my $key (sort keys %yal) {
	my $nb1=$yal{$key}->{'nb'};
	my $nb2=$blu{$key}->{'nb'};
	my $nb3=$ble{$key}->{'nb'};
	my $nb4=$com{$key}->{'nb'};

	my $fdr1=$yal{$key}->{'fdr'};
	my $fdr2=$blu{$key}->{'fdr'};
	my $fdr3=$ble{$key}->{'fdr'};
	my $fdr4=$com{$key}->{'fdr'};

	if ($nb1>0 && $nb2>0 && $$nb3>0 && $nb4>0) {
		if ($fdr1 <= 0.05 || $fdr2 <= 0.05 || $fdr3 <= 0.05 || $fdr4 <= 0.05) {
			print "$key\n"  if $yal{$key}->{'func'} eq "BIOLOGICAL_PROCESS";
		}
	}
}

sub Build_hash {
	my ($file)=@_; my %hash;
	open FILE, $file or die "can not open $file\n";
	while (<FILE>) {
		chomp;
		my @a=split /\t/;
		if (/^Tags/) {
			next;			
		} else {
			my $line;
			for (my $i = 0; $i < 10; $i++) {
				$line.=$a[$i]."\t";
			}
			$line=~s/\s+$//;
			$hash{$a[2]}={
				'info'=>$line,
				'fdr' =>$a[4],
				'nb'  =>\$a[6],
				'func'=>$a[3]
			};
		}
	}
	return %hash;
}
```

```bash
# V1 vs. C
# kangjingliang@KangdeMacBook-Pro-2 四  6 05 2025 16:29:20 ~/Documents/2023/WI/Enrichment_VsCs
perl functions_allspes.pl > functions_in_allspes.txt # 574 functions

# V2 vs. C
# kangjingliang@KangdeMacBook-Pro-2 四  6 05 2025 16:19:53 ~/Documents/2023/WI/Enrichment_VnCs
perl functions_allspes.pl > functions_in_allspes.txt # 348 functions

# kangjingliang@KangdeMacBook-Pro-2 一  6 09 2025 14:45:16 ~/Documents/2025/WI/SharedFunc
cp ~/Documents/2023/WI/Enrichment_VsCs/functions_in_allspes.txt functions_in_allspes_V1C.txt
cp ~/Documents/2023/WI/Enrichment_VnCs/functions_in_allspes.txt functions_in_allspes_V2C.txt
```

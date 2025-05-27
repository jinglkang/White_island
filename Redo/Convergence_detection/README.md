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

### 3.2. evolver for amino acid sequence simulation
#### Based on "mlc" (result file)
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

```bash
# create "MCaa.dat" in the current directory and run "evolver" which will use MCaa automatically
# (base) jlkang@hnu2024 Wed Apr 30 2025 20:24:53 /data2/jlkang/Nocturnal_fish/Orthologous/pep/OrthoFinder/Results_Jan15/Orthogroups/paml_input/OG0009073
evolver 7
ls -lt # list the files according the date (the lastest to earliest)
# the out put is: ancestral.txt; mc.txt; siterates.txt
# And then estimate the ratio of convergence and non-convergence in the selected site in the 1000 replicates
```

## only keep the position which was positively selected detected by hyphy and paml
```bash
# hyphy
# arthur25@hpc2021 Sun Aug 18 01:05:31 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000227
less Common_BUSTED.json|head -n 4|tail -n 1|perl -alne 's/\[//g;s/\,//g;@a=split;foreach $a (@a){$i++;print "$i" if $a>=10}'|wc -l

# paml
# kangjingliang@KangdeMacBook-Pro 日  8 18 00:57:32 ~/Documents/2023/WI/paml/WI_PSGs_plot
less 1.txt|perl -alne '@a=split /\;/;foreach my $a (@a){print $a}'

# the PSGs detected by both hyphy and paml
# arthur25@hpc2021 Sun Aug 18 11:58:24 /lustre1/g/sbs_schunter/Kang/orth16_new
vi Common_psg_id.txt
# arthur25@hpc2021 Sun Aug 18 12:29:11 /lustre1/g/sbs_schunter/Kang/orth16_new
perl temp6.pl > Common_psg_id_hyphy.txt
```

```temp6.pl
#!/usr/bin/perl
use strict;
use warnings;

my $PSGs="Common_psg_id.txt";
open PSG, $PSGs or die "can not open $PSGs\n";
while (<PSG>) {
	chomp;
	my $hyphy="paml_files/$_/Common_BUSTED.json";
	my $id=$_; my $i; 
	open HYPHY, $hyphy or die "can not open $hyphy\n";
	while (<HYPHY>) {
		chomp;
		$i++;
		if ($i==4) {
			s/\[//g;
			s/\]//g;
			s/\,//g;
			my @a=split;
			my $j; my $info;
			foreach my $a (@a) {
				$j++;
				$info.=$j.";" if $a>=10;
			}
			$info=~s/\;$//;
			$info=$id."\t".$info;
			print "$info\n";
			my $i=0;
			last;
		}
	}
}
```

## plot sequence alignments
```bash
# OG0002710: NMDE3
# arthur25@hpc2021 Sun Aug 18 16:30:19 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0002710
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0002710_NMDE3.fa
# kangjingliang@KangdeMacBook-Pro 日  8 18 16:34:45 ~/Documents/2023/WI/paml/WI_PSGs_plot_2
scp arthur25@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0002710/OG0002710_NMDE3.fa ./

# OG0017327: PDE6D
# arthur25@hpc2021 Sun Aug 18 17:05:16 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0017327
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0017327_PDE6D.fa
# kangjingliang@KangdeMacBook-Pro 日  8 18 16:34:45 ~/Documents/2023/WI/paml/WI_PSGs_plot_2
scp arthur25@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0017327/OG0017327_PDE6D.fa ./

# OG0026786: VATF
# arthur25@hpc2021 Sun Aug 18 18:37:41 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0026786
cp ../OG0000120_2/temp1.pl ./
perl temp1.pl > OG0026786_VATF.fa
# kangjingliang@KangdeMacBook-Pro 日  8 18 16:34:45 ~/Documents/2023/WI/paml/WI_PSGs_plot_2
scp arthur25@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0026786/OG0026786_VATF.fa ./

# OG0000117_2: AT2B2
# arthur25@hpc2021 Sun Aug 18 18:56:26 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000117_2
cp ../OG0000120_2/temp1.pl ./
# (base) kang1234@celia-PowerEdge-T640 Sun Aug 18 19:01:43 ~/software/bin
scp cds2pep.pl arthur25@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000117_2
./cds2pep.pl final_alignment.fa > final_alignment_pep.fa
perl temp1.pl > OG0000117_2_AT2B2.fa
# kangjingliang@KangdeMacBook-Pro 日  8 18 16:34:45 ~/Documents/2023/WI/paml/WI_PSGs_plot_2
scp arthur25@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000117_2/OG0000117_2_AT2B2.fa ./

# OG0000227: AT2A2
# arthur25@hpc2021 Sun Aug 18 19:28:26 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000227
cp ../OG0000120_2/temp1.pl ./
cp ../OG0000117_2/cds2pep.pl ./
./cds2pep.pl final_alignment.fa > final_alignment_pep.fa
perl temp1.pl > OG0000227_AT2A2.fa
# kangjingliang@KangdeMacBook-Pro 日  8 18 19:05:55 ~/Documents/2023/WI/paml/WI_PSGs_plot_2
scp arthur25@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0000227/OG0000227_AT2A2.fa ./

# OG0002651_2: HS12B
# arthur25@hpc2021 Sun Aug 18 20:01:57 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0002651_2
cp ../OG0000120_2/temp1.pl ./
cp ../OG0000117_2/cds2pep.pl ./
./cds2pep.pl final_alignment.fa > final_alignment_pep.fa
perl temp1.pl > OG0002651_2_HS12B.fa
# kangjingliang@KangdeMacBook-Pro 日  8 18 19:32:22 ~/Documents/2023/WI/paml/WI_PSGs_plot_2
scp arthur25@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0002651_2/OG0002651_2_HS12B.fa ./

# OG0003538: MICU2
# arthur25@hpc2021 Sun Aug 18 20:18:51 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0003538
cp ../OG0000120_2/temp1.pl ./
cp ../OG0000117_2/cds2pep.pl ./
./cds2pep.pl final_alignment.fa > final_alignment_pep.fa
perl temp1.pl > OG0003538_MICU2.fa
# kangjingliang@KangdeMacBook-Pro 日  8 18 19:32:22 ~/Documents/2023/WI/paml/WI_PSGs_plot_2
scp arthur25@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0003538/OG0003538_MICU2.fa ./

# OG0008436: IF2B2
# arthur25@hpc2021 Sun Aug 18 20:36:41 /lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0008436
cp ../OG0000120_2/temp1.pl ./
cp ../OG0000117_2/cds2pep.pl ./
./cds2pep.pl final_alignment.fa > final_alignment_pep.fa
perl temp1.pl > OG0008436_IF2B2.fa
scp arthur25@hpc2021-io1.hku.hk:/lustre1/g/sbs_schunter/Kang/orth16_new/paml_files/OG0008436/OG0008436_IF2B2.fa ./
```

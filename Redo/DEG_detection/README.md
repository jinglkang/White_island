# DEGs detection
## without consideration of individuals from North control
**Vs vs. Cs; Vn vs. Cs; Vn vs Vs**   
```bash
##############################################################################
# Blenny_Cs_Vn 156	Common_Cs_Vn 855	Yaldwyn_Cs_Vn 469	Blueeyed_Cs_Vn 88
# Blenny_Cs_Vs 108	Common_Cs_Vs 429	Yaldwyn_Cs_Vs 143	Blueeyed_Cs_Vs 29
# Blenny_Vn_Vs 89	Common_Vn_Vs 1143	Yaldwyn_Vn_Vs 159	Blueeyed_Vn_Vs 106
##############################################################################

########
# Common
########
# treat as up or down-regulated
# Cs_Vn: 855
# kangjingliang@kangjingliangdeMacBook-Pro 三  3 22 10:34:40 ~/Documents/2023/WI/DEGs/Common
perl up_or_dw.pl Common_Cs_Vn.DEGs.txt Common_Cs_Vn.csv Common_Cs_Vn
perl temp4.pl Common_relative_expression.txt Common_Cs_Vn_up.txt|perl -alne 'print unless /Cn/' > Common_relative_expression_plot_Cs_Vn_up.txt
perl temp4.pl Common_relative_expression.txt Common_Cs_Vn_dw.txt|perl -alne 'print unless /Cn/' > Common_relative_expression_plot_Cs_Vn_dw.txt

# Cs_Vs: 429
perl up_or_dw.pl Common_Cs_Vs.DEGs.txt Common_Cs_Vs.csv Common_Cs_Vs
perl temp4.pl Common_relative_expression.txt Common_Cs_Vs_up.txt|perl -alne 'print unless /Cn/' > Common_relative_expression_plot_Cs_Vs_up.txt
perl temp4.pl Common_relative_expression.txt Common_Cs_Vs_dw.txt|perl -alne 'print unless /Cn/' > Common_relative_expression_plot_Cs_Vs_dw.txt

# Vn_Vs: 1143
perl up_or_dw.pl Common_Vn_Vs.DEGs.txt Common_Vn_Vs.csv Common_Vn_Vs
perl temp4.pl Common_relative_expression.txt Common_Vn_Vs_up.txt|perl -alne 'print unless /Cn/' > Common_relative_expression_plot_Vn_Vs_up.txt
perl temp4.pl Common_relative_expression.txt Common_Vn_Vs_dw.txt|perl -alne 'print unless /Cn/' > Common_relative_expression_plot_Vn_Vs_dw.txt

########
# Yaldwyn
########
# treat as up or down-regulated
# Cs_Vn: 469
# kangjingliang@kangjingliangdeMacBook-Pro 三  3 22 17:32:52 ~/Documents/2023/WI/DEGs/Yaldwin
perl up_or_dw.pl Yaldwyn_Cs_Vn.DEGs.txt Yaldwyn_Cs_Vn.csv Yaldwyn_Cs_Vn
perl temp4.pl Yaldwyn_relative_expression.txt Yaldwyn_Cs_Vn_up.txt|perl -alne 'print unless /Cn/' > Yaldwyn_relative_expression_plot_Cs_Vn_up.txt
perl temp4.pl Yaldwyn_relative_expression.txt Yaldwyn_Cs_Vn_dw.txt|perl -alne 'print unless /Cn/' > Yaldwyn_relative_expression_plot_Cs_Vn_dw.txt

# Cs_Vs: 143
perl up_or_dw.pl Yaldwyn_Cs_Vs.DEGs.txt Yaldwyn_Cs_Vs.csv Yaldwyn_Cs_Vs
perl temp4.pl Yaldwyn_relative_expression.txt Yaldwyn_Cs_Vs_up.txt|perl -alne 'print unless /Cn/' > Yaldwyn_relative_expression_plot_Cs_Vs_up.txt
perl temp4.pl Yaldwyn_relative_expression.txt Yaldwyn_Cs_Vs_dw.txt|perl -alne 'print unless /Cn/' > Yaldwyn_relative_expression_plot_Cs_Vs_dw.txt

# Vn_Vs: 159
perl up_or_dw.pl Yaldwyn_Vn_Vs.DEGs.txt Yaldwyn_Vn_Vs.csv Yaldwyn_Vn_Vs
perl temp4.pl Yaldwyn_relative_expression.txt Yaldwyn_Vn_Vs_up.txt|perl -alne 'print unless /Cn/' > Yaldwyn_relative_expression_plot_Vn_Vs_up.txt
perl temp4.pl Yaldwyn_relative_expression.txt Yaldwyn_Vn_Vs_dw.txt|perl -alne 'print unless /Cn/' > Yaldwyn_relative_expression_plot_Vn_Vs_dw.txt

########
# Blenny
########
# treat as up or down-regulated
# Cs_Vn: 156
# kangjingliang@kangjingliangdeMacBook-Pro 三  3 22 17:32:52 ~/Documents/2023/WI/DEGs/Blenny
perl up_or_dw.pl Blenny_Cs_Vn.DEGs.txt Blenny_Cs_Vn.csv Blenny_Cs_Vn
perl temp4.pl Blenny_relative_expression.txt Blenny_Cs_Vn_up.txt|perl -alne 'print unless /Cn/' > Blenny_relative_expression_plot_Cs_Vn_up.txt
perl temp4.pl Blenny_relative_expression.txt Blenny_Cs_Vn_dw.txt|perl -alne 'print unless /Cn/' > Blenny_relative_expression_plot_Cs_Vn_dw.txt

# Cs_Vs: 108
perl up_or_dw.pl Blenny_Cs_Vs.DEGs.txt Blenny_Cs_Vs.csv Blenny_Cs_Vs
perl temp4.pl Blenny_relative_expression.txt Blenny_Cs_Vs_up.txt|perl -alne 'print unless /Cn/' > Blenny_relative_expression_plot_Cs_Vs_up.txt
perl temp4.pl Blenny_relative_expression.txt Blenny_Cs_Vs_dw.txt|perl -alne 'print unless /Cn/' > Blenny_relative_expression_plot_Cs_Vs_dw.txt

# Vn_Vs: 89
perl up_or_dw.pl Blenny_Vn_Vs.DEGs.txt Blenny_Vn_Vs.csv Blenny_Vn_Vs
perl temp4.pl Blenny_relative_expression.txt Blenny_Vn_Vs_up.txt|perl -alne 'print unless /Cn/' > Blenny_relative_expression_plot_Vn_Vs_up.txt
perl temp4.pl Blenny_relative_expression.txt Blenny_Vn_Vs_dw.txt|perl -alne 'print unless /Cn/' > Blenny_relative_expression_plot_Vn_Vs_dw.txt

########
# Blueeyed
########
# treat as up or down-regulated
# Cs_Vn: 88
# kangjingliang@kangjingliangdeMacBook-Pro 三  3 22 17:32:52 ~/Documents/2023/WI/DEGs/Blueeyed
perl up_or_dw.pl Blueeyed_Cs_Vn.DEGs.txt Blueeyed_Cs_Vn.csv Blueeyed_Cs_Vn
perl temp4.pl Blueeyed_relative_expression.txt Blueeyed_Cs_Vn_up.txt|perl -alne 'print unless /Cn/' > Blueeyed_relative_expression_plot_Cs_Vn_up.txt
perl temp4.pl Blueeyed_relative_expression.txt Blueeyed_Cs_Vn_dw.txt|perl -alne 'print unless /Cn/' > Blueeyed_relative_expression_plot_Cs_Vn_dw.txt

# Cs_Vs: 29
perl up_or_dw.pl Blueeyed_Cs_Vs.DEGs.txt Blueeyed_Cs_Vs.csv Blueeyed_Cs_Vs
perl temp4.pl Blueeyed_relative_expression.txt Blueeyed_Cs_Vs_up.txt|perl -alne 'print unless /Cn/' > Blueeyed_relative_expression_plot_Cs_Vs_up.txt
perl temp4.pl Blueeyed_relative_expression.txt Blueeyed_Cs_Vs_dw.txt|perl -alne 'print unless /Cn/' > Blueeyed_relative_expression_plot_Cs_Vs_dw.txt

# Vn_Vs: 106
perl up_or_dw.pl Blueeyed_Vn_Vs.DEGs.txt Blueeyed_Vn_Vs.csv Blueeyed_Vn_Vs
perl temp4.pl Blueeyed_relative_expression.txt Blueeyed_Vn_Vs_up.txt|perl -alne 'print unless /Cn/' > Blueeyed_relative_expression_plot_Vn_Vs_up.txt
perl temp4.pl Blueeyed_relative_expression.txt Blueeyed_Vn_Vs_dw.txt|perl -alne 'print unless /Cn/' > Blueeyed_relative_expression_plot_Vn_Vs_dw.txt
```
## Combine all DEGs for enrichment
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 四  3 23 05:25:48 ~/Documents/2023/WI/DEGs/Blenny
cat Blenny_Cs_Vn.DEGs.txt Blenny_Cs_Vs.DEGs.txt Blenny_Vn_Vs.DEGs.txt|sort -u >Blenny_all_DEGs.txt # 353 -> 313

# kangjingliang@kangjingliangdeMacBook-Pro 四  3 23 05:26:59 ~/Documents/2023/WI/DEGs/Blueeyed
cat Blueeyed_Cs_Vn.DEGs.txt Blueeyed_Cs_Vs.DEGs.txt Blueeyed_Vn_Vs.DEGs.txt|sort -u >Blueeyed_all_DEGs.txt # 223 -> 202

# kangjingliang@kangjingliangdeMacBook-Pro 四  3 23 05:29:37 ~/Documents/2023/WI/DEGs/Common
cat Common_Cs_Vn.DEGs.txt Common_Cs_Vs.DEGs.txt Common_Vn_Vs.DEGs.txt|sort -u >Common_all_DEGs.txt # 2427 -> 1909

# kangjingliang@kangjingliangdeMacBook-Pro 四  3 23 05:31:29 ~/Documents/2023/WI/DEGs/Yaldwin
cat Yaldwyn_Cs_Vn.DEGs.txt Yaldwyn_Cs_Vs.DEGs.txt Yaldwyn_Vn_Vs.DEGs.txt|sort -u >Yaldwyn_all_DEGs.txt # 771 -> 577

# kangjingliang@kangjingliangdeMacBook-Pro 五  3 24 05:22:39 ~/Documents/2023/WI/DEGs
mkdir Enrichment; cd Enrichment
cp ../Blenny/*enrichment* ../Blueeyed/*enrichment* ../Common/*enrichment* ../Yaldwin/*enrichment* ./
cp ~/Documents/2021/White_island/reads_number_matrix/South_enrichment/*_funcs.txt ./
cp ~/Documents/2021/White_island/reads_number_matrix/unprot_name_description_orthgroup.txt ./
extract_gene_functions -i *_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions pH_funcs.txt --output pH_DEGs
extract_gene_functions -i *_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions ion_funcs.txt --output ion_DEGs
extract_gene_functions -i *_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions immune_funcs.txt --output immune_DEGs
extract_gene_functions -i *_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions GABA_funcs.txt --output GABA_DEGs
extract_gene_functions -i *_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions CR_funcs.txt --output CR_DEGs
# up or down? which comparison
# kangjingliang@kangjingliangdeMacBook-Pro 五  3 24 06:24:57 ~/Documents/2023/WI/DEGs/Enrichment
perl Add_regu_info.pl CR_DEGs.txt > CR_DEGs_regulation.txt
```

```perl
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename; 

# Add_regu_info.pl
my %regu;

my @bleups=<../Blenny/*_up.txt>;
my @bledws=<../Blenny/*_dw.txt>;

my @bluups=<../Blueeyed/*_up.txt>;
my @bludws=<../Blueeyed/*_dw.txt>;

my @comups=<../Common/*_up.txt>;
my @comdws=<../Common/*_dw.txt>;

my @yalups=<../Yaldwin/*_up.txt>;
my @yaldws=<../Yaldwin/*_dw.txt>;

my @combis=(@bleups, @bledws, @bluups, @bludws, @comups, @comdws, @yalups, @yaldws);

foreach my $comb (@combis) {
	my ($spe, $reg)=$comb=~/\.\.\/(.*)\/(.*)\.txt/;
	open COMB, $comb or die "can not open $comb\n";
	while (<COMB>) {
		chomp;
		my $gene=$_;
		if ($regu{$gene}) {
			$regu{$gene}.=$reg.";";
		} else {
			$regu{$gene}=$reg.";";
		}
	}
	close COMB;
}

my $enri=$ARGV[0];
open ENRI, $enri or die "can not open $enri\n";
while (<ENRI>) {
	chomp;
	my @a=split /\t/;
	my $gene=$a[2];
	my $info;
	if ($regu{$gene}) {
		print "$_\t$regu{$gene}\n";
	}
}
```

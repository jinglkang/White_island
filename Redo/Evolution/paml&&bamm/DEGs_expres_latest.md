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

# kangjingliang@kangjingliangdeMacBook-Pro 三  4 19 18:09:02 ~/Documents/2023/WI/DEGs
perl combine_DEGs.pl -VsC Blenny/Blenny_Cs_Vs.DEGs.txt -VnC Blenny/Blenny_Cs_Vn.DEGs.txt -VnVs Blenny/Blenny_Vn_Vs.DEGs.txt -spe "Crested blenny" >Blenny_all_DEGs_ano.txt
perl combine_DEGs.pl -VsC Common/Common_Cs_Vs.DEGs.txt -VnC Common/Common_Cs_Vn.DEGs.txt -VnVs Common/Common_Vn_Vs.DEGs.txt -spe "Common triplefin" >Common_all_DEGs_ano.txt
perl combine_DEGs.pl -VsC Blueeyed/Blueeyed_Cs_Vs.DEGs.txt -VnC Blueeyed/Blueeyed_Cs_Vn.DEGs.txt -VnVs Blueeyed/Blueeyed_Vn_Vs.DEGs.txt -spe "Blue-eyed triplefin" >Blueeyed_all_DEGs_ano.txt
perl combine_DEGs.pl -VsC Yaldwyn/Yaldwyn_Cs_Vs.DEGs.txt -VnC Yaldwyn/Yaldwyn_Cs_Vn.DEGs.txt -VnVs Yaldwyn/Yaldwyn_Vn_Vs.DEGs.txt -spe "Yaldwin's triplefin" >Yaldwin_all_DEGs_ano.txt

# kangjingliang@kangjingliangdeMacBook-Pro 三  4 19 23:44:53 ~/Documents/2023/WI/DEGs
perl temp1.pl > Common_DEGs_2spe.txt

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
## 1. Core circadian rhythm genes: heatmap
```bash
# TPM normalized matrix
# kangjingliang@kangjingliangdeMacBook-Pro 一  3 27 15:43:57 ~/Documents/2023/WI/DEGs
scp Kang@147.8.76.229:/media/HDD/white_island/EVE_module/White_island.TPM.TMM.sqrt.rename.filtered.matrix ./

# kangjingliang@kangjingliangdeMacBook-Pro 一  3 27 14:41:09 ~/Documents/2023/WI/DEGs/Enrichment
perl Core_CR_select.pl >CR_DEGs_regulation_coreCR.txt
```
### Blenny
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 一  3 27 14:54:07 ~/Documents/2023/WI/DEGs/Blenny
less ../Enrichment/CR_DEGs_regulation_coreCR.txt|perl -alne '@a=split /\t/;my @coms=split /\;/,$a[-2];if (/^Blenny/){my $info;foreach my $com(@coms){$info.=$com.";" if $com=~/Blenny/};print "$a[2]\t$a[-1]\t$info"}' > Blenny_coreCR_DEGs.txt
less coldata_blenny.txt|perl -alne 'print unless $F[1] eq "Cn"' >coldata_blenny_remove_Cn.txt
# sort the codata as you need
# kangjingliang@kangjingliangdeMacBook-Pro 一  3 27 16:47:15 ~/Documents/2023/WI/DEGs/Blenny
cat coldata_blenny_remove_Cn.txt
#	Site
# Blenny_Vn_1	Vn
# Blenny_Vn_2	Vn
# Blenny_Vn_3	Vn
# Blenny_Vn_4	Vn
# Blenny_Vn_5	Vn
# Blenny_Vs_1	Vs
# Blenny_Vs_2	Vs
# Blenny_Vs_3	Vs
# Blenny_Vs_4	Vs
# Blenny_Vs_5	Vs
# Blenny_Cs_1     Cs # ==> change to "\t"
# Blenny_Cs_2     Cs # ==> change to "\t"
# Blenny_Cs_3	Cs
# Blenny_Cs_4	Cs
# Blenny_Cs_5	Cs
# Blenny_Cs_6	Cs

less coldata_blenny_remove_Cn.txt|perl -alne 's/\s+//g;s/\r//g;s/\s+/\t/g;print' >coldata_blenny_remove_Cn.txt.1
mv coldata_blenny_remove_Cn.txt.1 coldata_blenny_remove_Cn.txt

extract_reads_nb --matrix ../White_island.TPM.TMM.sqrt.rename.matrix --genes Blenny_coreCR_DEGs.txt --samples coldata_blenny_remove_Cn.txt >Blenny_coreCR_DEGs_norm.txt

# save "Blenny_coreCR_DEGs_norm.txt" as another one just keep the gene name
# "Blenny_coreCR_DEGs_norm_gm.txt"
# reorder
extract_reads_nb --matrix Blenny_coreCR_DEGs_norm_gm_relative.txt --samples coldata_blenny_remove_Cn.txt >Blenny_coreCR_DEGs_norm_gm_relative.txt.1
mv Blenny_coreCR_DEGs_norm_gm_relative.txt.1 Blenny_coreCR_DEGs_norm_gm_relative.txt
```
**R to get the relative gene expression**   
```R
setwd("~/Documents/2023/WI/DEGs/Blenny")
data<-read.table("Blenny_coreCR_DEGs_norm_gm.txt",header = T,row.names = 1)
for (i in 1:nrow(data)){
  data_1<-rescale(as.vector(as.matrix(data[i,])), to = c(-1, 1))
  if (i==1) {
    write.table(as.data.frame(t(data_1)),file="Blenny_coreCR_DEGs_norm_gm_relative.txt",sep = "\t")
  } else {
    write.table(as.data.frame(t(data_1)),file="Blenny_coreCR_DEGs_norm_gm_relative.txt",append = TRUE, col.names = FALSE,sep = "\t")
  }
}
```
### Blueeyed
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 一  3 27 17:04:39 ~/Documents/2023/WI/DEGs/Blueeyed
less ../Enrichment/CR_DEGs_regulation_coreCR.txt|perl -alne '@a=split /\t/;my @coms=split /\;/,$a[-2];if (/^Blueeyed/){my $info;foreach my $com(@coms){$info.=$com.";" if $com=~/Blueeyed/};print "$a[2]\t$a[-1]\t$info"}' > Blueeyed_coreCR_DEGs.txt
less coldata_blueeyed.txt|perl -alne 'print unless $F[1] eq "Cn"' >coldata_blueeyed_remove_Cn.txt
less coldata_blueeyed_remove_Cn.txt|perl -alne 's/\s+$//g;s/\r//g;s/\s+/\t/g;print' >coldata_blueeyed_remove_Cn.txt.1
mv coldata_blueeyed_remove_Cn.txt.1 coldata_blueeyed_remove_Cn.txt
extract_reads_nb --matrix ../White_island.TPM.TMM.sqrt.rename.matrix --genes Blueeyed_coreCR_DEGs.txt --samples coldata_blueeyed_remove_Cn.txt >Blueeyed_coreCR_DEGs_norm.txt
# reorder
extract_reads_nb --matrix Blueeyed_coreCR_DEGs_norm_gm_relative.txt --samples coldata_Blueeyed_remove_Cn.txt >Blueeyed_coreCR_DEGs_norm_gm_relative.txt.1
mv Blueeyed_coreCR_DEGs_norm_gm_relative.txt.1 Blueeyed_coreCR_DEGs_norm_gm_relative.txt
```
### Common
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 一  3 27 17:39:55 ~/Documents/2023/WI/DEGs/Yaldwin
less ../Enrichment/CR_DEGs_regulation_coreCR.txt|perl -alne '@a=split /\t/;my @coms=split /\;/,$a[-2];if (/^Yaldwyn/){my $info;foreach my $com(@coms){$info.=$com.";" if $com=~/Yaldwyn/};print "$a[2]\t$a[-1]\t$info"}' > Yaldwyn_coreCR_DEGs.txt
less coldata_Yaldwyn.txt|perl -alne 'print unless $F[1] eq "Cn"' >coldata_Yaldwyn_remove_Cn.txt
less coldata_Yaldwyn_remove_Cn.txt|perl -alne 's/\s+$//g;s/\r//g;s/\s+/\t/g;print' >coldata_Yaldwyn_remove_Cn.txt.1
mv coldata_Yaldwyn_remove_Cn.txt.1 coldata_Yaldwyn_remove_Cn.txt
extract_reads_nb --matrix ../White_island.TPM.TMM.sqrt.rename.matrix --genes Yaldwyn_coreCR_DEGs.txt --samples coldata_Yaldwyn_remove_Cn.txt >Yaldwyn_coreCR_DEGs_norm.txt
# reorder
extract_reads_nb --matrix Common_coreCR_DEGs_norm_gm_relative.txt --samples coldata_Common_remove_Cn.txt >Common_coreCR_DEGs_norm_gm_relative.txt.1
mv Common_coreCR_DEGs_norm_gm_relative.txt.1 Common_coreCR_DEGs_norm_gm_relative.txt
```
### Yaldwyn
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 一  3 27 17:41:30 ~/Documents/2023/WI/DEGs/Yaldwin
less ../Enrichment/CR_DEGs_regulation_coreCR.txt|perl -alne '@a=split /\t/;my @coms=split /\;/,$a[-2];if (/^Yaldwyn/){my $info;foreach my $com(@coms){$info.=$com.";" if $com=~/Yaldwyn/};print "$a[2]\t$a[-1]\t$info"}' > Yaldwyn_coreCR_DEGs.txt
less coldata_Yaldwyn.txt|perl -alne 'print unless $F[1] eq "Cn"' >coldata_Yaldwyn_remove_Cn.txt
less coldata_Yaldwyn_remove_Cn.txt|perl -alne 's/\s+$//g;s/\r//g;s/\s+/\t/g;print' >coldata_Yaldwyn_remove_Cn.txt.1
mv coldata_Yaldwyn_remove_Cn.txt.1 coldata_Yaldwyn_remove_Cn.txt
extract_reads_nb --matrix ../White_island.TPM.TMM.sqrt.rename.matrix --genes Yaldwyn_coreCR_DEGs.txt --samples coldata_Yaldwyn_remove_Cn.txt >Yaldwyn_coreCR_DEGs_norm.txt
# reorder
extract_reads_nb --matrix Yaldwyn_coreCR_DEGs_norm_gm_relative.txt --samples coldata_Yaldwyn_remove_Cn.txt >Yaldwyn_coreCR_DEGs_norm_gm_relative.txt.1
mv Yaldwyn_coreCR_DEGs_norm_gm_relative.txt.1 Yaldwyn_coreCR_DEGs_norm_gm_relative.txt
```
## 2. Core pH regulation genes: point line plot
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 二  3 28 09:57:23 ~/Documents/2023/WI/DEGs/Enrichment
perl Core_pH_select.pl >pH_DEGs_regulation_corepH.txt
```
### Blenny
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 二  3 28 14:17:06 ~/Documents/2023/WI/DEGs/Blenny
less ../Enrichment/pH_DEGs_regulation_corepH.txt|perl -alne '@a=split /\t/;my @coms=split /\;/,$a[-2];if (/^Blenny/){my $info;foreach my $com(@coms){$info.=$com.";" if $com=~/Blenny/};print "$a[2]\t$a[-1]\t$info"}' > Blenny_corepH_DEGs.txt
extract_reads_nb --matrix ../White_island.TPM.TMM.sqrt.rename.matrix --genes Blenny_corepH_DEGs.txt --samples coldata_blenny_remove_Cn.txt >Blenny_corepH_DEGs_norm.txt
cut -f 2 Blenny_corepH_DEGs.txt >Blenny_corepH_DEGs_nm.txt
# R
perl temp4.pl Blenny_corepH_DEGs_norm_gm_relative.txt Blenny_corepH_DEGs_nm.txt > Blenny_relative_expression_plot_pH.txt
perl temp5.pl >Blenny_relative_expression_plot_pH_category.txt

# R
perl temp6.pl Blenny_relative_expression_plot_pH_category_SE.txt >Blenny_relative_expression_plot_pH_category_SE_plot.txt
# plot again
perl temp7.pl Blenny_relative_expression_plot_pH_category_SE.txt >Blenny_relative_expression_plot_pH_category_SE_plot_2.txt
```
### Blueeyed
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 二  3 28 14:22:11 ~/Documents/2023/WI/DEGs/Blueeyed
# No pH regulation DEGs
```
### Yaldwyn
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 二  3 28 14:27:08 ~/Documents/2023/WI/DEGs/Yaldwyn
less ../Enrichment/pH_DEGs_regulation_corepH.txt|perl -alne '@a=split /\t/;my @coms=split /\;/,$a[-2];if (/^Yaldwyn/){my $info;foreach my $com(@coms){$info.=$com.";" if $com=~/Yaldwyn/};print "$a[2]\t$a[-1]\t$info"}' > Yaldwyn_corepH_DEGs.txt
extract_reads_nb --matrix ../White_island.TPM.TMM.sqrt.rename.matrix --genes Yaldwyn_corepH_DEGs.txt --samples coldata_Yaldwyn_remove_Cn.txt >Yaldwyn_corepH_DEGs_norm.txt
cut -f 2 Yaldwyn_corepH_DEGs.txt >Yaldwyn_corepH_DEGs_nm.txt
cp ../Blenny/temp5.pl ../Blenny/temp6.pl ./
# R
perl temp4.pl Yaldwyn_corepH_DEGs_norm_gm_relative.txt Yaldwyn_corepH_DEGs_nm.txt > Yaldwyn_relative_expression_plot_pH.txt
perl temp5.pl Yaldwyn_relative_expression_plot_pH.txt >Yaldwyn_relative_expression_plot_pH_category.txt
# R
perl temp6.pl Yaldwyn_relative_expression_plot_pH_category_SE.txt >Yaldwyn_relative_expression_plot_pH_category_SE_plot.txt
# plot again
perl temp7.pl Yaldwyn_relative_expression_plot_pH_category_SE.txt >Yaldwyn_relative_expression_plot_pH_category_SE_plot_2.txt
```
### Common
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 二  3 28 14:46:00 ~/Documents/2023/WI/DEGs/Common
less ../Enrichment/pH_DEGs_regulation_corepH.txt|perl -alne '@a=split /\t/;my @coms=split /\;/,$a[-2];if (/^Common/){my $info;foreach my $com(@coms){$info.=$com.";" if $com=~/Common/};print "$a[2]\t$a[-1]\t$info"}' > Common_corepH_DEGs.txt
extract_reads_nb --matrix ../White_island.TPM.TMM.sqrt.rename.matrix --genes Common_corepH_DEGs.txt --samples coldata_Common_remove_Cn.txt >Common_corepH_DEGs_norm.txt
cut -f 2 Common_corepH_DEGs.txt >Common_corepH_DEGs_nm.txt
cp ../Yaldwyn/temp5.pl ../Yaldwyn/temp6.pl ./
perl temp4.pl Common_corepH_DEGs_norm_gm_relative.txt Common_corepH_DEGs_nm.txt > Common_relative_expression_plot_pH.txt
perl temp5.pl Common_relative_expression_plot_pH.txt >Common_relative_expression_plot_pH_category.txt
# R
perl temp6.pl Common_relative_expression_plot_pH_category_SE.txt >Common_relative_expression_plot_pH_category_SE_plot.txt
# plot again
perl temp7.pl Common_relative_expression_plot_pH_category_SE.txt >Common_relative_expression_plot_pH_category_SE_plot_2.txt
```
## Plot KEGG pathway
### Transform Uniprot id to KEGG id
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 三  3 29 21:37:20 ~/Documents/2023/WI/DEGs/Enrichment
less unprot_name_description_orthgroup.txt|perl -alne 'next if /Orth/;my @a=split /\|/,$F[1];print "$a[2]"' >Orthgroup_uniprot_name.txt
# upload to "https://www.uniprot.org/id-mapping" to get kegg id (But it's not for human kegg id)
# transform to human id
less unprot_name_description_orthgroup.txt|perl -alne 'next if /Orth/;my @a=split /\|/,$F[1];(my $nm)=$a[-1]=~/(.*)\_/;my $b=$nm."_HUMAN";print "$b"' > Orthgroup_uniprot_name_human.txt # 12573 human id => 11543 hits
# Results: Orthgroup_uniprot_human_KEGG.tsv
# kangjingliang@kangjingliangdeMacBook-Pro 三  3 29 22:21:37 ~/Documents/2023/WI/DEGs/Enrichment
perl Extract_orth_DEGs_keggid.pl ../Blenny/Blenny_all_DEGs.txt|sort -u > Blenny_all_DEGs_keggid.txt
# submit to "https://pathview.uncc.edu/home"
perl Extract_orth_DEGs_keggid.pl ../Common/Common_all_DEGs.txt|sort -u > Common_all_DEGs_keggid.txt
perl Extract_orth_DEGs_keggid.pl ../Yaldwyn/Yaldwyn_all_DEGs.txt|sort -u > Yaldwyn_all_DEGs_keggid.txt
perl Extract_orth_DEGs_keggid.pl ../Blueeyed/Blueeyed_all_DEGs.txt|sort -u > Blueeyed_all_DEGs_keggid.txt

# extract gene from the selected pathway
# kangjingliang@kangjingliangdeMacBook-Pro 四  3 30 22:42:49 ~/Desktop/Blenny_KEGG
perl extract_gene_kegg_pathway.pl genedata.hsa04744.tsv

# extract vision related genes for a heatmap
# Acid-sensing ion channels (ASICs) are voltage-independent Na+ channels that are activated by extracellular acidification
# Mutations in the palm domain disrupt modulation of acid-sensing ion channel 1a currents by neuropeptides (scientific reports)
# Acid-Sensing Ion Channels in Zebrafish (animals: MDPI)
# kangjingliang@kangjingliangdeMacBook-Pro 五  3 31 10:12:52 ~/Documents/2023/WI/DEGs/Enrichment
vi Core_vision_DEGs_nm.txt # the name of selected vision DEGs
perl Add_regu_info.pl Vision_DEGs.txt >Vision_DEGs_regulation.txt
perl Get_unipro_nm_descrip.pl Core_vision_DEGs_nm.txt Vision_DEGs_regulation.txt|sort -u

perl Core_vision_select.pl > Vision_DEGs_regulation_coreVision.txt
# pheatmap: scale gene expression by row
# Blenny
less ../Enrichment/Vision_DEGs_regulation_coreVision.txt|perl -alne '@a=split /\t/;my @coms=split /\;/,$a[-2];if (/^Blenny/){my $info;foreach my $com(@coms){$info.=$com.";" if $com=~/Blenny/};print "$a[2]\t$a[-1]\t$info"}'|sort -u > Blenny_coreVision_DEGs.txt
extract_reads_nb --matrix ../WI_relative_expression_matrix.txt --genes Blenny_coreVision_DEGs.txt --samples coldata_Blenny_remove_Cn.txt >Blenny_coreVision_DEGs_norm.txt
# Blueeyed
less ../Enrichment/Vision_DEGs_regulation_coreVision.txt|perl -alne '@a=split /\t/;my @coms=split /\;/,$a[-2];if (/^Blueeyed/){my $info;foreach my $com(@coms){$info.=$com.";" if $com=~/Blueeyed/};print "$a[2]\t$a[-1]\t$info"}'|sort -u > Blueeyed_coreVision_DEGs.txt
extract_reads_nb --matrix ../WI_relative_expression_matrix.txt --genes Blueeyed_coreVision_DEGs.txt --samples coldata_Blueeyed_remove_Cn.txt >Blueeyed_coreVision_DEGs_norm.txt
# Common
less ../Enrichment/Vision_DEGs_regulation_coreVision.txt|perl -alne '@a=split /\t/;my @coms=split /\;/,$a[-2];if (/^Common/){my $info;foreach my $com(@coms){$info.=$com.";" if $com=~/Common/};print "$a[2]\t$a[-1]\t$info"}'|sort -u > Common_coreVision_DEGs.txt
extract_reads_nb --matrix ../WI_relative_expression_matrix.txt --genes Common_coreVision_DEGs.txt --samples coldata_Common_remove_Cn.txt >Common_coreVision_DEGs_norm.txt
# Yaldwyn
less ../Enrichment/Vision_DEGs_regulation_coreVision.txt|perl -alne '@a=split /\t/;my @coms=split /\;/,$a[-2];if (/^Yaldwyn/){my $info;foreach my $com(@coms){$info.=$com.";" if $com=~/Yaldwyn/};print "$a[2]\t$a[-1]\t$info"}'|sort -u > Yaldwyn_coreVision_DEGs.txt
extract_reads_nb --matrix ../WI_relative_expression_matrix.txt --genes Yaldwyn_coreVision_DEGs.txt --samples coldata_Yaldwyn_remove_Cn.txt >Yaldwyn_coreVision_DEGs_norm.txt

# Add V-type proton ATPase subunit S1 for pHi regulation in Blueeyed triplefin
# kangjingliang@kangjingliangdeMacBook-Pro 六  7 15 14:22:53 ~/Documents/2023/WI/DEGs/Blueeyed
vi Blueeyed_coreVision_DEGs.txt
extract_reads_nb --matrix ../WI_relative_expression_matrix.txt --genes Blueeyed_coreVision_DEGs.txt --samples coldata_Blueeyed_remove_Cn.txt >Blueeyed_coreVision_DEGs_norm.txt
```


## Select the top20 biological GO functions to extract the underlying DEGs
```bash
# ~/Documents/2023/WI/DEGs/Reduce_enrichment_top20_biological.csv
# kangjingliang@kangjingliangdeMacBook-Pro 四  4 27 15:27:30 ~/Documents/2023/WI/DEGs/Enrichment
extract_gene_functions -i *_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions top20_CR_funcs.txt --output top20_CR_funcs_DEGs
perl Add_regu_info.pl top20_CR_funcs_DEGs.txt > top20_CR_funcs_DEGs_regulation.txt
```

## Using mean expression to plot
```bash
# Blenny
# kangjingliang@kangjingliangdeMacBook-Pro 五  4 28 10:53:11 ~/Documents/2023/WI/DEGs/Blenny
perl estimate_mean_expression.pl Blenny_relative_expression_plot_Cs_Vs_up.txt > Blenny_relative_expression_plot_Cs_Vs_up_mean.txt
perl estimate_mean_expression.pl Blenny_relative_expression_plot_Cs_Vs_dw.txt > Blenny_relative_expression_plot_Cs_Vs_dw_mean.txt
perl estimate_mean_expression.pl Blenny_relative_expression_plot_Cs_Vn_dw.txt > Blenny_relative_expression_plot_Cs_Vn_dw_mean.txt
perl estimate_mean_expression.pl Blenny_relative_expression_plot_Cs_Vn_up.txt > Blenny_relative_expression_plot_Cs_Vn_up_mean.txt
perl estimate_mean_expression.pl Blenny_relative_expression_plot_Vn_Vs_dw.txt > Blenny_relative_expression_plot_Vn_Vs_dw_mean.txt
perl estimate_mean_expression.pl Blenny_relative_expression_plot_Vn_Vs_up.txt > Blenny_relative_expression_plot_Vn_Vs_up_mean.txt

# Blue-eyed
# kangjingliang@kangjingliangdeMacBook-Pro 五  4 28 11:04:51 ~/Documents/2023/WI/DEGs/Blueeyed
cp ../Blenny/estimate_mean_expression.pl ./
perl estimate_mean_expression.pl Blueeyed_relative_expression_plot_Cs_Vs_up.txt > Blueeyed_relative_expression_plot_Cs_Vs_up_mean.txt
perl estimate_mean_expression.pl Blueeyed_relative_expression_plot_Cs_Vs_dw.txt > Blueeyed_relative_expression_plot_Cs_Vs_dw_mean.txt
perl estimate_mean_expression.pl Blueeyed_relative_expression_plot_Cs_Vn_dw.txt > Blueeyed_relative_expression_plot_Cs_Vn_dw_mean.txt
perl estimate_mean_expression.pl Blueeyed_relative_expression_plot_Cs_Vn_up.txt > Blueeyed_relative_expression_plot_Cs_Vn_up_mean.txt
perl estimate_mean_expression.pl Blueeyed_relative_expression_plot_Vn_Vs_dw.txt > Blueeyed_relative_expression_plot_Vn_Vs_dw_mean.txt
perl estimate_mean_expression.pl Blueeyed_relative_expression_plot_Vn_Vs_up.txt > Blueeyed_relative_expression_plot_Vn_Vs_up_mean.txt

# Common
# kangjingliang@kangjingliangdeMacBook-Pro 五  4 28 11:22:19 ~/Documents/2023/WI/DEGs/Common
cp ../Blenny/estimate_mean_expression.pl ./
perl estimate_mean_expression.pl Common_relative_expression_plot_Cs_Vs_up.txt > Common_relative_expression_plot_Cs_Vs_up_mean.txt
perl estimate_mean_expression.pl Common_relative_expression_plot_Cs_Vs_dw.txt > Common_relative_expression_plot_Cs_Vs_dw_mean.txt
perl estimate_mean_expression.pl Common_relative_expression_plot_Cs_Vn_dw.txt > Common_relative_expression_plot_Cs_Vn_dw_mean.txt
perl estimate_mean_expression.pl Common_relative_expression_plot_Cs_Vn_up.txt > Common_relative_expression_plot_Cs_Vn_up_mean.txt
perl estimate_mean_expression.pl Common_relative_expression_plot_Vn_Vs_dw.txt > Common_relative_expression_plot_Vn_Vs_dw_mean.txt
perl estimate_mean_expression.pl Common_relative_expression_plot_Vn_Vs_up.txt > Common_relative_expression_plot_Vn_Vs_up_mean.txt

# Yaldwyn
# kangjingliang@kangjingliangdeMacBook-Pro 五  4 28 11:22:19 ~/Documents/2023/WI/DEGs/Yaldwyn
cp ../Blenny/estimate_mean_expression.pl ./
perl estimate_mean_expression.pl Yaldwyn_relative_expression_plot_Cs_Vs_up.txt > Yaldwyn_relative_expression_plot_Cs_Vs_up_mean.txt
perl estimate_mean_expression.pl Yaldwyn_relative_expression_plot_Cs_Vs_dw.txt > Yaldwyn_relative_expression_plot_Cs_Vs_dw_mean.txt
perl estimate_mean_expression.pl Yaldwyn_relative_expression_plot_Cs_Vn_dw.txt > Yaldwyn_relative_expression_plot_Cs_Vn_dw_mean.txt
perl estimate_mean_expression.pl Yaldwyn_relative_expression_plot_Cs_Vn_up.txt > Yaldwyn_relative_expression_plot_Cs_Vn_up_mean.txt
perl estimate_mean_expression.pl Yaldwyn_relative_expression_plot_Vn_Vs_dw.txt > Yaldwyn_relative_expression_plot_Vn_Vs_dw_mean.txt
perl estimate_mean_expression.pl Yaldwyn_relative_expression_plot_Vn_Vs_up.txt > Yaldwyn_relative_expression_plot_Vn_Vs_up_mean.txt

###################################################
# select the genes with specific expression pattern
###################################################
# kangjingliang@kangjingliangdeMacBook-Pro 五  4 28 16:05:09 ~/Documents/2023/WI/DEGs/Blenny
perl select_mean_expression.pl Blenny_relative_expression_plot_Cs_Vs_up_mean.txt > Blenny_relative_expression_plot_Cs_Vs_up_mean_select.txt
perl select_mean_expression.pl Blenny_relative_expression_plot_Cs_Vs_dw_mean.txt > Blenny_relative_expression_plot_Cs_Vs_dw_mean_select.txt
perl select_mean_expression.pl Blenny_relative_expression_plot_Cs_Vn_up_mean.txt > Blenny_relative_expression_plot_Cs_Vn_up_mean_select.txt
perl select_mean_expression.pl Blenny_relative_expression_plot_Cs_Vn_dw_mean.txt > Blenny_relative_expression_plot_Cs_Vn_dw_mean_select.txt
perl select_mean_expression.pl Blenny_relative_expression_plot_Vn_Vs_up_mean.txt > Blenny_relative_expression_plot_Vn_Vs_up_mean_select.txt
perl select_mean_expression.pl Blenny_relative_expression_plot_Vn_Vs_dw_mean.txt > Blenny_relative_expression_plot_Vn_Vs_dw_mean_select.txt

# Blue-eyed
# kangjingliang@kangjingliangdeMacBook-Pro 五  4 28 16:13:58 ~/Documents/2023/WI/DEGs/Blueeyed
cp ../Blenny/select_mean_expression.pl ./
perl select_mean_expression.pl Blueeyed_relative_expression_plot_Cs_Vs_up_mean.txt > Blueeyed_relative_expression_plot_Cs_Vs_up_mean_select.txt
perl select_mean_expression.pl Blueeyed_relative_expression_plot_Cs_Vs_dw_mean.txt > Blueeyed_relative_expression_plot_Cs_Vs_dw_mean_select.txt
perl select_mean_expression.pl Blueeyed_relative_expression_plot_Cs_Vn_up_mean.txt > Blueeyed_relative_expression_plot_Cs_Vn_up_mean_select.txt
perl select_mean_expression.pl Blueeyed_relative_expression_plot_Cs_Vn_dw_mean.txt > Blueeyed_relative_expression_plot_Cs_Vn_dw_mean_select.txt
perl select_mean_expression.pl Blueeyed_relative_expression_plot_Vn_Vs_up_mean.txt > Blueeyed_relative_expression_plot_Vn_Vs_up_mean_select.txt
perl select_mean_expression.pl Blueeyed_relative_expression_plot_Vn_Vs_dw_mean.txt > Blueeyed_relative_expression_plot_Vn_Vs_dw_mean_select.txt

# Common
# kangjingliang@kangjingliangdeMacBook-Pro 五  4 28 16:15:47 ~/Documents/2023/WI/DEGs/Common
cp ../Blenny/select_mean_expression.pl ./
perl select_mean_expression.pl Common_relative_expression_plot_Cs_Vs_up_mean.txt > Common_relative_expression_plot_Cs_Vs_up_mean_select.txt
perl select_mean_expression.pl Common_relative_expression_plot_Cs_Vs_dw_mean.txt > Common_relative_expression_plot_Cs_Vs_dw_mean_select.txt
perl select_mean_expression.pl Common_relative_expression_plot_Cs_Vn_up_mean.txt > Common_relative_expression_plot_Cs_Vn_up_mean_select.txt
perl select_mean_expression.pl Common_relative_expression_plot_Cs_Vn_dw_mean.txt > Common_relative_expression_plot_Cs_Vn_dw_mean_select.txt
perl select_mean_expression.pl Common_relative_expression_plot_Vn_Vs_up_mean.txt > Common_relative_expression_plot_Vn_Vs_up_mean_select.txt
perl select_mean_expression.pl Common_relative_expression_plot_Vn_Vs_dw_mean.txt > Common_relative_expression_plot_Vn_Vs_dw_mean_select.txt

# Yaldwyn
# kangjingliang@kangjingliangdeMacBook-Pro 五  4 28 16:17:07 ~/Documents/2023/WI/DEGs/Yaldwyn
cp ../Blenny/select_mean_expression.pl ./
perl select_mean_expression.pl Yaldwyn_relative_expression_plot_Cs_Vs_up_mean.txt > Yaldwyn_relative_expression_plot_Cs_Vs_up_mean_select.txt
perl select_mean_expression.pl Yaldwyn_relative_expression_plot_Cs_Vs_dw_mean.txt > Yaldwyn_relative_expression_plot_Cs_Vs_dw_mean_select.txt
perl select_mean_expression.pl Yaldwyn_relative_expression_plot_Cs_Vn_up_mean.txt > Yaldwyn_relative_expression_plot_Cs_Vn_up_mean_select.txt
perl select_mean_expression.pl Yaldwyn_relative_expression_plot_Cs_Vn_dw_mean.txt > Yaldwyn_relative_expression_plot_Cs_Vn_dw_mean_select.txt
perl select_mean_expression.pl Yaldwyn_relative_expression_plot_Vn_Vs_up_mean.txt > Yaldwyn_relative_expression_plot_Vn_Vs_up_mean_select.txt
perl select_mean_expression.pl Yaldwyn_relative_expression_plot_Vn_Vs_dw_mean.txt > Yaldwyn_relative_expression_plot_Vn_Vs_dw_mean_select.txt

# kangjingliang@kangjingliangdeMacBook-Pro 六  5 27 07:54:21 ~/Documents/2023/WI/DEGs/Yaldwyn
perl temp8.pl Yaldwyn_relative_expression_plot_Vn_Vs_dw_mean.txt |perl -alne '@a=split /\t/;$hash{$a[1]}++;END{foreach my $key (sort keys %hash){print "$key\t$hash{$key}"}}'
# Vn < Vs < Control	110 # 110 genes
# Vn < Vs > Control	23
```

### Divide into four patterns
```bash
# Blenny
# kangjingliang@kangjingliangdeMacBook-Pro 五  5 05 10:19:53 ~/Documents/2023/WI/DEGs/Blenny
cat Blenny_relative_expression_plot*_mean.txt|perl -alne 'next if /^Gene/;print' >Blenny_relative_expression_plot_total_mean.txt
perl divide_pattern.pl Blenny_relative_expression_plot_total_mean.txt

# Blueeyed
# kangjingliang@kangjingliangdeMacBook-Pro 五  5 05 10:27:03 ~/Documents/2023/WI/DEGs/Blueeyed
cp ../Blenny/divide_pattern.pl ./
cat Blueeyed_relative_expression_plot*_mean.txt|perl -alne 'next if /^Gene/;print' >Blueeyed_relative_expression_plot_total_mean.txt
perl divide_pattern.pl Blueeyed_relative_expression_plot_total_mean.txt

# Common
# kangjingliang@kangjingliangdeMacBook-Pro 五  5 05 10:27:03 ~/Documents/2023/WI/DEGs/Common
cp ../Blenny/divide_pattern.pl ./
cat Common_relative_expression_plot*_mean.txt|perl -alne 'next if /^Gene/;print' >Common_relative_expression_plot_total_mean.txt
perl divide_pattern.pl Common_relative_expression_plot_total_mean.txt

# Yaldwyn
# kangjingliang@kangjingliangdeMacBook-Pro 五  5 05 10:27:03 ~/Documents/2023/WI/DEGs/Yaldwyn
cp ../Blenny/divide_pattern.pl ./
cat Yaldwyn_relative_expression_plot*_mean.txt|perl -alne 'next if /^Gene/;print' >Yaldwyn_relative_expression_plot_total_mean.txt
perl divide_pattern.pl Yaldwyn_relative_expression_plot_total_mean.txt

# only 23/134 follow the pattern
perl divide_pattern.pl Yaldwyn_relative_expression_plot_Vn_Vs_dw_mean_select.txt
# two pattern: pattern_1b.txt (23 genes); pattern_2a.txt ()

```

### Check DEGs underlying the patterns
```bash
# Blenny
# Pattern1: V1 expression is alway highest or lowest among the individuals of three sites (V2->1, V1, Control)
# kangjingliang@kangjingliangdeMacBook-Pro 五  5 12 13:53:43 ~/Documents/2023/WI/DEGs/Blenny
# 169 DEGs
cat Blenny_relative_expression_plot_Cs_Vs_dw_mean_select.txt Blenny_relative_expression_plot_Cs_Vs_up_mean_select.txt Blenny_relative_expression_plot_Vn_Vs_dw_mean_select.txt Blenny_relative_expression_plot_Vn_Vs_up_mean.txt|perl -alne 'next if /^Gene/;$hash{$F[0]}++;print $F[0] if $hash{$F[0]}==1' > Blenny_pattern1_DEGs.txt
# Pattern2: V2->1 < V1 < Control; V2->1 > V1 > Control
# 123 DEGs
cat Blenny_relative_expression_plot_Cs_Vn_dw_mean_select.txt Blenny_relative_expression_plot_Cs_Vn_up_mean_select.txt|perl -alne 'next if /^Gene/;$hash{$F[0]}++;print $F[0] if $hash{$F[0]}==1' > Blenny_pattern2_DEGs.txt

# Blueeyed
# Pattern1: 98 DEGs
# kangjingliang@kangjingliangdeMacBook-Pro 五  5 12 14:46:48 ~/Documents/2023/WI/DEGs/Blueeyed
cat Blueeyed_relative_expression_plot_Cs_Vs_dw_mean_select.txt Blueeyed_relative_expression_plot_Cs_Vs_up_mean_select.txt Blueeyed_relative_expression_plot_Vn_Vs_dw_mean_select.txt Blueeyed_relative_expression_plot_Vn_Vs_up_mean_select.txt|perl -alne 'next if /^Gene/;$hash{$F[0]}++;print $F[0] if $hash{$F[0]}==1' > Blueeyed_pattern1_DEGs.txt
# Pattern2: 64 DEGs
cat Blueeyed_relative_expression_plot_Cs_Vn_dw_mean_select.txt Blueeyed_relative_expression_plot_Cs_Vn_up_mean_select.txt|perl -alne 'next if /^Gene/;$hash{$F[0]}++;print $F[0] if $hash{$F[0]}==1' > Blueeyed_pattern2_DEGs.txt

# Common
# Pattern1: 1198 DEGs
# kangjingliang@kangjingliangdeMacBook-Pro 五  5 12 14:53:33 ~/Documents/2023/WI/DEGs/Common
cat Common_relative_expression_plot_Cs_Vs_dw_mean_select.txt Common_relative_expression_plot_Cs_Vs_up_mean_select.txt Common_relative_expression_plot_Vn_Vs_dw_mean_select.txt Common_relative_expression_plot_Vn_Vs_up_mean_select.txt|perl -alne 'next if /^Gene/;$hash{$F[0]}++;print $F[0] if $hash{$F[0]}==1' > Common_pattern1_DEGs.txt
# Pattern2: 678 DEGs
cat Common_relative_expression_plot_Cs_Vn_dw_mean_select.txt Common_relative_expression_plot_Cs_Vn_up_mean_select.txt|perl -alne 'next if /^Gene/;$hash{$F[0]}++;print $F[0] if $hash{$F[0]}==1' > Common_pattern2_DEGs.txt

# Yaldwyn
# Pattern1: 118 DEGs
# kangjingliang@kangjingliangdeMacBook-Pro 五  5 12 15:07:23 ~/Documents/2023/WI/DEGs/Yaldwyn
cat Yaldwyn_relative_expression_plot_Cs_Vs_dw_mean_select.txt Yaldwyn_relative_expression_plot_Cs_Vs_up_mean_select.txt Yaldwyn_relative_expression_plot_Vn_Vs_dw_mean_select.txt Yaldwyn_relative_expression_plot_Vn_Vs_up_mean_select.txt|perl -alne 'next if /^Gene/;$hash{$F[0]}++;print $F[0] if $hash{$F[0]}==1' > Yaldwyn_pattern1_DEGs.txt
# Pattern2: 385 DEGs
cat Yaldwyn_relative_expression_plot_Cs_Vn_dw_mean_select.txt Yaldwyn_relative_expression_plot_Cs_Vn_up_mean_select.txt|perl -alne 'next if /^Gene/;$hash{$F[0]}++;print $F[0] if $hash{$F[0]}==1' > Yaldwyn_pattern2_DEGs.txt

extract_gene_functions -i *pattern*_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions pH_funcs.txt --output pH_DEGs_pattern
extract_gene_functions -i *pattern*_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions ion_funcs.txt --output ion_DEGs_pattern
extract_gene_functions -i *pattern*_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions immune_funcs.txt --output immune_DEGs_pattern
extract_gene_functions -i *pattern*_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions GABA_funcs.txt --output GABA_DEGs_pattern
extract_gene_functions -i *pattern*_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions CR_funcs.txt --output CR_DEGs_pattern
```
## Common DEGs after combination
```bash
# Common_DEGs_CR_funcs.txt; Common_DEGs_Vision_funcs.txt; Common_DEGs_immunity_funcs.txt
# kangjingliang@kangjingliangdeMacBook-Pro 三  5 17 10:36:45 ~/Documents/2023/WI/DEGs/Enrichment
extract_gene_functions -i Common_DEGs_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions Common_DEGs_CR_funcs.txt --output Common_DEGs_CR
extract_gene_functions -i Common_DEGs_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions Common_DEGs_Vision_funcs.txt --output Common_DEGs_Vision
extract_gene_functions -i Common_DEGs_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions Common_DEGs_immunity_funcs.txt --output Common_DEGs_immunity

# kangjingliang@kangjingliangdeMacBook-Pro 三  5 17 11:31:00 ~/Documents/2023/WI/DEGs/Enrichment
perl temp2.pl Common_DEGs_CR.txt|less
```

## Genes underlying GO functions
```bash
# Add the annotation information
# kangjingliang@kangjingliangdeMacBook-Pro 四  5 18 10:57:15 ~/Documents/2023/WI/DEGs
perl temp3.pl Blenny/Blenny_coreCR_DEGs.txt "Crested blenny" Blenny > Blenny/Blenny_coreCR_DEGs_info.txt
perl temp3.pl Blueeyed/Blueeyed_coreCR_DEGs.txt "Blue-eyed triplefin" Blueeyed >Blueeyed/Blueeyed_coreCR_DEGs_info.txt
perl temp3.pl Common/Common_coreCR_DEGs.txt "Common triplefin" Common > Common/Common_coreCR_DEGs_info.txt
perl temp3.pl Yaldwyn/Yaldwyn_coreCR_DEGs.txt "Yaldwin's triplefin" Yaldwyn >Yaldwyn/Yaldwyn_coreCR_DEGs_info.txt

# Combine the function for each gene
# kangjingliang@kangjingliangdeMacBook-Pro 四  5 18 16:34:15 ~/Documents/2023/WI/DEGs/Enrichment
perl temp3.pl Vision_DEGs_regulation.txt > Vision_DEGs_regulation_info.txt

# boxplot for vision KEGG genes
# kangjingliang@kangjingliangdeMacBook-Pro 四  5 18 19:38:49 ~/Documents/2023/WI/DEGs
perl Create_reads_nb_gene.pl --info Blenny/coldata_blenny_remove_Cn.txt --gene Blenny_vision_DEGs_KEGG.txt --spe Blenny > Blenny_vision_DEGs_KEGG_plot.txt
perl Create_reads_nb_gene.pl --info Blueeyed/coldata_blueeyed_remove_Cn.txt --gene Blueeyed_vision_DEGs_KEGG.txt --spe Blueeyed > Blueeyed_vision_DEGs_KEGG_plot.txt
perl Create_reads_nb_gene.pl --info Common/coldata_Common_remove_Cn.txt --gene Common_vision_DEGs_KEGG.txt --spe Common > Common_vision_DEGs_KEGG_plot.txt
perl Create_reads_nb_gene.pl --info Yaldwyn/coldata_Yaldwyn_remove_Cn.txt --gene Yaldwyn_vision_DEGs_KEGG.txt --spe Yaldwyn > Yaldwyn_vision_DEGs_KEGG_plot.txt
cat *_vision_DEGs_KEGG_plot.txt|perl -alne 'if (/^Gene/){$hash{$F[0]}++;print if $hash{$F[0]}==1}else{print}' > Vision_DEGs_KEGG_plot.txt

# boxplot for pH regulation DEGs
# kangjingliang@kangjingliangdeMacBook-Pro 日  5 21 11:36:12 ~/Documents/2023/WI/DEGs
perl Create_reads_nb_gene_pH.pl --info Blenny/coldata_blenny_remove_Cn.txt --gene Blenny/Blenny_corepH_DEGs.txt --spe Blenny > Blenny_corepH_DEGs_plot.txt
perl Create_reads_nb_gene_pH.pl --info Common/coldata_Common_remove_Cn.txt --gene Common/Common_corepH_DEGs.txt --spe Common > Common_corepH_DEGs_plot.txt
perl Create_reads_nb_gene_pH.pl --info Yaldwyn/coldata_Yaldwyn_remove_Cn.txt --gene Yaldwyn/Yaldwyn_corepH_DEGs.txt --spe Yaldwyn > Yaldwyn_corepH_DEGs_plot.txt
cat *_corepH_DEGs_plot.txt|perl -alne 'if (/^Gene/){$hash{$F[0]}++;print if $hash{$F[0]}==1}else{print}' > CorepH_DEGs_plot.txt
```

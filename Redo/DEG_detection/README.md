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
perl Add_regu_info_new.pl pH_DEGs_table.txt > pH_DEGs_table_regu.txt
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

## For some interesting functions
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 日  9 03 15:58:02 ~/Documents/2023/WI/DEGs/Enrichment
extract_gene_functions -i Common_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions interesting_func.txt --output Interesting_func_genes
less Interesting_func_genes.txt
```
## Investigate the specific DEGs in Common triplefin
```temp7.pl
#!/usr/bin/perl
use strict;
use warnings;

my %FUNC;
my @spes =qw(Common Blueeyed Yaldwin Blenny);
my @funcs=("behavioral fear response","aggressive behavior","larval feeding behavior","growth hormone secretion","glycogen biosynthetic process","regulation of glucose metabolic process","courtship behavior");
extract_info("Blenny", "Blenny_enrichment.txt");
extract_info("Blueeyed", "Blueeyed_enrichment.txt");
extract_info("Common", "Common_enrichment.txt");
extract_info("Yaldwin", "Yaldwyn_enrichment.txt");
foreach my $fun (@funcs) {
	foreach my $spe (@spes) {
		my ($num, $fdr);
		if ($FUNC{$spe}->{$fun}) {
			$num=$FUNC{$spe}->{$fun}->{'NUM'};
			$fdr=$FUNC{$spe}->{$fun}->{'FDR'};
		} else {
			$num=0;
			$fdr=1;
		}
		print "$fun\t$spe\t$num\t$fdr\n";
	}
}

sub extract_info {
	my ($spe, $enrich)=@_;
	open ENRICH, $enrich or die "can not open $enrich\n";
	while (<ENRICH>) {
		chomp;
		next if /^Tags/;
		my @a=split /\t/;
		my ($name, $fdr, $num)=($a[2], $a[4], $a[6]);
		$FUNC{$spe}->{$name}={
			FDR => $fdr,
			NUM => $num
		};
	}
}
```
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 一  9 04 10:35:31 ~/Documents/2023/WI/DEGs/Enrichment
perl temp7.pl > Interesting_func_compare.txt
# check the DEGs of each species underlying these interesting functions
extract_gene_functions -i *_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions Interesting_func_plot.txt --output Interesting_func_plot_DEGs
# rank the genes according to the functions
perl temp8.pl > Interesting_func_plot_DEGs_summany.txt

# make a heatmap of ion channel and glutamate receptors genes for all the four species
# Blenny
# kangjingliang@kangjingliangdeMacBook-Pro 一  9 04 16:43:09 ~/Documents/2023/WI/DEGs/Blenny
vi ion_glutamate_DEGs.txt
extract_reads_nb --matrix ../White_island.TPM.TMM.sqrt.rename.matrix --genes ion_glutamate_DEGs.txt --samples coldata_blenny_remove_Cn.txt >Blenny_ion_glutamate_DEGs_norm.txt

# Blueeyed
cd ../Blueeyed
# kangjingliang@kangjingliangdeMacBook-Pro 一  9 04 16:54:13 ~/Documents/2023/WI/DEGs/Blueeyed
cp ../Blenny/ion_glutamate_DEGs.txt ./
extract_reads_nb --matrix ../White_island.TPM.TMM.sqrt.rename.matrix --genes ion_glutamate_DEGs.txt --samples coldata_blueeyed_remove_Cn.txt >Blueeyed_ion_glutamate_DEGs_norm.txt

# Common
# kangjingliang@kangjingliangdeMacBook-Pro 一  9 04 16:58:58 ~/Documents/2023/WI/DEGs/Common
cp ../Blenny/ion_glutamate_DEGs.txt ./
extract_reads_nb --matrix ../White_island.TPM.TMM.sqrt.rename.matrix --genes ion_glutamate_DEGs.txt --samples coldata_Common_remove_Cn.txt >Common_ion_glutamate_DEGs_norm.txt

# Yaldwyn
cd ../Yaldwyn
# kangjingliang@kangjingliangdeMacBook-Pro 一  9 04 17:17:18 ~/Documents/2023/WI/DEGs/Yaldwyn
cp ../Blenny/ion_glutamate_DEGs.txt ./
extract_reads_nb --matrix ../White_island.TPM.TMM.sqrt.rename.matrix --genes ion_glutamate_DEGs.txt --samples coldata_Yaldwyn_remove_Cn.txt >Yaldwyn_ion_glutamate_DEGs_norm.txt

# make a whole heatmap based on all species
# kangjingliang@kangjingliangdeMacBook-Pro 一  9 04 20:41:23 ~/Documents/2023/WI/DEGs/Common
cat coldata_Common_remove_Cn.txt ../Blueeyed/coldata_blueeyed_remove_Cn.txt ../Yaldwyn/coldata_Yaldwyn_remove_Cn.txt ../Blenny/coldata_blenny_remove_Cn.txt > coldata_all_spe__remove_Cn.txt
less coldata_all_spe__remove_Cn.txt|perl -alne 's/\r//g;print' > coldata_all_spe__remove_Cn.txt.1
mv coldata_all_spe__remove_Cn.txt.1 coldata_all_spe__remove_Cn.txt
mv coldata_all_spe__remove_Cn.txt coldata_all_spe_remove_Cn.txt
extract_reads_nb --matrix ../White_island.TPM.TMM.sqrt.rename.matrix --genes ion_glutamate_DEGs.txt --samples coldata_all_spe_remove_Cn.txt >Allspe_ion_glutamate_DEGs_norm.txt

# Compare the big functions
# kangjingliang@kangjingliangdeMacBook-Pro 二  9 05 22:07:06 ~/Documents/2023/WI/DEGs/Enrichment
perl temp9.pl > Interesting_Bigfunc_compare.txt

# Combine both the big and specific functions
# kangjingliang@kangjingliangdeMacBook-Pro 二  9 05 22:27:13 ~/Documents/2023/WI/DEGs/Enrichment
cat Interesting_Bigfunc_compare.txt Interesting_func_compare.txt > Interesting_Totalfunc_compare.txt
vi Interesting_Totalfunc_compare.txt # remove "response to stress"

# ComplexUpset: display the overlapped genes between functions
# kangjingliang@kangjingliangdeMacBook-Pro 三  9 13 02:00:07 ~/Documents/2023/WI/DEGs/Enrichment
perl temp10.pl > ComplexUpset_input.txt

```

## Common genes in at least two spcies
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 日  9 10 15:36:51 ~/Documents/2023/WI/WI_enric
# CR
cp ~/Documents/2023/WI/DEGs/Enrichment/unprot_name_description_orthgroup.txt ./
extract_gene_functions -i Common2spe_DEG_reduce_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions Common2spe_CR_funcs.txt --output Common2spe_CRfunc_DEGs
perl ../DEGs/Enrichment/temp2.pl Common2spe_CRfunc_DEGs.txt > Common2spe_CRfunc_DEGs_info.txt

# Vision
extract_gene_functions -i Common2spe_DEG_reduce_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions Common2spe_Vision_funcs.txt --output Common2spe_Visionfunc_DEGs
perl ../DEGs/Enrichment/temp2.pl Common2spe_Visionfunc_DEGs.txt > Common2spe_Visionfunc_DEGs_info.txt

# plot the intersection by ComplexUpset

# kangjingliang@kangjingliangdeMacBook-Pro 三  9 13 10:56:38 ~/Documents/2023/WI/DEGs/Enrichment
vi Interesting_func_plot_Final.txt
extract_gene_functions -i *_enrichment.txt -a unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions Interesting_func_plot_Final.txt --output Interesting_func_plot_Final_DEGs
perl temp8.pl > Interesting_func_plot_Final_DEGs_summany.txt
```

## Check the DEGs underlying KEGG pathway
```bash
# kangjingliang@kangjingliangdeMacBook-Pro 四  9 14 11:47:16 ~/Documents/2023/WI/DEGs/Enrichment/Common_KEGG
perl extract_gene_kegg_pathway_2.pl genedata.hsa04713.tsv > Circadian_entrainment_common.txt
# 131 DEGs of the common triplefin involved in Circadian entrainment

# kangjingliang@kangjingliangdeMacBook-Pro 四  9 14 11:56:00 ~/Documents/2023/WI/DEGs/Enrichment/Blenny_KEGG
perl extract_gene_kegg_pathway_2.pl genedata.hsa04713.tsv > Circadian_entrainment_blenny.txt
# 23 DEGs of the crest blenny involved in Circadian entrainment

# None in the blue-eyed triplefin

# kangjingliang@kangjingliangdeMacBook-Pro 四  9 14 12:01:55 ~/Documents/2023/WI/DEGs/Enrichment/Yaldwyn_KEGG
perl extract_gene_kegg_pathway_2.pl genedata.hsa04713.tsv > Circadian_entrainment_Yaldwin.txt
# 21 DEGs of the crest blenny involved in Circadian entrainment
```

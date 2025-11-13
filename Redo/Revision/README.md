# Molecular Ecology revision
## The overlap between Non-DEGs in V1/V2 and V1/C or V2/C
```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;

my $sepc=$ARGV[0];
my $unip="unprot_name_description_orthgroup.txt";
my %hash;
open UNIP, $unip or die "can not open $unip\n";
while (<UNIP>) {
	chomp;
	s/\r+//g;
	next if /^Orth_id/;
	my @a=split /\t/;
	my $info=$a[1]."\t".$a[2];
	$hash{$a[0]}=$info;
}

open SPEC, $sepc or die "can not open $sepc\n";
while (<SPEC>) {
	chomp;
	my @a=split /\t/;
	s/\r+//g;
	if (/^Orth/) {
		print "$_\tUni_id\tGene_des\n";
	} else {
		print "$_\t$hash{$a[0]}\n";
	}
}
```

```bash
# kangjingliang@KangdeMBP-2 日 10 05 2025 15:17:34 ~/Documents/2025/WI/ME/Revison/NonV1V2/
perl temp1.pl Blenny.txt > Blenny_ano.txt
perl temp1.pl Blueeyed.txt > Blueeyed_ano.txt
perl temp1.pl Commontriplefin.txt > Commontriplefin_ano.txt
perl temp1.pl Yaldwin.txt > Yaldwin_ano.txt
```

## Compare log2FC for each species
```temp1.pl
#!/usr/bin/perl
use strict;
use warnings;

my @DEGs; my %hash;
my $num=$ARGV[0];
my $speDEG=$ARGV[1];
&Detect_DEGs($speDEG);
&Build_hash($ARGV[2]);
&Build_hash($ARGV[3]);
&Build_hash($ARGV[4]);
my @spes=qw(Blenny Blueeyed Common Yaldwyn);
print "OrthID\tSpecies\tLog2FC\n";

foreach my $gene (@DEGs) {
	# my $info=$gene."\t";
	# my ($info1, $info2);
	foreach my $spe (sort {$hash{$gene}->{$b} <=> $hash{$gene}->{$a}} @spes) {
		my $logFC;
		($hash{$gene}->{$spe})?($logFC=$hash{$gene}->{$spe}):($logFC="NA");
		#my $logFC=$hash{$gene}->{$spe};
		#$info1.=$spe."($logFC);";
		#$info2.=$spe.">";
		print "$gene\t$spe\t$logFC\n";
	}
	#$info1=~s/\;$//;
	#$info2=~s/\>$//;
	#print "$gene\t$info1\t$info2\n";
}


sub Detect_DEGs {
	my $spe=$_[0];
	my $csv=$spe."_V".$num."C.csv";
	open CSV, $csv or die "can not open $csv\n";
	while (<CSV>) {
		chomp;
		s/\"//g;
		next if /baseMean/ || /NA/;
		my @a=split /\,/;
		if ($a[1]>=10 && abs($a[2])>=0.3 && $a[-1]<=0.05) {
			$hash{$a[0]}->{$spe}=abs($a[2]);
			push @DEGs, $a[0];
		}
	}
}

sub Build_hash {
	my $spe=$_[0];
	my $csv=$spe."_V".$num."C.csv";
	open CSV, $csv or die "can not open $csv\n";
	while (<CSV>) {
		chomp;
		s/\"//g;
		my @a=split /\,/;
		if (/baseMean/ || /NA/) {
			next;
		} else {
			$hash{$a[0]}->{$spe}=abs($a[2]);			
		}
		# next if /baseMean/; # || /NA/;		
		# $hash{$a[0]}->{$spe}=abs($a[2]);
	}	
}
```

```bash
# V1 vs. control
# kangjingliang@KangdeMBP-2 日 10 12 2025 22:46:12 ~/Documents/2025/WI/ME/Revison/Log2FC/V1C
# Blenny
perl temp1.pl 1 Blenny Blueeyed Common Yaldwyn > Blenny_DEGs_log2FC.txt
# R test the differences between the DEG log2FC of Blenny across species

# Blueeyed
perl temp1.pl 1 Blueeyed Blenny Common Yaldwyn > Blueeyed_DEGs_log2FC.txt
# R test the differences between the DEG log2FC of Blueeyed across species

# Common
perl temp1.pl 1 Common Blenny Blueeyed Yaldwyn > Common_DEGs_log2FC.txt
# R test the differences between the DEG log2FC of Common across species

# Yaldwyn
perl temp1.pl 1 Yaldwyn Blenny Blueeyed Common > Yaldwyn_DEGs_log2FC.txt
# R test the differences between the DEG log2FC of Yaldwyn across species
```


```V1C_log2FC_compare.R
library(tidyverse)
library(ggpubr)
library(rstatix)
setwd("~/Documents/2025/WI/ME/Revison/Log2FC/V1C")

# Blenny
# Log2FC of Blennty DEGs across species
data <- read.table("Blenny_DEGs_log2FC.txt", header = T, sep="\t")
names(data)
head(data)
oneway <- aov(Log2FC ~ Species, data = data)
summary(oneway)
data %>%  pairwise_t_test(
  Log2FC ~ Species, 
  p.adjust.method = "bonferroni"
)

p1<-ggplot(data, aes(x=Species, y=Log2FC)) + 
  geom_boxplot()+ theme_bw() +
  xlab("") + theme_classic()+
  geom_hline(yintercept=0.3, linetype="dashed")+
  theme(axis.text.x=element_text(colour="black",family="Times",size=15), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=15,face="plain",colour = "black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 15,face="plain"), #设置y轴标题的字体属性
        axis.ticks.length=unit(.2, "cm"), # set tick length
        # panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Times", colour="black", size=15),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Times", colour="black", size=15), #设置图例的总标题的字体属性
        panel.grid.major = element_blank(), legend.position = c(0.15,0.88),  #不显示网格线
        panel.grid.minor = element_blank())+ylab("Log2FC of crested Blenny DEGs across species")

# Blueeyed
# Log2FC of Blueeyed DEGs across species
data <- read.table("Blueeyed_DEGs_log2FC.txt", header = T, sep="\t")
names(data)
head(data)
oneway <- aov(Log2FC ~ Species, data = data)
summary(oneway)
data %>%  pairwise_t_test(
  Log2FC ~ Species, 
  p.adjust.method = "bonferroni"
)
p2<-ggplot(data, aes(x=Species, y=Log2FC)) + 
  geom_boxplot()+ theme_bw() +
  xlab("") + theme_classic()+
  geom_hline(yintercept=0.3, linetype="dashed")+
  theme(axis.text.x=element_text(colour="black",family="Times",size=15), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=15,face="plain",colour = "black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 15,face="plain"), #设置y轴标题的字体属性
        axis.ticks.length=unit(.2, "cm"), # set tick length
        # panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Times", colour="black", size=15),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Times", colour="black", size=15), #设置图例的总标题的字体属性
        panel.grid.major = element_blank(), legend.position = c(0.15,0.88),  #不显示网格线
        panel.grid.minor = element_blank())+ylab("Log2FC of Blue-eyed triplefin DEGs across species")

# Common
# Log2FC of Common DEGs across species
data <- read.table("Common_DEGs_log2FC.txt", header = T, sep="\t")
names(data)
head(data)
oneway <- aov(Log2FC ~ Species, data = data)
summary(oneway)
data %>%  pairwise_t_test(
  Log2FC ~ Species, 
  p.adjust.method = "bonferroni"
)
p3<-ggplot(data, aes(x=Species, y=Log2FC)) + 
  geom_boxplot()+ theme_bw() +
  xlab("") + theme_classic()+
  geom_hline(yintercept=0.3, linetype="dashed")+
  theme(axis.text.x=element_text(colour="black",family="Times",size=15), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=15,face="plain",colour = "black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 15,face="plain"), #设置y轴标题的字体属性
        axis.ticks.length=unit(.2, "cm"), # set tick length
        # panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Times", colour="black", size=15),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Times", colour="black", size=15), #设置图例的总标题的字体属性
        panel.grid.major = element_blank(), legend.position = c(0.15,0.88),  #不显示网格线
        panel.grid.minor = element_blank())+ylab("Log2FC of common triplefin DEGs across species")

# Yaldwyn
# Log2FC of Yaldwyn DEGs across species
data <- read.table("Yaldwyn_DEGs_log2FC.txt", header = T, sep="\t")
names(data)
head(data)
oneway <- aov(Log2FC ~ Species, data = data)
summary(oneway)
data %>%  pairwise_t_test(
  Log2FC ~ Species, 
  p.adjust.method = "bonferroni"
)
p4<-ggplot(data, aes(x=Species, y=Log2FC)) + 
  geom_boxplot()+ theme_bw() +
  xlab("") + theme_classic()+
  geom_hline(yintercept=0.3, linetype="dashed")+
  theme(axis.text.x=element_text(colour="black",family="Times",size=15), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=15,face="plain",colour = "black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 15,face="plain"), #设置y轴标题的字体属性
        axis.ticks.length=unit(.2, "cm"), # set tick length
        # panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Times", colour="black", size=15),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Times", colour="black", size=15), #设置图例的总标题的字体属性
        panel.grid.major = element_blank(), legend.position = c(0.15,0.88),  #不显示网格线
        panel.grid.minor = element_blank())+ylab("Log2FC of Yaldwin's triplefin DEGs across species")

ggarrange(p1,p2,p3,p4,ncol = 2, nrow = 2, widths = c(1,1), heights = c(1,1),
          common.legend = T,align ="hv", 
          font.label=list(size = 14, color = "black", face = "bold",family="Times"))
```

```bash
# V2 vs. control
# kangjingliang@KangdeMBP-2 日 10 12 2025 22:46:12 ~/Documents/2025/WI/ME/Revison/Log2FC/V2C
# Blenny
perl temp1.pl 2 Blenny Blueeyed Common Yaldwyn > Blenny_DEGs_log2FC.txt
# R test the differences between the DEG log2FC of Blenny across species

# Blueeyed
perl temp1.pl 2 Blueeyed Blenny Common Yaldwyn > Blueeyed_DEGs_log2FC.txt
# R test the differences between the DEG log2FC of Blueeyed across species

# Common
perl temp1.pl 2 Common Blenny Blueeyed Yaldwyn > Common_DEGs_log2FC.txt
# R test the differences between the DEG log2FC of Common across species

# Yaldwyn
perl temp1.pl 2 Yaldwyn Blenny Blueeyed Common > Yaldwyn_DEGs_log2FC.txt
# R test the differences between the DEG log2FC of Yaldwyn across species
```

### The reasons that the DEGs of V2C is greater than V1C
```bash
# kangjingliang@KangdeMacBook-Pro-2 三 10 22 2025 20:31:05 ~/Documents/2025/WI/ME/Revison/Enrichment_V1V2
perl temp2.pl Common
# V1>>V2	2118
# V1<<V2	7748
# V1==V2	906

perl temp2.pl Yaldwyn
# V1>>V2	490
# V1<<V2	6846
# V1==V2	646

perl temp2.pl Blueeyed
# V1>>V2	447
# V1<<V2	3787
# V1==V2	276

perl temp2.pl Blenny
# V1>>V2	2508
# V1<<V2	3566
# V1==V2	678

# For the crested blenny and Yaldwin's triplefin
# extract the significantly enriched functions (BIOLOGICAL_PROCESS) for V1C and V2C
# the number of genes in the significantly enriched functions is more than 5
# V2C-V1C > 5 || V1C-V2C > 5
perl Blenny_Yaldwin.pl Blenny
# V1>>V2	0	>5
# V1<<V2	7	>5
# V1==V2	0

perl Blenny_Yaldwin.pl Yaldwyn
# V1>>V2	0	>5
# V1<<V2	33	>5
# V1==V2	2

# For the common triplefin
# extract the significantly enriched functions (BIOLOGICAL_PROCESS) for V1C and V2C
# the number of genes in the significantly enriched functions is more than 10
# V2C-V1C > 80 || V1C-V2C > 80
perl Common_triplefin.pl Common
# V1>>V2	0	(V1-V2)>80
# V1<<V2	126	(V2-V1)>80
# V1==V2	16	V1==V2

# only use the significant enriched functions (reduced) in V1C and V2C for plot
# kangjingliang@KangdeMacBook-Pro-2 四 10 23 2025 20:02:22 ~/Documents/2025/WI/ME/Revison/Enrichment_V1V2
perl temp3.pl Blenny > Blenny_go_plot.txt
perl temp3.pl Yaldwyn > Yaldwyn_go_plot.txt
perl temp3.pl Common > Common_go_plot.txt
```

## Extract the DEGs related to behavior
```bash
# V1 vs. control
# 1. crested blenny
# kangjingliang@KangdeMacBook-Pro-2 日 10 26 2025 10:13:52 ~/Documents/2025/WI/ME/Revison/Enrichment_V1V2/Enrichment_V1C
less BlennyNofilter_enrichment.txt|perl -alne 'next if /^Tags/;my @a=split /\t/;if ($a[3] eq "BIOLOGICAL_PROCESS" && $a[2]=~/behavior/i && $a[6] > 0){print $a[2]}' > Blenny_behaviorFuncs.txt
extract_gene_functions -i BlennyNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions Blenny_behaviorFuncs.txt --output Blenny_behaviorFuncs_DEGs
less Blenny_behaviorFuncs_DEGs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[3]\t$a[4]"'|sort -u > Blenny_behaviorFuncs_DEGs_uniq.txt # 21 DEGs

# 2. Blueeyed triplefin
less BlueeyedNofilter_enrichment.txt|perl -alne 'next if /^Tags/;my @a=split /\t/;if ($a[3] eq "BIOLOGICAL_PROCESS" && $a[2]=~/behavior/i && $a[6] > 0){print $a[2]}' > Blueeyed_behaviorFuncs.txt
extract_gene_functions -i BlueeyedNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions Blueeyed_behaviorFuncs.txt --output Blueeyed_behaviorFuncs_DEGs
less Blueeyed_behaviorFuncs_DEGs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[3]\t$a[4]"'|sort -u > Blueeyed_behaviorFuncs_DEGs_uniq.txt # 5 DEGs

# 3. Common triplefin
less CommonNofilter_enrichment.txt|perl -alne 'next if /^Tags/;my @a=split /\t/;if ($a[3] eq "BIOLOGICAL_PROCESS" && $a[2]=~/behavior/i && $a[6] > 0){print $a[2]}' > Common_behaviorFuncs.txt
extract_gene_functions -i CommonNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions Common_behaviorFuncs.txt --output Common_behaviorFuncs_DEGs
less Common_behaviorFuncs_DEGs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[3]\t$a[4]"'|sort -u > Common_behaviorFuncs_DEGs_uniq.txt # 127 DEGs

# 4. Yaldwyn's triplefin
less YaldwynNofilter_enrichment.txt|perl -alne 'next if /^Tags/;my @a=split /\t/;if ($a[3] eq "BIOLOGICAL_PROCESS" && $a[2]=~/behavior/i && $a[6] > 0){print $a[2]}' > Yaldwyn_behaviorFuncs.txt
extract_gene_functions -i YaldwynNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions Yaldwyn_behaviorFuncs.txt --output Yaldwyn_behaviorFuncs_DEGs
less Yaldwyn_behaviorFuncs_DEGs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[3]\t$a[4]"'|sort -u > Yaldwyn_behaviorFuncs_DEGs_uniq.txt # 44 DEGs


###########################
# V2 vs. control
# 1. crested blenny
# kangjingliang@KangdeMacBook-Pro-2 日 10 26 2025 10:13:52 ~/Documents/2025/WI/ME/Revison/Enrichment_V1V2/Enrichment_V2C
less BlennyNofilter_enrichment.txt|perl -alne 'next if /^Tags/;my @a=split /\t/;if ($a[3] eq "BIOLOGICAL_PROCESS" && $a[2]=~/behavior/i && $a[6] > 0){print $a[2]}' > Blenny_behaviorFuncs.txt
extract_gene_functions -i BlennyNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions Blenny_behaviorFuncs.txt --output Blenny_behaviorFuncs_DEGs
less Blenny_behaviorFuncs_DEGs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[3]\t$a[4]"'|sort -u > Blenny_behaviorFuncs_DEGs_uniq.txt # 33 DEGs

# 2. Blueeyed triplefin
less BlueeyedNofilter_enrichment.txt|perl -alne 'next if /^Tags/;my @a=split /\t/;if ($a[3] eq "BIOLOGICAL_PROCESS" && $a[2]=~/behavior/i && $a[6] > 0){print $a[2]}' > Blueeyed_behaviorFuncs.txt
extract_gene_functions -i BlueeyedNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions Blueeyed_behaviorFuncs.txt --output Blueeyed_behaviorFuncs_DEGs
less Blueeyed_behaviorFuncs_DEGs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[3]\t$a[4]"'|sort -u > Blueeyed_behaviorFuncs_DEGs_uniq.txt # 15 DEGs

# 3. Common triplefin
less CommonNofilter_enrichment.txt|perl -alne 'next if /^Tags/;my @a=split /\t/;if ($a[3] eq "BIOLOGICAL_PROCESS" && $a[2]=~/behavior/i && $a[6] > 0){print $a[2]}' > Common_behaviorFuncs.txt
extract_gene_functions -i CommonNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions Common_behaviorFuncs.txt --output Common_behaviorFuncs_DEGs
less Common_behaviorFuncs_DEGs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[3]\t$a[4]"'|sort -u > Common_behaviorFuncs_DEGs_uniq.txt # 206 DEGs

# 4. Yaldwyn's triplefin
less YaldwynNofilter_enrichment.txt|perl -alne 'next if /^Tags/;my @a=split /\t/;if ($a[3] eq "BIOLOGICAL_PROCESS" && $a[2]=~/behavior/i && $a[6] > 0){print $a[2]}' > Yaldwyn_behaviorFuncs.txt
extract_gene_functions -i YaldwynNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions Yaldwyn_behaviorFuncs.txt --output Yaldwyn_behaviorFuncs_DEGs
less Yaldwyn_behaviorFuncs_DEGs.txt|perl -alne 'my @a=split /\t/;print "$a[2]\t$a[3]\t$a[4]"'|sort -u > Yaldwyn_behaviorFuncs_DEGs_uniq.txt # 99 DEGs

# Target two GOs: 1. behavior; 2. synaptic signaling
# kangjingliang@KangdeMBP-2 一 11 10 2025 17:15:35 ~/Documents/2025/WI/ME/Revison/Enrichment_V1V2/Enrichment_V1C
vi target_two_GOs.txt
extract_gene_functions -i BlennyNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions target_two_GOs.txt --output Blenny_target_two_GOs_DEGs
extract_gene_functions -i BlueeyedNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions target_two_GOs.txt --output Blueeyed_target_two_GOs_DEGs
extract_gene_functions -i CommonNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions target_two_GOs.txt --output Common_target_two_GOs_DEGs
extract_gene_functions -i YaldwynNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions target_two_GOs.txt --output Yaldwyn_target_two_GOs_DEGs

less Common_target_two_GOs_DEGs.txt |perl -alne '@a=split /\t/;$hash1{$a[2]}="$a[3]\t$a[4]";$hash2{$a[2]}.=$a[1]." & ";END{foreach my $key (keys %hash1){my $des=$hash1{$key}; my $go=$hash2{$key};$go=~s/\s+&\s+$//;print "common triplefin\tV1vsC\t$go\t$key\t$des"}}' > Common_target_two_GOs_DEGs_uniq_V1C.txt
less Yaldwyn_target_two_GOs_DEGs.txt |perl -alne '@a=split /\t/;$hash1{$a[2]}="$a[3]\t$a[4]";$hash2{$a[2]}.=$a[1]." & ";END{foreach my $key (keys %hash1){my $des=$hash1{$key}; my $go=$hash2{$key};$go=~s/\s+&\s+$//;print "Yaldwin triplefin\tV1vsC\t$go\t$key\t$des"}}' > Yaldwyn_target_two_GOs_DEGs_uniq_V1C.txt
less Blenny_target_two_GOs_DEGs.txt |perl -alne '@a=split /\t/;$hash1{$a[2]}="$a[3]\t$a[4]";$hash2{$a[2]}.=$a[1]." & ";END{foreach my $key (keys %hash1){my $des=$hash1{$key}; my $go=$hash2{$key};$go=~s/\s+&\s+$//;print "crested blenny\tV1vsC\t$go\t$key\t$des"}}' > Blenny_target_two_GOs_DEGs_uniq_V1C.txt
less Blueeyed_target_two_GOs_DEGs.txt |perl -alne '@a=split /\t/;$hash1{$a[2]}="$a[3]\t$a[4]";$hash2{$a[2]}.=$a[1]." & ";END{foreach my $key (keys %hash1){my $des=$hash1{$key}; my $go=$hash2{$key};$go=~s/\s+&\s+$//;print "blueeyed triplefin\tV1vsC\t$go\t$key\t$des"}}' > Blueeyed_target_two_GOs_DEGs_uniq_V1C.txt

# kangjingliang@KangdeMBP-2 一 11 10 2025 17:39:23 ~/Documents/2025/WI/ME/Revison/Enrichment_V1V2/Enrichment_V2C
extract_gene_functions -i BlennyNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions target_two_GOs.txt --output Blenny_target_two_GOs_DEGs
extract_gene_functions -i BlueeyedNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions target_two_GOs.txt --output Blueeyed_target_two_GOs_DEGs
extract_gene_functions -i CommonNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions target_two_GOs.txt --output Common_target_two_GOs_DEGs
extract_gene_functions -i YaldwynNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions target_two_GOs.txt --output Yaldwyn_target_two_GOs_DEGs

less Common_target_two_GOs_DEGs.txt |perl -alne '@a=split /\t/;$hash1{$a[2]}="$a[3]\t$a[4]";$hash2{$a[2]}.=$a[1]." & ";END{foreach my $key (keys %hash1){my $des=$hash1{$key}; my $go=$hash2{$key};$go=~s/\s+&\s+$//;print "common triplefin\tV2vsC\t$go\t$key\t$des"}}' > Common_target_two_GOs_DEGs_uniq_V2C.txt
less Yaldwyn_target_two_GOs_DEGs.txt |perl -alne '@a=split /\t/;$hash1{$a[2]}="$a[3]\t$a[4]";$hash2{$a[2]}.=$a[1]." & ";END{foreach my $key (keys %hash1){my $des=$hash1{$key}; my $go=$hash2{$key};$go=~s/\s+&\s+$//;print "Yaldwin triplefin\tV2vsC\t$go\t$key\t$des"}}' > Yaldwyn_target_two_GOs_DEGs_uniq_V2C.txt
less Blenny_target_two_GOs_DEGs.txt |perl -alne '@a=split /\t/;$hash1{$a[2]}="$a[3]\t$a[4]";$hash2{$a[2]}.=$a[1]." & ";END{foreach my $key (keys %hash1){my $des=$hash1{$key}; my $go=$hash2{$key};$go=~s/\s+&\s+$//;print "crested blenny\tV2vsC\t$go\t$key\t$des"}}' > Blenny_target_two_GOs_DEGs_uniq_V2C.txt
less Blueeyed_target_two_GOs_DEGs.txt |perl -alne '@a=split /\t/;$hash1{$a[2]}="$a[3]\t$a[4]";$hash2{$a[2]}.=$a[1]." & ";END{foreach my $key (keys %hash1){my $des=$hash1{$key}; my $go=$hash2{$key};$go=~s/\s+&\s+$//;print "blueeyed triplefin\tV2vsC\t$go\t$key\t$des"}}' > Blueeyed_target_two_GOs_DEGs_uniq_V2C.txt
```

### only focus on the target behavioral functions
```bash
# kangjingliang@KangdeMacBook-Pro-2 日 10 26 2025 13:16:47 ~/Documents/2025/WI/ME/Revison/Enrichment_V1V2/Enrichment_V1C
vi target_behavior_funcs.txt
# visual behavior
# olfactory behavior
# auditory behavior
# exploration behavoir

extract_gene_functions -i BlennyNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions target_behavior_funcs.txt --output Blenny_target_behavior_funcs_DEGs
extract_gene_functions -i BlueeyedNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions target_behavior_funcs.txt --output Blueeyed_target_behavior_funcs_DEGs
extract_gene_functions -i CommonNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions target_behavior_funcs.txt --output Common_target_behavior_funcs_DEGs
extract_gene_functions -i YaldwynNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions target_behavior_funcs.txt --output Yaldwyn_target_behavior_funcs_DEGs

# kangjingliang@KangdeMacBook-Pro-2 日 10 26 2025 13:27:47 ~/Documents/2025/WI/ME/Revison/Enrichment_V1V2/Enrichment_V2C
extract_gene_functions -i BlennyNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions target_behavior_funcs.txt --output Blenny_target_behavior_funcs_DEGs
extract_gene_functions -i BlueeyedNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions target_behavior_funcs.txt --output Blueeyed_target_behavior_funcs_DEGs
extract_gene_functions -i CommonNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions target_behavior_funcs.txt --output Common_target_behavior_funcs_DEGs
extract_gene_functions -i YaldwynNofilter_enrichment.txt -a ../unprot_name_description_orthgroup.txt --gene_column 1 --func_column 3 --functions target_behavior_funcs.txt --output Yaldwyn_target_behavior_funcs_DEGs
```

### list the significantly enriched functions across species
```bash
# kangjingliang@KangdeMacBook-Pro-2 一 11 03 2025 21:31:14 ~/Documents/2025/WI/ME/Revison/Enrichment_V1V2/Enrichment_V1C
less Blenny_reduced_enrichment.txt|perl -alne 's/\r//g;my @a=split /\t/;if (/^Tags/){print "$_\tSpecies\tCMP"}else{print "$_\tblenny\tV1vs.C" if $a[4] <= 0.05}' > 1.txt
less Blueeyed_reduced_enrichment.txt|perl -alne 's/\r//g;my @a=split /\t/;if (/^Tags/){next}else{print "$_\tblue-eyed\tV1vs.C" if $a[4] <= 0.05}' >> 1.txt
less Common_reduced_enrichment.txt|perl -alne 's/\r//g;my @a=split /\t/;if (/^Tags/){next}else{print "$_\tcommon\tV1vs.C" if $a[4] <= 0.05}' >> 1.txt
less Yaldwyn_reduced_enrichment.txt|perl -alne 's/\r//g;my @a=split /\t/;if (/^Tags/){next}else{print "$_\tYaldwin\tV1vs.C" if $a[4] <= 0.05}' >> 1.txt
cp 1.txt ../Enrichment_V2C/
# kangjingliang@KangdeMacBook-Pro-2 一 11 03 2025 21:35:37 ~/Documents/2025/WI/ME/Revison/Enrichment_V1V2/Enrichment_V2C
less Blenny_reduced_enrichment.txt|perl -alne 's/\r//g;my @a=split /\t/;if (/^Tags/){next}else{print "$_\tblenny\tV2vs.C" if $a[4] <= 0.05}' >> 1.txt
less Blueeyed_reduced_enrichment.txt|perl -alne 's/\r//g;my @a=split /\t/;if (/^Tags/){next}else{print "$_\tblue-eyed\tV2vs.C" if $a[4] <= 0.05}' >> 1.txt
less Common_reduced_enrichment.txt|perl -alne 's/\r//g;my @a=split /\t/;if (/^Tags/){next}else{print "$_\tcommon\tV2vs.C" if $a[4] <= 0.05}' >> 1.txt
less Yaldwyn_reduced_enrichment.txt|perl -alne 's/\r//g;my @a=split /\t/;if (/^Tags/){next}else{print "$_\tYaldwin\tV2vs.C" if $a[4] <= 0.05}' >> 1.txt
```
### hyphy: positive selective analyse
```bash
# (base) jlkang@hnu2024 Tue Nov 04 2025 11:06:35 ~/software
conda create -n hyphy_env -c bioconda -c conda-forge hyphy
conda activate hyphy_env
# (hyphy_env) jlkang@hnu2024 Tue Nov 04 2025 11:15:28 ~/software
cp BUSTED-MH.bf ~/WI_convergence
# (hyphy_env) jlkang@hnu2024 Tue Nov 04 2025 13:27:35 ~/WI_convergence
ll|perl -alne 'my @a=split;$a[-1]=~s/\/$//;print $a[-1] if /^d/ && /OG/' > CEGs_list.txt
cat spe_Blenny_Common.tre
# (((Ocomp,((Stickleback,Fugu),(((Daru,((Acura,Apoly),(Pmol,Padel))),(Blenny{Foreground},((Yaldwyn,Blueeyed),Common{Foreground}))),(Platyfish,Medaka)))),Zebrafish),Spottedgar);
nohup perl run_hyphy.pl CEGs_list.txt spe_Blenny_Common.tre Blenny_Common > hyphy_process.txt 2>&1 &
# [1] 873254
# json result: final_alignment.fa.BUSTED.json
# summary the results
# (base) jlkang@hnu2024 Tue Nov 04 2025 16:55:23 ~/WI_convergence/Hyphy
perl temp1.pl > Hyphy_pValue.txt
R
# p_apoly<-read.table(file="Hyphy_pValue.txt")
# p_apoly$fdr<- p.adjust(p_apoly$V2,method="fdr",length(p_apoly$V2))
# write.table(p_apoly, file="Hyphy_pValue_fdr.txt",row.names=F,col.names=F,quote=F,sep="\t")
 
# obtain the positive selection site
# (base) jlkang@hnu2024 Tue Nov 04 2025 21:57:39 ~/WI_convergence/OG0006045
less final_alignment.fa.BUSTED.json|head -n 4|tail -n 1|perl -alne 's/\[|]//g;s/\,//g;my @a=split;for(my $i=0;$i<@a;$i++){my $j=$i+1;print "$j\t$a[$i]" if $a[$i]>=10}'
```

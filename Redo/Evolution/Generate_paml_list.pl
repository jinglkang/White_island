#!/usr/bin/perl
use strict;
use warnings;
use File::Basename; 

my $ano="all_swissprot_diamond_ano.txt"; # annotation file
my $seq="sub_orth_id_final.txt"; # the name of seqs per sub_orth

my %ANO;
open ANO, $ano or die "can not open $ano\n";
while (<ANO>) {
	chomp;
	my @a=split;
	my $gene;
	if ($a[0]=~/ENSTRUG/) {
		$gene="Fugu_".$a[0];
	} elsif ($a[0]=~/ENSORLG/) {
		$gene="Medaka_".$a[0];
	} elsif ($a[0]=~/ENSXMAG/) {
		$gene="Platyfish_".$a[0];
	} elsif ($a[0]=~/ENSLOCG/) {
		$gene="Spottedgar_".$a[0];
	} elsif ($a[0]=~/ENSGACG/) {
		$gene="Stickleback_".$a[0];
	} elsif ($a[0]=~/ENSDARG/) {
		$gene="Zebrafish_".$a[0];
	} else {
		$gene=$a[0];
	}
	my @b=split /\|/, $a[1];
	(my $name)=$b[-1]=~/(.*)\_.*/;
	my $score=$a[-1];
	$ANO{$gene}={
		'name'  => $name,
		'score' => $score
	};
}

open SEQ, $seq or die "can not open $seq\n";
while (<SEQ>) {
	chomp;
	my @a=split;
	if (/Suborth/) {
		print "$_\n";
	} else {
		my $info=$a[0]."\t";
		for (my $i = 1; $i < @a-1; $i++) {
			my @genes=split /\,/, $a[$i];
			my ($targe, $score); # the final selected gene
			foreach my $gene (@genes) {
				if ($score) {
					if ($ANO{$gene}->{'score'} >= $score) {
						$targe=$gene;
						$score=$ANO{$gene}->{'score'};
					}
				} else {
					$targe=$gene;
					$score=$ANO{$gene}->{'score'};
				}
			}
			$info.=$targe."\t";
		}
		$info=~s/\s+$//;
		print "$info\n";
	}
}

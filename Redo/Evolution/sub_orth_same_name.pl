#!/usr/bin/perl
use strict;
use warnings;
use File::Basename; 

my @spes=qw(Acura Apoly Blenny Blue_eyed Common Daru Fugu Medaka Ocomp Padel Platyfish Pmol Spottedgar Stickleback Yaldwyn Zebrafish);
my $ano="all_swissprot_diamond_ano.txt"; # annotation file
my $gnm="sub_orth_genecount.txt"; # the seqs number per species per sub_orth
my $seq="sub_orth_id.txt"; # the name of seqs per sub_orth

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
	my $score=$a[0];
	$ANO{$gene}={
		'name'  => $name,
		'score' => $score
	};
}

my %ORT;
open GNM, $gnm or die "can not open $gnm\n";
while (<GNM>) {
	chomp;
	my @a=split;
	next if /^Orthogroup/;
	my $j;
	for (my $i = 1; $i < @a; $i++) {
		if ($a[$i]>0) {
			$j++;
		}
	}
	if ($j && $j==16) {
		$ORT{$a[0]}++;
	}
}

my $header="Suborth\t";
foreach my $spe (@spes) {
	$header.=$spe."\t";
}
$header.="\tname";

my $matrix="sub_orth_genecount_final.txt";
my $matrid="sub_orth_id_final.txt";

open MATRIX, ">$matrix" or die "can not create $matrix\n";
print MATRIX "$header\n";

open MATRID, ">$matrid" or die "can not create $matrid\n";
print MATRID "$header\n";

open SEQ, $seq or die "can not open $seq\n";
while (<SEQ>) {
	chomp;
	my @a=split;
	my $orth=$a[0];
	if ($ORT{$orth}) {
		my %hash;
		for (my $i = 1; $i < @a; $i++) {
			$a[$i]=~s/\,$//;
			my $id;
			my $spe;
			if ($a[$i]=~/eyed/) {
				$id="Blue_".$a[$i];
				$spe="Blue_eyed";
			} else {
				$id=$a[$i];
				($spe)=$a[$i]=~/(.*)\_/;
			}
#			print "$id\t";
			if ($ANO{$id}) {
				my $name=$ANO{$id}->{'name'};
#				print "$name\t";
				$hash{$name}->{$spe}=[] unless exists $hash{$name}->{$spe};
				push @{$hash{$name}->{$spe}}, $id;
			}
		}
		foreach my $name (sort keys %hash) {
			my $k;
			my ($info1, $info2);
			$info1=$orth."\t";
			$info2=$orth."\t";
			foreach my $spe (@spes) {
				$k++ if $hash{$name}->{$spe} && @{$hash{$name}->{$spe}}>0;
				if ($hash{$name}->{$spe} && @{$hash{$name}->{$spe}}>0) {
					my $h=@{$hash{$name}->{$spe}};
#					print "$h\t";
				}
			}
#			print "$k\t";
			if ($k==16) {
				foreach my $spe (@spes) {
					my $nb=@{$hash{$name}->{$spe}};
					$info1.=$nb."\t";
					my $info3;
					foreach my $gene (@{$hash{$name}->{$spe}}) {
						$info3.=$gene.",";
					}
					$info3=~s/\,$//;
					$info2.=$info3."\t";
				}
				$info1.="$name";
				print MATRIX "$info1\n";
				$info2.=$info2."\t$name";
				print MATRID "$info2\n";
			}
		}
	}
}

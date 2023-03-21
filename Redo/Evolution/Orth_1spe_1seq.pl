#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

my $orth="Orthogroups/Orthogroups.GeneCount.tsv";
my %hash1;
open ORTH, $orth or die "can not open $orth\n";
while (<ORTH>) {
	chomp;
	if (/^Orthogroup/) {
		next;
	} else {
		my @a=split;
		my $j;
		for (my $i = 1; $i < @a; $i++) {
			$j++ if $a[$i]>=1;
			last if $a[$i]==0;
		}
		$hash1{$a[0]}++ if $j && $j==@a-1;
	}
}
close ORTH;

my $dir="Ortho_1spe_1seq";
unless (-e $dir) {
	system("mkdir $dir");
}
# less Orthogroup_Sequences/OG0000000.fa|grep '>'|perl -alne 's/\>//g;if (/\_/){s/_\d+//;print}else{s/\d+//g;print}'|sort -u
my @fastas=<Orthogroup_Sequences/*.fa>;
foreach my $fasta (@fastas) {
	my $orthid;
	my %hash2;
	($orthid)=$fasta=~/Orthogroup_Sequences\/(.*)\.fa/;
	if ($hash1{$orthid}) {
		my ($spe, $id);
		open ORTH, $fasta or die "can not open $fasta\n";
		while (<ORTH>) {
			chomp;
			if (/>/) {
				s/>//;
				if (/\_/) {
					$id =$_;
					s/_\d+$//;
					$spe=$_;
				} else {
					$id =$_;
					s/\d+$//g;
					$spe=$_;
				}
			} else {
				my ($seq, $len);
				$seq=$_;
				$len=length($_);
				if ($hash2{$spe}) {
					my $oldlen=$hash2{$spe}->{'len'};
					if ($len>=$oldlen) {
						$hash2{$spe}={
							'id' =>$id,
							'seq'=>$seq,
							'len'=>$len
						};
					}
				} else {
					$hash2{$spe}={
						'id' =>$id,
						'seq'=>$seq,
						'len'=>$len
					};
				}
			}
		}
		my $name=basename($fasta);
		my $newf="$dir/$name";
		open NEW, ">$newf" or die "can not create $newf\n";
		foreach my $key (sort keys %hash2) {
			my ($id, $seq);
			$id =$hash2{$key}->{'id'};
			$seq=$hash2{$key}->{'seq'};
			print NEW ">$id\n$seq\n";
		}
		close NEW;
	}
}

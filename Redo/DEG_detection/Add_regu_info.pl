#!/usr/bin/perl
use strict;
use warnings;
use File::Basename; 

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

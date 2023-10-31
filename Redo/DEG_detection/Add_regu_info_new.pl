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

my @yalups=<../Yaldwyn/*_up.txt>;
my @yaldws=<../Yaldwyn/*_dw.txt>;

my %hash1=(
        'Cs'=>'C',
        'Vs'=>'V1',
        'Vn'=>'V2->1'
        );

my %hash2=(
        'Blenny'=>'Crested blenny',
        'Blueeyed'=>'Blue-eyed triplefin',
        'Common'=>'Common triplefin',
        'Yaldwyn'=>'Yaldwin\'s triplefin'
        );

my @combis=(@bleups, @bledws, @bluups, @bludws, @comups, @comdws, @yalups, @yaldws);

foreach my $comb (@combis) {
        my ($spe, $reg)=$comb=~/\.\.\/(.*)\/(.*)\.txt/;
        next if $reg=~/plot/i;
        my @a=split /\_/, $reg;
        #my ($site1, $site2, $cha)=$reg=~/.*\_(.*?)\_(.*?)\_(.*?)/;
        my ($site1, $site2, $cha)=($a[1], $a[2], $a[3]);
        my $regNew=$hash1{$site1}."_".$hash1{$site2}."_".$cha;
        my $speNew=$hash2{$spe};
        open COMB, $comb or die "can not open $comb\n";
        while (<COMB>) {
                chomp;
                my $gene=$_;
                if ($regu{$speNew}->{$gene}) {
                        $regu{$speNew}->{$gene}.=$regNew.";";
                } else {
                        $regu{$speNew}->{$gene}=$regNew.";";
                }
        }
        close COMB;
}

my $enri=$ARGV[0];
open ENRI, $enri or die "can not open $enri\n";
while (<ENRI>) {
        chomp;
        my @a=split /\t/;
        my $spe=$a[0];
        my $gene=$a[1];
        my $info;
        if ($regu{$spe}->{$gene}) {
                $regu{$spe}->{$gene}=~s/\;$//;
                print "$_\t$regu{$spe}->{$gene}\n";
        }
}

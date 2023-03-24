#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

my $deg=$ARGV[0];
my $csv=$ARGV[1];
my $nm =$ARGV[2];
my $up =$nm."_up.txt";
my $dw =$nm."_dw.txt";

open UP, ">$up" or die "can not create $up\n";
open DW, ">$dw" or die "can not create $dw\n";

my %hash;
open DEG, $deg or die "can not open $deg\n";
while (<DEG>) {
        chomp;
        $hash{$_}++;
}

open CSV, $csv or die "can not open $csv\n";
while (<CSV>) {
        chomp;
        my @a=split /\,/;
        my $gene=$a[0];
        $gene=~s/\"//g;
        if ($hash{$gene}) {
                if ($a[2]>0) {
                        print UP "$gene\n";
                } elsif ($a[2]<0) {
                        print DW "$gene\n";
                }
        }
}

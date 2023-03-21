#!/usr/bin/perl
use strict;
use warnings;

system("mkdir Orth_Seq_changeHeader") unless -e "Orth_Seq_changeHeader";
my @fas=<Orthogroup_Sequences/*.fa>;
foreach my $fa (@fas) {
	(my $name)=$fa=~/Orthogroup_Sequences\/(.*\.fa)/;
    open FIL1, "$fa" or die "can not open $fa";
    open FIL2, ">Orth_Seq_changeHeader/$name" or die "can not create Orth_Seq_changeHeader/$name";
    while (<FIL1>) {
    	chomp;
    	if (/>/) {
    		s/>//;
            if (/ENSTRUG/) {
            	my $header="Fugu_".$_;
                print FIL2 ">$header\n";
            } elsif (/ENSORLG/) {
                my $header="Medaka_".$_;
                print FIL2 ">$header\n";
            } elsif (/ENSXMAG/) {
            	my $header="Platyfish_".$_;
                print FIL2 ">$header\n";
            } elsif (/ENSLOCG/) {
                my $header="Spottedgar_".$_;
                print FIL2 ">$header\n";
            } elsif (/ENSGACG/) {
                my $header="Stickleback_".$_;
                print FIL2 ">$header\n";
            } elsif (/ENSDARG/) {
                my $header="Zebrafish_".$_;
                print FIL2 ">$header\n";
            } elsif (/Blue_eyed/) {
                s/\_//;
                print FIL2 ">$_\n";
            }else {
                my $header=$_;
                print FIL2 ">$header\n";
            }
        } else {
            print FIL2 "$_\n";
        }
    }
    close FIL1;
    close FIL2;
};

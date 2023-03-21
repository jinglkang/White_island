#!/usr/bin/perl
use strict;
use warnings;

system("mkdir single_copy") unless -e "single_copy";
my @fas=<Ortho_1spe_1seq/*.fa>;
foreach my $fa (@fas) {
        (my $name)=$fa=~/Ortho_1spe_1seq\/(.*\.fa)/;
        open FIL1, "$fa" or die "can not open $fa";
        open FIL2, ">single_copy/$name" or die "can not create single_copy/$name";
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
                        } else {
                                my $header=$_;
                                print FIL2 ">$header\n";
                        }
                } else {
                        print FIL2 "$_\n";
                }
        }
        close FIL1;
        close FIL2;
}

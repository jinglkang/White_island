#!/usr/bin/perl
use strict;
use warnings;

my $matr=$ARGV[0];
my $DEGs=$ARGV[1];

my @heads; my %hash;
my %rcod=(
        'Cn'=>'1Cn',
        'Vn'=>'2Vn',
        'Vs'=>'3Vs',
        'Cs'=>'4Cs',
        );
open MATR, $matr or die "can not open $matr\n";
while (<MATR>) {
        chomp;
        s/\r//g;
        if (/^\s+/) {
                @heads=split /\t/;
#               print "$heads[0]\n";
        } else {
                my @a=split /\t/;
                for (my $i = 1; $i < @a; $i++) {
                        my ($spe, $sou, $nb)=split /\_/, $heads[$i];
                        $hash{$a[0]}->{$heads[$i]}={
                                'SPE'=>$spe,
                                'SOU1'=>$sou,
                                'SOU2'=>$rcod{$sou}.$nb,
#                               'NB' =>$nb,
                                'EXP'=>$a[$i]
                        };
                }
        }
}

print "Gene\tExpression\tSpecies\tsource1\tsource2\n";
open DEG, $DEGs or die "can not open $DEGs\n";
while (<DEG>) {
        chomp;
        my $gene=$_;
        my ($spe, $sou1, $sou2, $nb, $exp);
        my $i=0;
        foreach my $head (@heads) {
                $i++;
                if ($i>1 && $hash{$gene}) {
                        $spe=$hash{$gene}->{$head}->{'SPE'};
                        $sou1=$hash{$gene}->{$head}->{'SOU1'};
                        $sou2=$hash{$gene}->{$head}->{'SOU2'};
#                       $nb =$hash{$gene}->{$head}->{'NB'};
                        $exp=$hash{$gene}->{$head}->{'EXP'};
                        print "$gene\t$exp\t$spe\t$sou1\t$sou2\n";
                }
        }
}

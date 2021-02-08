#!/usr/bin/perl -w

# Copyright (c) Goro Terai. All rights reserved.
#
# This source code or any portion thereof must not be
# reproduced or used in any manner whatsoever.

use strict;

unless($#ARGV == 0){
    die "usage: $0 [*.gz]\n";
}

my $Eth = 1e-50;

#my $counter = 0;

my %Best;
open(FH, "zcat $ARGV[0]|") || die;
while(<FH>){
    chomp;
    next if(/^\#/);
    
    my ($score, $name1, $start1, $alnSize1, $strand1, $seqSize1,
	$name2,	$start2, $alnSize2, $strand2, $seqSize2, $blocks, $EG2, $E) = split /\t/;
#    $EG2 =~ s/^EG2\=//;
    $E =~ s/^E\=//;
    
    next if($E > $Eth);
#    $counter++;
    
    if(!defined $Best{$name2}){
	$Best{$name2}{SC} = $score;
	$Best{$name2}{DATA} = $_;
    }
    else{
	if($Best{$name2}{SC} < $score){
	    $Best{$name2}{SC} = $score;
	    $Best{$name2}{DATA} = $_;
	    $Best{$name2}{C}++;
	}
	elsif($Best{$name2}{SC} == $score){
	    $Best{$name2}{C}++;
	    my $r = rand(1);
	    if($r < 1/$Best{$name2}{C}){ #ランダムでヒットを更新する。
#	    if($counter % 2 == 0){ # ランダムでヒットを更新する。
		$Best{$name2}{SC} = $score;
		$Best{$name2}{DATA} = $_;
	    }
	}
    }

}
close(FH);

foreach my $id (keys %Best){
    print "$Best{$id}{DATA}\n";
}

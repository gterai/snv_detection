#!/usr/bin/perl -w

# Copyright (c) Goro Terai. All rights reserved.
#
# This source code or any portion thereof must not be
# reproduced or used in any manner whatsoever.

use strict;
use FileHandle;

unless($#ARGV == 1){
    die "usage: $0 [*_mutation.txt] [*_mutation.txt]\n";
}

my %Mut1 = &readMutation($ARGV[0]);
my %Mut2 = &readMutation($ARGV[1]);

foreach my $chr (sort keys %Mut1){
    foreach my $p (sort {$a <=> $b} keys %{$Mut1{$chr}}){
	if(defined $Mut2{$chr}{$p}){
	    my ($type1, $type2) = ($Mut1{$chr}{$p}{TYPE}, $Mut2{$chr}{$p}{TYPE});
	    if($type1 eq $type2){
		print "common\t$Mut1{$chr}{$p}{DATA}\n";
	    }
	    else{
		print "common(?)\t$Mut1{$chr}{$p}{DATA}\n";
	    }
	}
        else{
	    print "uniq\t$Mut1{$chr}{$p}{DATA}\n";
	}
    }
}

sub readMutation($){
    my %h;
    
    my $fh = new FileHandle($_[0]) || die;
    my $chr;
    while(<$fh>){
	chomp;
	if(/^>(\S+)/){
	    $chr = $1;
	}
	else{
	    my @a = split /\t/;
	    my ($pos, $type) = ($a[0], $a[3]);
	    if($type ne "OK"){
		$h{$chr}{$pos}{DATA} = $_;
		$h{$chr}{$pos}{TYPE} = $type;
	    }
	}
    }

    return %h;
}

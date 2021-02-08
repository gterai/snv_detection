#!/usr/bin/perl -w

# Copyright (c) Goro Terai. All rights reserved.
#
# This source code or any portion thereof must not be
# reproduced or used in any manner whatsoever.

use strict;
use FileHandle;

unless($#ARGV >= 1){
    die "usage: $0 [genome.mfa] [read1.fastq\@read2.best] (read2.fastq\@read2.best)\n";
}

my %Genome = &readSequence($ARGV[0]); # �Ǹ�˥����å��˻Ȥ���

my %Map; # 0-based
for(my $f = 1; $f <= $#ARGV; $f++){
    
    my ($read_file, $bestmap_file) = split /\@/, $ARGV[$f];
    my $read_href = &readFastQ_ref($read_file);
    
    ## test
    #foreach my $id (keys %{$read_href}){
    #print "$id\t$read_href->{$id}\n";
    #}
    #exit(0);
    
    my $p = 0;
    if($bestmap_file =~ /\.gz$/){
	open(FH, "zcat $bestmap_file|") || die;
    }
    else{
	open(FH, $bestmap_file) || die;
    }
    while(<FH>){
	chomp;

	$p++;
	print STDERR "process: $p $bestmap_file\n" if($p % 1000 == 0);
	
	my ($score, $chr, $g_fm, $g_aln, $g_str, $g_len,
	    $r_id, $r_fm, $r_aln, $r_str, $r_len, $blocks, $EG2, $E) = split /\t/;
	
	my @read; # ���饤�󤵤줿��ʬ�Τߤ����
	if($r_str eq "+"){
	    @read = split "", substr($read_href->{$r_id}, $r_fm, $r_aln); 
	}
	else{
	    @read = split "", substr(&revComp($read_href->{$r_id}), $r_fm, $r_aln);
	}
	
	# blocks�ν���
	{
	    my $local_g_fm = $g_fm; # ���Υ����а���
	    my $rel_r_fm = 0; # ���饤�󤵤줿�꡼�����Ҿ�Τ����а���
	    my $g_ins_total; # �����å���
	    if($blocks =~ /\,/){
		# ��������»����
		my @b = split /\,/, $blocks;
		for(my $i = 0; $i <= $#b; $i++){
		    if($b[$i] =~ /:/){
			# ��������»�Υ��Υơ������
			my ($g_ins, $r_ins) = split ':', $b[$i];
			$g_ins_total += $g_ins;
			
			if($g_ins != 0 && $r_ins != 0){ # �ɤ��餫������0�Ǥʤ��Ȥ����ʤ���
			    die "Douple insertion found in $blocks\n";
			}
			else{
			    if($g_ins != 0){ # ���Υ�¦������(Read¦�˷�»��
				# g_ins >= 1, r_ins = 0
				for(my $i = 0; $i < $g_ins; $i++){
				    splice(@read, $rel_r_fm, 0, '-');
				}
				
				$local_g_fm += $g_ins;
				$rel_r_fm += $g_ins; # "-"��$g_ins�ο���������Ƥ��롣
				
			    }
			    else{ # ���Υ�¦�˷�»(Read¦������)
				# g_ins = 0, r_ins >= 1
				my $r_ins_seq = join "", @read[$rel_r_fm..$rel_r_fm+$r_ins-1];
				#print "TEST:($r_ins) $rel_r_fm..$rel_r_fm+$r_ins-1 $r_ins_seq\n";
				$Map{$chr}[$local_g_fm]{"I$r_ins_seq"}++;
				# �������줿��ʬ��Read�����롣
				splice(@read, $rel_r_fm, $r_ins);
				
				#$rel_r_fm -= $r_ins; # ���Τǡ����ʤ��Ƥ褤��
				
			    }
			}
			
#			$local_g_fm += $g_ins;
#			$local_r_fm += $r_ins;

		    }
		    elsif($b[$i] =~ /^\d+$/){
			$local_g_fm += $b[$i];
			$rel_r_fm += $b[$i];
		    }
		    else{
			die "Unexpected block annotation in $blocks\n";
		    }
		}

		# check
		my $array_len = scalar @read;
		#print "$r_aln+$g_ins_total\t$array_len\t$rel_r_fm\n";
	    }
	    else{
		# read�ϲ��⤤����ʤ���
	    }
	    
	}
	#print @read, "\n";
	
	my $g = 0;
	for(my $i = 0; $i <= $#read; $i++){
	    $Map{$chr}[$g_fm + $g]{$read[$i]}++;
	    $g++;
	}
#	print STDERR "$p\n";
#	last if($p == 1000000);
#	last if($p == 100);
    }
    close(FH);
    
}

#my $DEPTH_TH = 5; # for test
my $DEPTH_TH = 10;

foreach my $chr (keys %Map){
    print ">$chr\n";
    my @chr = split "", $Genome{$chr};
    for(my $i = 0; $i <= $#chr; $i++){
	my $i_disp = $i;
	if(defined $Map{$chr}[$i]){
	    my $depth = 0;
	    my $max_ins = 0;
	    my $max_ins_char = 0;
	    my $max_nuc = 0;       # 20201102�ɵ��������ݥ������˺Ǥ�¿�����饤�󤵤�����ο���read����
	    my $max_nuc_char = ""; # 20201102�ɵ��������ݥ������˺Ǥ�¿�����饤�󤵤�����
	    my $total_nuc = 0; # 20201106�������ݥ������˥��饤�󤵤������
	    foreach my $n (keys %{$Map{$chr}[$i]}){
		if($n !~ /^I/){ # �����ǤϤʤ�
		    $depth += $Map{$chr}[$i]{$n};
		    if($max_nuc < $Map{$chr}[$i]{$n}){
			$max_nuc = $Map{$chr}[$i]{$n};
			$max_nuc_char = $n;
		    }
		    $total_nuc += $Map{$chr}[$i]{$n} if($n ne "-"); # 20201106�������ݥ������˥��饤�󤵤������
		    
		}
		else{
		    #$sum_ins += $Map{$chr}[$i]{$n};
		    if($max_ins < $Map{$chr}[$i]{$n}){
			$max_ins = $Map{$chr}[$i]{$n};
			$max_ins_char = $n;
		    }
		}
	    }
	    
	    # ����ΰ㤤
	    my $type = "OK";
	    my $mut = "";
	    if($depth > $DEPTH_TH){
		if($max_nuc_char ne $chr[$i]){ # 20201102�ɵ����Ǥ�¿�����饤�󤵤줿���𤬡���ե���󥹤Ȱۤʤ롣
		    if($max_nuc_char eq "-"){
			$type = "Del";
			$mut = "$chr[$i]->$max_nuc_char";
		    }
		    else{
			if($total_nuc > 0 && $max_nuc/$total_nuc >= 0.9){  # 20201106 �����ɲ�
			    $type = "SNP";  # 20201106 ���������
			}
			else{
			    $type = "snp";  # 20201106 �ˤ����
			}
			$mut = "$chr[$i]->$max_nuc_char";  
		    }
		}
		elsif($depth / 2 < $max_ins){ # 20201102�ɵ����ޥåפ��줿�꡼�ɤ�Ⱦ���ʾ����������
		    $type = "Ins";
		    my $disp = $max_ins_char;
		    $disp =~ s/^I//;
		    $mut = "$disp before $chr[$i]"
		}
	    }
	    

	    print "$i_disp\t($chr[$i])";
	    
	    if($type eq "OK"){
		print "\t-\t$type\t$mut";
	    }
	    else{
		print "\t$chr\t$type\t$mut";
	    }
	    # ��������
	    my $seg = "";
	    my $AT_annot = "";
	    if($type ne "OK"){
		$seg .= join "", @chr[$i-10..$i-1];
		my $n = $chr[$i]; $n =~tr /A-Z/a-z/;
		$seg .= $n;
		$seg .= join "", @chr[$i+1..$i+10];
		# AT-rich���Υơ������
		my $gc_ratio = &calcGC(split "", $seg);
		if($gc_ratio < 0.3){
		    $AT_annot = "AT-rich";
		}
	    }
	    print "\t$seg\t$AT_annot";
	    
	    
	    
	    # detail
	    foreach my $n (keys %{$Map{$chr}[$i]}){
		print "\t$n:$Map{$chr}[$i]{$n}";
	    }
#	    print "\n";
	    print "($total_nuc\t$depth)\n";
	    
	}
    }
}

sub readSequence(){
    #my $href = readSequence_ref($_[0],$_[1]);
    my $href;
    if(defined $_[1]){
        $href = &readSequence_ref($_[0], $_[1]);
    }
    else{
        $href = &readSequence_ref($_[0]);
    }

    return %{$href};
}

sub readSequence_ref(){
    my ($delim, $id, %hash);

    if($#_ == 0){
        $delim = ' ';
    }
    else{
        $delim = "$_[1]";
    }

    my $fh;
    if($_[0] =~ /\.gz$/){
        $fh = new FileHandle("zcat $_[0]|") || die "Can't open $_[0] ($!)\n";;
    }
    else{
        $fh = new FileHandle($_[0]) || die "Can't open $_[0] ($!)\n";;
    }
    while(<$fh>){
        chomp;
        next if(/^\#/);
        if(/^>(.+)$/){
            #print "=$_\n";
            #print "=$1\n";
            #print STDERR "$delim\n";
            ($id) = split /$delim/, $1;
        }
        else{
            $hash{$id} .= $_;
        }
    }
    $fh->close();

    return \%hash;

}


sub readFastQ_ref(){
    my ($id, %hash);

    my $fh;
    if($_[0] =~ /\.gz$/){
        $fh = new FileHandle("zcat $_[0]|") || die "Can't open $_[0] ($!)\n";;
    }
    else{
        $fh = new FileHandle($_[0]) || die "Can't open $_[0] ($!)\n";;
    }
    while(<$fh>){
        chomp;
        if(/^\@(\S+)/){
            $id = $1;
            chomp($_ = <$fh>);
            $hash{$id} = $_;
            chomp($_ = <$fh>); # + line
            if(!/^\+/){
                die "Unexpected line $_ for $id\n";
            }
            chomp($_ = <$fh>); # quality score
        }
        else{
            die "Unexpected line $_\n";
        }
    }
    $fh->close();

    return \%hash;

}

sub revComp($){
    my $seq = $_[0];
    $seq = reverse $seq;
    $seq =~ tr/ATGCUatgcu/TACGAtacga/;

    return $seq;
}

sub calcGC{

    my $gc = 0;
    for(my $i = 0; $i <= $#_; $i++){
        $gc++ if($_[$i] =~ /^[GCgc]$/);
    }

    return $gc/scalar(@_);
}

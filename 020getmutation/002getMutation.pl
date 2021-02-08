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

my %Genome = &readSequence($ARGV[0]); # 最後にチェックに使う。

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
	
	my @read; # アラインされた部分のみを取得
	if($r_str eq "+"){
	    @read = split "", substr($read_href->{$r_id}, $r_fm, $r_aln); 
	}
	else{
	    @read = split "", substr(&revComp($read_href->{$r_id}), $r_fm, $r_aln);
	}
	
	# blocksの処理
	{
	    my $local_g_fm = $g_fm; # ゲノム絶対位置
	    my $rel_r_fm = 0; # アラインされたリード断片上のの相対位置
	    my $g_ins_total; # チェック用
	    if($blocks =~ /\,/){
		# 挿入、欠損あり
		my @b = split /\,/, $blocks;
		for(my $i = 0; $i <= $#b; $i++){
		    if($b[$i] =~ /:/){
			# 挿入、欠損のアノテーション
			my ($g_ins, $r_ins) = split ':', $b[$i];
			$g_ins_total += $g_ins;
			
			if($g_ins != 0 && $r_ins != 0){ # どちらか片方は0でないといけない。
			    die "Douple insertion found in $blocks\n";
			}
			else{
			    if($g_ins != 0){ # ゲノム側に挿入(Read側に欠損）
				# g_ins >= 1, r_ins = 0
				for(my $i = 0; $i < $g_ins; $i++){
				    splice(@read, $rel_r_fm, 0, '-');
				}
				
				$local_g_fm += $g_ins;
				$rel_r_fm += $g_ins; # "-"を$g_insの数だけ入れている。
				
			    }
			    else{ # ゲノム側に欠損(Read側に挿入)
				# g_ins = 0, r_ins >= 1
				my $r_ins_seq = join "", @read[$rel_r_fm..$rel_r_fm+$r_ins-1];
				#print "TEST:($r_ins) $rel_r_fm..$rel_r_fm+$r_ins-1 $r_ins_seq\n";
				$Map{$chr}[$local_g_fm]{"I$r_ins_seq"}++;
				# 挿入された部分をReadから削る。
				splice(@read, $rel_r_fm, $r_ins);
				
				#$rel_r_fm -= $r_ins; # 削るので、やんなくてよい。
				
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
		# readは何もいじらない。
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
	    my $max_nuc = 0;       # 20201102追記、下記ポジションに最も多くアラインされる塩基の数（read数）
	    my $max_nuc_char = ""; # 20201102追記、下記ポジションに最も多くアラインされる塩基
	    my $total_nuc = 0; # 20201106　下記ポジションにアラインされる塩基数
	    foreach my $n (keys %{$Map{$chr}[$i]}){
		if($n !~ /^I/){ # 挿入ではない
		    $depth += $Map{$chr}[$i]{$n};
		    if($max_nuc < $Map{$chr}[$i]{$n}){
			$max_nuc = $Map{$chr}[$i]{$n};
			$max_nuc_char = $n;
		    }
		    $total_nuc += $Map{$chr}[$i]{$n} if($n ne "-"); # 20201106　下記ポジションにアラインされる塩基数
		    
		}
		else{
		    #$sum_ins += $Map{$chr}[$i]{$n};
		    if($max_ins < $Map{$chr}[$i]{$n}){
			$max_ins = $Map{$chr}[$i]{$n};
			$max_ins_char = $n;
		    }
		}
	    }
	    
	    # 塩基の違い
	    my $type = "OK";
	    my $mut = "";
	    if($depth > $DEPTH_TH){
		if($max_nuc_char ne $chr[$i]){ # 20201102追記、最も多くアラインされた塩基が、リファレンスと異なる。
		    if($max_nuc_char eq "-"){
			$type = "Del";
			$mut = "$chr[$i]->$max_nuc_char";
		    }
		    else{
			if($total_nuc > 0 && $max_nuc/$total_nuc >= 0.9){  # 20201106 条件の追加
			    $type = "SNP";  # 20201106 厳しい基準
			}
			else{
			    $type = "snp";  # 20201106 緩い基準
			}
			$mut = "$chr[$i]->$max_nuc_char";  
		    }
		}
		elsif($depth / 2 < $max_ins){ # 20201102追記、マップされたリードの半数以上に挿入あり
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
	    # 周辺配列
	    my $seg = "";
	    my $AT_annot = "";
	    if($type ne "OK"){
		$seg .= join "", @chr[$i-10..$i-1];
		my $n = $chr[$i]; $n =~tr /A-Z/a-z/;
		$seg .= $n;
		$seg .= join "", @chr[$i+1..$i+10];
		# AT-richアノテーション
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

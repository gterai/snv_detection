#!/usr/bin/perl -w

my $FDIR = "../fastq/";
#my $FDIR = "/home/terai/TRAHED/PPins/from_ito/";
opendir(D, $FDIR) || die;
my @fastq = grep /\.fastq$/, readdir(D);
closedir(D);

my %S2F;
foreach my $file (@fastq){
    $file =~ /^(.+)_R[12]_001\.fastq$/ or die "Can't parse strain from $file.\n";
    my $strain = $1;
    
    $S2F{$strain}{$file} = 1;
    
}

my $BDIR = "../010besthit/";
#my $BDIR = "/home/terai/TRAHED/PPins/script/bestmap_by_score/";
foreach my $strain (sort keys %S2F){

    my $line;
    foreach my $fastq (sort keys %{$S2F{$strain}}){
        $line .= " $FDIR/$fastq\@$BDIR/$fastq.gmap.gz.best.gz";
    }

    print "perl 002getMutation.pl ../index/Merged.mfa $line 1> $strain\_mutation.txt 2> $strain\_mutation.stderr\n";
}

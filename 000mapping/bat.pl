

foreach my $fastq ("78G10_S8_L001_R1_001.fastq", "78G10_S8_L001_R2_001.fastq", "W166_S7_L001_R1_001.fastq", "W166_S7_L001_R2_001.fastq"){
    
    print "lastal -Q sanger ../index/LastDB ../fastq/$fastq -f TAB -P 10 | gzip -c > $fastq.gmap.gz\n";

}


#!/usr/bin/perl

package SanityCheck;

sub check_sample_table_format{
    my $sample_metadata_table = $_[0];
    open(my $metadata, "<$sample_metadata_table") or die ("sample metadata table file, $sample_metadata_table, does not exist in work directory\n");
    my %groups;
    my %samples;
    my $die=0;
    while (my $line = <$metadata>) {
	chomp $line;
	$line =~ s/\"//g;
	if ($line eq ''){ next;}
	if ($line =~ /^#/) {next;}
	my @coln = split(/\,/, $line);
	my $group_id = $coln[3].'_'.$coln[4];
	my $sample_id = $coln[1].'_'.$coln[0];
	my $label = $coln[5];
	#count colns;
	my $coln_no = 0;
	foreach my $i (@coln){
	    $i =~ s/^\s+|\s+$//g; # remove white space
	    if (length($i)>0) {
		$coln_no++;
	    }
	}
	$groups{$group_id}{'id'}{$sample_id}++;
	$groups{$group_id}{'label'}{$label}++;
	$groups{$group_id}{'coln'}{$coln_no}++;
	#check the format of sample types
	my @seq_type = ('illumina.pe', 'illumina.se', 'nanopore');
	if (! grep {/$coln[4]/} @seq_type) {
	    print("sample, $sample_id, has an unknown seq_type, $coln[4]!\n");
	    $die=1;
	}
	#check illegal letters
	foreach my $i (@coln){
            chomp($i);
            if ($i ne ''){
                if ($i =~ /[^a-zA-Z0-9\.]/) {
		    print("$i in the row of sample $sample_id contains illegal letters!\n");
		    $die=1;
		}
            }
		
	}
	#check total columns
	if ((scalar@coln) < 5) {
	    print("sample $sample_id has at least five columns separated by \',\'\n");
	    $die=1;
	}
    }
    foreach my $i (keys%groups){
	my %id = %{$groups{$i}{'id'}};
	my %label = %{$groups{$i}{'label'}};
	my @coln_counts = keys(%{$groups{$i}{'coln'}});
	foreach my $j (keys%id){
	    if ($id{$j}>1) {
		print ("Group $i, sample $j, is not unique!\n");
		$die=1;
	    }
	    
	}
	foreach my $j (keys%label){
	    if ($label{$j}>1) {
		print("Group $i, label $j, is not unique!\n");
		$die=1;
	    }
	    
	}
	if (scalar(@coln_counts)>1) {
	    print ("Group $i has unequal columns!\n");
	    $die=1;
	}
	
    }
    close $metadata;
    
    if ($die) {
	die("Sample info table fails sanity check\n");
    }
}

sub check_rawdata_files{
    my $sample_metadata_table = $_[0];
    my $rawdata = $_[1];
    open(my $metadata, "<$sample_metadata_table")
    or die ("sample metadata table file, $sample_metadata_table, does not exist in work directory\n");
    my %samples;
    my %labels;
    my $die=0;
    while (my $line = <$metadata>) {
	chomp $line;
	$line =~ s/\"//g;
	if ($line eq ''){
            next;
        }
	if ($line =~ /^#/) {
	    next;
	}else{
	    my @coln = split(/\,/, $line);
	    my $sample = $coln[1].'_'.$coln[0];
	    my $seq_type = $coln[4];
	    if ($seq_type =~ /\.pe$/) {
		my $R1 = $sample.'_R1.fastq.gz';
		if (! -e "$rawdata/$R1") {
		    print("$rawdata/$R1 does not exist!\n");
		    $die=1;
		}
		my $R2 = $sample.'_R2.fastq.gz';
		if (! -e "$rawdata/$R2") {
		    print("$rawdata/$R2 does not exist!\n");
		    $die=1;
		}	
		
	    }else{
		my $SE = $sample.'.fastq.gz';
		if (! -e "$rawdata/$SE") {
		    print("$rawdata/$SE does not exist!\n");
		    $die=1;
		}
	    }
	    
	    	
	}
    }
    if ($die) {
	die("Some raw data files are missing!\n");
    }
    
    close $metadata;
}


sub check_abundance_table{
    my $sample_metadata_table = $_[0];
    my $abun_table = $_[1];
    open(my $metadata, "<$sample_metadata_table")
    or die ("sample metadata table file, $sample_metadata_table, does not exist in work directory\n");
    open(my $abun_file, "<$abun_table") or die;
    my $title = <$abun_file>;
    chomp $title;
    my @title = split(/\t/, $title);
    shift(@title);
    my %samples;
    foreach my $i (@title){
	$samples{$i}++;
    }
    while (my $line = <$metadata>) {
	chomp $line;
	$line =~ s/\"//g;
	if ($line eq ''){
            next;
        }
	if ($line =~ /^#/) {
	    next;
	}else{
	    my @coln = split(/\,/, $line);
	    my $internal_id = $coln[1].'_'.$coln[0]; 
	    if (!exists($samples{$internal_id})) {
		print("Warnings: sample $internal_id is not included in the all_abun_table.tsv.\n");
	    }
	    
	}
    }
    
    close $metadata;
}

1;

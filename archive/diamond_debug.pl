#!/usr/bin/perl -I /home/stdell/Desktop/shotgun_pipeline/ubuntu/scripts/bin
use strict;
use warnings;
use NCBITaxonomy;
package DiamondAnalysis;

sub AMR_Diamond {
    my $input_dir = $_[0];
    my $output_dir = $_[1];
    my $sample_info = $_[2];
    my $bin_path = $_[3];
    my $threads = $_[4];
    my $database_folder = $_[5];
    mkdir("$output_dir/raw_results");
    mkdir("$output_dir/summary");
    foreach my $i (keys%{$sample_info}){
	my $seq_type = $sample_info -> {$i}->{'seq_type'};
	my $R1;
	if ($seq_type =~ /\.pe/) {
	    $R1 = "$input_dir/$i/R1_bk.fastq";
	    #my $R1_SE = "$input_dir/$i.R1.unpaired.fastq.gz";
	}else{
	    $R1 = "$input_dir/$i.fastq";
	}
	my $command = "diamond blastx -d $database_folder/diamond/NCBIATB.dmnd -q $R1 --outfmt 6 stitle qseqid slen -o $output_dir/raw_results/$i.results.txt ".
	       "--max-target-seqs 10 -e 0.000000000000001 --id 80 --query-cover 90 --threads $threads --quiet >/dev/null 2>&1";
	system($command);
	summarize_results("$output_dir/raw_results/$i.results.txt", "$input_dir/$i/read2hit.tsv", 'atb', "$output_dir/summary/$i.summary.tsv");
    }
    
    
}

sub viru_Diamond {
    my $input_dir = $_[0];
    my $output_dir = $_[1];
    my $sample_info = $_[2];
    my $bin_path = $_[3];
    my $threads = $_[4];
    my $database_folder = $_[5];
    mkdir("$output_dir/raw_results");
    mkdir("$output_dir/summary");
    foreach my $i (keys%{$sample_info}){
	my $seq_type = $sample_info -> {$i}->{'seq_type'};
	my $R1;
	if ($seq_type =~ /\.pe/) {
	    $R1 = "$input_dir/$i/R1_bk.fastq";
	    #my $R1_SE = "$input_dir/$i.R1.unpaired.fastq.gz";
	}else{
	    $R1 = "$input_dir/$i.fastq";
	    
	}
	my $command= "diamond blastx -d $database_folder/diamond/VFDBS.dmnd -q $R1 --outfmt 6 stitle qseqid slen -o $output_dir/raw_results/$i.results.txt ".
	       "--max-target-seqs 10 -e 0.000000000000001 --id 80 --query-cover 90 --threads $threads --quiet >/dev/null 2>&1";
	system($command);
	summarize_results("$output_dir/raw_results/$i.results.txt", "$input_dir/$i/read2hit.tsv", 'vf', "$output_dir/summary/$i.summary.tsv");

    }
    
}


sub summarize_results{
    my $diamond_file = shift;
    my $read2hit_file = shift;
    my $tag=shift;
    my $output_file=shift;
    my %read2taxa;
    open(my $f1, "<$read2hit") or die;
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	$read2taxa{$coln[0]}=$coln[3];
    }
    my %gene_abun;
    close $f1;
    open(my $f2, "<$diamond_file") or die;
    my @records;
    my $last_read = '';
    my %taxa_abun;
    while (my line = <$line>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	if (length($last_read)>0 && $last_read ne $coln[1]) {
	    my $gene=get_concensus(\@records, $tag);
	    my $taxa='NA';
	    if (exists($reads2taxa{$last_read})) {
		$taxa = $reads2taxa{$last_read};
	    }
	    $gene_abun{$taxa}{$gene}{'abun'}++;
	    $taxa_abun{$taxa}++;
	    @records=();
	    
	}
	push(@records, $line);
	$last_read = $coln[1];
    }
    
    my $gene=get_concensus(\@records, $tag);
    my $taxa='NA';
    if (exists($reads2taxa{$last_read})) {
	$taxa = $reads2taxa{$last_read};
    }
    $gene_abun{$taxa}{$gene}{'abun'}++;
    
    open(my $f2, ">$output_file") or die;
    print($f2 "species\tprotein\tread_counts\n");
    my @taxa = sort{$taxa_abun{$b}<=>$taxa_abun{$a}}keys%gene_abun;
    foreach my $i (@taxa){
	foreach my $j (keys%{$gene_abun{$i}}){
	    my $abun = $gene_abun{$i}{$j}{'abun'};
	    print($f2 "$i\t$j\t$abun\n");
	    
	}
    }
    close $f2;
    
    
}

sub get_concensus{
    my $records = $_[0];
    my $tag = $_[1];
    my %names;
    if ($tag eq 'vf') {
	foreach my $i (@{$records}){
	    my @coln = split(/   /,$i);
	    $names{$coln[1]}++;
	}
	
    }elsif($tag eq 'atb'){
	foreach my $i (@{$records}){
	    my @coln = split(/ \[/,$i);
	    @coln = split(/\ /, $coln[0]);
	    if ($i =~ /MULTISPECIES/) {
		shift(@coln);
		shift(@coln);
	    }else{
		shift(@coln);
	    }
	    my $new = join(" ", @coln);
	    $names{$new}++;
	}
    }
    my @name = sort{$names{$b}<=>$names{$a}}keys%names;
    return($name[0]);
}


1;

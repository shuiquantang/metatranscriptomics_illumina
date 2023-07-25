#!/usr/bin/perl -I /home/stdell/Desktop/shotgun_pipeline/ubuntu/scripts-mf4/bin
use strict;
use warnings;
package MetaPhlan;
use Parallel::Loops;

#MetaPhlan::run_mf4($trim_folder, $output_dir, $sample_info, $threads);
sub run_mf4{
    my $trim_folder = shift;
    my $output_folder = shift;
    my $sample_info = shift;
    my $threads = shift;
    my @samples = keys%{$sample_info};
    my $default_db_folder = "/home/ubuntu/shotgun_db/mf4";
    my $cpu = `nproc`;
    chomp $cpu;
    my %outputs;
    my $pl = Parallel::Loops->new($threads);
    $pl -> foreach (\@samples, sub{
    #foreach my $i (@samples){
	my $i = $_;
	my $output_file = "$output_folder/$i.txt";
	my $cmd = "metaphlan $trim_folder/$i.fastq.gz --bowtie2out $output_folder/$i.bt2.bz2 --nproc $cpu --input_type fastq -o $output_file > /dev/null 2>&1";
	#my $cmd = "metaphlan $trim_folder/$i.fastq.gz --bowtie2out $output_folder/$i.bt2.bz2 --nproc $cpu --input_type fastq -o $output_file --bowtie2db $default_db_folder > /dev/null 2>&1";	
	system($cmd);
	system("rm $output_folder/$i.bt2.bz2");
    #}
    });
    foreach my $i (@samples){
	my $output_file = "$output_folder/$i.txt";
	if (-e $output_file) {
	    my @id = split(/\//, $output_file);
	    my @name = split(/\.|\_/, $id[-1]);
	    $outputs{$output_file}=$name[1];
	}
    }
    
    my @output = sort {$a cmp $b}keys(%outputs);
    my $list = join(" ", @output);
    my $cmd = "merge_metaphlan_tables.py $list > $output_folder/metaphlan_output.tsv";
    system($cmd);
    $cmd = "rm $list";
    system($cmd);
    break_metaphlan_outputs($output_folder, "metaphlan_output.tsv");
    my $input_table = "$output_folder/8.mf4_strain.tsv";
    my $abun_table_dir = "abun_table";
    system("mkdir $abun_table_dir");
    prepare_inputs($input_table, $abun_table_dir);
}

sub prepare_inputs{
    my $input =  shift;
    my $output_dir = shift;
    open(my $f1, "<$input") or die;
    open(my $all, ">$output_dir/all_abun_table.tsv") or die;
    open(my $euk, ">$output_dir/eukaryote_abun_table.tsv") or die;
    open(my $prok, ">$output_dir/prokaryote_abun_table.tsv") or die;
    open(my $virus, ">$output_dir/virus_abun_table.tsv") or die;
    my $header = <$f1>; chomp $header;
    my @header = split(/\t/, $header);
    $header[0] = '#taxon_id';
    my $new_header = join("\t", @header);
    print($all "$new_header\n");
    print($euk "$new_header\n");
    print($prok "$new_header\n");
    print($virus "$new_header\n");
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	my $tax = $coln[0];
	my @ranks = split(/\|/, $tax);
	my @new_ranks;
	my $kingdom = 'prokaryote';
	foreach  my $i (@ranks){
	    my $rank = substr($i, 3, length($i)-3);
	    $rank =~ s/\_/\ /g;
	    $i = substr($i, 0, 3).$rank;
	    
	    if ($i =~ /^c__/) {
		next;
	    }
	    if ($i =~ /^t__/) {
		$i = 'z'.substr($i, 1, length($i)-1);
	    }
	    if ($i =~ /^k__/) {
		if ($i =~ /^k__Eukaryota/) {
		    $kingdom = 'eukaryote';
		}elsif($i =~ /^k__Viruses/){
		    $kingdom = 'virus';
		}
		
	    }
	    push(@new_ranks, $i);
	    
	}
	my $new_tax = join(";", @new_ranks);
	$coln[0] = $new_tax;
	my $new_line = join("\t",@coln);
	print($all "$new_line\n");
	if ($kingdom eq 'eukaryote') {
	    print($euk "$new_line\n");
	}elsif($kingdom eq 'virus'){
	    print($virus "$new_line\n");
	}else{
	    print($prok "$new_line\n");
	}
    }
    close $all;
    close $euk;
    close $prok;
    close $virus;
    
}

sub break_metaphlan_outputs{
    my $dir = shift;
    my $input = shift;
    my %fhs;
    open(my $kingdom, ">$dir/1.mf4_kingdom.tsv") or die;
    $fhs{'k'}=$kingdom;
    open(my $phylum, ">$dir/2.mf4_pylum.tsv") or die;
    $fhs{'p'}=$phylum;
    open(my $class, ">$dir/3.mf4_class.tsv") or die;
    $fhs{'c'}=$class;
    open(my $order, ">$dir/4.mf4_order.tsv") or die;
    $fhs{'o'}=$order;
    open(my $family, ">$dir/5.mf4_family.tsv") or die;
    $fhs{'f'}=$family;
    open(my $genus, ">$dir/6.mf4_genus.tsv") or die;
    $fhs{'g'}=$genus;
    open(my $species, ">$dir/7.mf4_species.tsv") or die;
    $fhs{'s'}=$species;
    open(my $strain, ">$dir/8.mf4_strain.tsv") or die;
    $fhs{'t'}=$strain;
    open(my $f1, "<$dir/$input") or die;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line =~ /^#/) {
	    next;
	}
	my @coln = split(/\t/, $line);
	if (scalar(@coln)>1) {
	    if ($coln[0] eq 'clade_name') {
		foreach my $i (values(%fhs)){
		    print($i "$line\n");
		}
	    }else{
		my $tax = $coln[0];
		my @rank = split(/\|/, $tax);
		my $last = pop@rank;
		if ($last =~ /^[a-z]__/) {
		    my $tag = substr($last, 0, 1);
		    my $fh = $fhs{$tag};
		    print($fh "$line\n");
		}
	    }
	}else{
	    next;
	}
	
	
    }
    foreach my $i (values(%fhs)){
	close $i;
    }
    
}

1;

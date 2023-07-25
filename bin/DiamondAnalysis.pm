#!/usr/bin/perl -I /home/stdell/Desktop/shotgun_pipeline/ubuntu/scripts/bin
use strict;
use warnings;
use NCBITaxonomy;
package DiamondAnalysis;
use Parallel::Loops;

sub AMR_Diamond {
    my $input_dir = $_[0];
    my $output_dir = $_[1];
    my $sample_info = $_[2];
    my $species_abun = $_[3];
    my $threads = $_[4];
    my $cpu = `nproc`;
    chomp $cpu;
    mkdir("$output_dir/raw_results");
    mkdir("$output_dir/summary");
    my $pl = Parallel::Loops->new($threads);
    my @samples = keys%{$sample_info};
    $pl -> foreach (\@samples, sub{
    #foreach my $i (keys%{$sample_info}){
	my $i = $_;
	my $R1 = "$input_dir/$i.fastq.gz";
	my $command = "diamond blastx -d /diamond_db/NCBIATB.dmnd -q $R1 --outfmt 6 stitle qseqid slen -o $output_dir/raw_results/$i.results.txt ".
	       "--max-target-seqs 10 -e 0.000000000000001 --id 80 --query-cover 90 --threads $cpu --quiet >/dev/null 2>&1";
	system($command);
    #}
    });
    foreach my $i (keys%{$sample_info}){
	my %abun;
	if (exists($species_abun->{$i})) {
	    %abun = %{$species_abun->{$i}};
	}
	summarize_results("$output_dir/raw_results/$i.results.txt", "$output_dir/summary/$i.summary.txt", \%abun);
    }
    system("rm -r $output_dir/raw_results/");
}

sub viru_Diamond {
    my $input_dir = $_[0];
    my $output_dir = $_[1];
    my $sample_info = $_[2];
    my $species_abun = $_[3];
    my $threads = $_[4];
    my $cpu = `nproc`;
    chomp $cpu;
    mkdir("$output_dir/raw_results");
    mkdir("$output_dir/summary");
    my $pl = Parallel::Loops->new($threads);
    my @samples = keys%{$sample_info};
    $pl -> foreach (\@samples, sub{
	my $i = $_;
    #foreach my $i (keys%{$sample_info}){
	my $R1 = "$input_dir/$i.fastq.gz";
	my $command= "diamond blastx -d /diamond_db/VFDBS.dmnd -q $R1 --outfmt 6 stitle qseqid slen -o $output_dir/raw_results/$i.results.txt ".
	       "--max-target-seqs 10 -e 0.000000000000001 --id 80 --query-cover 90 --threads $cpu --quiet >/dev/null 2>&1";
	system($command);
    #}
    });
    foreach my $i (keys%{$sample_info}){
	my %abun;
	if (exists($species_abun->{$i})) {
	    %abun = %{$species_abun->{$i}};
	}
	summarize_results("$output_dir/raw_results/$i.results.txt", "$output_dir/summary/$i.summary.txt", \%abun);
    }
    system("rm -r $output_dir/raw_results/");
    
}


sub summarize_results{
    my $diamond_file = shift;
    my $output_file=shift;
    my $sp_abun = shift;
    
    my %hit_counts;
    open(my $f1, "<$diamond_file") or die;
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	$hit_counts{$coln[0]}++;
    }
    close $f1;
    my %gene_abun;
    open(my $f2, "<$diamond_file") or die;
    my @records;
    my $last_read = '';
    my %species_hits;
    while (my $line = <$f2>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	if (length($last_read)>0 && ($last_read ne $coln[1])) {
	    my ($species, $gene) = get_best_hit(\@records, \%hit_counts, $sp_abun);
	    $gene_abun{$species}{$gene}{'abun'}++;
	    $species_hits{$species}++;
	    @records=();
	    
	}
	push(@records, $line);
	$last_read = $coln[1];
    }
    close $f2;
    if (scalar(@records)>0) {
	my ($species, $gene) = get_best_hit(\@records, \%hit_counts, $sp_abun);
	$gene_abun{$species}{$gene}{'abun'}++;
	$species_hits{$species}++;
    }
    
    open(my $f3, ">$output_file") or die;
    print($f3 "protein\torigin_species\tfound_in_taxa_profiling\tread_counts\n");
    my @taxa = sort{$species_hits{$b}<=>$species_hits{$a}}keys%gene_abun;
    foreach my $i (@taxa){
	
	foreach my $j (keys%{$gene_abun{$i}}){
	    my $abun = $gene_abun{$i}{$j}{'abun'};
	    my $present=0;
	    if (length($i)>0) {
		foreach my $k (keys(%{$sp_abun})){
		    if ($k =~ /$i/) {
			$present=1;
		    }
		}
	    }
	    if ($present==1) {
		print($f3 "$j\t$i\tYes\t$abun\n");
	    }else{
		print($f3 "$j\t$i\tNo\t$abun\n");
	    }
	    
	    
	}
    }
    close $f3;
    
    
}

sub get_best_hit{
    my $records = $_[0];
    my $hit_counts = $_[1];
    my $species_abun = $_[2];
    my %names;
    my %counts_hash;
    my %good_hits;
    my %better_hits;
    
    foreach my $i (@{$records}){
	my @arr = split(/\t/, $i);
	my $hit = $arr[0];
	my $counts = $hit_counts->{$hit};
	if ($hit =~ /   /) {
	    if ($hit =~ /\|/) {
		my @coln = split(/   /, $hit);
		    if (scalar(@coln)>2) {
			my $species = $coln[2];
			my @name = split(/\[|\]|\ /, $species);
			if (scalar(@name)<4) {
			    $species = $name[1];
			}else{
			    $species = "$name[1] $name[2]";
			}
			my $gene = $coln[1];
			$good_hits{$counts}{$species}{$gene}++;
			if (exists($species_abun->{$species})) {
			    $better_hits{$counts}{$species}{$gene}++;
			}
			
		    }else{
			my $species = $coln[1];
			my @name = split(/\[|\]|\ /, $species);
			if (scalar(@name)<4) {
			    $species = $name[1];
			}else{
			    $species = "$name[1] $name[2]";
			}
			my @id = split(/\ /, $coln[0]);
			shift(@id);
			my $gene = join(" ", @id);
			$good_hits{$counts}{$species}{$gene}++;
			if (exists($species_abun->{$species})) {
			    $better_hits{$counts}{$species}{$gene}++;
			}
		    }
		
            }
        }else{
	    my @coln = split(/\[|\]/, $hit);
	    my $gene = substr($coln[0],0,length($coln[0])-1);
	    my @id = split(/\ /, $gene);
	    shift(@id);
	    $gene = join(" ", @id);
	    my $species = $coln[-1];
	    $good_hits{$counts}{$species}{$gene}++;
	    if (exists($species_abun->{$species})) {
		$better_hits{$counts}{$species}{$gene}++;
	    }
	    
	}

    }
    my $species;
    my $gene;
    if (scalar(keys(%better_hits))>0) {
	my @counts = sort {$b <=> $a} keys(%better_hits);
	my @species = keys%{$better_hits{$counts[0]}};
	@species = sort{$species_abun ->{$b} <=> $species_abun ->{$b}}@species;
	$species = $species[0];
	my @genes = keys%{$better_hits{$counts[0]}{$species}};
	@genes = sort @genes;
	$gene = $genes[0];
    }else{
	my @counts = sort {$b <=> $a} keys(%good_hits);
	my @species = keys%{$good_hits{$counts[0]}};
	$species = $species[0];
	my @genes = keys%{$good_hits{$counts[0]}{$species}};
	@genes = sort @genes;
	$gene = $genes[0];
    }
    return($species, $gene);
}


1;

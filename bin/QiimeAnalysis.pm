#!/usr/bin/perl
use strict;
use warnings;
package QiimeAnalysis;

# rename_sample_in_abun_table($abun_table, $sample_info, $analysis_group_info, $new_abun_table, $i)
sub rename_sample_in_abun_table{
    # return a subset of abundance table for samples in this group and change the sample_id from internal_id to customer labels.
    my $abun_table = $_[0];
    my $group = $_[1];
    my $analysis_group_info = $_[2];
    my $new_abun_table = $_[3];
    my $group_id = $_[4];
    my @sample_labels = @{$analysis_group_info -> {$group_id}->{'sample_order'}};
    my @index;
    open(my $f1, "<$abun_table") or die;
    my $taxa_no=0;
    my %sample2count_no;
    my %internal_id2taxa_abun;
    my @internal_id;
    my @taxa;
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	if (scalar(@coln)<2) {
	    next;
	}
	if ($line =~ /^#/) {
	    shift(@coln);
	    @internal_id = @coln;
	}else{
	    my $taxa = shift(@coln);
	    push(@taxa, $taxa);
	    for (my $i=0;$i<scalar(@coln);$i++){
		$internal_id2taxa_abun{$taxa}{$internal_id[$i]}= $coln[$i];
		$sample2count_no{$internal_id[$i]}+=$coln[$i];
	    }
	    
	}
	
    }
    close $f1;
    
    
    # remove empty samples otherwise the taxaplot won't work
    my @new_samples;
    foreach my $i (@sample_labels){
	my $internal_id = $group ->{$group_id}->{$i}->{'internal_id'};
	if (exists($sample2count_no{$internal_id})) {
	    if ($sample2count_no{$internal_id}==0) {
		next;
	    }else{
	        push(@new_samples, $i);
	    }
	}
	
	
    }
    
    open(my $f2, ">$new_abun_table") or die;
    my $title_line = join("\t", @new_samples);
    print ($f2 "#OTU ID\t$title_line\n");
    foreach my $i (@taxa){
	my @coln = ($i);
	my $sum = 0;
	foreach my $j (@new_samples){
	    my $internal_id = $group ->{$group_id}->{$j}->{'internal_id'};
	    my $abun = $internal_id2taxa_abun{$i}{$internal_id};
	    if (length($abun)>0) {
		$sum+=$abun;
	    }
	    
	    
	    push(@coln, $abun);
	}
	if ($sum>0) { # only record the taxa that has abundance >0;
	    my $line = join("\t", @coln);
	    print($f2 "$line\n");
	    $taxa_no++; # count existing taxa
	}
	
	
    }

    close $f2;
    return ($taxa_no);
}

sub biom_conversion{
    my $abun_table = shift;
    my $dir = shift;
    my $phylogeny = shift;
    mkdir($dir);
    my @ranks= @{$phylogeny->{'ranks'}};
    my @rank_prefix=@{$phylogeny->{'prefix'}};
    my %prefix2rank;
    for (my $i=0;$i<scalar(@rank_prefix);$i++) {
	$prefix2rank{$rank_prefix[$i]}=$ranks[$i];
    }
    my @index = @{$phylogeny -> {'ranks_to_use'}};
    my @ranks_to_use;
    foreach my $i (@index){
	push (@ranks_to_use, $ranks[$i-1]);
    }
    my %abun_table_at_ranks;
    my %total_counts_at_ranks; # the total number of reads assigned to a sample
    #->{'genus'->}->{'k__bacteria;p__firmicutes;..g__Dehalobacter'}->(0.1,0,0.2,0,0.5,...)
    my @sample_order;
    my $head_line;
    open(my $f1, "<$abun_table") or die;
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	if (scalar(@coln)<2) {
	    next;
	}
	if ($line =~ /^#/) {
	    $head_line = $line;
	    shift(@coln);
	    @sample_order = @coln;
	}else{
	    my $taxa = $coln[0];
	    shift(@coln);
	    my @abun = @coln;
	    my @taxa = split(/\;/, $taxa);
	    for (my $i=scalar(@taxa); $i>0;$i--){
		my $index = $i-1;
		my $prefix = substr($taxa[$index],0,3);
		my $rank = $prefix2rank{$prefix};
		my @new_taxa = @taxa[0..$index];
		my $new_taxa = join(";", @new_taxa);
		for(my $j=0;$j<scalar(@sample_order);$j++){
		    $abun_table_at_ranks{$rank}{$new_taxa}[$j]+=$abun[$j];
		    $total_counts_at_ranks{$rank}[$j]+=$abun[$j];
		}
		
	    }
	}
	
    }
    # output abun_table of read counts.
    
    foreach my $i (@ranks_to_use){
	mkdir("$dir/$i");
	# output the abundance table
	my $abun_table_file = "$dir/$i/$i.tsv";
	my $biom_file = "$dir/$i/abun_table.biom";
	open(my $f3, ">$dir/$i/abun_table_read_counts.tsv") or die;
	open(my $f4, ">$abun_table_file") or die;
	print ($f3 "# Constructed from biom file\n");
	print ($f3 "$head_line\n");
	print ($f4 "# Constructed from biom file\n");
	print ($f4 "$head_line\n");
	my %abun = %{$abun_table_at_ranks{$i}};
	my @taxa = sort {$a cmp $b} keys(%abun);
	my @total_counts = @{$total_counts_at_ranks{$i}};
	foreach my $t (@taxa){
	    my @abun_counts = @{$abun{$t}};
	    my @abun;
	    for (my $i=0; $i<scalar(@abun_counts);$i++){
		if ($total_counts[$i]==0) {
		    $abun[$i]=0;
		}else{
		    $abun[$i]= $abun_counts[$i]/$total_counts[$i];
		    $abun[$i]= sprintf("%.6f", $abun[$i]);
		}
	    }
	    my @coln_counts = ($t);
	    push(@coln_counts, @abun_counts);
	    my $line_counts = join("\t", @coln_counts);
	    print($f3 "$line_counts\n");
	    
	    my @coln = ($t);
	    push(@coln, @abun);
	    my $line = join("\t", @coln);
	    print($f4 "$line\n");
	}
	close $f3;
	close $f4;
	# convert the abundance table to biom format
	my $command ="biom convert -i $abun_table_file -o $biom_file --table-type=\"OTU table\" --to-json";
	
	system ($command);
	system ("cp $abun_table_file $dir/$i/abun_table.tsv");
    }
    
}

#QiimeAnalysis::taxa_composition_barplot($abun_table_folder, $barplots);
sub taxa_composition_barplot{
    my $abun_table_folder = shift;
    my $barplot_folder = shift;
    my $phylogeny = shift;
    my @ranks= @{$phylogeny->{'ranks'}};
    my @index = @{$phylogeny -> {'ranks_to_use'}};
    my $size = scalar(@index);
    my @ranks_to_use;
    my @abun_tables;
    foreach my $i (@index){
	push (@ranks_to_use, $ranks[$i-1]);
    }
    foreach my $i (@ranks_to_use){
	my $abun_table = "$abun_table_folder/$i/$i.tsv";
	push (@abun_tables, $abun_table);
    }
    shift(@abun_tables);
    my $str = join(",", @abun_tables);
    my $command = "plot_taxa_summary.py -i $str -o $barplot_folder -l Phylum,Class,Order,Family,Genus,Species --chart_type bar >/dev/null 2>&1";
    system($command);
}


sub count_host_DNA{
    my $abun_file=shift;
    open(my $f1, "<$abun_file") or die;
    open(my $f2, ">host_dna_abun.csv") or die;
    my %abun;
    my @list=('Bacteria', 'Archaea', 'Virus', 'Fungi', 'Homo sapiens', 'Mus', 'Other eukaryotes');  
    my $title = <$f1>;
    print ($f2 "$title");
    my @samples = split(/\t/, $title);
    shift(@samples);
    for (my $i=0;$i<scalar(@list);$i++){
	for (my $j=0;$j<scalar(@samples);$j++){
	    $abun{$list[$i]}[$j]=0;
	}
	
    }
    
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	if (scalar(@coln)<1) {
	    next;
	}
	
	my $taxa = shift(@coln);
	my @ranks = split(/\;/, $taxa);
	my $key ='';
	foreach my $i (@ranks){
	    if ($i eq 'k__Bacteria') {
		$key = $list[0];
	    }elsif($i eq 'k__Archaea'){
		$key = $list[1];
	    }elsif($i eq 'k__Viruses'){
		$key = $list[2];
	    }elsif(($i eq 'p__Ascomycota') || ($i eq 'p__Basidiomycota') ){
		$key = $list[3];
	    }elsif($i eq 's__Homo sapiens'){
		$key = $list[4];
	    }elsif($i eq 'g_Mus'){
		$key = $list[5];
	    }elsif ($i eq 'k__Eukaryota'){
		$key = $list[6];
	    }
	    
	}
	foreach (my $i=0;$i<scalar(@coln);$i++){
	    $abun{$key}[$i]+=$coln[$i];
	}
	
    }
    
    foreach my $i (@list){
	my @abun = @{$abun{$i}};
	my $sum=0;
	foreach my $i (@abun){
	    $sum+=$i;
	}
	if ($sum>0) {
	    unshift(@abun, $i);
	    my $line = join("\t", @abun);
	    print($f2 "$line\n");
	}
    }
    close $f1;
    close $f2;
}

1;
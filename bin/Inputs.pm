#!/usr/bin/perl
use strict;
use warnings;
package Inputs;

#Inputs::read_abun_table($taxa_abun_table, \%species_abun);
sub print_and_execute{
    my $cmd = shift;
    my $log = shift;
    # Get the current local time
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime();
    # Format the year correctly (localtime returns the year as the number of years since 1900)
    $year += 1900;
    # Format the month correctly (localtime returns months as 0-11)
    $mon += 1;
    # Print the formatted date and time
    my $time = "Current date and time: $year-$mon-$mday $hour:$min:$sec";
    open(my $f1, ">>$log") or die;
    print($f1 "\n#-----------------------------------------------------\n$time\n$cmd\n#-----------------------------------------------------\n\n");
    #print("\n#-----------------------------------------------------\n$cmd\n#-----------------------------------------------------\n\n");
    close $f1;
    system($cmd);
}

sub log_only{
    my $cmd = shift;
    my $log = shift;
    # Get the current local time
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime();
    # Format the year correctly (localtime returns the year as the number of years since 1900)
    $year += 1900;
    # Format the month correctly (localtime returns months as 0-11)
    $mon += 1;
    # Print the formatted date and time
    my $time = "Current date and time: $year-$mon-$mday $hour:$min:$sec";
    open(my $f1, ">>$log") or die;
    print($f1 "\n#-----------------------------------------------------\n$time\n$cmd\n#-----------------------------------------------------\n\n");
    #print("\n#-----------------------------------------------------\n$cmd\n#-----------------------------------------------------\n\n");
    close $f1;
}

sub read_abun_table{
    my $file = shift;
    my $species_abun = shift;
    open(my $f1, "<$file") or die;
    my @titles;
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	if ($line =~ /^#/) {
	    @titles = @coln;
	}else{
	    my $taxa = $coln[0];
	    my @ranks = split(/\;/, $taxa);
	    my $sp='';
	    foreach my $i (@ranks){
		if ($i=~ /^s__/) {
		    $sp = substr($i,3, length($i)-3);
		}
	    }
	    for(my $i=1; $i<scalar(@titles); $i++){
		if ($coln[$i]>0) {
		    $species_abun->{$titles[$i]}{$sp} = $coln[$i];
		}
		
		
	    }
	}
	
    }
    close $f1;
    
    
}

sub read_sample_table{
    my $table = shift;
    open(my $sample_table, "<$table") or die ("failed to open MasterTable");
    my %groups;
    my %group_info;
    my %sample_info;
    my %error_model_groups;
    my %error_models;
    my @cat_titles;
    while (my $line = <$sample_table>) {
        chomp $line;
        $line =~ s/\"//g;
	$line =~ s/\r//g;
        if (length($line) == 0) { #skip empty lines
            next;
        }
        if ($line =~/^#/) {
            my @line = split(/\,|\t/, $line);
	    @cat_titles = splice(@line, 6);
        }else{
            my @line = split(/\,|\t/, $line);
            my $group_id = "$line[3]\.$line[4]";
            my $run_id = $line[2];
            my $internal_id = "$line[1]\_$line[0]";
            my $sample_label = '';
            if (length($line[5])>0) {
                $sample_label = $line[5];
            }else{
                $sample_label = $internal_id;
            }
	    push (@{$group_info{$group_id}{'sample_order'}}, $sample_label);
            $groups{$group_id}{$sample_label}{'internal_id'}=$internal_id;
	    $sample_info{$internal_id}{'label'}=$sample_label;
            $group_info{$group_id}{'seq_type'} = $line[4];
	    $sample_info{$internal_id}{'seq_type'}=$line[4];
	    
	    my @category = splice(@line, 6);
	    my $category_size=0;
	    foreach my $i (@category){
		if (length($i)>0) {
		    $category_size++;
		}
		
	    }
	    $group_info{$group_id}{'category_size'} = $category_size;
            $groups{$group_id}{$sample_label}{'categories'}=\@category;
            $sample_info{$internal_id}{'run_id'}=$run_id;
	    
	    # paired-end reads or not
	    if ($line[4] =~ /\.pe$/) {
		my $R1=$internal_id.'_R1.fastq.gz';
		$sample_info{$internal_id}{'R1'}=$R1;
		my $R2=$internal_id.'_R2.fastq.gz';
		$sample_info{$internal_id}{'R2'}=$R2;
	    }else{
		my $R1=$internal_id.'.fastq.gz';
		$sample_info{$internal_id}{'R1'}=$R1;
	    }
        }
    }
    return (\%groups, \%group_info, \%sample_info, \@cat_titles);
}

# determine the number of parallel processes for each step.

sub parallel_process_allocation{
    # determine the number of parallel processes for each step.
    my $max_memory_consumption_per_process = shift;
    my $max_memory = `free -m | awk 'FNR == 2 {print \$2}'`; # max memory available in the system in Mb
    chomp($max_memory);
    my $max_core = `nproc`;
    chomp($max_core);
    if ($max_core==1) { # if the system only has one cpu
	return 1;
    }
    
    my $max_core_to_use;
    my $a = $max_core;
    my $b = int($max_memory/$max_memory_consumption_per_process);
    if ($a>=$b) {
	$max_core_to_use=$b;
    }else{
	$max_core_to_use=$a;
    }
    if ($max_core_to_use>1) { # save at least two cpu cores for system
	return $max_core_to_use;
    }else{
	return 1; # only one can use
    }
}

# rank a parameter file to determine what phylogeny ranks to use in the analysis;
sub ranks_to_use {
    my $script_dir = shift;
    my %phylogeny;
    open(my $f1, "<$script_dir/ranks_to_use.csv") or die;
    while (my $line = <$f1>) {
	chomp $line;
	$line =~ s/\"//g;
	$line =~ s/\r//g;
	my @coln = split(/\,/, $line);
	if (scalar(@coln)<1) {
	    next;
	}else{
	    my $key = shift(@coln);
	    $phylogeny{$key}=\@coln;
	}
    }
    return(\%phylogeny);
}


sub create_mapping_file{
    my $group_id = shift;
    my $group_info = shift;
    my $group = shift;
    my $cat_titles=shift;
    my $min_category_no=scalar(@{$cat_titles});
    open(my $map_file, ">qiime.mapping.file.txt") or die;
    open(my $sample_list, ">sample.list.txt") or die;
    my @sample_order = @{$group_info -> {$group_id} -> {'sample_order'}};
    my @title = ('#SampleID', 'BarcodeSequence', 'LinkerPrimerSequence');
    push (@title, @{$cat_titles});
    my $title = join("\t", @title);
    print ($map_file "$title\n");
    foreach my $i (@sample_order){
	my @categories = @{$group->{$group_id}->{$i}->{'categories'}};
	my $size=0;
	# checking how many categories
	foreach my $i (@categories){
	    $i =~ s/\s+//g;
	    if (length($i)>0) {
		$size++;
	    }
	    
	}
	if ($size<$min_category_no) {
	    $min_category_no = $size;
	}
	my @line = ($i, 'NA', 'GTGCCAGCMGCCGCGGTAA');
	push (@line, @categories);
	my $line = join("\t", @line);
	print ($map_file "$line\n");
	print ($sample_list "$i\n");
    }
    close $map_file;
    close $sample_list;
    #prepare category mapping files for heatmap analysis
    for (my $i=0; $i<$min_category_no; $i++){
	my $cat = $cat_titles->[$i];
	open(my $map_file, ">$cat.mapping.file.txt") or die;
	foreach my $j (@sample_order){
	    my @categories = @{$group->{$group_id}->{$j}->{'categories'}};
	    my @line = ($j, 'NA', 'GTGCCAGCMGCCGCGGTAA');
	    push (@line, $categories[$i]);
	    my $line = join("\t", @line);
	    print ($map_file "$line\n");
	}
    }
    return ($min_category_no);
}

1;

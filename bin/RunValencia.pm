#!/usr/bin/perl -I /home/stdell/Desktop/Dada2.pipeline/Bin/
use strict;
use warnings;
use Inputs;

package RunValencia;

sub run_valencia{
    my $species_file = $_[0];
    my $output_dir = $_[1];
    my $sample2species;
    my $sample2species_bac;
    my %tot_reads;
    my %lineage_from_sp;    
    open(my $f1, "<$species_file") or die;
    my @species_order;
    my @samples;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line eq '# Constructed from biom file') {
	    next;
	}
	my @coln = split(/\t/, $line);
	if (scalar(@coln)>1) {
	    if ($line=~ /^#OTU ID/) {
		@samples = @coln;
		$samples[0]='Taxon_ID';
	    }else{
		my $lineage = $coln[0];
		my $species = extract_species_name($lineage);
		push(@species_order, $species);
		my $sp_lineage  = extract_species_lineage($lineage);
		$lineage_from_sp{$species} = $sp_lineage;
		for (my $i=1; $i<scalar(@coln); $i++){
		    $sample2species->{$samples[$i]}{$species}+=$coln[$i];
		    if ($lineage =~ /k__Bacteria|k__Archaea/) { #only bacteria species are considered for valencia analysis
			$tot_reads{$samples[$i]}+=$coln[$i];
			$sample2species_bac->{$samples[$i]}{$species}+=$coln[$i];
		    }
		}
	    }
	}
    }
    # output abundance table for each sample
    my @species;
    system("mkdir $output_dir/abun_tables");
    foreach my $i (keys(%{$sample2species})){
	open(my $f2, ">$output_dir/abun_tables/$i.txt") or die;
	print($f2 "#species\tlineage\tabun_by_read_count\n");
	my %abun = %{$sample2species->{$i}};
	my $tot_reads = $tot_reads{$i};
	@species = keys(%abun);
	my @sp = sort{$abun{$b}<=>$abun{$a}} @species;
	foreach my $i (@sp){
	    my $abun = $abun{$i}/$tot_reads;
	    if ($abun>0) {
		$abun = sprintf("%.5f", $abun);
		my $lineage = $lineage_from_sp{$i};
		print($f2 "$i\t$lineage\t$abun\n");
	    }
	}
	close $f2;
    }
    #prepare valencia inputs, bacteria species only
    system("mkdir $output_dir/valencia");
    open(my $f3, ">$output_dir/valencia/valencia_input.csv") or die;
    @samples = keys(%{$sample2species_bac});
    my @bac_sp = keys(%{$sample2species_bac->{$samples[0]}});
    my @valencia_species = valencia_format(\@bac_sp);
    my $headline = join(",", @valencia_species);
    print($f3 "sampleID,read_count,$headline\n");
    foreach my $i (@samples){
	my $tot_read_count=$tot_reads{$i};
	my @coln = ($i, $tot_read_count);
	for (my $j=0;$j<scalar(@bac_sp);$j++){
	    my $sp = $bac_sp[$j];
	    my $counts = $sample2species_bac->{$i}->{$sp};
	    push(@coln, $counts);
	}
	my $line = join(",",@coln);
	print($f3 "$line\n");
    }
    close $f3;
    
    #run valencia
    my $valencia_input_file = "$output_dir/valencia/valencia_input.csv";
   
    my $cmd = "python3 /VALENCIA-master/Valencia.py -ref /VALENCIA-master/CST_centroids_012920.csv -i $valencia_input_file -o $output_dir/valencia/valencia_output";
    #print("$cmd\n");
    system($cmd);
}

sub valencia_format{
    my @species = @{$_[0]};
    my @valencia_sp;
    foreach my $i (@species){
	my $new_name;
	if ($i =~ /BVAB1|BVAB-1/) {
	    $new_name = "g_BVAB1";
	    push(@valencia_sp, $new_name);
	    next;
	}
	if ($i =~ /BVAV2|BVAB-2/) {
	    $new_name = "g_BVAB2";
	    push(@valencia_sp, $new_name);
	    next;
	}
	if ($i =~ /BVAV3|BVAB-3/) {
	    $new_name = "g_BVAB3";
	    push(@valencia_sp, $new_name);
	    next;
	}
	my @ids = split(/\(|\)/, $i);
	my $rank = "g";
	my $name = $i;
	if (scalar(@ids)>1) {
	    $rank = $ids[1];
	    $name=$ids[2];
	}
	
        
	my @sp = split(/\ /, $name);
	my $species = pop(@sp);
	my $genus = join(" ", @sp);
	if ($genus =~ /^(Lactobacillus|Gardnerella|Prevotella|Atopobium|Sneathia)$/) {
	    if ($species eq 'sp.') {
		$new_name = "$rank\_$genus";
	    }else{
		$new_name = "$genus\_$species";
	    }
	    
	    
	}else{
	    $new_name = "$rank\_$genus";
	}
	
	push(@valencia_sp, $new_name);
    }
    return (@valencia_sp);
}

sub extract_species_lineage{
    my $lineage = shift;
    my @ranks = split(/;/, $lineage);
    my @new_ranks;
    foreach my $i (@ranks){
	if ($i =~ /^(x__|y__|z__)/) {
	    next;
	}else{
	    push(@new_ranks, $i)
	}
	
    }
    my $new_lineage = join(";", @new_ranks);
    return ($new_lineage);
}

sub extract_species_name{
    my $lineage = shift;
    my @ranks = split(/;/, $lineage);
    my $genus='';
    my $species='';
    foreach my $i (reverse(@ranks)){
	if ($i =~ /^z__/) {
	    my $strain = substr($i, 3, length($i)-3);
	    if ($strain =~ /Human papillomavirus/) {
		return($strain);
	    }
	    next;
	    
	}
	
	if ($i =~ /^(x__|y__)/) {
	    next;
	}
	
	
	if ($i =~ /^s__/) {
	    $species = substr($i, 3, length($i)-3);
	    my @names = split(/\ /, $species);
	    if ($lineage =~ /^k__Viruses/) {
		return($species);
	    }
	    
	    if ($species =~ /sp\.|unculture|unknown| bacterium/) {
		$species = 'sp.';
	    }elsif($species =~ /^Candidatus/){
		shift(@names);
		shift(@names);
		$species = join(" ", @names);
	    }else{
		shift(@names);
		$species = join(" ", @names);
	    }
	    
	    
	}elsif ($i =~ /^g__/) {
	    $genus = substr($i, 3, length($i)-3);
	    if ($genus eq 'unknown') {
		$genus='';
	    }elsif($genus =~ /\-/){# if there is ambiguity in the genus assignment
		$genus='';
	    }else{
		return("$genus $species");
	    }
	    
	    
	}else{
	    $genus = substr($i, 3, length($i)-3);
	    my $rank = substr($i,0,1);
	    if ($genus eq 'unknown') {
		$genus = '';
	    }else{
		return("($rank)$genus $species");
	    }    
	}
    }   
}


sub copy_abun_files{
    my $input_dir = shift;
    my $output_dir = shift;
    my $rank = shift;
    
    my $rel_abun_file = "$input_dir/$rank.tsv";
    my $new_rel_abun_file = "$output_dir/$rank\_rel_abun.tsv";
    copy_file ($rel_abun_file, $new_rel_abun_file);
    
    my $abs_abun_file = "$input_dir/abun_table_read_counts.tsv";
    my $new_abs_abun_file = "$output_dir/$rank\_read_counts.tsv";
    copy_file ($abs_abun_file, $new_abs_abun_file);
}

sub copy_file{
    my $input_file = shift;
    my $output_file = shift;
    open(my $f1, "<$input_file") or die;
    open(my $f2, ">$output_file") or die;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line eq '# Constructed from biom file') {
	    next;
	}elsif($line =~ /^#OTU ID/){
	    my @coln = split(/\t/, $line);
	    $coln[0] = 'taxa';
	    my $new_line = join("\t", @coln);
	    print($f2 "$new_line\n");
	}else{
	    print($f2 "$line\n");
	}
	
    }
    close $f1;
    close $f2;
    
}

1;

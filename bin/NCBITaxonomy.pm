#!/usr/bin/perl
use strict;
use warnings;
use Parallel::Loops;
use Inputs;
package NCBITaxonomy;

sub interpret_centrifuge_results{
    my $hit_table=shift;
    my $nodes=shift;
    my $names=shift;
    my $phylogeny = shift;
    my $seq_type = shift;
    my %unique_hits;
    my @ambi_reads;
    open(my $f1, "<$hit_table") or die;
    my $title = <$f1>;
    my $id='';
    my %multiple_hits;
    my %lineage;
    while (my $line = <$f1>) {
        chomp $line;
        my @coln = split(/\t/, $line);
	if (scalar@coln==0) {next;}
	my $taxid=$coln[2];
	if ($taxid==0) {next;}
	my $hit_length_cutoff=80;# filter reads that have short_hit_length
	if ($seq_type =~ /illumina\.pe/) {# illumina paired
	    $hit_length_cutoff = 80;
	}elsif($seq_type =~ /illumina\.se/){#nanopore
	    $hit_length_cutoff = 80;
	}elsif($seq_type =~ /nanopore/){#nanopore
	    $hit_length_cutoff = 300;
	}elsif($seq_type =~ /pacbio/){#nanopore
	    $hit_length_cutoff = 300;
	}
	if ($coln[5]<$hit_length_cutoff) {# filter hits that have a hit_length
	    next;
	}
	# record the lineage information
	my $phylo_lineage;
	if (exists($lineage{$taxid})) {
	    $phylo_lineage = $lineage{$taxid};
	}else{
	    $phylo_lineage = get_lineage($taxid, $nodes, $names, $phylogeny);
	    if (scalar(@{$phylo_lineage})>0) { # skip the hits to return no lineage due to errors in the database
		$lineage{$taxid}=$phylo_lineage;
	    }else{
		next;
	    }
	    
	    if ($phylo_lineage->[5] =~/^g__Shigella$/) { # merge Shigella genus into Escherichia
		$phylo_lineage->[5] = 'g__Escherichia';
	    }
	    
	}
	# count useful reads, # not counting human reads as useful reads

	# enrol previous read infor if this line belongs to a new read.
	if ($id ne $coln[0]) {
	    my @hits = keys(%multiple_hits);
	    %multiple_hits=();
	    if (scalar(@hits)>1) {
		push (@ambi_reads, \@hits);
	    }elsif(scalar(@hits)==1){
		$unique_hits{$hits[0]}++;
	    }
	}
	    
	# iterate
	$id = $coln[0];
	$multiple_hits{$taxid}++;
    }
    my @hits = keys(%multiple_hits);
    if (scalar(@hits)>1) {
	push (@ambi_reads, \@hits);
    }elsif(scalar(@hits)==1){
	$unique_hits{$hits[0]}++;
    }
    # reassign reads that have multiple hits if the hits have survived so far.
    foreach my $i (@ambi_reads){
        my %unique_hits_shared;
        foreach my $j (@{$i}){
            if (exists($unique_hits{$j})) {
                $unique_hits_shared{$j}=$unique_hits{$j};
            }
        }
	if (keys%unique_hits_shared==0) { # if no hits can be found in the pool of unique hits, assign by common ancestor
	    assign_ambi_reads_by_common_ancestor($i,\%unique_hits, $nodes, $names, $phylogeny, \%lineage);
	}else{ # if hits can be found in the the pool of unique hits, assign randomly according to the abundance of unique hits.
	    assign_ambi_reads_by_abundance(\%unique_hits_shared,\%unique_hits, \%lineage);
	}        
    }
    
    #count useful microbial reads
    my $total_useful_reads=count_microbial_reads(\%unique_hits, \%lineage);
    
    
    #trim closely related strains by absolute counts and relative counts
    my $rel_abun_cutoff = 0.0001; # set this value to 0.00001 for QC of the standard for purity check.
    my $abs_cutoff = $total_useful_reads*$rel_abun_cutoff;
    if ($abs_cutoff<10) {$abs_cutoff=10;}
    my $abs_cutoff_at_strains = $abs_cutoff; # at strain level only keep the most abundance one
    my $rel_cutoff_between_species = 0.1;  # at species level, keep the one that has an abundance >10% abundance of the most abundant one.
    my $abs_cutoff_at_genus = $abs_cutoff;
    
    my %strain_size;
    my %taxa2taxid; # $taxa2taxid{$genus}{$species}{$strain}=();
    foreach my $i (keys%unique_hits){
	my $lineage = $lineage{$i};
	my $strain = $lineage->[9];
	my $species = $lineage->[6];
	my $genus = $lineage->[5];
	my $size = $unique_hits{$i};
	my $rel_abun = $size/$total_useful_reads;
	if ($genus =~ /unknown/) {
	    if ($rel_abun<0.002) {# remove low abundance unknown genus.
		next;
	    }
	}
	
	push(@{$taxa2taxid{$genus}{$species}{$strain}},$i);
	$strain_size{$genus}{$species}{$strain}+=$size;
    }
    
    #trim closely related strains
    my %species_size;
    foreach my $g (keys%strain_size){
	foreach my $s (keys%{$strain_size{$g}}){	    
	    my $strain_info = $strain_size{$g}{$s};
	    trim_strains($strain_info, $abs_cutoff_at_strains);
	    #from here, each species only contain one known strain with or without an 'unknown' strain
	    my @counts = values(%{$strain_size{$g}{$s}});
	    my $sum=0;
	    foreach my $i (@counts){
		$sum+= $i;
	    }
	    if ($sum==0) {
		delete($strain_size{$g}{$s});
	    }else{
		$species_size{$g}{$s}=$sum;
	    }

	}
    
    }
    #trim clostly related species
    foreach my $g (keys%strain_size){
	if ($g =~ /unknown/) { # skip the unknown genus because they don't belong to the same genus
	    next;
	}
        my $species_info = $species_size{$g};
	my $strain_info = $strain_size{$g};
	trim_species($species_info, $strain_info, $rel_cutoff_between_species);
        my @counts = values(%{$species_size{$g}});
        my $sum=0;
        foreach my $i (@counts){
            $sum+= $i;
        }
	if ($sum==0) {
	    delete($species_size{$g});
	    delete($strain_size{$g});
	}
	
        
    }
    my %survived_hits;
    foreach my $g(keys%strain_size){
	
        foreach my $s (keys%{$strain_size{$g}}){
            foreach my $z (keys%{$strain_size{$g}{$s}}){
                my $size = $strain_size{$g}{$s}{$z};
		my @taxid = sort{$unique_hits{$b}<=>$unique_hits{$a}}@{$taxa2taxid{$g}{$s}{$z}};
		$survived_hits{$taxid[0]}=$size;
            }
        }
    }
    
    my%tax_abun;
    foreach my $i(keys%survived_hits){
        my $abun = $survived_hits{$i};
        my $lineage = $lineage{$i};
        $tax_abun{$i}{'abun'}=$abun;
        $tax_abun{$i}{'lineage'}=$lineage;
    }

    my $taxa2abun = transform_abun_table(\%tax_abun, $phylogeny);

}

sub assign_ambi_reads_by_common_ancestor{
    my $ambi_reads=shift;
    my $unique_hits=shift;
    my $nodes=shift;
    my $names=shift;
    my $phylogeny=shift;
    my $rank_lineage=shift;
    my %node_lineage;
    foreach my $i (@{$ambi_reads}){
	my $lineage = get_lineage_nodes($i, $nodes);
	$node_lineage{$i}=$lineage;
    }
    my @hits = keys(%node_lineage);
    my $hit1 = pop@hits;
    my $lineage1 = $node_lineage{$hit1};
    #iterate through all nodes in the lineage to find the shared lowest nodes.
    for (my $i=0; $i<scalar(@{$lineage1}); $i++) {
	my $found=0;
	foreach my $j (@hits){
	    my $lineage2 = $node_lineage{$j};
	    foreach my $x (@{$lineage2}){
		if ($lineage1->[$i] eq $x) {
		    $found++;
		    last;
		}
		
	    }
	    
	}
	if ($found==scalar(@hits)) {# if the node_id is found in all of the hits, this means it is a common ancestor.
	    my $node_id = $lineage1->[$i];
	    my $taxa_lineage;
	    if (exists($rank_lineage->{$node_id})) {
		$taxa_lineage = $rank_lineage->{$node_id};
	    }else{
		$taxa_lineage = get_lineage($node_id, $nodes, $names, $phylogeny);
		$rank_lineage->{$node_id}=$taxa_lineage;
	    }
	    $unique_hits->{$node_id}++;
	    return;
	}
	
    }
    
    
}

sub assign_ambi_reads_by_abundance{
    my $ambi_hits=shift;
    my $unique_hits = shift;
    my @hits = keys%{$ambi_hits};
    if (scalar(@hits)==0) {
	return;
    }elsif(scalar(@hits)==1){
	$unique_hits->{$hits[0]}++;
    }else{
	#when there are more than one ambiguous hits, assign the read to one of them based on their existing abundance
	@hits = sort{$ambi_hits->{$a}<=>$ambi_hits->{$b}} @hits;
	my %levels;
	my $sum=0;
	foreach my $i (@hits){
	    $sum+=$ambi_hits->{$i};
	    $levels{$i}=$sum;
	}
	my $random = rand($sum);
	foreach my $i (@hits){
	    my $level = $levels{$i};
	    if ($random<=$level) {
		$unique_hits->{$i}++;
		last;
	    }
	    
	}
    }
    
    
}


sub trim_species{
    #at strain level, only keep the taxa at highest abundance.
    my $species_info = shift;
    my $strain_info = shift;
    my $rel_cutoff_between_species=shift;
    my $reads_deleted=0;
    # if there is more than one strain and one of them is unknown, delete the unknow
    if (scalar(keys%{$species_info})>1 && exists($species_info->{'s__unknown'})) {
	$reads_deleted+=$species_info->{'s__unknown'};
	delete($species_info->{'s__unknown'});
	delete($strain_info->{'s__unknown'});
	
    }
    
    my @species = sort{$species_info->{$b}<=>$species_info->{$a}} keys%{$species_info};
    if (scalar@species > 0) {
	my $species_of_highest_abun = shift(@species);
	my $abun_cutoff=$species_info->{$species_of_highest_abun}*$rel_cutoff_between_species;
	foreach my $i (@species){
	    my $abun = $species_info->{$i};
	    if ($abun<$abun_cutoff) {
		$reads_deleted+=$species_info->{$i};
	        delete($species_info->{$i});
		delete($strain_info->{$i});
	    }
    
	}
    }
    
    #add the deleted reads to existing species and strain
    my $reads_remained=0;
    foreach my $s (keys(%{$species_info})){
	    $reads_remained+=$species_info->{$s};
    }
    foreach my $s (keys%{$strain_info}){
	$species_info->{$s}+=int($reads_deleted*($species_info->{$s})/$reads_remained);
	foreach my $z (keys%{$strain_info->{$s}}){
	    $strain_info->{$s}->{$z} += int($reads_deleted*($strain_info->{$s}->{$z})/$reads_remained)
	}
    }

}

sub trim_strains{
    #at strain level, only keep the taxa at highest abundance.
    my $strain_info = shift;
    my $abs_cutoff = shift;
    my $reads_deleted=0;
    foreach my $i (keys%{$strain_info}){
	my $abun = $strain_info->{$i};
	if ($abun<$abs_cutoff) {
	    $reads_deleted+=$abun;
	    delete($strain_info->{$i});
	}
    }
    # if there is more than one strain and one of them is unknown, delete the unknow
    if (scalar(keys%{$strain_info})>1 and exists($strain_info->{'z__unknown'})) {
	$reads_deleted+=$strain_info->{'z__unknown'};
	delete($strain_info->{'z__unknown'});
    }
    my @taxa = sort{$strain_info->{$b}<=>$strain_info->{$a}} keys%{$strain_info};
    my $taxa_of_highest_abun='';
    #find the most abundant taxa if existing and delete others
    if (scalar@taxa>0) {
	$taxa_of_highest_abun = shift(@taxa);
	foreach my $i (@taxa){
	    $reads_deleted+=$strain_info->{$i};
	    delete($strain_info->{$i});
	}
    }
    # if the most abundant taxa exists, add the deleted reads to its abundance.
    if (length($taxa_of_highest_abun)>0) {
	$strain_info->{$taxa_of_highest_abun}+=$reads_deleted
    }
}


sub get_lineage_nodes{
    my $node_id =shift;
    my $nodes=shift;
    my @lineage_nodes;
    my $num=0;
    while ($node_id != 1 and $num<50) {
	push(@lineage_nodes, $node_id);
	my $parent_id = $nodes -> {$node_id} -> {'parent'};
	$node_id = $parent_id;
	$num++;
    }
    return (\@lineage_nodes);
}

sub get_lineage{
    my $node_id =shift;
    my $nodes=shift;
    my $names = shift;
    my $phylogeny = shift;
    my $base_name = $names->{$node_id}{'name'};
    my %rank2name;
    my @lineage;
    my @rank_lvl= @{$phylogeny -> {'ranks'}};
    my @rank_prefix= @{$phylogeny -> {'prefix'}};
    my $num=0;
    my $rank='';
    
    if (length($node_id)==0) {
	return(\@lineage);
    }
    
    while ($node_id != 1 and $num<50) {
	my $name='';
	if (exists($nodes->{$node_id}->{'rank'})) {
	    $rank=$nodes->{$node_id}->{'rank'};
	    $name=$names->{$node_id}->{'name'};
	}else{
	    return(\@lineage);
	}
	
	$rank2name{$rank}=$name;
	my $parent_id = $nodes -> {$node_id} -> {'parent'};
	$node_id = $parent_id;
	$num++;
    }
    
    if (exists($rank2name{'species'})) {
	extract_species_info($base_name, \%rank2name);
    }
    
    for (my $i=0; $i<scalar(@rank_prefix);$i++){
	my $id = $rank_lvl[$i];
	my $rank_id ="unknown";
	if (exists($rank2name{$id})) {
	    $rank_id=$rank2name{$id};
	}
	my $prefix = $rank_prefix[$i];
	my $new_id = $prefix.$rank_id;
	push(@lineage, $new_id);
    }
    
    
    return \@lineage;
    
}

sub transform_abun_table{
    my $abun_table=shift;
    my $phylogeny = shift;
    my $total_counts=0;
    my %taxa2abun;
    foreach my $node_id (keys%{$abun_table}){
	my $lineage = $abun_table -> {$node_id} -> {'lineage'};
	my $abun = $abun_table -> {$node_id} -> {'abun'};
	$total_counts += $abun;
	my @ranks = @{$lineage};
	my @rank_index=@{$phylogeny -> {'ranks_to_use'}};
	my @new_ranks;
	foreach my $i (@rank_index){
	    push (@new_ranks, $ranks[$i-1]);
	}
	my $new_lineage = join(";", @new_ranks);
	$abun_table -> {$node_id} -> {'lineage'} = $new_lineage;
	
    }
    foreach my $node_id (keys%{$abun_table}){
	my $lineage = $abun_table -> {$node_id} -> {'lineage'};
	my $abun = $abun_table -> {$node_id} -> {'abun'};
	#$abun = $abun/$total_counts;
	#$abun = sprintf("%.6f", $abun);
	$taxa2abun{$lineage}=$abun;
    }
    return (\%taxa2abun);
}

sub extract_species_info{
    my $base_name = shift;
    my $rank2name = shift;
    if (exists($rank2name->{'subspecies'})) {
	if ($base_name =~ / serotype | serovar /) {
	    my @name = split(/ serotype | serovar | str. | str | strain /, $base_name);
	    $rank2name->{'serotype'} = $name[1];
	}else{
	    $rank2name->{'serotype'} = 'unknown'
	}
	if ($base_name =~ / str. | str | strain /) {
	    my @name = split(/ str. | str | strain /, $base_name);
	    $rank2name->{'strain'} = $name[1];
	}else{
	    $rank2name->{'strain'} = 'unknown';
	}
	
    }elsif ($rank2name->{'species'} ne $base_name) {
	$rank2name->{'subspecies'} ='unknown'; # if subspecies are not recorded before hand it is unknown.
	if ($base_name =~ / serotype | serovar /) {
	    my @name = split(/ serotype | serovar | str. | str | strain /, $base_name);
	    $rank2name->{'serotype'} = $name[1];
	}else{
	    $rank2name->{'serotype'} = 'unknown'
	}
	if ($base_name =~ / str. | str | strain /) {
	    my @name = split(/ str. | str | strain /, $base_name);
	    $rank2name->{'strain'} = $name[1];
	}else{
	    $rank2name->{'strain'} = $base_name;
	}
    }else{
	$rank2name->{'subspecies'} ='unknown';
	$rank2name->{'serotype'} = 'unknown';
	if ($base_name =~ / sp. /) {
	    $rank2name->{'species'} = 'unknown'; 
	    my @name = split(/ sp. | str. | str | strain /, $base_name); # Dehalobacter sp. CF, Dehalobacter sp. str. CF, Limnodrilus sp. strain 12
	    $rank2name->{'strain'} = pop(@name);
	}else{
	    $rank2name->{'strain'} = 'unknown';
	}
	
    }
}
sub count_microbial_reads{
    my $unique_hits = shift;
    my $lineage = shift;
    my $sum=0;
    foreach my $taxid (keys%{$unique_hits}){
	my $taxa_lineage=$lineage->{$taxid};
	my $abun = $unique_hits->{$taxid};
	my $genus = $taxa_lineage->[5];
	if ($genus =~ /^g__[Hh]omo$|^g__[Mm]us$/) { #skip human hits
	    next;
	}else{
	    $sum+=$abun;
	}
	
    }
    return $sum;
}

1;
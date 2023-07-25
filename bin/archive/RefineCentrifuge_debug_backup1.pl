#!/usr/bin/perl
use strict;
use warnings;
#use Inputs;
use List::MoreUtils qw(any);
#package RefineCentrifuge;
#----------------------------
my $input_dir = 'Trimmomatic';
my $output_dir = 'Centrifuge';
#my $bin_path = $_[3];
my $threads = 4;
my $database_folder = '/home/stdell/Desktop/shotgun_pipeline/ubuntu/database/';
my $script_dir = '/home/stdell/Desktop/shotgun_pipeline/ubuntu/scripts/';
my $phylogeny = ranks_to_use($script_dir);
my $database = "$database_folder/zymo_centrifuge";
my $ncbi_tax_database = "$database_folder/NCBI_taxonomy";
my $min_hit_length = 80;
my ($nodes, $names)=read_ncbi_taxa_database($ncbi_tax_database);
my $R1 = "$input_dir/in709_1_R1.paired.fastq.gz";
my $R2 = "$input_dir/in709_1_R2.paired.fastq.gz";
my $seq_type = "illumina.pe";
my ($taxa2abun, $hits_info, $hit_fate) = refine_centrifuge_results("$output_dir", "in709_1_class.tsv", $nodes, $names, $phylogeny, $seq_type, $R1, $R2, $database, $threads);

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


sub read_ncbi_taxa_database{
    my $ncbi_database_folder=shift;
    my %nodes;
    open(my $NCBI_nodes, "<$ncbi_database_folder/zymo.tree") or die;
    while (my $line = <$NCBI_nodes>) {
        chomp $line;
        my @coln = split(/\t/, $line);
        if (scalar@coln>1) {
            $nodes{$coln[0]}{'parent'}=$coln[2];
	    $nodes{$coln[0]}{'rank'}=$coln[4];
        }
    }
    close $NCBI_nodes;

    my %names;
    open(my $NCBI_names, "<$ncbi_database_folder/zymo.name") or die;
    while (my $line = <$NCBI_names>) {
        chomp $line;
        my @coln = split(/\t/, $line);
        my $name_class = $coln[6];
        if (exists($names{$coln[0]}{'name'})) {
	    if ($name_class eq 'scientific name') {
		$names{$coln[0]}{'name_class'}=$coln[6];
		$names{$coln[0]}{'note'}=$coln[4];
		$names{$coln[0]}{'name'}=$coln[2];
	    }
	}else{
	    $names{$coln[0]}{'name_class'}=$coln[6];
	    $names{$coln[0]}{'note'}=$coln[4];
	    $names{$coln[0]}{'name'}=$coln[2];
	}
    }
    close $NCBI_names;
    return(\%nodes, \%names);
}

#--------------------
sub refine_centrifuge_results{
    my $output_dir=shift;
    my $hit_table=shift;
    my $nodes=shift;
    my $names=shift;
    my $phylogeny = shift;
    my $seq_type = shift;
    my $R1 = shift;
    my $R2 = shift;
    my $database = shift;
    my $threads = shift;
    my ($hit_info, $read_info, $lineage) = get_unique_hits("$output_dir/$hit_table", $seq_type, $nodes, $names, $phylogeny, $database); # hit1->{$read1}=1
    my $hit_fate;
    my %opt_hit_info=();
    
    # make a copy of R1 and R2 for subsequent analysis
    my $new_R1 = 'R1.fastq';
    my $new_R2 = 'R2.fastq';
    system("gunzip -c $R1 > $output_dir/$new_R1");
    system("gunzip -c $R2 > $output_dir/$new_R2");
    
    while (scalar(keys%{$hit_info})>0) {
	my $best_hit = get_best_hit($hit_info); # assuming the best hit has the most reads;
	print ("$best_hit\n");
	$hit_fate->{$best_hit}='best_hits';
	move_best_hit($best_hit, $hit_info, \%opt_hit_info, $read_info);
	my $query_reads = get_all_reads_from_hit_hash($hit_info);
	if (scalar(keys%{$query_reads})==0) {
	    last;
	}
	
	download_ref_genome($output_dir, $best_hit, $database);
	update_raw_reads($output_dir, $new_R1, $new_R2, $query_reads); #update the same file and only keep the reads for the subsequent read mapping
	read_mapping_with_bwa($output_dir, $new_R1, $new_R2, $threads, $database);
	my $reassigned_reads=parse_bam_file($output_dir);
	update_hit_info($reassigned_reads, $hit_info, $read_info, $best_hit, $hit_fate, \%opt_hit_info);
    }
    # redo mapping of all reads to the optimzed hit list.
    #match_by_read_mapping2(\%opt_hit_info, $new_R1, $new_R2);
    #filter_by_genome_coverage(\%opt_hit_info, $hit_fate);
    my $taxa2abun;
    foreach my $i (keys%opt_hit_info){
	my %reads = %{$opt_hit_info{$i}{'reads'}};
	my $taxid = $opt_hit_info{$i}{'taxid'};
	$taxa2abun->{$taxid}=scalar(keys%reads);
    }
    return($taxa2abun, \%opt_hit_info, $hit_fate);
}

sub update_hit_info{
    my $reassigned_reads=$_[0];
    my $hit_info = $_[1];
    my $read_info = $_[2];
    my $best_hit = $_[3];
    my $hit_fate = $_[4];
    my $opt_hit_info = $_[5];
    foreach my $i (keys%{$reassigned_reads}){
	my %old_hits = %{$read_info->{$i} ->{'hits'}};
	$opt_hit_info->{$best_hit}->{'reads'}->{$i}='strain';
	foreach my $j (keys%old_hits){
	    if ($j ne $best_hit) {
		delete($hit_info->{$j}->{'reads'}->{$i});
	    }

	}
	
	$read_info->{$i} ->{'hits'}->{$best_hit}='reassigned';
    }
    foreach my $i (keys%{$hit_info}){
	my %reads = %{$hit_info -> {$i} -> {'reads'}};
	if (scalar(keys%reads)==0) {
	    $hit_fate->{$i} = 'removed';
	    delete($hit_info -> {$i});
	}
	
    }
    
}



sub parse_bam_file{
    my $dir = shift;
    open(my $f1, "<$dir/results.sam") or die;
    my %reads;
    #my @match_codes = (77, 141);
    while (my $line = <$f1>) {
	chomp $line;
	if ($line !~ /^@/) {
	    my @coln = split(/\t/, $line);
	    if ($coln[1] != 77 && $coln[1] != 141) {
		$reads{$coln[0]}++;
	    }
	    
	    
	}
	
    }
   return(\%reads);
}

sub read_mapping_with_bwa{
    my $dir = shift;
    my $R1 = shift;
    my $R2 = shift;
    my $threads = shift;
    my $database = shift;
    system("bwa index -a bwtsw $dir/ref.fasta");
    system("bwa mem -P -t $threads $dir/ref.fasta $dir/$R1 $dir/$R2 > $dir/results.sam");

}


sub update_raw_reads{
    my $dir = shift;
    my $R1 = shift;
    my $R2 = shift;
    my $read_list = shift;
    my $R1_s = 'R1_s.fastq';
    my $R2_s = 'R2_s.fastq';
    open(my $f1, "<$dir/$R1") or die;
    open(my $f2, "<$dir/$R2") or die;
    open(my $f3, ">$dir/$R1_s") or die;
    open(my $f4, ">$dir/$R2_s") or die;
    while (my $R1_line1 = <$f1>) {
	my $R1_line2 = <$f1>; my $R1_line3 = <$f1>; my $R1_line4 = <$f1>;
	my $R2_line1 = <$f2>; my $R2_line2 = <$f2>; my $R2_line3 = <$f2>;my $R2_line4 = <$f2>;
	my @coln = split(/\ /, $R1_line1);
	my $read_id = substr($coln[0], 1, length($coln[0])-1);
	if (exists($read_list->{$read_id})) {
	    print($f3 "$R1_line1$R1_line2$R1_line3$R1_line4");
	    print($f4 "$R2_line1$R2_line2$R2_line3$R2_line4");
	}
	
    }
    close $f1;
    close $f2;
    close $f3;
    close $f4;
    system("rm -f $dir/$R1");
    system("rm -f $dir/$R2");
    system("mv $dir/$R1_s $dir/$R1");
    system("mv $dir/$R2_s $dir/$R2");
    
}

sub download_ref_genome{
    my $dir = shift;
    my $best_hit = shift;
    my $database = shift;
    if (length($best_hit)>0) {
	#system("aws s3 cp s3://zymo-files/WGS_Pipeline/shotgun_database/20200327_genomes/$best_hit $dir/ref.fasta.gz");
	#debug scripts
	system("cp /home/stdell/Desktop/shotgun_pipeline/ubuntu/database/zymo_centrifuge/fasta/$best_hit $dir/ref.fasta.gz");
	system("gunzip -c $dir/ref.fasta.gz > $dir/ref.fasta");
    }
}

sub get_best_hit{
    my $hit_info = shift;
    my %hit_no;
    foreach my $i (keys%{$hit_info}){
	$hit_no{$i} = scalar(keys(%{$hit_info->{$i}->{'reads'}}));
    }
    my @hits = sort{$hit_no{$b}<=> $hit_no{$a}}keys%hit_no;
    my $best_hit = shift(@hits);
    return($best_hit);
}

sub move_best_hit{
    my $best_hit=$_[0];
    my $hit_info=$_[1];
    my $opt_hit_info = $_[2];
    my $read_info=$_[3];
    %{$opt_hit_info->{$best_hit}}=%{$hit_info->{$best_hit}};
    delete($hit_info->{$best_hit});
    foreach my $i (keys%{$opt_hit_info->{$best_hit}->{'reads'}}){
	my %hits=%{$read_info->{$i}->{'hits'}};
	foreach my $j (keys%hits){
	    if ($j ne $best_hit) {
		delete($hit_info->{$j}->{'reads'}->{$i});
	    }
	    delete($read_info->{$i}->{'hits'}->{$j});
	}
    }
    
}

sub get_all_reads_from_hit_hash{
    my $hit_info = shift;
    my %reads;
    foreach my $i (keys%{$hit_info}){
	my %read_list = %{$hit_info->{$i}->{'reads'}};
	%reads = (%reads, %read_list);
    }
    return (\%reads);
}


sub get_unique_hits{
    my $hit_table = shift;
    my $seq_type = shift;
    my $nodes = shift;
    my $names = shift;
    my $phylogeny = shift;
    my $database = shift;
    my @ambi_reads;
    open(my $f1, "<$hit_table") or die;
    my $title = <$f1>;
    my $last_read='';
    my $read_index='';
    my %read_info;# read info
    my %lineage;
    my %hit_info; # list of reads uniquely assigned to each hit
    my @read_list;
    my %reads_above_strain;
    
    my $acc2genome_file = read_acc2genome_table($database);
    
    while (my $line = <$f1>) {
        chomp $line;
        my @coln = split(/\t/, $line);
	if (scalar@coln==0) {next;}
	my $taxid=$coln[2];
	my $acc_no = $coln[1];
	if ($taxid eq '703612') {
	    next;
	}
	if ($taxid==0) {next;}
	my $hit_length_cutoff=80;# filter reads that have short_hit_length
	if ($seq_type =~ /illumina\.pe/) {# illumina paired
	    $hit_length_cutoff = 80;
	}elsif($seq_type =~ /illumina\.se/){#nanopore
	    $hit_length_cutoff = 80;
	}elsif($seq_type =~ /nanopore/){#nanopore
	    $hit_length_cutoff = 80;
	}elsif($seq_type =~ /pacbio/){#nanopore
	    $hit_length_cutoff = 80;
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
	# enrol previous read infor if this line belongs to a new read.
	my $read_id = $coln[0];#use the array index as the read id
	my $genome_id = '';
	if ($coln[1] =~ /^zymo/) {
	    my @name = split(/s[0-9]/, $coln[1]);
	    $genome_id = $name[0].'.fasta.gz';
	}else{
	    if (exists($acc2genome_file->{$coln[1]})) {
		$genome_id = $acc2genome_file->{$coln[1]};
	    }else{
		print("error: the genome of $coln[1] was not found!\n");
	    }
	    
	}
	
	if ($coln[1] =~ /\.|^zymo/) {
	    $hit_info{$genome_id}{'reads'}{$read_id}='strain';
	    $hit_info{$genome_id}{'taxid'}=$coln[2];
	    $read_info{$read_id}{'size'}=$coln[6];
	    $read_info{$read_id}{'hits'}{$genome_id}='Centrifuge.unique';
	}else{
	    @{$reads_above_strain{$coln[2]}{$read_id}} = ($genome_id,$coln[6]);
	}
	$last_read = $coln[0];
	
    }
    
    #assign reads without strain-level resolution to existing strains.
    my @taxid_reassigned;
    foreach my $i (keys%hit_info){
	my $taxid = $hit_info{$i}{'taxid'};
	my $node_lineage = get_lineage_nodes($taxid, $nodes);
	foreach my $j (@{$node_lineage}){
	    if (exists($reads_above_strain{$j})) {
		push(@taxid_reassigned, $j);
		my %reads = %{$reads_above_strain{$j}};
		foreach my $m(keys%reads){
		    $hit_info{$i}{'reads'}{$m}=$reads{$m}->[0];
		    $read_info{$m}{'size'}=$reads{$m}->[1];
		    $read_info{$m}{'hits'}{$i}='Centrifuge_above_strains_reassigned';
		}
		
	    }
	    
	}
    }
    # delete the ones that have been reassigned
    foreach my $i (@taxid_reassigned){
	delete($reads_above_strain{$i});
    }
    # for reads of taxid that cannot be assigned to a existing strain.
    foreach my $i (keys%reads_above_strain){
	my %reads = %{$reads_above_strain{$i}};
	foreach my $m(keys%reads){
	    $hit_info{$i}{'reads'}{$m}=$reads{$m}->[0];
	    $hit_info{$i}{'taxid'}=$i;
	    $read_info{$m}{'size'}=$reads{$m}->[1];
	    $read_info{$m}{'hits'}{$i}='Centrifuge_cannot_assigned_to_known_strains';
	}
    }
       
    return (\%hit_info, \%read_info, \%lineage);
}

sub read_acc2genome_table {
    my $database = shift;
    my %acc2genome_file;
    open(my $f1, "<$database/acc2genome.tsv") or die;
    while (my $line = <$f1>) {
	my @coln = split(/\t/, $line);
	$acc2genome_file{$coln[0]}=$coln[1];
    }
    return \%acc2genome_file;
}

sub get_genome_id{
    my $acc_no = shift;
    my $database = shift;
    my $line = `grep $acc_no $database/acc2genome.tsv`;
    chomp $line;
    my @coln =split(/\t/, $line);
    my $genome_id = $coln[1];
    return $genome_id;
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


#1;
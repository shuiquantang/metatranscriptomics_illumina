#!/usr/bin/perl -I /home/stdell/Desktop/shotgun_pipeline/ubuntu/scripts/bin
use strict;
use warnings;
package StrainScan;
use Parallel::Loops;
use Data::Dumper;

#StrainScan::select_species($sample_info, $subset_size, $abun_table_folder, $output_dir);
sub select_species{
    my $subset_size = shift;
    my $abun_table_dir = shift;
    my $output_dir = shift;
    my $top_species_number = shift;
    my $reads_cutoff = shift;
    my $chosen_species = shift;
    my $trim_folder = shift;
    my $sp_wishlist = shift;
    my $StrainScan_ref_dict_file =shift;
    my %species2taxaID;
    my $strainscan_tmp_dir = 'StrainScan_tmp';
    if (-d $strainscan_tmp_dir) {
	
    }else{
	system("mkdir $strainscan_tmp_dir");
    }
    read_strain_scan_dict($StrainScan_ref_dict_file, \%species2taxaID, $strainscan_tmp_dir);
    
    my %rel_species_abun;
    read_abun_table("$abun_table_dir/all_abun_table.tsv", \%rel_species_abun);
    my %rel_domain_abun;
    read_host_DNA_table("$abun_table_dir/host_dna_abun.csv", \%rel_domain_abun);
    
    my %tot_reads;
    get_tot_reads(\%tot_reads, "$trim_folder/summary.tsv");
    foreach my $i (keys(%tot_reads)){
	my $tot_reads = $tot_reads{$i};
	if ($tot_reads>$subset_size) {
	    $tot_reads = $subset_size;
	}
	my $microbial_read_perc = $rel_domain_abun{$i}{'Microbial_reads'};
	my @species = sort {$rel_species_abun{$i}{$b} <=> $rel_species_abun{$i}{$a}} keys(%{$rel_species_abun{$i}});
	my $num=0;
	foreach my $j (@species){
	    if ($j !~ /^k__Bacteria/) {
		next;
	    }
	    if ($num == $top_species_number) {
		last;
	    }
	    my $rel_abun = $rel_species_abun{$i}{$j};
	    my $reads = $tot_reads*$microbial_read_perc*$rel_abun;
	    my $species_name = retrieve_species_name($j);
	    my $present = 0;
	    foreach my $m (@{$sp_wishlist}){
		if ($species_name =~ /^$m/) { # it can be search by species or genus, Escherichia coli or Escherichia
		    $present = 1;
		    last;
		}
		
	    }
	    if (exists($species2taxaID{$species_name})) { # the species is a bacteria or in the wishlist, and its reference genome in strainscan exists, it has enough reads and the chosen species is less than 10
		my $taxaID = $species2taxaID{$species_name};
		if ($present == 1) {
		    $chosen_species->{"$i\@$species_name\@$taxaID"}=$rel_abun;
		}elsif($reads>=$reads_cutoff){
		    $chosen_species->{"$i\@$species_name\@$taxaID"}=$rel_abun;
		    $num++;
		}
	    }
	    
	    
	    
	}
    }
}

sub get_tot_reads{
    my $tot_reads = shift;
    my $summary = shift;
    open(my $f1, "<$summary") or die;
    my $header = <$f1>;
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	$tot_reads->{$coln[0]}=$coln[1];
    }
    close $f1;
    
}

sub retrieve_species_name{
    my $lineage = shift;
    my @rank = split(/\;/, $lineage);
    my $species = 'NA';
    foreach my $i (@rank){
	if ($i =~ /^s__/) {
	    $species = substr($i, 3, length($i)-3);
	}
	
    }
    return($species);
}

sub read_host_DNA_table{
    my $file = shift;
    my $rel_domain = shift;
    open(my $f1, "<$file") or die;
    my @header;
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	
	if ($line =~ /^#/) {
	    @header = @coln;
	}else{
	    if (scalar(@coln) != scalar(@header)) {
		next;
	    }
	    for(my $i=1;$i<scalar(@header);$i++){
		$rel_domain->{$header[$i]}->{$coln[0]}+=$coln[$i];
	    }
	}
	
    }
}

sub read_abun_table {
    my $file = shift;
    my $rel_abun = shift;
    open(my $f1, "<$file") or die;
    my @header;
    my %tot_reads;
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	
	if ($line =~ /^#/) {
	    @header = @coln;
	}else{
	    if (scalar(@coln) != scalar(@header)) {
		next;
	    }
	    for(my $i=1;$i<scalar(@header);$i++){
		my $species = extract_species($coln[0]);
		$rel_abun->{$header[$i]}->{$species}+=$coln[$i];
		$tot_reads{$header[$i]}+=$coln[$i];
	    }
	}
	
    }
    foreach my $i (keys(%tot_reads)){
	my $tot = $tot_reads{$i};
	if ($tot>0) {
	    foreach my $j (keys(%{$rel_abun->{$i}})){
		my $rel = ($rel_abun->{$i}->{$j})/$tot;
		$rel = sprintf("%.5f", $rel);
		$rel_abun->{$i}->{$j} = $rel;
	    }
	}
	
	
    }
    
}

sub extract_species{
    my $lineage = shift;
    my @ranks = split(/;/, $lineage);
    my @new_ranks;
    foreach my $i (@ranks){
	push(@new_ranks, $i);
	if ($i =~ /^s__/) {
	    last;
	}
	
    }
    my $new_lineage = join(";", @new_ranks);
    return $new_lineage;
}

#StrainScan::scan_strains($trim_folder, $output_dir, \%chosen_species, $StrainScan_ref_dict_file);

sub scan_strains{
    my $fastq_dir = shift;
    my $output_dir = shift;
    my $chosen_species = shift;
    my $StrainScan_dict = shift;
    my $s3_ref = shift;
    my $strainscan_tmp_dir = 'StrainScan_tmp';
    
    my %unique_sp;
    my %samples;
    my @chosen_species = keys(%{$chosen_species});
    my %species2taxaID;
    foreach my $i (@chosen_species){
	my @ids = split(/\@/,$i);
	$unique_sp{$ids[1]}++;
	$samples{$ids[0]}++;
	$species2taxaID{$ids[1]}=$ids[2];
    }
    
    my @unique_sp = keys(%unique_sp);
    my $cpu = `nproc`;
    chomp $cpu;
    my $threads = int($cpu/3);
    my $pl = Parallel::Loops->new($threads);
    $pl -> foreach (\@unique_sp, sub{
    #foreach my $i (@unique_sp){
	my $i = $_;
	my $taxaID = $species2taxaID{$i};
	system("aws s3 sync $s3_ref/taxid_$taxaID $strainscan_tmp_dir/$taxaID");
    });
    #}
    my @samples = keys(%samples);
    $pl = Parallel::Loops->new($cpu);
    $pl -> foreach (\@samples, sub{
    #foreach my $i (@samples){
	my $i = $_;
	my $R1 = "$fastq_dir/$i\_R1.fastq";
	my $R2 = "$fastq_dir/$i\_R2.fastq";
	system("cat $R1 $R2 >$fastq_dir/$i.fastq");
    });
    #}
    $pl = Parallel::Loops->new($threads);
    my %strain_abun;
    my %all_strains;
    $pl -> share(\%all_strains, \%strain_abun);
    $pl -> foreach (\@chosen_species, sub{
    #foreach my $i (keys(%{$chosen_species})){
	my $i = $_;
	my @ids = split(/\@/, $i);
	my $sample_id = $ids[0];
	my $sp = $ids[1];
	my $rel_abun = $chosen_species->{$i};
	my $taxaID = $species2taxaID{$sp};
	# run StrainScan
	$sp =~ s/\ /\_/g;
	system("mkdir -p \"$output_dir/$sample_id/$sp\"");
	system("python /strainscan/StrainScan.py -i $fastq_dir/$sample_id.fastq -d $strainscan_tmp_dir/$taxaID -o \"$output_dir/$sample_id/$sp\"");
	record_strain_abun($sample_id, $rel_abun, "$output_dir/$sample_id/$sp/final_report.txt", \%strain_abun, $sp, \%all_strains);
    });
    #}
    #print Dumper(%strain_abun);
    #print Dumper(%all_strains);
    my %abun;
    foreach my $i (keys(%strain_abun)){
	my @ids = split(/\t/, $i);
	$abun{$ids[0]}{$ids[1]} = $strain_abun{$i};
    }
    %strain_abun = %abun;
    foreach my $i (@samples){
	system("rm $fastq_dir/$i.fastq");
	#system("rm $fastq_dir/$i\_*.fastq");
    }
    @samples = sort{$a cmp $b}keys(%strain_abun);
    my @strains = sort {$a cmp $b}keys(%all_strains);
    open(my $f1, ">$output_dir/strain_summary.tsv") or die;
    my @coln = @samples;
    unshift(@coln, "strain_id");
    my $header = join("\t", @coln);
    print($f1 "$header\n");
    foreach my $i (@strains){
	@coln = ($i);
	foreach my $j (@samples){
	    my $abun = 0;
	    if (exists($strain_abun{$j}{$i})) {
		$abun = $strain_abun{$j}{$i};
	    }else{
		print ("$j\t$i\t$abun\n");
	    }
	    push(@coln, $abun);
	    
	}
	my $line = join("\t", @coln);
	print($f1 "$line\n");
    }
    close $f1;
    
}

sub record_strain_abun{
    my $sample_id = shift;
    my $rel_abun = shift;
    my $report = shift;
    my $strain_abun = shift;
    my $species = shift;
    my $all_strains = shift;
    open(my $f1, "<$report") or do {print("$report does not exist!\n"); return};
    my $header = <$f1>;
    while (my $line = <$f1>) {
	chomp $line;
	my @coln  = split(/\t/, $line);
	my $strain = "$species;$coln[1]";
	my $abun = $rel_abun*$coln[3];
	#$abun = printf("%.6f", $abun);
	$strain_abun->{"$sample_id\t$strain"} = "$abun";
	$all_strains ->{$strain}++;
	#print("$strain\t$abun\n");
    }
    close $f1;
    
}


sub read_strain_scan_dict{
    my $file = shift;
    my $sp2taxaID = shift;
    my $tmp_dir = shift;
    system("aws s3 cp $file $tmp_dir/StrainScan_ref_dict.tsv");
    open(my $f1, "<$tmp_dir/StrainScan_ref_dict.tsv") or die;
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split(/\t|\,/, $line);
	if ($line =~ /^#/) {
	    next;
	}else{
	    if (scalar(@coln)>1) {
		$sp2taxaID ->{$coln[0]} = $coln[1];
	    }
	}
	
    }

    close $f1;
    
}

1;

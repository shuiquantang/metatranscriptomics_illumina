#!/usr/bin/perl -I /media/sf_VirtualBox/scripts/metaillu/scripts/bin

package AnalysisPreparation;


sub rename_sample_and_prepare_inputs{
    #$work_dir, $i, $analysis_groups, $analysis_group_info, $sample_info, \@taxa_domains
    my $work_dir = $_[0];
    my $output_dir = $_[1];
    my $group_id = $_[2];
    my $groups = $_[3];
    my $group_info = $_[4];
    my $sample_info = $_[5];
    my $domain_category = $_[6];
    my @sample_labels = @{$group_info -> {$group_id}->{'sample_order'}};
    my @internal_ids;
    my %internal_id_dict;
    foreach my $i (@sample_labels){
	my $internal_id = $groups->{$group_id}->{$i}->{'internal_id'};
	push (@internal_ids, $internal_id);
	$internal_id_dict{$internal_id}++;
    }
    # AMR
    my $AMR_parent_dir = "$work_dir/diamond/AMR";
    my $AMR_child_dir = "$group_id/AMR_diamond";
    #virulence factor
    move_diamond_outputs(\@internal_ids, $sample_info, $AMR_parent_dir, $AMR_child_dir);
    my $virulence_parent_dir = "$work_dir/diamond/virulence_genes";
    my $virulence_child_dir = "$group_id/virulence_diamond";
    move_diamond_outputs(\@internal_ids, $sample_info, $virulence_parent_dir, $virulence_child_dir);
    #centrifuge
    #my $centrifuge_parent_dir = "$work_dir/Centrifuge";
    #my $centrifuge_child_dir = "$group_id/Centrifuge";
    #move_folder_by_sample_id(\%internal_id_dict, $centrifuge_parent_dir, $centrifuge_child_dir);
    #FastQC
    my $fastqc_parent_dir = "$work_dir/read_processing/FastQC";
    my $fastqc_chld_dir = "$group_id/FastQC";
    move_folder_by_sample_id(\%internal_id_dict, $fastqc_parent_dir, $fastqc_chld_dir);
    #Trimmomatic
    my $trimmomatic_parent_dir = "$work_dir/read_processing/Trimmomatic";
    my $trimmomatic_child_dir = "$group_id/Trimmomatic";
    rename_sample_in_trim_summary(\@internal_ids, $sample_info, $trimmomatic_parent_dir, $trimmomatic_child_dir);
    #Humann
    if (-d "$work_dir/humann/") {
	my $humann_parent_dir = "$work_dir/humann/summary";
	my $humann_child_dir = "$group_id/humann";
	rename_sample_in_humann_outputs(\@internal_ids, $sample_info, $humann_parent_dir, $humann_child_dir);
    }
    
    
}


sub rename_sample_in_trim_summary{
    my $internal_ids = $_[0];
    my $sample_info = $_[1];
    my $parent_dir = $_[2];
    my $child_dir = $_[3];
    mkdir($child_dir);
    my %trim_info;
    unless (-e "$parent_dir/summary.tsv") {
	return;
    }
    
    open(my $f1, "<$parent_dir/summary.tsv") or die;
    open(my $f2, ">$child_dir/summary.tsv") or die;
    my $title = <$f1>;
    my @title = split(/\t/, $title);
    shift(@title);
    unshift (@title, "sample_lable");
    unshift (@title, "internal_id");
    $title = join("\t", @title);
    print($f2 "$title");
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	my $internal_id = shift(@coln);
	$trim_info{$internal_id}=\@coln;
	
    }
    foreach my $i (@{$internal_ids}){
	my $sample_label = $sample_info->{$i}->{'label'};
	if (exists($trim_info{$i})) {
	    my @coln = @{$trim_info{$i}};
	    unshift(@coln, $sample_label);
	    unshift(@coln, $i);
	    my $line = join("\t", @coln);
	    print($f2 "$line\n");
	}
	
    }
    close $f1;
    close $f2;
    
}

sub move_folder_by_sample_id{
    my $internal_id_dict = $_[0];
    my $parent_dir = $_[1];
    my $child_dir = $_[2];
    mkdir($child_dir);
    opendir(my $dir, $parent_dir) || die "can't opendir $parent_dir: $!";
    my @files = grep { !/^\./ } readdir($dir);
    foreach my $i (@files){
	if ($i =~ /_class.tsv/) { #skip some centrifuge files that are big
	    next;
	}
	if ($i =~ /_class_filtered.tsv/) {
	    next;
	}
	my @name = split(/\_/,$i);
	if (scalar(@name)>2) {
	    my $internal_id = $name[0].'_'.$name[1];
	    if ($internal_id_dict->{$internal_id}) {
		system("cp $parent_dir/$i $child_dir/$i");
	    }
	}
    }
}

sub move_diamond_outputs{
    my $internal_ids = $_[0];
    my $sample_info = $_[1];
    my $parent_dir = $_[2];
    my $child_dir = $_[3];
    mkdir($child_dir);
    #mkdir("$child_dir/raw_results");
    mkdir("$child_dir/summary");
    foreach my $i (@{$internal_ids}){
	my $sample_label = $sample_info->{$i}->{'label'};
	my $parent_file = "$i.summary.txt";
	my $child_file = "$sample_label.summary.txt";
	if (-e "$parent_dir/summary/$parent_file") {
	    system("cp $parent_dir/summary/$parent_file $child_dir/summary/$child_file");
	}
	
    }
}


sub rename_sample_in_humann_outputs{
    my $internal_ids = $_[0];
    my $sample_info = $_[1];
    my $parent_dir = $_[2];
    my $child_dir = $_[3];
    mkdir($child_dir);
    my @file_types = ("gene_fam_cpm", "pathway_abun_cpm", "pathway_cov",
		      "combined_pathway_abun_filt", "pathway_abun_filt",
		      "species_pathway_cov_filt", "gene_fam_cpm_filt", 
		      "species_gene_fam_cpm_filt", "pathway_cov_filt",
		      "species_pathway_abun_filt");
    
    foreach my $i (@file_types){
	my %abun_info;
	unless (-e "$parent_dir/$i.tsv") {
	    next;
	}
	open(my $f1, "<$parent_dir/$i.tsv") or die;
	open(my $f2, ">$child_dir/$i.tsv") or die;
	my $title = <$f1>;
	my @title = split(/\t/, $title);
	shift(@title);
	for (my $i=0; $i<scalar(@title); $i++){
	    my @ids = split(/\_/, $title[$i]);
	    $title[$i] =$ids[0].'_'.$ids[1];
	}
	my @tax_order;
	my %used_samples;
	while (my $line = <$f1>) {
	    chomp $line;
	    my @coln = split(/\t/, $line);
	    my $tax = shift(@coln);
	    push(@tax_order, $tax);
	    for (my $i=0; $i<scalar(@coln); $i++){
		my $internal_id = $title[$i];
		my $abun = $coln[$i];
		$abun_info{$tax}{$internal_id}=$abun;
		$used_samples{$internal_id}++;
	    }
	    
	}
	
	my @new_internal_id;
	foreach my $i (@{$internal_ids}){
	    if (exists($used_samples{$i})) {
		push(@new_internal_id, $i);
	    }
	    
	}
	
	my @new_title=('#OTU ID');
	foreach my $i (@new_internal_id){
	    my $sample_label = $sample_info->{$i}->{'label'};
	    push(@new_title, $sample_label);
	}
	my $new_header = join("\t",@new_title);
	print($f2 "# for formatting consistency\n");
	print($f2 "$new_header\n");
	foreach my $x (@tax_order){
	    my @coln = ($x);
	    foreach my $y (@new_internal_id){
		my $abun = $abun_info{$x}{$y};
		push(@coln, $abun);
	    }
	    my $line = join("\t", @coln);
	    print($f2 "$line\n");
	}
	close $f1;
	close $f2;
    }

}




1;
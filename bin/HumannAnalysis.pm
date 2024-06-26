#!/usr/bin/perl -I /home/stdell/Desktop/shotgun_pipeline/ubuntu/scripts/bin
use strict;
use warnings;
#use NCBITaxonomy;
use Parallel::Loops;
package HumannAnalysis;
use Inputs;

sub humann3 {
    #$trim_folder, $humann_folder, $sample_info, $threads, $shotgun_database_folder
    my $input_dir = $_[0];
    my $output_dir = $_[1];
    my $sample_info = $_[2];
    my $threads = $_[3];
    my $subset_size = $_[4]*2*4; #convert paired-end read number to single reads
    my $cpu = `nproc`;
    chomp $cpu;
    # run humann
    my @samples = sort(keys%{$sample_info});
    my $pl = Parallel::Loops->new($cpu);
    $pl -> foreach (\@samples, sub{
        my $i = $_;
    #$foreach my $i (@samples){
        system("mkdir -p $input_dir/subset");
	my $input_file = "$input_dir/$i.fastq.gz";
	my $output_file = "$input_dir/subset/$i.fastq.gz";
	system("gunzip -c $input_file | head -$subset_size | gzip > $output_file");
	#system("rm $input_file");
    #}
    });
    
    
    $pl = Parallel::Loops->new($threads);
    $pl -> foreach (\@samples, sub{
        my $i = $_;
    #foreach my $i (@samples){
	my $tag='--bypass-translated-search';
	system("mkdir $output_dir/$i");
	system("touch $output_dir/$i.log");
	#my $databases = "--nucleotide-database /home/ubuntu/metaillu_database/humann3/chocophlan --protein-database /home/ubuntu/metaillu_database/humann3/uniref90";
	my $cmd = "humann -i $input_dir/subset/$i.fastq.gz -o $output_dir/$i --threads $cpu $tag --input-format fastq.gz 2>&1 >>$output_dir/$i.log";
	Inputs::print_and_execute($cmd, "$output_dir/$i.log"); 
	$cmd = "humann_renorm_table --input $output_dir/$i/$i\_pathabundance.tsv --output $output_dir/$i\_pathway_abun_cpm.tsv --units cpm --update-snames 2>&1 >>$output_dir/$i.log";
	Inputs::print_and_execute($cmd, "$output_dir/$i.log"); 
	$cmd = "humann_renorm_table --input $output_dir/$i/$i\_genefamilies.tsv --output $output_dir/$i\_gene_fam_cpm.tsv --units cpm --update-snames 2>&1 >> $output_dir/$i.log";
	Inputs::print_and_execute($cmd, "$output_dir/$i.log"); 
	$cmd = "mv $output_dir/$i/$i\_pathcoverage.tsv $output_dir/$i\_pathway_cov.tsv";
	Inputs::print_and_execute($cmd, "$output_dir/$i.log");
	system("rm -r $output_dir/$i");
    #}
    });
    system("rm $input_dir/subset");
    my @file_types = ("gene_fam_cpm", "pathway_abun_cpm", "pathway_cov");
    system("mkdir $output_dir/summary");
    foreach my $i (@file_types){
	my $cmd = "humann_join_tables --file_name $i -i $output_dir -o $output_dir/summary/$i.tsv";#  > /dev/null 2>&1";
	system($cmd);
    }
    system("rm -r $output_dir/*.tsv"); # remove intermediate files
    #fitler to remove unmapped and unintegrated abundance
    filter_break_gene_family_table("$output_dir/summary/gene_fam_cpm.tsv", "$output_dir/summary/gene_fam_cpm_filt.tsv", "$output_dir/summary/species_gene_fam_cpm_filt.tsv");
    #filter and separate the pathways from pathways+microbes
    my @pathway_file_names = ("pathway_abun_filt.tsv", "species_pathway_abun_filt.tsv", "pathway_cov_filt.tsv", "species_pathway_cov_filt.tsv", "combined_pathway_abun_filt.tsv");
    my $dir = "$output_dir/summary";
    filter_break_pathway_table($dir, "pathway_abun_cpm.tsv", "pathway_cov.tsv", \@pathway_file_names);
}

sub filter_break_gene_family_table{
    my $old_table = shift;
    my $new_table = shift;
    my $species_table = shift;
    open(my $f1, "<$old_table") or return();
    open(my $f2, ">$new_table") or die;
    open(my $f3, ">$species_table") or die;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line =~ /^#/) {
	    print ($f2 "$line\n");
	    print ($f3 "$line\n");
	}else{
	    my @coln = split(/\t/, $line);
	    my @ids = split(/\|/, $coln[0]);
	    if ($ids[0] =~ /\_unknown$|^UNMAPPED$/) {
		
	    }else{
		if (scalar(@ids)==1) {
		    print($f2 "$line\n");
		}elsif(scalar(@ids)>1){
		    if ($ids[1] =~ /^g__|^unclassified$/) {
			print($f3 "$line\n");
		    }else{
			print ("Formatting errors in gene family output files:\n$line\n");
		    }
		    
		    
		}
		
		
	    }
	    
	}
	
	
    }
    close $f1;
    close $f2;
    close $f3;
    
}

sub filter_break_pathway_table{
    my $dir = $_[0];
    my $pathway_abun_file = $_[1];
    my $pathway_cov_file = $_[2];
    my @new_files =@{$_[3]};
    open(my $abun_file, "<$dir/$pathway_abun_file") or return();
    open(my $cov_file, "<$dir/$pathway_cov_file") or return();
    open(my $f1, ">$dir/$new_files[0]") or die;
    open(my $f2, ">$dir/$new_files[1]") or die;
    open(my $f3, ">$dir/$new_files[2]") or die;
    open(my $f4, ">$dir/$new_files[3]") or die;
    open(my $f5, ">$dir/$new_files[4]") or die;
    while (my $line1 = <$abun_file>) {
	my $line2 = <$cov_file>;
	chomp $line1;
	chomp $line2;
	if ($line1 =~ /^#/) {
	    print ($f1 "$line1\n");
	    print ($f2 "$line1\n");
	    print ($f3 "$line2\n");
	    print ($f4 "$line2\n");
	    print ($f5 "$line2\n");
	}else{
	    my @coln = split(/\t/, $line1);
	    my @ids = split(/\|/, $coln[0]);
	    if ($ids[0] =~ /^UNMAPPED$|^UNINTEGRATED$/) {
		next;
	    }
	    
	    print ($f5 "$line1\n");
	    
	    if (scalar(@ids)==1){
		print ($f1 "$line1\n");
		print ($f3 "$line2\n");
	    }elsif(scalar(@ids)>1){
		if ($ids[1]=~/^g__|^unclassified$/) {
		    print ($f2 "$line1\n");
		    print ($f4 "$line2\n");
		}else{
		    print ("Formatting errors in pathway output files:\n$line1\n");
		}
		
	    }
	    
	}
	
	
    }
    close $f1;
    close $f2;
    close $f3;
    close $f4;
    close $f5;
    close $abun_file;
    close $cov_file;
    
}

1;


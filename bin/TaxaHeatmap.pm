#!/usr/bin/perl -I /home/stdell/Desktop/Dada2.pipeline/Bin/
use strict;
use warnings;
use Inputs;

package TaxaHeatmap;

sub plot_heatmap_by_category{
    my $abun_table_folder = $_[0];
    my $category_number = $_[1];
    my $cat_titles = $_[2];
    my $bin_path = $_[3];
    my $heatmap_folder = $_[4];
    my $phylogeny = $_[5];
    
    if ($category_number>0) {
        for (my $j = 0; $j<$category_number; $j++){
            my $cat = $cat_titles->[$j];
            my $map_file = "$cat.mapping.file.txt";
            my $heatmap_dir = "$heatmap_folder\_$cat";
	    mkdir($heatmap_dir);
            plot_taxa_heatmap($abun_table_folder, $bin_path, $map_file, $heatmap_dir, $phylogeny);
        }
    }else{
        my $map_file = "qiime.mapping.file.txt";
        my $heatmap_dir = $heatmap_folder;
	mkdir($heatmap_dir);
        plot_taxa_heatmap($abun_table_folder, $bin_path, $map_file, $heatmap_dir, $phylogeny);
    }
}


sub plot_taxa_heatmap{
    my $abun_table_folder= $_[0];
    my $bin_path = $_[1];
    my $map_file = $_[2];
    my $heatmap_dir = $_[3];
    my $phylogeny = $_[4];
    my @ranks= @{$phylogeny->{'ranks'}};
    my @rank_prefix=@{$phylogeny->{'prefix'}};
    my @index = @{$phylogeny -> {'ranks_to_use'}};
    my @ranks_to_use;
    foreach my $i (@index){
	push (@ranks_to_use, $ranks[$i-1]);
    }
    shift(@ranks_to_use);
    foreach my $i (@ranks_to_use){
	system("mkdir $heatmap_dir/$i");
	my $heatmap_script = "$bin_path/draw_heatmap.py";
	my $output1 = "$heatmap_dir/$i/heatmap_with_sample_clustering.pdf";
	my $output2 = "$heatmap_dir/$i/heatmap_without_sample_clustering.pdf";
	my $old_file = "$abun_table_folder/$i/$i.tsv";
	my $new_file = "$heatmap_dir/$i/new_abun_table.tsv";
	open(my $old, "<$old_file") or die ("$old_file does not exist\n");
	open(my $new, ">$new_file") or die;
	my $title = <$old>;
	$title = <$old>;
	my @title = split(/\t/, $title);
	$title[0] = 'Taxon';
	$title = join("\t", @title);
	print ($new "$title");
	while (my $line = <$old>) {
	    chomp $line;
	    my @line = split(/\t/, $line);
	    my @lineage = split(/;/, $line[0]);
	    my @new_lineage;
	    for (my $j = scalar@lineage; $j>0; $j--){
		my $taxa = $lineage[$j-1];
		if ($taxa =~ /__$/) {
		        
		}elsif($taxa =~ /__unknown$/){
		        
		}elsif (lc($taxa) =~ /s__/){
		    unshift(@new_lineage, $taxa);
		}else{
		    unshift(@new_lineage, $taxa);
		    last;
		}
	    }
	    my $label=join(";", @new_lineage);
	    if ($label eq 'None') {
		$label = 'unknown'
	    }
	    
	    $line[0] = $label;
	    $line = join("\t", @line);
	    print ($new "$line\n");
	    
	}
	close $new;
	close $old;
	my $command = "python $heatmap_script -i $new_file -m $map_file -a $output1 -b $output2 -t 50 -r True >/dev/null 2>&1";
	system($command); 
    }
}


1;

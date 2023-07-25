#!/usr/bin/perl -I /home/stdell/Desktop/Dada2.pipeline/Bin/
use strict;
use warnings;
use Inputs;

package FuncHeatmap;

sub plot_heatmap_by_category{
    my $abun_file = $_[0];
    my $category_number = $_[1];
    my $cat_titles = $_[2];
    my $bin_path = $_[3];
    my $tag = $_[4];
   
    if ($category_number>0) {
        for (my $j = 0; $j<$category_number; $j++){
            my $cat = $cat_titles->[$j];
            my $map_file = "$cat.mapping.file.txt";
            my $heatmap_dir = "heatmap_$tag\_$cat";
	    mkdir($heatmap_dir);
            plot_taxa_heatmap($abun_file, $bin_path, $map_file, $heatmap_dir);
        }
    }else{
        my $map_file = "qiime.mapping.file.txt";
        my $heatmap_dir = "heatmap_$tag";
	mkdir($heatmap_dir);
        plot_taxa_heatmap($abun_file, $bin_path, $map_file, $heatmap_dir);
    }
}


sub plot_taxa_heatmap{
    my $old_file= $_[0];
    my $bin_path = $_[1];
    my $map_file = $_[2];
    my $heatmap_dir = $_[3];
    my $top_rows = 75;
    system("mkdir $heatmap_dir");
    my $heatmap_script = "$bin_path/draw_heatmap.py";
    my $output1 = "$heatmap_dir/heatmap_with_sample_clustering.pdf";
    my $output2 = "$heatmap_dir/heatmap_without_sample_clustering.pdf";
    my $new_file = "$heatmap_dir/new_abun_table.tsv";
    open(my $old, "<$old_file") or die ("$old_file does not exist\n");
    open(my $new, ">$new_file") or die;
    my $title = <$old>;
    $title = <$old>;
    my @title = split(/\t/, $title);
    my %abun;
    while (my $line = <$old>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	if (scalar@coln>1) {
	    my $sum=0;
	    for (my $i=1; $i<scalar(@title); $i++){
		$abun{$coln[0]}{$title[$i]}=$coln[$i];		
		$sum+=$coln[$i];
	    }
	    $abun{$coln[0]}{'#sum#'}=$sum;
	}
	
    }
    
    #get the top 50 taxa
    my @taxa = sort{$abun{$b}{'#sum#'}<=>$abun{$a}{'#sum#'}}keys%abun;
    if (scalar(@taxa)>=$top_rows) {
	@taxa = @taxa[0..$top_rows-1];
    }
    
    #calculate the sum of the abundance of all taxa in a sample
    my %sample_abun_sum;
    shift(@title);
    foreach my $i (@title){
	foreach my $j (@taxa){
	    if (exists($abun{$j}{$i})) {
		if ($abun{$j}{$i}) {
		    $sample_abun_sum{$i}+=$abun{$j}{$i};
		}else{
		    $sample_abun_sum{$i}+=0;
		}
	    }else{
		$sample_abun_sum{$i}+=0;
	    }
	}
    }
    
    # if the total abundance of a sample equals to 0, assign a small non-zero value to it.
    foreach my $i (@title){
	if ($sample_abun_sum{$i}==0) {
	    foreach my $j (@taxa){
		$abun{$j}{$i}=0.000000001;
	    }
	}
	
    }
    
    unshift(@title, 'Taxon');
    my $new_title = join("\t", @title);
    print($new "$new_title\n");
    shift(@title);
    foreach my $i (@taxa){
	my @row = ($i);
	foreach my $j (@title){
	    push(@row, $abun{$i}{$j})
	}
	my $new_line = join("\t", @row);
	print($new "$new_line\n");
    }    
    close $new;
    close $old;
    my $command = "python $heatmap_script -i $new_file -m $map_file -a $output1 -b $output2 -t $top_rows -r True >/dev/null 2>&1";
    system($command); 
}


1;

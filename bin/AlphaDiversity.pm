#!/usr/bin/perl -I /home/stdell/Desktop/Dada2.pipeline/Bin/
use strict;
use warnings;
use Inputs;

package AlphaDiversity;

sub alpha_diversity{
    my $abun_table_folder = $_[0];
    my $groups = $_[1];
    my $category_number = $_[2];
    my $cat_titles = $_[3];
    my $bin_path = $_[4];
    my $alpha_diversity_folder = $_[5];
    my $phylogeny = $_[6];
    my $group_id = $_[7];
    
    if ($category_number>0) {
        for (my $j = 0; $j<$category_number; $j++){
            my $cat = $cat_titles->[$j];
            my $alpha_div_dir = "$alpha_diversity_folder\_$cat";
	    mkdir($alpha_div_dir);
            alpha_div_plot_by_cat($abun_table_folder, $bin_path, $groups, $j, $alpha_div_dir, $phylogeny, $group_id);
        }
    }else{
        my $alpha_div_dir = $alpha_diversity_folder;
	mkdir($alpha_div_dir);
        alpha_div_plot_by_cat($abun_table_folder, $bin_path, $groups, '', $alpha_div_dir, $phylogeny, $group_id);
    }
}

sub alpha_div_plot_by_cat{
    my $abun_table_folder=shift;
    my $bin_path = shift;
    my $groups = shift;
    my $category_index = shift;
    my $alpha_div_dir = shift;
    my $phylogeny = shift;
    my $group_id = shift;
    my @ranks= @{$phylogeny->{'ranks'}};
    my @rank_prefix=@{$phylogeny->{'prefix'}};
    my @index = (6,7,10);
    my @ranks_to_use;
    foreach my $i (@index){
	push (@ranks_to_use, $ranks[$i-1]);
    }
    foreach my $i (@ranks_to_use){
	system("mkdir $alpha_div_dir/$i");
	my $abun_table_file = "$abun_table_folder/$i/$i.tsv";
	my ($abun_table, $total_abun, $sample_order)=read_abun_table($abun_table_file);
	my $alpha_div_rawdata_file = "$alpha_div_dir/$i/cat2shannon.tsv";
	my $alpha_div_cat_pdf = "$alpha_div_dir/$i/group_shannon.pdf";
	my $alpha_div_cat_png = "$alpha_div_dir/$i/group_shannon.png";
	my $alpha_div_pdf = "$alpha_div_dir/$i/sample_shannon.pdf";
	my $alpha_div_png = "$alpha_div_dir/$i/sample_shannon.png";
	open(my $f1, ">$alpha_div_rawdata_file") or die;
	print ($f1 "sample_id\tshannon_index\tcategory\n");
	foreach my $j (@{$sample_order}){
	    my $total_counts = $total_abun->{$j};
	    my $shannon_index = shannon_index_calculator(\%{$abun_table->{$j}}, $total_counts);
	    if (length($category_index)>0) {
		my $cat = $groups -> {$group_id} -> {$j} -> {'categories'} -> [$category_index];
		print ($f1 "$j\t$shannon_index\t$cat\n");
	    }else{
		print ($f1 "$j\t$shannon_index\t\n");
	    }
	    
	}
	close $f1;
	my $command = "Rscript $bin_path/alpha.diversity.plot.r -i $alpha_div_rawdata_file -a $alpha_div_cat_pdf -b $alpha_div_cat_png -c $alpha_div_pdf -d $alpha_div_png -t shannon_index_$i >/dev/null 2>&1";
	system($command);
	
    }
    
    
}

sub read_abun_table{
    my $abun_table_file = shift;
    open(my $f1, "<$abun_table_file") or die;
    my @sample_order;
    my %abun_table;
    my %total_abun;
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	if (scalar(@coln)<2) {
	    next;
	}
	if ($line=~ /^#/) {
	    @sample_order = @coln;
	}else{
	    for(my $i=1;$i<scalar(@coln);$i++){
		$abun_table{$sample_order[$i]}{$coln[0]}=$coln[$i];
		$total_abun{$sample_order[$i]}+=$coln[$i];
	    }
	}
	
    }
    shift(@sample_order);
    return(\%abun_table, , \%total_abun, \@sample_order);
    
}

sub shannon_index_calculator{
    my $taxa_abun = shift;
    my $total_counts = shift;
    my $shannon_index = 0;
    if ($total_counts == 0) {
	return(0);
    }
    
    foreach my $i (values(%{$taxa_abun})){
	my $rel_abun = $i/$total_counts;
	if ($rel_abun == 0) {
	    return (0);
	    next;
	}
	my $x = -$rel_abun*log($rel_abun);
	$shannon_index+=$x;
    }
    return $shannon_index;
}

1;

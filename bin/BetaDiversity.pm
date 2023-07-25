#!/usr/bin/perl -I /home/stdell/Desktop/Dada2.pipeline/Bin/
use strict;
use warnings;

package BetaDiversity;

sub beta_diversity{
    my $abun_table_folder = $_[0];
    my $groups = $_[1];
    my $category_number = $_[2];
    my $cat_titles = $_[3];
    my $bin_path = $_[4];
    my $beta_diversity_folder = $_[5];
    my $phylogeny=$_[6];
    my $group_info = $_[7];
    my $group_id = $_[8];
    my @ranks= @{$phylogeny->{'ranks'}};
    my @rank_prefix=@{$phylogeny->{'prefix'}};
    my @index = (6,7,10);
    my @ranks_to_use;
    my $sample_no = scalar(@{$group_info -> {$group_id} -> {'sample_order'}});
    # if there are less than 6 samples, skip beta-diversity analysis
    if ($sample_no<5) {
	print("Less than 5 samples -> skipped the beta-diversity analysis\n");
	return;
    }
    
        
    foreach my $i (@index){
	push (@ranks_to_use, $ranks[$i-1]);
    }
    mkdir($beta_diversity_folder);
    foreach my $i (@ranks_to_use){
	
	my $biom_file = "$abun_table_folder/$i/abun_table.biom";
	my $abun_table = "$abun_table_folder/$i/$i.tsv";
	my $taxa_no = get_taxa_no($abun_table);
	if ($taxa_no < 5) {
	    print("Less than 5 taxa -> skipped the beta-diversity analysis\n");
	    return;
	}
	system("mkdir $beta_diversity_folder/$i");
	my $map_file = "qiime.mapping.file.txt";
	my $command = "beta_diversity.py -i $biom_file -o $beta_diversity_folder/$i/dist --metrics binary_jaccard,bray_curtis >/dev/null 2>&1";
	system($command);
	$command = "principal_coordinates.py -i $beta_diversity_folder/$i/dist -o $beta_diversity_folder/$i/cordinates >/dev/null 2>&1";
	system($command);
	# might need to modify the format of tha taxa abundance table
	$command = "make_emperor.py -i $beta_diversity_folder/$i/cordinates/pcoa_bray_curtis_abun_table.txt -o $beta_diversity_folder/$i/biplot -m $map_file -t $abun_table >/dev/null 2>&1";
	system($command);
	
    }
}

sub get_taxa_no{
    my $abun_table = shift;
    open(my $f1, "<$abun_table") or die;
    my $taxa_no=0;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line =~ /^#/) {
	    next;
	}else{
	    $taxa_no++;
	}
	
    }
    return ($taxa_no);
}

1;
